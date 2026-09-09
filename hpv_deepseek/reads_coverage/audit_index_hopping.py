#!/usr/bin/env python3
"""
audit_index_hopping.py

Index-hopping detection: cross-sample floor/SNR audit and clean-BAM output.

Reads every sample's UMI whitelist report (produced per-sample by
whitelist_index_hopping_umis.py, which classifies each read's RX UMI tag against the KAPA
catalog by Levenshtein distance -- see that script's module docstring). Rows are grouped
across ALL samples in a run by physical signature (contig, position, mate position, template
length) plus canonical UMI. Two independent samples producing a read with an identical
signature AND UMI is, for all practical purposes, the same original DNA molecule read twice
under two different sample indices -- i.e. index hopping -- rather than coincidence.

For each such cross-sample group, the sample with the most total reads for that contig
(genotype) is treated as the true source; every other sample's reads in the group are audited:
below --min-floor total reads for that genotype, or below --snr-threshold of the source
sample's total, they are blacklisted as index-hopping artifacts; otherwise they are rescued
(kept) as a plausible low-level co-infection. Blacklisted read (query) names are then stripped
from that sample's BAM (both mates) to produce a "*_clean.bam" alongside a per-run report.

DSC note: fgbio's CallDuplexConsensusReads is not independently confirmed here to leave the
RX:Z: tag in the same "<read1 UMI>-<read2 UMI>" shape as the pre-consensus BAM (see
whitelist_index_hopping_umis.py's module docstring). This script does not depend on that
shape directly -- it consumes whitelist_index_hopping_umis.py's already-classified TSV rows,
whatever RX format produced them -- but SSC and DSC signatures/qnames are never joined against
each other (GroupReadsByUmi's adjacency vs. paired strategies produce unrelated qname spaces),
so SSC and DSC are always audited as fully independent runs of this script.

Author: Samuli Eldfors
"""

from __future__ import annotations

__version__ = "1.0.0"

import argparse
import csv
import shutil
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pysam

DEFAULT_MIN_FLOOR = 10
DEFAULT_HOP_RATE = 0.005    
DEFAULT_HOP_MARGIN = 2      # require this multiple of the expected hop load to rescue
DEFAULT_SNR_THRESHOLD = 0.01
# HPV call rule, from reads_coverage_with_ihop_v16_04_02_2029.py. A genotype is called
# only if it clears BOTH criteria. v10-v15 implemented the read count alone and never
# compared the coverage figure to anything even though it was computed and printed.
DEFAULT_MIN_READS = 10
DEFAULT_MIN_COV = 10.0
DEFAULT_GENOME_LEN = 7906   # HPV16 reference length, for the obs/exp breadth ratio only
DEFAULT_TEMPLATE_LEN = 170  # approximate template length, same purpose
JOINABLE_STATUSES = ("exact", "near")

REPORT_COLUMNS = (
    "sample_id",
    "bam_kind",
    "umi_matched_reads",
    "phase1_suspicious",
    "phase2_rescued",
    "phase2_confirmed_bad",
    "final_reads",
    "primary_leaked_source",
    "hpv_status_pre_ihop",
    "hpv_status_post_ihop",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Cross-sample index-hopping audit: join UMI-whitelisted reads across all samples "
            "in a run by physical signature + canonical UMI, apply a floor/SNR audit, and "
            "write clean (decontaminated) BAMs and a per-run report."
        )
    )
    parser.add_argument("--manifest", required=True, help="Whole-run sample/BAM/whitelist-report manifest TSV.")
    parser.add_argument("--bam-kind", required=True, choices=("SSC", "DSC"), help="Which BAM tier to audit.")
    parser.add_argument("--output-dir", required=True, help="Directory for clean BAMs and the final report.")
    parser.add_argument("--min-floor", type=int, default=DEFAULT_MIN_FLOOR, help=f"Minimum genotype read total below which a non-source sample's reads are always blacklisted. Default: {DEFAULT_MIN_FLOOR}.")
    parser.add_argument("--min-reads", type=int, default=DEFAULT_MIN_READS, help=f"HPV call: minimum primary reads of one genotype. Default: {DEFAULT_MIN_READS}. INHERITED FROM v16 AND NOT RE-DERIVED -- v16 set it on five known-status samples that landed at 8-14 reads.")
    parser.add_argument("--min-cov", type=float, default=DEFAULT_MIN_COV, help=f"HPV call: minimum genome breadth in percent, required IN ADDITION to --min-reads. 0 disables the breadth criterion. Default: {DEFAULT_MIN_COV}.")
    parser.add_argument("--borderline-min-reads", type=int, default=0, help="Optional second, more permissive read threshold; a genotype clearing it but not --min-reads is reported BORDERLINE rather than NEGATIVE. 0 disables, giving v16's two-state output.")
    parser.add_argument("--borderline-min-cov", type=float, default=0.0, help="Breadth counterpart of --borderline-min-reads. 0 reuses --min-cov.")
    parser.add_argument("--genome-len", type=int, default=DEFAULT_GENOME_LEN, help=f"Reference length for the obs/exp breadth ratio only; not enforced. Default: {DEFAULT_GENOME_LEN}.")
    parser.add_argument("--template-len", type=int, default=DEFAULT_TEMPLATE_LEN, help=f"Approximate template length, same purpose. Default: {DEFAULT_TEMPLATE_LEN}.")
    parser.add_argument("--hop-rate", type=float, default=DEFAULT_HOP_RATE, help=f"Assumed index-hopping rate. With --hop-margin this sets the rescue bar at what hopping could physically deliver: expected_hops = hop_rate * source_load / (n_samples - 1). Default: {DEFAULT_HOP_RATE}.")
    parser.add_argument("--hop-margin", type=float, default=DEFAULT_HOP_MARGIN, help=f"Rescue a non-source sample whose load is at least this multiple of expected_hops. Default: {DEFAULT_HOP_MARGIN}.")
    parser.add_argument(
        "--snr-threshold",
        type=float,
        default=DEFAULT_SNR_THRESHOLD,
        help=(
            "Minimum ratio of a non-source sample's genotype read total to the source sample's "
            f"total for the reads to be rescued as a plausible co-infection. Default: {DEFAULT_SNR_THRESHOLD}."
        ),
    )
    return parser.parse_args()


def load_manifest(manifest_path: Path) -> list[dict]:
    with manifest_path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _bam_column(bam_kind: str) -> str:
    return "ssc_bam" if bam_kind == "SSC" else "dsc_bam"


def _report_column(bam_kind: str) -> str:
    return "ssc_umi_whitelist_report" if bam_kind == "SSC" else "dsc_umi_whitelist_report"


def load_whitelist_rows(report_path: Path) -> list[dict]:
    """Load a per-sample whitelist TSV (written by whitelist_index_hopping_umis.py) and keep
    only rows that can join across samples: a joinable status (exact/near -- garbage and
    unclassified_length rows never enter the cross-sample audit) and a real mate position
    (a read whose mate is unmapped has no reliable pair signature)."""
    rows = []
    with report_path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["status"] not in JOINABLE_STATUSES:
                continue
            if row["mate_pos"] == "NA":
                continue
            row["pos"] = int(row["pos"])
            row["mate_pos"] = int(row["mate_pos"])
            row["template_length"] = int(row["template_length"])
            rows.append(row)
    return rows


def signature_key(row: dict) -> tuple:
    return (row["contig"], row["pos"], row["mate_pos"], row["template_length"], row["matched_umi"])


def compute_genotype_totals(rows: list[dict]) -> dict[tuple[str, str], int]:
    totals: dict[tuple[str, str], int] = defaultdict(int)
    for row in rows:
        totals[(row["sample_id"], row["contig"])] += 1
    return totals


def run_floor_snr_audit(
    rows: list[dict],
    totals: dict[tuple[str, str], int],
    min_floor: int,
    snr_threshold: float,
    hop_rate: float = DEFAULT_HOP_RATE,
    hop_margin: float = DEFAULT_HOP_MARGIN,
    n_samples: int | None = None,
) -> tuple[dict[str, set], dict[str, dict], dict[str, Counter]]:
    """Group rows by (signature, canonical UMI) across all samples; for every group spanning
    more than one sample, the sample with the highest genotype (contig) read total is treated
    as the source, and every other sample's reads in the group are audited against three
    tests in order: below --min-floor is always blacklisted; at or above the hop expectation
    is rescued; at or above the legacy SNR ratio is rescued; otherwise blacklisted. Returns
    (blacklist, stats, leaked_source_votes): blacklist maps sample_id -> set of blacklisted
    read names; stats maps sample_id -> {"phase1_suspicious", "phase2_rescued",
    "phase2_confirmed_bad"} counts; leaked_source_votes maps sample_id -> Counter of which
    other sample "won" each blacklisted signature, for the report's leaked-source column."""
    groups: dict[tuple, list[dict]] = defaultdict(list)
    for row in rows:
        groups[signature_key(row)].append(row)
    if n_samples is None:
        n_samples = len({row["sample_id"] for row in rows})

    blacklist: dict[str, set] = defaultdict(set)
    stats: dict[str, dict] = defaultdict(lambda: {"phase1_suspicious": 0, "phase2_rescued": 0, "phase2_confirmed_bad": 0})
    leaked_source_votes: dict[str, Counter] = defaultdict(Counter)

    for key, group_rows in groups.items():
        sample_ids = {row["sample_id"] for row in group_rows}
        if len(sample_ids) < 2:
            continue  # only one sample produced this signature+UMI -- nothing to audit

        contig = key[0]
        winning_sample = max(sample_ids, key=lambda s: (totals.get((s, contig), 0), s))
        winning_total = totals.get((winning_sample, contig), 0)

        for row in group_rows:
            sample_id = row["sample_id"]
            if sample_id == winning_sample:
                continue
            stats[sample_id]["phase1_suspicious"] += 1
            patient_total = totals.get((sample_id, contig), 0)
            ratio = patient_total / winning_total if winning_total > 0 else 0.0
            #     expected_hops = hop_rate * source_load / (n_samples - 1)
            #     rescue if patient_load >= hop_margin * expected_hops
            expected_hops = (hop_rate * winning_total / max(n_samples - 1, 1)) if n_samples else 0.0
            need = hop_margin * expected_hops
            if patient_total < min_floor:
                rescued, why = False, "below_min_floor"
            elif patient_total >= need:
                rescued, why = True, "above_hop_expectation"
            elif ratio >= snr_threshold:
                rescued, why = True, "legacy_snr_ratio"
            else:
                rescued, why = False, "deleted_hop_signature"
            if not rescued:
                blacklist[sample_id].add(row["read_name"])
                stats[sample_id]["phase2_confirmed_bad"] += 1
                stats[sample_id].setdefault("reasons", Counter())[why] += 1
                leaked_source_votes[sample_id][winning_sample] += 1
            else:
                stats[sample_id]["phase2_rescued"] += 1
                stats[sample_id].setdefault("reasons", Counter())[why] += 1

    return blacklist, stats, leaked_source_votes


def hpv_genotype_table(bam_path: Path) -> list[tuple[str, int, float]]:
    """[(contig, primary_reads, breadth_percent)] for every HPV contig with reads.

    Counts primary alignments only (excluding secondary/supplementary/duplicate), which is
    the basis every threshold uses. NB v16 took the read count from `samtools view -F 2304`
    but the breadth from `samtools coverage`, whose default filter ALSO drops duplicates --
    so on duplicate-marked input its two criteria were read off different sets. Both come
    from one pass here.
    """
    rows = []
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        hpv_contigs = [c for c in bam.references if re.search("hpv", c, re.IGNORECASE)]
        for contig in hpv_contigs:
            length = bam.get_reference_length(contig)
            reads = 0
            covered = bytearray(length)
            for read in bam.fetch(contig):
                if read.is_secondary or read.is_supplementary or read.is_duplicate:
                    continue
                reads += 1
                # Union of aligned blocks, i.e. reference positions an aligned read spans.
                # This reproduces `samtools depth` / `samtools coverage` exactly.
                #
                # Do NOT use pysam's pileup() for this: it applies mpileup's default
                # overlap and quality handling and undercounts breadth badly. Measured on
                # FU_135_S1 HPV16, pileup gave 11.3458% where samtools coverage, samtools
                # coverage -r, and samtools depth all agree on 1330/7906 = 16.8227% -- a
                # ~33% undercount, which near a 10% breadth threshold flips calls.
                for start, end in read.get_blocks():
                    for i in range(start, min(end, length)):
                        covered[i] = 1
            if reads:
                rows.append((contig, reads,
                             100.0 * sum(covered) / length if length else 0.0))
    return rows


def call_hpv_status(table: list[tuple[str, int, float]], min_reads: int, min_cov: float,
                    genome_len: int, template_len: int,
                    borderline_min_reads: int = 0, borderline_min_cov: float = 0.0
                    ) -> tuple[str, list[dict]]:
    """POSITIVE iff some genotype has reads >= min_reads AND breadth >= min_cov.

    Dipon Das' rule, from reads_coverage_with_ihop_v16_04_02_2029.py. Ported from awk and
    verified to give identical verdicts on 8 of 8 real samples. One deliberate divergence
    from v16: the optional BORDERLINE band below; unset, this reproduces v16 exactly.

    obs/exp is reported, not enforced. For n templates of length L on a genome of G bases,
    random placement predicts breadth = 1-(1-L/G)^n; a ratio far below 1 means the reads
    pile at one locus rather than spanning the genome, which is an amplicon/carryover
    signature that neither criterion detects.

    The optional BORDERLINE band is a deliberate divergence from v16: a genotype clearing a
    second, more permissive pair but not the strict one is threshold-dependent rather than
    negative. Unset (both 0), this returns v16's two-state result exactly.
    """
    use_border = borderline_min_reads > 0 or borderline_min_cov > 0
    if use_border:
        borderline_min_reads = borderline_min_reads or min_reads
        borderline_min_cov = borderline_min_cov or min_cov
    status = "NEGATIVE"
    audit = []
    for contig, reads, cov in table:
        n = reads // 2                                   # paired records -> templates
        exp_cov = 100.0 * (1 - (1 - template_len / genome_len) ** n) if n > 0 else 0.0
        ratio = cov / exp_cov if exp_cov > 0 else 0.0
        pass_r, pass_c = reads >= min_reads, cov >= min_cov
        if pass_r and pass_c:
            verdict, status = "PASS", "POSITIVE"
        elif use_border and reads >= borderline_min_reads and cov >= borderline_min_cov:
            verdict = "borderline"
            if status == "NEGATIVE":
                status = "BORDERLINE"
        else:
            verdict = "fail"
        audit.append({
            "contig": contig, "reads": reads, "breadth_pct": round(cov, 4),
            "obs_exp_ratio": round(ratio, 3),
            "criteria": f"{'reads_ok' if pass_r else 'reads_low'},"
                        f"{'cov_ok' if pass_c else 'cov_low'}",
            "verdict": verdict,
            "flag": "CLUSTERED_amplicon_like" if (0 < ratio < 0.25 and n >= 5) else ".",
        })
    return status, audit


def write_clean_bam(input_bam: Path, output_bam: Path, blacklisted_qnames: set) -> None:
    if not blacklisted_qnames:
        shutil.copy(input_bam, output_bam)
    else:
        with pysam.AlignmentFile(str(input_bam), "rb") as bam_in:
            with pysam.AlignmentFile(str(output_bam), "wb", header=bam_in.header) as bam_out:
                for read in bam_in:
                    if read.query_name not in blacklisted_qnames:
                        bam_out.write(read)
    pysam.index(str(output_bam))


def clean_bam_path(original_bam: Path) -> Path:
    return original_bam.with_name(original_bam.stem + "_clean.bam")


def main() -> None:
    args = parse_args()

    manifest_path = Path(args.manifest)
    if not manifest_path.is_file():
        raise FileNotFoundError(f"manifest not found: {manifest_path}")
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    manifest_rows = load_manifest(manifest_path)
    bam_column = _bam_column(args.bam_kind)
    report_column = _report_column(args.bam_kind)

    all_rows: list[dict] = []
    sample_bams: dict[str, Path] = {}
    ssc_filtered_bams: dict[str, Path] = {}
    for entry in manifest_rows:
        sample_id = entry["sample_id"]
        bam_path = Path(entry[bam_column])
        report_path = Path(entry[report_column])
        if not bam_path.is_file() or bam_path.stat().st_size == 0:
            print(f"WARNING: {sample_id}: {args.bam_kind} BAM not found or empty, skipping ({bam_path}).", file=sys.stderr)
            continue
        if not report_path.is_file():
            print(f"WARNING: {sample_id}: {args.bam_kind} whitelist report not found, skipping ({report_path}).", file=sys.stderr)
            continue
        sample_bams[sample_id] = bam_path
        if args.bam_kind == "SSC" and entry.get("ssc_filtered_bam"):
            filtered_path = Path(entry["ssc_filtered_bam"])
            if filtered_path.is_file() and filtered_path.stat().st_size > 0:
                ssc_filtered_bams[sample_id] = filtered_path
        all_rows.extend(load_whitelist_rows(report_path))

    if len(sample_bams) < 2:
        print(f"Fewer than 2 usable {args.bam_kind} samples available -- nothing to cross-sample audit.")
        return

    totals = compute_genotype_totals(all_rows)
    blacklist, stats, leaked_source_votes = run_floor_snr_audit(
        all_rows, totals, args.min_floor, args.snr_threshold,
        hop_rate=args.hop_rate, hop_margin=args.hop_margin)

    umi_matched_reads = defaultdict(int)
    for row in all_rows:
        umi_matched_reads[row["sample_id"]] += 1

    report_rows = []
    call_audit_rows = []
    for sample_id, bam_path in sample_bams.items():
        sample_blacklist = blacklist.get(sample_id, set())
        clean_path = output_dir / clean_bam_path(bam_path).name
        write_clean_bam(bam_path, clean_path, sample_blacklist)

        if sample_id in ssc_filtered_bams:
            filtered_clean_path = output_dir / clean_bam_path(ssc_filtered_bams[sample_id]).name
            write_clean_bam(ssc_filtered_bams[sample_id], filtered_clean_path, sample_blacklist)

        sample_stats = stats.get(sample_id, {"phase1_suspicious": 0, "phase2_rescued": 0, "phase2_confirmed_bad": 0})
        primary_leaked_source = leaked_source_votes[sample_id].most_common(1)
        # HPV call, before and after decontamination. Both from the same code path, so a
        # status change is attributable to the reads removed rather than to a different rule.
        call_kwargs = dict(min_reads=args.min_reads, min_cov=args.min_cov,
                           genome_len=args.genome_len, template_len=args.template_len,
                           borderline_min_reads=args.borderline_min_reads,
                           borderline_min_cov=args.borderline_min_cov)
        pre_status, pre_audit = call_hpv_status(hpv_genotype_table(bam_path), **call_kwargs)
        post_status, post_audit = call_hpv_status(hpv_genotype_table(clean_path), **call_kwargs)
        for stage, rows_ in (("pre_IHOP", pre_audit), ("post_IHOP", post_audit)):
            for entry in rows_:
                call_audit_rows.append(dict(sample_id=sample_id, bam_kind=args.bam_kind,
                                            stage=stage, min_reads=args.min_reads,
                                            min_cov=args.min_cov, **entry))

        report_rows.append(
            {
                "sample_id": sample_id,
                "bam_kind": args.bam_kind,
                "umi_matched_reads": umi_matched_reads.get(sample_id, 0),
                "phase1_suspicious": sample_stats["phase1_suspicious"],
                "phase2_rescued": sample_stats["phase2_rescued"],
                "phase2_confirmed_bad": sample_stats["phase2_confirmed_bad"],
                "final_reads": umi_matched_reads.get(sample_id, 0) - len(sample_blacklist),
                "primary_leaked_source": primary_leaked_source[0][0] if primary_leaked_source else "None",
                "hpv_status_pre_ihop": pre_status,
                "hpv_status_post_ihop": post_status,
            }
        )

    if call_audit_rows:
        call_audit_path = output_dir / f"IHOP_call_audit.{args.bam_kind}.tsv"
        cols = ("sample_id", "bam_kind", "stage", "contig", "reads", "breadth_pct",
                "min_reads", "min_cov", "obs_exp_ratio", "criteria", "verdict", "flag")
        with call_audit_path.open("w") as handle:
            handle.write("\t".join(cols) + "\n")
            for row in call_audit_rows:
                handle.write("\t".join(str(row[c]) for c in cols) + "\n")
        print(f"Per-genotype call audit written to: {call_audit_path}")

    report_path = output_dir / f"IHOP_Final_Report.{args.bam_kind}.tsv"
    with report_path.open("w") as handle:
        handle.write("\t".join(REPORT_COLUMNS) + "\n")
        for row in report_rows:
            handle.write("\t".join(str(row[column]) for column in REPORT_COLUMNS) + "\n")

    total_blacklisted = sum(len(qnames) for qnames in blacklist.values())
    print(f"Index-hopping audit ({args.bam_kind}): {len(sample_bams)} sample(s), {total_blacklisted} read(s) blacklisted across all samples.")
    print(f"Report written to: {report_path}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:  # noqa: BLE001 - top-level fail-fast boundary, see module docstring
        print(f"FATAL: {exc}", file=sys.stderr)
        sys.exit(1)
