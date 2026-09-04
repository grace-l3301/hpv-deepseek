#!/usr/bin/env python3
"""
whitelist_index_hopping_umis.py

Index-hopping detection: UMI whitelisting step.

For a final SSC or DSC consensus BAM, classify each templates's RX:Z: UMI tag against the
KAPA Universal UMI catalog using the same Levenshtein-distance matcher that
extract_umis_by_species.py uses to classify raw reads before extraction (umi_matching.py),
instead of a Hamming/positional-mismatch whitelist. Levenshtein distance is length-aware
(each observed UMI segment is checked against only the length-appropriate catalog) and
correctly scores insertion/deletion errors that a positional mismatch count cannot detect
at all.

This script performs the whitelisting/classification step only: it reports, per read, the
genomic signature (contig, position, mate position, template length) and UMI match status
(exact / near / garbage / unclassified_length), which is the foundation a later cross-sample
index-hopping audit (matching identical signature+UMI across different samples) would join
on. It does not itself perform that cross-sample join, apply any accept/reject threshold
across samples, or modify the input BAM.

RX tag format: fgbio's ExtractUmisFromBam writes a single RX:Z: tag shared by both mates of
a template, formatted as "<read1 UMI>-<read2 UMI>" when both mates carry a UMI. This has been
confirmed for the pre-consensus, pre-GroupReadsByUmi BAM that extract_umis_by_species.py
produces. Whether GroupReadsByUmi / CallMolecularConsensusReads / CallDuplexConsensusReads
alter this value before it reaches the final SSC/DSC BAM RX tag has not been independently
verified here -- this script handles a bare single-segment RX value as well as a "-"-joined
pair, and reports any segment whose length isn't 3 or 5 as "unclassified_length" rather than
guessing, so an unexpected tag shape shows up in the summary counts instead of silently
mis-classifying reads.

Author: Samuli Eldfors
"""

from __future__ import annotations

__version__ = "1.0.0"

import argparse
import re
import sys
from pathlib import Path

import pysam

from umi_matching import DEFAULT_MAX_EDIT_DISTANCE, UmiWhitelistResult, classify_umi_for_whitelist, load_umi_catalogs

DEFAULT_CONTIG_FILTER_REGEX = r"(?i)hpv"

REPORT_COLUMNS = (
    "sample_id",
    "bam_kind",
    "read_name",
    "contig",
    "pos",
    "mate_pos",
    "template_length",
    "raw_rx",
    "matched_umi",
    "edit_distance",
    "status",
)

# Severity order used to pick one overall status when an RX tag has more than one segment
# (e.g. "<read1 UMI>-<read2 UMI>"): the worst segment status wins for the read as a whole.
STATUS_SEVERITY = {"exact": 0, "near": 1, "garbage": 2, "unclassified_length": 3}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Whitelist RX:Z: UMI tags in a final SSC/DSC consensus BAM against the KAPA "
            "Universal UMI catalog via Levenshtein distance, as the whitelisting step of "
            "index-hopping detection."
        )
    )
    parser.add_argument("--input-bam", required=True, help="Final, coordinate-sorted, indexed SSC or DSC consensus BAM.")
    parser.add_argument("--sample-id", required=True, help="Sample ID, recorded in the output report.")
    parser.add_argument("--bam-kind", required=True, help="Free-text label for this BAM tier, e.g. SSC, SSC_filtered, or DSC.")
    parser.add_argument("--umi-catalog-file", required=True, help="Path to the KAPA Universal UMI catalog text file.")
    parser.add_argument("--output-report", required=True, help="Path for the per-read whitelist TSV report.")
    parser.add_argument(
        "--max-edit-distance",
        type=int,
        default=DEFAULT_MAX_EDIT_DISTANCE,
        help=f"Maximum Levenshtein distance accepted as a 'near' catalog match. Default: {DEFAULT_MAX_EDIT_DISTANCE}.",
    )
    parser.add_argument(
        "--contig-filter-regex",
        default=DEFAULT_CONTIG_FILTER_REGEX,
        help=(
            "Only reads mapped to a contig matching this regex are classified; others are "
            f"skipped entirely (avoids scanning the bulk human-genome-mapped reads that are "
            f"irrelevant to viral cross-sample contamination). Default: {DEFAULT_CONTIG_FILTER_REGEX!r}. "
            "Pass an empty string to disable filtering and classify every contig."
        ),
    )
    return parser.parse_args()


def classify_rx_tag(
    rx: str,
    three_nt_catalog: frozenset[str],
    five_nt_catalog: frozenset[str],
    max_edit_distance: int,
    cache: dict[str, UmiWhitelistResult],
) -> tuple[str, ...]:
    """Split an RX tag into its per-mate segments (on '-' if present, else treat the whole
    value as one segment) and classify each independently, memoizing by raw segment string
    so a segment observed on many reads is only run through Levenshtein once."""
    segments = rx.split("-") if "-" in rx else [rx]
    results = []
    for segment in segments:
        if segment not in cache:
            cache[segment] = classify_umi_for_whitelist(segment, three_nt_catalog, five_nt_catalog, max_edit_distance)
        results.append(cache[segment])
    return tuple(results)


def overall_status(segment_results: tuple[UmiWhitelistResult, ...]) -> str:
    return max((r.status for r in segment_results), key=lambda status: STATUS_SEVERITY[status])


def whitelist_bam(
    input_bam: Path,
    sample_id: str,
    bam_kind: str,
    three_nt_catalog: frozenset[str],
    five_nt_catalog: frozenset[str],
    max_edit_distance: int,
    contig_filter_regex: str,
) -> tuple[list[dict], dict[str, int]]:
    contig_pattern = re.compile(contig_filter_regex) if contig_filter_regex else None
    cache: dict[str, UmiWhitelistResult] = {}
    rows: list[dict] = []
    counts = {"exact": 0, "near": 0, "garbage": 0, "unclassified_length": 0}

    with pysam.AlignmentFile(str(input_bam), "rb") as bam_in:
        for read in bam_in:
            if read.is_unmapped or read.is_secondary or read.is_supplementary or not read.is_read1:
                continue
            if not read.has_tag("RX"):
                continue
            contig = read.reference_name
            if contig_pattern is not None and not contig_pattern.search(contig or ""):
                continue

            rx = read.get_tag("RX")
            segment_results = classify_rx_tag(rx, three_nt_catalog, five_nt_catalog, max_edit_distance, cache)
            status = overall_status(segment_results)
            counts[status] += 1

            rows.append(
                {
                    "sample_id": sample_id,
                    "bam_kind": bam_kind,
                    "read_name": read.query_name,
                    "contig": contig,
                    "pos": read.reference_start + 1,
                    "mate_pos": (read.next_reference_start + 1) if not read.mate_is_unmapped else "NA",
                    "template_length": read.template_length,
                    "raw_rx": rx,
                    "matched_umi": "-".join(r.matched_umi or "NA" for r in segment_results),
                    "edit_distance": ",".join("NA" if r.edit_distance is None else str(r.edit_distance) for r in segment_results),
                    "status": status,
                }
            )

    return rows, counts


def write_report(output_report: Path, rows: list[dict]) -> None:
    output_report.parent.mkdir(parents=True, exist_ok=True)
    with output_report.open("w") as handle:
        handle.write("\t".join(REPORT_COLUMNS) + "\n")
        for row in rows:
            handle.write("\t".join(str(row[column]) for column in REPORT_COLUMNS) + "\n")


def main() -> None:
    args = parse_args()

    input_bam = Path(args.input_bam)
    if not input_bam.is_file():
        raise FileNotFoundError(f"input BAM not found: {input_bam}")
    catalog_path = Path(args.umi_catalog_file)
    if not catalog_path.is_file():
        raise FileNotFoundError(f"UMI catalog file not found: {catalog_path}")

    three_nt_catalog, five_nt_catalog = load_umi_catalogs(catalog_path)

    rows, counts = whitelist_bam(
        input_bam,
        args.sample_id,
        args.bam_kind,
        three_nt_catalog,
        five_nt_catalog,
        args.max_edit_distance,
        args.contig_filter_regex,
    )

    write_report(Path(args.output_report), rows)

    total = sum(counts.values())
    print(f"UMI whitelist summary for {args.sample_id} ({args.bam_kind}):")
    for status in ("exact", "near", "garbage", "unclassified_length"):
        print(f"  {status}: {counts[status]} read(s)")
    print(f"  total: {total} read(s)")
    if counts["unclassified_length"] > 0:
        print(
            f"  WARNING: {counts['unclassified_length']} read(s) had an RX segment whose length "
            "was neither 3 nor 5 -- verify the RX tag format on this BAM tier matches the "
            "assumption documented in this script's module docstring.",
            file=sys.stderr,
        )
    print(f"Report written to: {args.output_report}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:  # noqa: BLE001 - top-level fail-fast boundary, see module docstring
        print(f"FATAL: {exc}", file=sys.stderr)
        sys.exit(1)
