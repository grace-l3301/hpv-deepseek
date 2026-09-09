#!/usr/bin/env python3
"""
extract_umis_by_species.py

Classify each read pair's KAPA Universal UMI Adapter species (3-nt vs 5-nt) directly from
read sequence content, then run fgbio ExtractUmisFromBam once per species-pair combination
with the read structure that is exactly correct for that combination, and merge the results
back into a single unmapped, UMI-tagged BAM at the path the caller expects.

Invoked by ssc_and_dsc.py in place of a single blind, fixed-length ExtractUmisFromBam call.

Catalog loading and Levenshtein matching live in umi_matching.py, shared with
whitelist_index_hopping_umis.py, which whitelists already-extracted RX tags against the
same catalog as the index-hopping detection whitelisting step.

Author: Samuli Eldfors
"""

from __future__ import annotations

__version__ = "1.1.0"

import argparse
import subprocess
import sys
from pathlib import Path

import pysam

from umi_matching import (DEFAULT_MAX_EDIT_DISTANCE, OVERHANG_BASE, UNCLASSIFIED,
                          UmiMatch, load_umi_catalogs, match_umi_with_overhang)

DEFAULT_THREE_NT_TRIM_LENGTH = 4  # 3-nt UMI + 1-base T-overhang; no padding assumed
DEFAULT_FIVE_NT_TRIM_LENGTH = 6  # 5-nt UMI + 1-base T-overhang

SPECIES_LENGTH = {"3nt": 3, "5nt": 5}
CLASSIFIED_BUCKETS = ("3nt_3nt", "3nt_5nt", "5nt_3nt", "5nt_5nt")
REJECT_BUCKET = "unclassified"
ALL_BUCKETS = CLASSIFIED_BUCKETS + (REJECT_BUCKET,)


def classify_leading_sequence(
    seq: str,
    three_nt_catalog: frozenset[str],
    five_nt_catalog: frozenset[str],
    max_edit_distance: int = DEFAULT_MAX_EDIT_DISTANCE,
) -> UmiMatch:
    """Classify a read's UMI species from its own leading bases: exact match first, then
    Levenshtein fallback (<=max_edit_distance) against each catalog. 5-nt matching is
    skipped entirely for reads shorter than 5 bases, to avoid a spurious distance-1 match
    via pure deletion on a too-short string."""
    if len(seq) < 3:
        return UNCLASSIFIED

    prefix3 = seq[:3]
    if prefix3 in three_nt_catalog:
        return UmiMatch("3nt", prefix3, 0)

    if len(seq) >= 5:
        prefix5 = seq[:5]
        if prefix5 in five_nt_catalog:
            return UmiMatch("5nt", prefix5, 0)

    # Needs one base beyond the UMI, so a read too short to carry it cannot be fuzzy-matched at that length.
    if len(seq) >= 4:
        match3 = match_umi_with_overhang(seq, three_nt_catalog, max_edit_distance)
        if match3.accepted:
            return UmiMatch("3nt", match3.umi, match3.edit_distance)

    if len(seq) >= 6:
        match5 = match_umi_with_overhang(seq, five_nt_catalog, max_edit_distance)
        if match5.accepted:
            return UmiMatch("5nt", match5.umi, match5.edit_distance)

    return UNCLASSIFIED


def apply_umi_correction(read, match: UmiMatch) -> bool:
    """Rewrite leading UMI+overhang in place to the matched catalog entry."""
    if match.umi is None or match.species == UNCLASSIFIED.species:
        return False
    entry = match.umi + OVERHANG_BASE
    seq = read.query_sequence
    if seq is None or len(seq) < len(entry) or seq[:len(entry)] == entry:
        return False
    qual = read.query_qualities            # save: the next line discards it
    read.query_sequence = entry + seq[len(entry):]
    read.query_qualities = qual            # restore
    return True


def classify_pair(r1_match: UmiMatch, r2_match: UmiMatch) -> str:
    """Returns a CLASSIFIED_BUCKETS entry, or REJECT_BUCKET if either mate is unclassified.
    Never guesses a bucket on partial evidence."""
    if r1_match.species == "unclassified" or r2_match.species == "unclassified":
        return REJECT_BUCKET
    return f"{r1_match.species}_{r2_match.species}"


def build_read_structure(species: str, trim_length: int) -> str:
    """Build an fgbio read-structure string that captures the full UMI into M and skips
    exactly the remaining non-genomic bases (the T-overhang, and any padding beyond it)
    before the template segment."""
    umi_length = SPECIES_LENGTH[species]
    skip_length = trim_length - umi_length
    if skip_length < 1:
        raise ValueError(
            f"trim_length {trim_length} leaves no room for the T-overhang skip for "
            f"species {species!r} (UMI length {umi_length})"
        )
    return f"{umi_length}M{skip_length}S+T"


def read_structures_for_bucket(bucket: str, three_nt_trim: int, five_nt_trim: int) -> tuple[str, str]:
    r1_species, r2_species = bucket.split("_")
    trim_by_species = {"3nt": three_nt_trim, "5nt": five_nt_trim}
    r1_structure = build_read_structure(r1_species, trim_by_species[r1_species])
    r2_structure = build_read_structure(r2_species, trim_by_species[r2_species])
    return r1_structure, r2_structure


def read_bam_header(bam_path: Path) -> pysam.AlignmentHeader:
    with pysam.AlignmentFile(str(bam_path), "rb", check_sq=False) as bam_in:
        return bam_in.header


def iter_read_pairs(bam_path: Path):
    """Yield (read1, read2) for each interleaved pair in an unmapped, paired BAM. Raises
    RuntimeError on any pairing irregularity (odd record count, mismatched query names, or
    a pair that isn't exactly one read1 + one read2) rather than guessing."""
    with pysam.AlignmentFile(str(bam_path), "rb", check_sq=False) as bam_in:
        iterator = iter(bam_in)
        for first in iterator:
            try:
                second = next(iterator)
            except StopIteration:
                raise RuntimeError(f"{bam_path}: trailing unpaired record {first.query_name!r}")
            if first.query_name != second.query_name:
                raise RuntimeError(
                    f"{bam_path}: adjacent records have mismatched query names "
                    f"({first.query_name!r} != {second.query_name!r}); expected interleaved "
                    "read1/read2 pairs"
                )
            if first.is_read1 and second.is_read2:
                yield first, second
            elif first.is_read2 and second.is_read1:
                yield second, first
            else:
                raise RuntimeError(f"{bam_path}: pair {first.query_name!r} is not one read1 + one read2")


def split_bam_into_buckets(
    input_bam: Path,
    work_dir: Path,
    sample_id: str,
    three_nt_catalog: frozenset[str],
    five_nt_catalog: frozenset[str],
    max_edit_distance: int,
) -> tuple[dict, dict, int]:
    """Classify every pair in input_bam and write it to its bucket's BAM under work_dir.
    Returns (counts, bucket_paths, corrected): counts maps every bucket in ALL_BUCKETS to
    its pair count (0 if empty) and NOTHING else, because the caller sums its values to get
    the pair total; bucket_paths maps only non-empty buckets to their Path; corrected is the
    number of reads whose UMI bases were rewritten."""
    header = read_bam_header(input_bam)
    counts = {bucket: 0 for bucket in ALL_BUCKETS}
    corrected = 0
    writers: dict[str, pysam.AlignmentFile] = {}
    bucket_paths: dict[str, Path] = {}

    try:
        for r1, r2 in iter_read_pairs(input_bam):
            r1_match = classify_leading_sequence(r1.query_sequence, three_nt_catalog, five_nt_catalog, max_edit_distance)
            r2_match = classify_leading_sequence(r2.query_sequence, three_nt_catalog, five_nt_catalog, max_edit_distance)
            bucket = classify_pair(r1_match, r2_match)
            counts[bucket] += 1

            # Correct the UMI bases before extraction, so a 1-error UMI joins its true family instead of forming a singleton.
            if bucket != REJECT_BUCKET:
                # Tallied separately, NOT into `counts`: the caller reports
                # sum(counts.values()) as the pair total, so an extra key there would
                # silently inflate it by the correction count (~2% of reads).
                corrected += int(apply_umi_correction(r1, r1_match))
                corrected += int(apply_umi_correction(r2, r2_match))

            if bucket not in writers:
                bucket_path = work_dir / f"{sample_id}_unmapped_bucket_{bucket}.bam"
                writers[bucket] = pysam.AlignmentFile(str(bucket_path), "wb", header=header)
                bucket_paths[bucket] = bucket_path

            writers[bucket].write(r1)
            writers[bucket].write(r2)
    finally:
        for writer in writers.values():
            writer.close()

    return counts, bucket_paths, corrected


def run_extract_umis(
    fgbio_jar: Path,
    jvm_options: str,
    input_bam: Path,
    output_bam: Path,
    r1_structure: str,
    r2_structure: str,
) -> None:
    cmd = [
        "java",
        *jvm_options.split(),
        "-jar",
        str(fgbio_jar),
        "ExtractUmisFromBam",
        "-i",
        str(input_bam),
        "-o",
        str(output_bam),
        "-r",
        r1_structure,
        r2_structure,
        "-t",
        "RX",
        "-a",
        "true",
    ]
    subprocess.run(cmd, check=True)


def merge_bam_files(extracted_paths: list[Path], output_bam: Path, header: pysam.AlignmentHeader) -> None:
    if not extracted_paths:
        with pysam.AlignmentFile(str(output_bam), "wb", header=header):
            pass
        return
    cmd = ["samtools", "cat", "-o", str(output_bam), *[str(p) for p in extracted_paths]]
    subprocess.run(cmd, check=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Classify KAPA Universal UMI Adapter species from read content and run fgbio "
            "ExtractUmisFromBam once per species-pair combination with the correct read "
            "structure, then merge the results into a single UMI-tagged unmapped BAM."
        )
    )
    parser.add_argument("--input-bam", required=True, help="Unmapped, paired BAM to classify and extract.")
    parser.add_argument("--output-bam", required=True, help="Path for the merged, UMI-extracted unmapped BAM.")
    parser.add_argument("--work-dir", required=True, help="Directory for per-bucket intermediate BAMs.")
    parser.add_argument("--sample-id", required=True, help="Sample ID, used as the intermediate-file prefix.")
    parser.add_argument("--fgbio-jar", required=True, help="Full path to the fgbio jar.")
    parser.add_argument("--jvm-options", default="-Xmx4g", help="Java JVM memory option(s) for fgbio.")
    parser.add_argument("--umi-catalog-file", required=True, help="Path to the KAPA UMI catalog text file.")
    parser.add_argument(
        "--three-nt-trim-length",
        type=int,
        default=DEFAULT_THREE_NT_TRIM_LENGTH,
        help=f"Total non-genomic prefix length to trim for the 3-nt UMI species. Default: {DEFAULT_THREE_NT_TRIM_LENGTH}.",
    )
    parser.add_argument(
        "--five-nt-trim-length",
        type=int,
        default=DEFAULT_FIVE_NT_TRIM_LENGTH,
        help=f"Total non-genomic prefix length to trim for the 5-nt UMI species. Default: {DEFAULT_FIVE_NT_TRIM_LENGTH}.",
    )
    parser.add_argument(
        "--max-edit-distance",
        type=int,
        default=DEFAULT_MAX_EDIT_DISTANCE,
        help=f"Maximum Levenshtein distance allowed for a fuzzy catalog match. Default: {DEFAULT_MAX_EDIT_DISTANCE}.",
    )
    parser.add_argument(
        "--keep-intermediate-bams",
        action="store_true",
        help="Keep per-bucket intermediate BAMs in --work-dir instead of deleting them after the merge.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    fgbio_jar = Path(args.fgbio_jar)
    if not fgbio_jar.is_file():
        raise FileNotFoundError(f"fgbio jar not found: {fgbio_jar}")
    catalog_path = Path(args.umi_catalog_file)
    if not catalog_path.is_file():
        raise FileNotFoundError(f"UMI catalog file not found: {catalog_path}")
    input_bam = Path(args.input_bam)
    if not input_bam.is_file():
        raise FileNotFoundError(f"input BAM not found: {input_bam}")

    work_dir = Path(args.work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    output_bam = Path(args.output_bam)

    three_nt_catalog, five_nt_catalog = load_umi_catalogs(catalog_path)

    counts, bucket_paths, corrected = split_bam_into_buckets(
        input_bam, work_dir, args.sample_id, three_nt_catalog, five_nt_catalog, args.max_edit_distance
    )

    extracted_paths = []
    for bucket in CLASSIFIED_BUCKETS:
        if counts[bucket] == 0:
            continue
        r1_structure, r2_structure = read_structures_for_bucket(bucket, args.three_nt_trim_length, args.five_nt_trim_length)
        extracted_path = work_dir / f"{args.sample_id}_unmapped_umi_extracted_bucket_{bucket}.bam"
        run_extract_umis(fgbio_jar, args.jvm_options, bucket_paths[bucket], extracted_path, r1_structure, r2_structure)
        extracted_paths.append(extracted_path)

    header = read_bam_header(input_bam)
    merge_bam_files(extracted_paths, output_bam, header)

    if not args.keep_intermediate_bams:
        for path in bucket_paths.values():
            path.unlink(missing_ok=True)
        for path in extracted_paths:
            path.unlink(missing_ok=True)

    total = sum(counts.values())
    print("UMI classification summary:")
    for bucket in ALL_BUCKETS:
        print(f"  {bucket}: {counts[bucket]} read pair(s)")
    print(f"  total: {total} read pair(s)")
    print(f"  UMI bases corrected: {corrected} read(s)"
          f"{f' ({100 * corrected / (2 * total):.2f}% of reads)' if total else ''}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:  # noqa: BLE001 - top-level fail-fast boundary, see module docstring
        print(f"FATAL: {exc}", file=sys.stderr)
        sys.exit(1)
