#!/usr/bin/env python3
"""
umi_matching.py

Shared KAPA Universal UMI catalog loading and Levenshtein-distance fuzzy matching.
Used by extract_umis_by_species.py to classify a read's UMI species from its own
leading bases before extraction, and by whitelist_index_hopping_umis.py to classify
an already-extracted UMI (a BAM RX tag value) against the same catalog after
extraction, as the whitelisting step of index-hopping detection.

Author: Samuli Eldfors
"""

from __future__ import annotations

__version__ = "1.0.0"

from dataclasses import dataclass
from pathlib import Path

import Levenshtein

DEFAULT_MAX_EDIT_DISTANCE = 1


@dataclass(frozen=True)
class UmiMatch:
    species: str  # "3nt", "5nt", or "unclassified"
    umi: str | None
    edit_distance: int | None


UNCLASSIFIED = UmiMatch("unclassified", None, None)


@dataclass(frozen=True)
class CatalogMatch:
    umi: str | None  # nearest catalog entry; None only if the catalog itself is empty
    edit_distance: int | None  # true nearest Levenshtein distance, even when it exceeds max_edit_distance
    accepted: bool  # True iff edit_distance is within the caller's max_edit_distance


@dataclass(frozen=True)
class UmiWhitelistResult:
    species: str | None  # "3nt" or "5nt" (from the observed UMI's own length), or None
    matched_umi: str | None
    edit_distance: int | None
    status: str  # "exact", "near", "garbage", or "unclassified_length"


def load_umi_catalogs(catalog_path: Path) -> tuple[frozenset[str], frozenset[str]]:
    """Parse the KAPA UMI catalog (tab-separated: index, bare UMI, full sequence with
    T-overhang; '#'-prefixed lines are comments). Returns (three_nt_catalog, five_nt_catalog)
    of bare UMI strings, partitioned by length."""
    three_nt = set()
    five_nt = set()
    with open(catalog_path) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 2:
                continue
            umi = fields[1].strip()
            if len(umi) == 3:
                three_nt.add(umi)
            elif len(umi) == 5:
                five_nt.add(umi)
            else:
                raise ValueError(f"{catalog_path}: UMI {umi!r} has length {len(umi)}, expected 3 or 5")
    if not three_nt or not five_nt:
        raise ValueError(f"{catalog_path}: expected both 3-nt and 5-nt UMI catalog entries")
    return frozenset(three_nt), frozenset(five_nt)


OVERHANG_BASE = "T"


def match_umi_with_overhang(leading: str, catalog: frozenset[str],
                            max_edit_distance: int) -> CatalogMatch:
    """Match `leading` against `catalog` scoring the UMI and the T-overhang together."""
    best_umi = None
    best_distance = None
    for umi in sorted(catalog):
        full = umi + OVERHANG_BASE
        distance = Levenshtein.distance(leading[:len(full)], full)
        if best_distance is None or distance < best_distance:
            best_distance = distance
            best_umi = umi
    accepted = best_distance is not None and best_distance <= max_edit_distance
    return CatalogMatch(best_umi, best_distance, accepted)


def match_umi(query: str, catalog: frozenset[str], max_edit_distance: int) -> CatalogMatch:
    """Return the catalog entry closest to query by Levenshtein distance, always including
    the true nearest distance (even when it exceeds max_edit_distance) so callers that need
    the real magnitude of a rejected match -- not just accept/reject -- don't lose it.
    Iterates sorted(catalog) so a genuine tie between two equidistant entries resolves
    deterministically to the lexicographically-first one, rather than depending on set
    iteration order."""
    best_umi = None
    best_distance = None
    for umi in sorted(catalog):
        distance = Levenshtein.distance(query, umi)
        if best_distance is None or distance < best_distance:
            best_distance = distance
            best_umi = umi
    accepted = best_distance is not None and best_distance <= max_edit_distance
    return CatalogMatch(best_umi, best_distance, accepted)


def classify_umi_for_whitelist(
    observed_umi: str,
    three_nt_catalog: frozenset[str],
    five_nt_catalog: frozenset[str],
    max_edit_distance: int = DEFAULT_MAX_EDIT_DISTANCE,
) -> UmiWhitelistResult:
    """Whitelist an already-extracted UMI (e.g. one segment of a BAM RX tag, already
    isolated from any genomic sequence or T-overhang) against the length-appropriate
    catalog. Unlike classify_leading_sequence in extract_umis_by_species.py, there is no
    genomic tail to fall back on and no length ambiguity to resolve: the observed UMI's
    own length unambiguously selects the one catalog to check, which avoids the
    mixed-length whitelist ambiguity of comparing against a flat combined catalog."""
    if len(observed_umi) == 3:
        species, catalog = "3nt", three_nt_catalog
    elif len(observed_umi) == 5:
        species, catalog = "5nt", five_nt_catalog
    else:
        return UmiWhitelistResult(None, None, None, "unclassified_length")

    if observed_umi in catalog:
        return UmiWhitelistResult(species, observed_umi, 0, "exact")

    match = match_umi(observed_umi, catalog, max_edit_distance)
    status = "near" if match.accepted else "garbage"
    return UmiWhitelistResult(species, match.umi, match.edit_distance, status)
