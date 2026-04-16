#!/usr/bin/env python3

import argparse
import csv
import sys
from typing import Iterable, List


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Filter merged MAF records for non-UMI panel data. By default, only SNP and "
            "indel records are evaluated. A record is filtered if COMMON == 1, if Mean_VAF "
            "is < 0.01 or > 0.99, if mean DEPTH is < 500, or if an indel length is > 20 bp. "
            "Single-caller records are removed. Records supported only by LoFreq, Pisces, "
            "and/or VarDict are removed. If Mutect2 supports a record, the corresponding "
            "Mutect2 FILTER entry must be PASS. Indels overlapping RepeatMasker "
            "(RMSK != NA) are removed."
        )
    )
    parser.add_argument("maf", help="Input MAF/TSV file")
    parser.add_argument(
        "--min-vaf",
        type=float,
        default=0.01,
        help="Minimum VAF required for SNP/indel records (default: 0.01)",
    )
    parser.add_argument(
        "--max-vaf",
        type=float,
        default=0.99,
        help="Maximum VAF allowed for SNP/indel records (default: 0.99)",
    )
    parser.add_argument(
        "--min-depth",
        type=float,
        default=500,
        help="Minimum depth required for SNP records (default: 500)",
    )
    parser.add_argument(
        "--keep-filtered",
        action="store_true",
        help="Keep filtered rows and append FILTER_C instead of dropping them",
    )
    parser.add_argument(
        "--max-indel-len",
        type=int,
        default=20,
        help="Maximum allowed indel length in bp (default: 20)",
    )
    parser.add_argument(
        "--single-caller-min-vaf",
        "--lpv-only-min-vaf",
        dest="lpv_only_min_vaf",
        type=float,
        default=0.05,
        help=(
            "Minimum Mean_VAF for variants supported only by LoFreq/Pisces/VarDict "
            "(default: 0.05)"
        ),
    )
    parser.add_argument(
        "--keep-lpv-only",
        dest="remove_lpv_only",
        action="store_false",
        default=True,
        help="Do not remove records supported only by LoFreq/Pisces/VarDict",
    )
    parser.add_argument(
        "--keep-single-caller",
        dest="remove_single_caller",
        action="store_false",
        default=True,
        help="Do not remove records supported by only one caller",
    )
    parser.add_argument(
        "--no-require-mutect2-pass",
        dest="require_mutect2_pass",
        action="store_false",
        default=True,
        help="Do not require the Mutect2-specific FILTER value to be PASS",
    )
    parser.add_argument(
        "--keep-rmsk-indel",
        dest="remove_rmsk_indel",
        action="store_false",
        default=True,
        help="Do not remove indels with RMSK value not equal to NA",
    )
    return parser.parse_args()


def split_values(raw: str) -> List[str]:
    return [item.strip() for item in str(raw).split("|") if item.strip() not in {"", "."}]


def numeric_mean(raw: str, default: float = 0.0) -> float:
    values = []
    for item in split_values(raw):
        try:
            values.append(float(item))
        except ValueError:
            continue
    if not values:
        return default
    return sum(values) / len(values)


def any_common(raw: str) -> bool:
    for item in split_values(raw):
        if item == "1":
            return True
        try:
            if float(item) >= 1:
                return True
        except ValueError:
            continue
    return False


def is_snp(row: dict) -> bool:
    variant_type = row.get("Variant_Type", "").upper()
    if variant_type == "SNP":
        return True

    ref = row.get("Reference_Allele", "")
    alt = row.get("Tumor_Seq_Allele2", "")
    return len(ref) == 1 and len(alt) == 1


def is_indel(row: dict) -> bool:
    variant_type = row.get("Variant_Type", "").upper()
    if variant_type in {"DEL", "INS", "INDEL"}:
        return True

    ref = row.get("Reference_Allele", "")
    alt = row.get("Tumor_Seq_Allele2", "")
    return len(ref) != len(alt)


def indel_length(row: dict) -> int:
    ref = row.get("Reference_Allele", "")
    alt = row.get("Tumor_Seq_Allele2", "")
    return abs(len(ref) - len(alt))


def split_callers(row: dict) -> List[str]:
    return split_values(row.get("Mutation_Status", ""))


def caller_names(row: dict) -> List[str]:
    return [caller.lower() for caller in split_callers(row)]


def numeric_int(raw: str, default: int = 0) -> int:
    try:
        return int(float(str(raw)))
    except (TypeError, ValueError):
        return default


def filter_values(row: dict) -> List[str]:
    return split_values(row.get("FILTER", ""))


def all_filters_pass(row: dict) -> bool:
    filters = filter_values(row)
    return bool(filters) and all(filter_value.upper() == "PASS" for filter_value in filters)


def is_mutect2_only(callers: List[str]) -> bool:
    return len(callers) == 1 and callers[0] == "mutect2"


def is_lpv_only(callers: List[str]) -> bool:
    allowed = {"lofreq", "pisces", "vardict"}
    return bool(callers) and set(callers).issubset(allowed)


def get_mutect2_filter_values(row: dict) -> List[str]:
    callers = split_callers(row)
    filters = filter_values(row)
    return [
        filters[idx]
        for idx, caller in enumerate(callers)
        if idx < len(filters) and caller.lower() == "mutect2"
    ]


def evaluate_filters(
    row: dict,
    min_vaf: float,
    max_vaf: float,
    min_depth: float,
    max_indel_len: int,
    lpv_only_min_vaf: float,
    remove_lpv_only: bool,
    remove_single_caller: bool,
    require_mutect2_pass: bool,
    remove_rmsk_indel: bool,
) -> List[str]:
    reasons: List[str] = []

    if not (is_snp(row) or is_indel(row)):
        return reasons

    if any_common(row.get("COMMON", "0")):
        reasons.append("COMMON")

    vaf = numeric_mean(row.get("Mean_VAF", "0"))
    if vaf < min_vaf:
        reasons.append("LVAF")
    if vaf > max_vaf:
        reasons.append("HVAF")

    depth = numeric_mean(row.get("DEPTH", "0"))
    if depth < min_depth:
        reasons.append("LDEP")

    if is_indel(row) and indel_length(row) > max_indel_len:
        reasons.append("LINDEL")
    if remove_rmsk_indel and is_indel(row) and row.get("RMSK", "NA") != "NA":
        reasons.append("RMSK_INDEL")

    callers = caller_names(row)
    if remove_single_caller and len(callers) == 1:
        reasons.append("SINGLE")
    if remove_lpv_only and is_lpv_only(callers):
        reasons.append("LPV_ONLY")
    if require_mutect2_pass and "mutect2" in callers:
        mutect2_filters = get_mutect2_filter_values(row)
        if not mutect2_filters or any(filter_value.upper() != "PASS" for filter_value in mutect2_filters):
            reasons.append("M2NPASS")
    elif is_mutect2_only(callers):
        if not all_filters_pass(row):
            reasons.append("NPASS")
    elif is_lpv_only(callers):
        if not all_filters_pass(row):
            reasons.append("NPASS")
        if vaf <= lpv_only_min_vaf:
            reasons.append("LPVVAF")

    return reasons


def main() -> int:
    args = parse_args()

    with open(args.maf, "r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise SystemExit("Input file is missing a header")

        fieldnames: Iterable[str]
        if args.keep_filtered and "FILTER_C" not in reader.fieldnames:
            fieldnames = [*reader.fieldnames, "FILTER_C"]
        else:
            fieldnames = reader.fieldnames

        writer = csv.DictWriter(
            sys.stdout,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()

        for row in reader:
            reasons = evaluate_filters(
                row,
                min_vaf=args.min_vaf,
                max_vaf=args.max_vaf,
                min_depth=args.min_depth,
                max_indel_len=args.max_indel_len,
                lpv_only_min_vaf=args.lpv_only_min_vaf,
                remove_lpv_only=args.remove_lpv_only,
                remove_single_caller=args.remove_single_caller,
                require_mutect2_pass=args.require_mutect2_pass,
                remove_rmsk_indel=args.remove_rmsk_indel,
            )
            if args.keep_filtered:
                row["FILTER_C"] = "PASS" if not reasons else ";".join(reasons)
                writer.writerow(row)
            elif not reasons:
                writer.writerow(row)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
