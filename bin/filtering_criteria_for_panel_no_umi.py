#!/usr/bin/env python3

import argparse
import csv
import sys
from typing import Iterable, List


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Filter MAF records for non-UMI panel data. By default, filters SNP/indel "
            "records with COMMON == 1, VAF < 0.01, VAF > 0.99, DEPTH < 500, or indel "
            "length > 20 bp."
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
        type=float,
        default=0.05,
        help="Minimum Mean_VAF for non-Mutect2 single-caller records (default: 0.05)",
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


def single_caller_name(row: dict) -> str:
    callers = split_callers(row)
    if len(callers) == 1:
        return callers[0]
    return ""


def numeric_int(raw: str, default: int = 0) -> int:
    try:
        return int(float(str(raw)))
    except (TypeError, ValueError):
        return default


def filter_values(row: dict) -> List[str]:
    return split_values(row.get("FILTER", ""))


def is_pass_only(row: dict) -> bool:
    filters = filter_values(row)
    return len(filters) == 1 and filters[0].upper() == "PASS"
def evaluate_filters(
    row: dict,
    min_vaf: float,
    max_vaf: float,
    min_depth: float,
    max_indel_len: int,
    single_caller_min_vaf: float,
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

    n_callers = numeric_int(row.get("N_Callers", "0"))
    filters = filter_values(row)
    caller = single_caller_name(row).lower()

    if n_callers == 1:
        if not is_pass_only(row):
            reasons.append("NPASS")
        if caller == "mutect2":
            pass
        else:
            if vaf <= single_caller_min_vaf:
                reasons.append("SCVAF")

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
                single_caller_min_vaf=args.single_caller_min_vaf,
            )
            if args.keep_filtered:
                row["FILTER_C"] = "PASS" if not reasons else ";".join(reasons)
                writer.writerow(row)
            elif not reasons:
                writer.writerow(row)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
