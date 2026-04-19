#!/usr/bin/env python3

import argparse
import csv
import gzip
import math
import re
import sys
from collections import defaultdict
from statistics import NormalDist
from typing import DefaultDict, Dict, Iterable, List, Optional, TextIO, Tuple, TypeVar


VariantKey = Tuple[str, int, str, str]
PositionKey = Tuple[str, int]
PositionRecord = Tuple[str, str, str]
MafVariantKey = Tuple[str, int, int, str, str]
T = TypeVar("T")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description=(
            "Annotate a MAF/TSV file with INFO values from a VCF, and optionally derive\n"
            "case-vs-ChinaMAP association statistics from AC/AN and Variant_Count.\n\n"
            "Matching strategy:\n"
            "1. Raw VCF-style exact match when the MAF contains usable REF/ALT columns.\n"
            "2. MAF/ANNOVAR-style exact match using Chromosome, Start_Position,\n"
            "   End_Position, Reference_Allele, and Tumor_Seq_Allele2.\n"
            "   This supports indels represented with '-' and shifted coordinates.\n"
            "3. Position-only fallback is used only when the MAF row does not carry usable\n"
            "   alleles and the VCF has exactly one record at that position.\n\n"
            "Indel normalization examples:\n"
            "VCF  chr1 200 A   AT   -> MAF chr1 201 201 -  T\n"
            "VCF  chr1 300 ATC A    -> MAF chr1 301 302 TC -\n\n"
            "Statistics column:\n"
            "If the variant is matched and the INFO string contains AC and AN, the script\n"
            "builds a 2x2 table using:\n"
            "  Sample cohort   ALT = Variant_Count, REF = total_allele_count - Variant_Count\n"
            "  ChinaMAP        ALT = AC,            REF = AN - AC\n"
            "and reports Odds Ratio, 95% CI, and two-sided Fisher exact P value."
        ),
        epilog=(
            "Examples:\n"
            "  python3 annotate_maf_with_vcf_info.py input.maf ChinaMAP.vcf.gz -o output.maf\n"
            "  python3 annotate_maf_with_vcf_info.py input.maf ChinaMAP.vcf.gz \\\n"
            "      --info-column ChinaMAP --stats-column ChinaMAP_Assoc \\\n"
            "      --variant-count-column Variant_Count --total-allele-count 1100"
        ),
    )
    parser.add_argument("maf", help="Input MAF/TSV file")
    parser.add_argument("vcf", help="Input VCF/VCF.GZ file used as annotation source")
    parser.add_argument(
        "-o",
        "--output",
        default="-",
        help="Output MAF path (default: stdout)",
    )
    parser.add_argument(
        "--info-column",
        "--column-name",
        dest="column_name",
        default="ChinaMAP",
        help="Output column name for the matched VCF INFO string (default: ChinaMAP)",
    )
    parser.add_argument(
        "--missing-value",
        default="NA",
        help="Value written when no INFO annotation is found (default: NA)",
    )
    parser.add_argument(
        "--variant-count-column",
        default="Variant_Count",
        help=(
            "MAF column name containing the sample/cohort ALT allele count used in the\n"
            "association test (default: Variant_Count)"
        ),
    )
    parser.add_argument(
        "--total-allele-count",
        type=int,
        default=1100,
        help=(
            "Total allele count in the sample/cohort used together with Variant_Count.\n"
            "For diploid 550 samples this would be 1100. (default: 1100)"
        ),
    )
    parser.add_argument(
        "--stats-column",
        default=None,
        help=(
            "Output column name for the association statistics string formatted as\n"
            "OR=...;CI95=low-high;P=... (default: <info-column>_Stats)"
        ),
    )
    return parser.parse_args()


def open_maybe_gzip(path: str, mode: str = "rt") -> TextIO:
    if path.endswith(".gz"):
        return gzip.open(path, mode, newline="")
    return open(path, mode, newline="")


def detect_delimiter(path: str) -> str:
    with open_maybe_gzip(path, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith("##"):
                continue
            return "\t" if "\t" in line else ","
    return "\t"


def normalize_chr(value: str) -> str:
    chrom = str(value).strip()
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]
    if chrom.upper() in {"M", "MT"}:
        return "MT"
    return chrom


def parse_int(value: str) -> int:
    return int(float(str(value).strip()))


def clean_value(value: Optional[str]) -> Optional[str]:
    if value is None:
        return None
    cleaned = str(value).strip()
    if cleaned in {"", ".", "NA", "N/A", "NULL", "None"}:
        return None
    return cleaned


def deduplicate(values: Iterable[T]) -> List[T]:
    seen = set()
    ordered: List[T] = []
    for value in values:
        if value not in seen:
            seen.add(value)
            ordered.append(value)
    return ordered


def parse_info_header_number(header_line: str) -> Optional[Tuple[str, str]]:
    match = re.match(r"##INFO=<ID=([^,>]+),Number=([^,>]+)", header_line.strip())
    if not match:
        return None
    return match.group(1), match.group(2)


def per_alt_info_string(
    raw_info: str,
    alt_index: int,
    alt_count: int,
    info_number_by_id: Dict[str, str],
) -> str:
    if raw_info in {"", "."}:
        return raw_info or "."

    tokens: List[str] = []
    for entry in raw_info.split(";"):
        if not entry:
            continue
        if "=" not in entry:
            tokens.append(entry)
            continue

        key, value = entry.split("=", 1)
        number = info_number_by_id.get(key)
        values = value.split(",")

        if number == "A" and len(values) == alt_count:
            tokens.append(f"{key}={values[alt_index]}")
        elif number == "R" and len(values) == alt_count + 1:
            tokens.append(f"{key}={values[0]},{values[alt_index + 1]}")
        else:
            tokens.append(entry)

    return ";".join(tokens) if tokens else "."


def trim_common_suffix(ref: str, alt: str) -> Tuple[str, str]:
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref = ref[:-1]
        alt = alt[:-1]
    return ref, alt


def trim_common_prefix(pos: int, ref: str, alt: str) -> Tuple[int, str, str]:
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        pos += 1
        ref = ref[1:]
        alt = alt[1:]
    return pos, ref, alt


def vcf_alt_to_maf_key(chrom: str, pos: int, ref: str, alt: str) -> MafVariantKey:
    ref = ref.strip()
    alt = alt.strip()
    pos, ref, alt = trim_common_prefix(pos, *trim_common_suffix(ref, alt))

    if len(ref) == 1 and len(alt) == 1:
        return chrom, pos, pos, ref, alt

    if len(ref) == 1 and len(alt) > 1 and alt.startswith(ref):
        return chrom, pos + 1, pos + 1, "-", alt[1:]

    if len(ref) > 1 and len(alt) == 1 and ref.startswith(alt):
        return chrom, pos + 1, pos + len(ref) - 1, ref[1:], "-"

    return chrom, pos, pos + len(ref) - 1, ref, alt


def load_vcf_annotations(
    vcf_path: str,
) -> Tuple[
    Dict[VariantKey, List[str]],
    Dict[PositionKey, List[PositionRecord]],
    Dict[MafVariantKey, List[str]],
    Dict[Tuple[str, int, int], List[PositionRecord]],
]:
    exact_matches: DefaultDict[VariantKey, List[str]] = defaultdict(list)
    position_matches: DefaultDict[PositionKey, List[PositionRecord]] = defaultdict(list)
    maf_exact_matches: DefaultDict[MafVariantKey, List[str]] = defaultdict(list)
    maf_position_matches: DefaultDict[Tuple[str, int, int], List[PositionRecord]] = defaultdict(list)
    info_number_by_id: Dict[str, str] = {}

    with open_maybe_gzip(vcf_path, "rt") as handle:
        for raw_line in handle:
            if not raw_line.strip():
                continue
            if raw_line.startswith("##INFO="):
                parsed = parse_info_header_number(raw_line)
                if parsed is not None:
                    info_id, number = parsed
                    info_number_by_id[info_id] = number
                continue
            if raw_line.startswith("##"):
                continue
            if raw_line.startswith("#CHROM"):
                continue

            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise SystemExit(f"Invalid VCF record with fewer than 8 columns: {raw_line.strip()}")

            chrom = normalize_chr(fields[0])
            pos = parse_int(fields[1])
            ref = fields[3].strip()
            alt_values = [alt.strip() for alt in fields[4].split(",") if alt.strip()]
            raw_info = fields[7].strip() or "."

            for alt_index, alt in enumerate(alt_values):
                info = per_alt_info_string(
                    raw_info=raw_info,
                    alt_index=alt_index,
                    alt_count=len(alt_values),
                    info_number_by_id=info_number_by_id,
                )
                exact_matches[(chrom, pos, ref, alt)].append(info)
                position_matches[(chrom, pos)].append((ref, alt, info))
                maf_key = vcf_alt_to_maf_key(chrom, pos, ref, alt)
                maf_exact_matches[maf_key].append(info)
                maf_position_matches[(maf_key[0], maf_key[1], maf_key[2])].append((maf_key[3], maf_key[4], info))

    return (
        dict(exact_matches),
        dict(position_matches),
        dict(maf_exact_matches),
        dict(maf_position_matches),
    )


def maf_exact_keys(row: dict) -> List[VariantKey]:
    chrom = normalize_chr(row["Chromosome"])
    pos = parse_int(row["Start_Position"])
    candidate_keys: List[VariantKey] = []

    ref = clean_value(row.get("REF"))
    alt = clean_value(row.get("ALT"))
    if ref and alt:
        for alt_allele in alt.split(","):
            alt_allele = alt_allele.strip()
            if alt_allele:
                candidate_keys.append((chrom, pos, ref, alt_allele))

    maf_ref = clean_value(row.get("Reference_Allele"))
    maf_alt = clean_value(row.get("Tumor_Seq_Allele2"))
    if maf_ref and maf_alt and maf_ref != "-" and maf_alt != "-":
        candidate_keys.append((chrom, pos, maf_ref, maf_alt))

    return deduplicate(candidate_keys)


def maf_normalized_keys(row: dict) -> List[MafVariantKey]:
    chrom = normalize_chr(row["Chromosome"])
    start = parse_int(row["Start_Position"])
    end_value = clean_value(row.get("End_Position"))
    if end_value is None:
        return []
    end = parse_int(end_value)
    candidate_keys: List[MafVariantKey] = []

    maf_ref = clean_value(row.get("Reference_Allele"))
    maf_alt = clean_value(row.get("Tumor_Seq_Allele2"))
    if maf_ref and maf_alt:
        candidate_keys.append((chrom, start, end, maf_ref, maf_alt))

    return deduplicate(candidate_keys)


def has_usable_alleles(row: dict) -> bool:
    ref = clean_value(row.get("Reference_Allele"))
    alt = clean_value(row.get("Tumor_Seq_Allele2"))
    return ref is not None and alt is not None


def parse_info_string(raw_info: str) -> Dict[str, str]:
    info_map: Dict[str, str] = {}
    if raw_info in {"", ".", "NA", "N/A"}:
        return info_map
    for entry in raw_info.split(";"):
        if not entry:
            continue
        if "=" in entry:
            key, value = entry.split("=", 1)
            info_map[key] = value
        else:
            info_map[entry] = "1"
    return info_map


def parse_nonnegative_int(value: Optional[str]) -> Optional[int]:
    cleaned = clean_value(value)
    if cleaned is None:
        return None
    try:
        parsed = int(float(cleaned))
    except ValueError:
        return None
    if parsed < 0:
        return None
    return parsed


def format_stat_value(value: float) -> str:
    if math.isinf(value):
        return "inf"
    if value == 0:
        return "0"
    if abs(value) >= 1000 or abs(value) < 0.001:
        return f"{value:.3e}"
    return f"{value:.6f}".rstrip("0").rstrip(".")


def fisher_exact_two_sided(a: int, b: int, c: int, d: int) -> float:
    row1 = a + b
    row2 = c + d
    col1 = a + c
    total = row1 + row2

    min_x = max(0, col1 - row2)
    max_x = min(row1, col1)

    def log_choose(n: int, k: int) -> float:
        return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)

    def log_prob(x: int) -> float:
        return log_choose(col1, x) + log_choose(total - col1, row1 - x) - log_choose(total, row1)

    observed = log_prob(a)
    probabilities: List[float] = []
    for x in range(min_x, max_x + 1):
        current = log_prob(x)
        if current <= observed + 1e-12:
            probabilities.append(math.exp(current))

    p_value = sum(probabilities)
    return min(p_value, 1.0)


def association_stats_string(
    info_value: str,
    variant_count_value: Optional[str],
    total_allele_count: int,
    missing_value: str,
) -> str:
    info_map = parse_info_string(info_value)
    ac = parse_nonnegative_int(info_map.get("AC"))
    an = parse_nonnegative_int(info_map.get("AN"))
    variant_count = parse_nonnegative_int(variant_count_value)

    if ac is None or an is None or variant_count is None:
        return missing_value
    if total_allele_count <= 0 or an <= 0:
        return missing_value
    if variant_count > total_allele_count or ac > an:
        return missing_value

    a = variant_count
    b = total_allele_count - variant_count
    c = ac
    d = an - ac

    if min(a, b, c, d) < 0:
        return missing_value

    p_value = fisher_exact_two_sided(a, b, c, d)

    a_ci = a + 0.5 if 0 in {a, b, c, d} else float(a)
    b_ci = b + 0.5 if 0 in {a, b, c, d} else float(b)
    c_ci = c + 0.5 if 0 in {a, b, c, d} else float(c)
    d_ci = d + 0.5 if 0 in {a, b, c, d} else float(d)

    odds_ratio = (a_ci * d_ci) / (b_ci * c_ci)
    se = math.sqrt((1.0 / a_ci) + (1.0 / b_ci) + (1.0 / c_ci) + (1.0 / d_ci))
    z_value = NormalDist().inv_cdf(0.975)
    log_or = math.log(odds_ratio)
    ci_low = math.exp(log_or - z_value * se)
    ci_high = math.exp(log_or + z_value * se)

    return (
        f"OR={format_stat_value(odds_ratio)};"
        f"CI95={format_stat_value(ci_low)}-{format_stat_value(ci_high)};"
        f"P={format_stat_value(p_value)}"
    )


def annotate_row(
    row: dict,
    exact_matches: Dict[VariantKey, List[str]],
    position_matches: Dict[PositionKey, List[PositionRecord]],
    maf_exact_matches: Dict[MafVariantKey, List[str]],
    maf_position_matches: Dict[Tuple[str, int, int], List[PositionRecord]],
    missing_value: str,
) -> Tuple[str, str]:
    for key in maf_exact_keys(row):
        if key in exact_matches:
            return "|".join(deduplicate(exact_matches[key])), "exact_raw"

    for key in maf_normalized_keys(row):
        if key in maf_exact_matches:
            return "|".join(deduplicate(maf_exact_matches[key])), "exact_maf"

    # If MAF already carries usable REF/ALT alleles, a failed allele-aware match
    # should stay unmatched rather than falling back to position-only annotation.
    if has_usable_alleles(row):
        return missing_value, "no_match"

    pos_key = (normalize_chr(row["Chromosome"]), parse_int(row["Start_Position"]))
    records = deduplicate(position_matches.get(pos_key, []))
    if len(records) == 1:
        return records[0][2], "position_only_raw"
    if len(records) > 1:
        return missing_value, "ambiguous_position"

    end_value = clean_value(row.get("End_Position"))
    if end_value is not None:
        maf_pos_key = (
            normalize_chr(row["Chromosome"]),
            parse_int(row["Start_Position"]),
            parse_int(end_value),
        )
        maf_records = deduplicate(maf_position_matches.get(maf_pos_key, []))
        if len(maf_records) == 1:
            return maf_records[0][2], "position_only_maf"
        if len(maf_records) > 1:
            return missing_value, "ambiguous_position"
    return missing_value, "no_match"


def annotate_maf(
    maf_path: str,
    maf_delimiter: str,
    output_path: str,
    exact_matches: Dict[VariantKey, List[str]],
    position_matches: Dict[PositionKey, List[PositionRecord]],
    maf_exact_matches: Dict[MafVariantKey, List[str]],
    maf_position_matches: Dict[Tuple[str, int, int], List[PositionRecord]],
    column_name: str,
    stats_column: str,
    variant_count_column: str,
    total_allele_count: int,
    missing_value: str,
) -> None:
    with open_maybe_gzip(maf_path, "rt") as in_handle:
        reader = csv.DictReader(in_handle, delimiter=maf_delimiter)
        required = {"Chromosome", "Start_Position"}
        if reader.fieldnames is None or not required.issubset(set(reader.fieldnames)):
            missing = sorted(required.difference(set(reader.fieldnames or [])))
            raise SystemExit(f"MAF file missing columns: {', '.join(missing)}")

        fieldnames = list(reader.fieldnames)
        if column_name not in fieldnames:
            fieldnames.append(column_name)
        if stats_column not in fieldnames:
            fieldnames.append(stats_column)

        out_handle: TextIO
        if output_path == "-":
            out_handle = sys.stdout
        else:
            out_handle = open_maybe_gzip(output_path, "wt")

        try:
            writer = csv.DictWriter(
                out_handle,
                fieldnames=fieldnames,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()

            for row in reader:
                info_value, _match_status = annotate_row(
                    row=row,
                    exact_matches=exact_matches,
                    position_matches=position_matches,
                    maf_exact_matches=maf_exact_matches,
                    maf_position_matches=maf_position_matches,
                    missing_value=missing_value,
                )
                row[column_name] = info_value
                row[stats_column] = association_stats_string(
                    info_value=info_value,
                    variant_count_value=row.get(variant_count_column),
                    total_allele_count=total_allele_count,
                    missing_value=missing_value,
                )
                writer.writerow(row)
        finally:
            if output_path != "-":
                out_handle.close()


def main() -> int:
    args = parse_args()
    maf_delimiter = detect_delimiter(args.maf)
    exact_matches, position_matches, maf_exact_matches, maf_position_matches = load_vcf_annotations(args.vcf)
    stats_column = args.stats_column or f"{args.column_name}_Stats"

    annotate_maf(
        maf_path=args.maf,
        maf_delimiter=maf_delimiter,
        output_path=args.output,
        exact_matches=exact_matches,
        position_matches=position_matches,
        maf_exact_matches=maf_exact_matches,
        maf_position_matches=maf_position_matches,
        column_name=args.column_name,
        stats_column=stats_column,
        variant_count_column=args.variant_count_column,
        total_allele_count=args.total_allele_count,
        missing_value=args.missing_value,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
