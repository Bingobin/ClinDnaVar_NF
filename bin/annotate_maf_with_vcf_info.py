#!/usr/bin/env python3

import argparse
import csv
import gzip
import sys
from collections import defaultdict
from typing import DefaultDict, Dict, Iterable, List, Optional, TextIO, Tuple, TypeVar


VariantKey = Tuple[str, int, str, str]
PositionKey = Tuple[str, int]
PositionRecord = Tuple[str, str, str]
MafVariantKey = Tuple[str, int, int, str, str]
T = TypeVar("T")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Annotate a MAF/TSV file with raw INFO values from a VCF. "
            "The script first attempts exact matching by raw VCF-style REF/ALT when the "
            "MAF contains REF/ALT columns. It also converts VCF variants to MAF/ANNOVAR-"
            "style coordinates so indels represented with '-' and shifted positions can "
            "be matched through Chromosome, Start_Position, End_Position, "
            "Reference_Allele, and Tumor_Seq_Allele2. If allele-aware matching fails, it "
            "falls back to position-only matching when the VCF has exactly one record at "
            "that position."
        )
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
        help="Output column name for raw VCF INFO values (default: ChinaMAP)",
    )
    parser.add_argument(
        "--missing-value",
        default="NA",
        help="Value written when no INFO annotation is found (default: NA)",
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
        return chrom, pos + 1, pos, "-", alt[1:]

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

    with open_maybe_gzip(vcf_path, "rt") as handle:
        for raw_line in handle:
            if not raw_line.strip() or raw_line.startswith("##"):
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
            info = fields[7].strip() or "."

            for alt in alt_values:
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
                writer.writerow(row)
        finally:
            if output_path != "-":
                out_handle.close()


def main() -> int:
    args = parse_args()
    maf_delimiter = detect_delimiter(args.maf)
    exact_matches, position_matches, maf_exact_matches, maf_position_matches = load_vcf_annotations(args.vcf)

    annotate_maf(
        maf_path=args.maf,
        maf_delimiter=maf_delimiter,
        output_path=args.output,
        exact_matches=exact_matches,
        position_matches=position_matches,
        maf_exact_matches=maf_exact_matches,
        maf_position_matches=maf_position_matches,
        column_name=args.column_name,
        missing_value=args.missing_value,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
