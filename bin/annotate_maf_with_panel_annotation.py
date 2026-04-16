#!/usr/bin/env python3

import argparse
import csv
import gzip
import sys
from collections import Counter, defaultdict
from typing import Dict, Iterable, List, TextIO


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Annotate a MAF/TSV file with panel annotation and per-variant occurrence count. "
            "Panel annotation is assigned by overlap between MAF Chromosome/Start_Position/"
            "End_Position and panel intervals. Variant count is computed using "
            "Chromosome, Start_Position, End_Position, Reference_Allele, and "
            "Tumor_Seq_Allele2."
        )
    )
    parser.add_argument("maf", help="Input MAF/TSV file")
    parser.add_argument("panel", help="Panel annotation CSV/TSV file")
    parser.add_argument(
        "--panel-chr-column",
        default="panel_chr",
        help="Chromosome column name in panel annotation file (default: panel_chr)",
    )
    parser.add_argument(
        "--panel-start-column",
        default="panel_start",
        help="Start column name in panel annotation file (default: panel_start)",
    )
    parser.add_argument(
        "--panel-end-column",
        default="panel_end",
        help="End column name in panel annotation file (default: panel_end)",
    )
    parser.add_argument(
        "--panel-annotation-column",
        default="annotation",
        help="Annotation column name in panel annotation file (default: annotation)",
    )
    parser.add_argument(
        "--annotation-output-column",
        default="Panel_Annotation",
        help="Output column name for panel annotation (default: Panel_Annotation)",
    )
    parser.add_argument(
        "--count-output-column",
        default="Variant_Count",
        help="Output column name for variant occurrence count (default: Variant_Count)",
    )
    return parser.parse_args()


def open_maybe_gzip(path: str, mode: str = "rt") -> TextIO:
    if path.endswith(".gz"):
        return gzip.open(path, mode, newline="")
    return open(path, mode, newline="")


def detect_delimiter(path: str) -> str:
    with open_maybe_gzip(path, "rt") as handle:
        first_line = handle.readline()
    return "\t" if "\t" in first_line else ","


def normalize_chr(value: str) -> str:
    return str(value).strip()


def parse_int(value: str) -> int:
    return int(float(str(value).strip()))


def load_panel_annotations(
    panel_path: str,
    delimiter: str,
    chr_col: str,
    start_col: str,
    end_col: str,
    annotation_col: str,
) -> Dict[str, List[dict]]:
    panel_by_chr: Dict[str, List[dict]] = defaultdict(list)
    with open_maybe_gzip(panel_path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        required = {chr_col, start_col, end_col, annotation_col}
        if reader.fieldnames is None or not required.issubset(set(reader.fieldnames)):
            missing = sorted(required.difference(set(reader.fieldnames or [])))
            raise SystemExit(f"Panel annotation file missing columns: {', '.join(missing)}")

        for row in reader:
            chrom = normalize_chr(row[chr_col])
            start = parse_int(row[start_col])
            end = parse_int(row[end_col])
            annotation = row[annotation_col].strip() or "NA"
            panel_by_chr[chrom].append(
                {"start": start, "end": end, "annotation": annotation}
            )

    for chrom in panel_by_chr:
        panel_by_chr[chrom].sort(key=lambda item: (item["start"], item["end"]))
    return panel_by_chr


def variant_key(row: dict) -> tuple:
    return (
        row["Chromosome"],
        row["Start_Position"],
        row["End_Position"],
        row["Reference_Allele"],
        row["Tumor_Seq_Allele2"],
    )


def count_variants(maf_path: str, delimiter: str) -> Counter:
    counts: Counter = Counter()
    with open_maybe_gzip(maf_path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        required = {
            "Chromosome",
            "Start_Position",
            "End_Position",
            "Reference_Allele",
            "Tumor_Seq_Allele2",
        }
        if reader.fieldnames is None or not required.issubset(set(reader.fieldnames)):
            missing = sorted(required.difference(set(reader.fieldnames or [])))
            raise SystemExit(f"MAF file missing columns: {', '.join(missing)}")
        for row in reader:
            counts[variant_key(row)] += 1
    return counts


def find_panel_annotation(panel_by_chr: Dict[str, List[dict]], chrom: str, start: int, end: int) -> str:
    matches: List[str] = []
    for interval in panel_by_chr.get(chrom, []):
        if interval["start"] > end:
            break
        if interval["end"] >= start:
            matches.append(interval["annotation"])

    if not matches:
        return "NA"
    return "|".join(dict.fromkeys(matches))


def annotate_maf(
    maf_path: str,
    maf_delimiter: str,
    panel_by_chr: Dict[str, List[dict]],
    variant_counts: Counter,
    annotation_output_column: str,
    count_output_column: str,
) -> None:
    with open_maybe_gzip(maf_path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter=maf_delimiter)
        if reader.fieldnames is None:
            raise SystemExit("Input MAF file is missing a header")

        fieldnames: List[str] = list(reader.fieldnames)
        if annotation_output_column not in fieldnames:
            fieldnames.append(annotation_output_column)
        if count_output_column not in fieldnames:
            fieldnames.append(count_output_column)

        writer = csv.DictWriter(
            sys.stdout,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()

        for row in reader:
            chrom = normalize_chr(row["Chromosome"])
            start = parse_int(row["Start_Position"])
            end = parse_int(row["End_Position"])

            row[annotation_output_column] = find_panel_annotation(
                panel_by_chr, chrom=chrom, start=start, end=end
            )
            row[count_output_column] = str(variant_counts[variant_key(row)])
            writer.writerow(row)


def main() -> int:
    args = parse_args()
    maf_delimiter = detect_delimiter(args.maf)
    panel_delimiter = detect_delimiter(args.panel)

    panel_by_chr = load_panel_annotations(
        panel_path=args.panel,
        delimiter=panel_delimiter,
        chr_col=args.panel_chr_column,
        start_col=args.panel_start_column,
        end_col=args.panel_end_column,
        annotation_col=args.panel_annotation_column,
    )
    variant_counts = count_variants(args.maf, delimiter=maf_delimiter)

    annotate_maf(
        maf_path=args.maf,
        maf_delimiter=maf_delimiter,
        panel_by_chr=panel_by_chr,
        variant_counts=variant_counts,
        annotation_output_column=args.annotation_output_column,
        count_output_column=args.count_output_column,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
