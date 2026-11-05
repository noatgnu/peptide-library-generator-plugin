#!/usr/bin/env python3
"""Generate a peptide library from exon/splicing data, or a novel tryptic digest from any FASTA.

Three library types, selected with --library_type:

- orf_fasta: translates rMATS-style exon inclusion/exclusion nucleotide sequences in all
  three forward reading frames into a full-ORF FASTA, with a tracking table recording each
  entry's gene, splice type, source row, event, frame, genomic location, and any other
  entry sharing its exact peptide sequence.
- translated_aa: builds a full-length AA-sequence FASTA from an already-translated
  splice-junction CSV (Transcript, AA_seq, genes, novel_junction_type, novel_junction_dpsi,
  NMD, TSL columns), with a matching tracking table.
- novel_tryptic: digests a source FASTA into tryptic peptides (fixed missed-cleavage count,
  length range), keeps only peptides that are not a substring of any sequence in a reference
  FASTA, and writes a target FASTA plus a reversed-sequence decoy FASTA.
"""

import argparse
import sys
from collections import Counter
from csv import DictWriter
from pathlib import Path

import pandas as pd

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

READ_KWARGS = {"dtype": str, "keep_default_na": False, "na_filter": False}

ORF_TRACKING_COLUMNS = [
    "Header", "Unique.ID", "Gene", "Gene.Name", "Cell.Tag", "Group.Label", "Source.File",
    "Row.ID", "Event", "Frame", "Chromosome", "Strand", "Location", "Comparison",
    "Peptide.Length", "Sequence", "Duplicate.Of",
]

TRANSLATED_AA_TRACKING_COLUMNS = [
    "Header", "Unique.ID", "Transcript", "Gene", "Cell.Tag", "Junction.Type", "dPSI",
    "NMD", "TSL", "Sequence.Length", "Sequence",
]

VALID_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWYXBZJUO*")


def translate_frame(nucleotide_seq: str, frame_offset: int) -> str:
    """Translate nucleotide_seq from frame_offset, stopping at the first in-frame stop codon."""
    seq = nucleotide_seq.upper()
    peptide = []
    for i in range(frame_offset, len(seq) - 2, 3):
        amino_acid = CODON_TABLE.get(seq[i:i + 3], "X")
        if amino_acid == "*":
            break
        peptide.append(amino_acid)
    return "".join(peptide)


def read_rows(path: str) -> list:
    """Read an AS_results CSV or XLSX file into a list of row dicts, every value kept as a string."""
    reader = pd.read_excel if path.lower().endswith(".xlsx") else pd.read_csv
    return reader(path, **READ_KWARGS).to_dict("records")


def find_comparison_column(column_names) -> str:
    """Return the column name ending in '.l2fc' (e.g. 'KD_vs_Control.l2fc' or 'KO_vs_WT.l2fc')."""
    for name in column_names:
        if name.endswith(".l2fc"):
            return name
    raise ValueError("No '.l2fc' column found; set --comparison_column explicitly")


def build_orf_entries(source_path: str, columns: dict, cell_tag: str, comparison_column: str, group_label_column: str, min_length: int):
    """Yield (header, peptide, tracking_row) for every ORF translation in source_path at least min_length long.

    columns maps the logical field name (e.g. "chromosome") to the actual column name in
    source_path, so any tool's naming can be used without renaming the source file.
    """
    rows = read_rows(source_path)
    if not rows:
        return
    comparison_col = comparison_column or find_comparison_column(rows[0].keys())
    comparison = comparison_col[: -len(".l2fc")] if comparison_col.endswith(".l2fc") else comparison_col

    for row in rows:
        gene, gene_name, row_id = row[columns["gene"]], row[columns["gene_name"]], row[columns["id"]]
        chromosome, strand = row[columns["chromosome"]], row[columns["strand"]]
        location = f"{chromosome}:{row[columns['exon_start']]}-{row[columns['exon_end']]}({strand})"
        group_label = row.get(group_label_column, "") if group_label_column else ""

        for event, seq_col in (("INC", "sequence_inclusion"), ("EXC", "sequence_exclusion")):
            nucleotide_seq = row[columns[seq_col]]
            for frame in (1, 2, 3):
                peptide = translate_frame(nucleotide_seq, frame - 1)
                if len(peptide) < min_length:
                    continue
                unique_id = f"{gene_name}_{cell_tag}_{row_id}_{event}_frame{frame}"
                header = f">cp|{unique_id}|{gene_name}_{comparison}"
                tracking_row = {
                    "Header": header,
                    "Unique.ID": unique_id,
                    "Gene": gene,
                    "Gene.Name": gene_name,
                    "Cell.Tag": cell_tag,
                    "Group.Label": group_label,
                    "Source.File": source_path,
                    "Row.ID": row_id,
                    "Event": event,
                    "Frame": frame,
                    "Chromosome": chromosome,
                    "Strand": strand,
                    "Location": location,
                    "Comparison": comparison,
                    "Peptide.Length": len(peptide),
                    "Sequence": peptide,
                }
                yield header, peptide, tracking_row


def annotate_duplicates(entries: list) -> None:
    """Set each tracking_row's "Duplicate.Of" to the semicolon-separated Unique.IDs of every
    other entry sharing its exact peptide sequence (blank if it is unique)."""
    ids_by_sequence = {}
    for _, peptide, tracking_row in entries:
        ids_by_sequence.setdefault(peptide, []).append(tracking_row["Unique.ID"])

    for _, peptide, tracking_row in entries:
        others = [uid for uid in ids_by_sequence[peptide] if uid != tracking_row["Unique.ID"]]
        tracking_row["Duplicate.Of"] = ";".join(others)


def generate_orf_fasta(args, output_dir: Path):
    """Translate exon inclusion/exclusion nucleotide sequences into a full-ORF FASTA.

    The source file's column names are mapped via the --*_column arguments (defaulted to
    rMATS's own AS_results column names); no particular splice-type vocabulary or filename
    convention is assumed.
    """
    if not args.exon_data_file:
        raise SystemExit("--exon_data_file is required for library_type=orf_fasta")
    if not args.cell_tag:
        raise SystemExit("--cell_tag is required for library_type=orf_fasta")

    columns = {
        "chromosome": args.chromosome_column,
        "strand": args.strand_column,
        "gene": args.gene_column,
        "gene_name": args.gene_name_column,
        "id": args.id_column,
        "sequence_inclusion": args.sequence_inclusion_column,
        "sequence_exclusion": args.sequence_exclusion_column,
        "exon_start": args.exon_start_column,
        "exon_end": args.exon_end_column,
    }

    entries = list(build_orf_entries(
        args.exon_data_file, columns, args.cell_tag, args.comparison_column, args.group_label_column, args.min_length,
    ))
    annotate_duplicates(entries)

    output_file = output_dir / "orf_library.fasta"
    tracking_file = output_dir / "orf_tracking.tsv"
    with open(output_file, "w") as out, open(tracking_file, "w", newline="") as track:
        tracking_writer = DictWriter(track, fieldnames=ORF_TRACKING_COLUMNS, delimiter="\t")
        tracking_writer.writeheader()
        for header, peptide, tracking_row in entries:
            out.write(f"{header}\n{peptide}\n")
            tracking_writer.writerow(tracking_row)

    duplicate_count = sum(1 for _, _, row in entries if row["Duplicate.Of"])
    print(
        f"Wrote {len(entries)} ORF entries from {args.exon_data_file} to {output_file} and "
        f"{tracking_file} ({duplicate_count} entries share their sequence with at least one other entry)"
    )


def build_translated_aa_entries(csv_path: str, cell_tag: str, min_length: int):
    """Yield (header, sequence, tracking_row) for every row in csv_path whose AA_seq is at least min_length long."""
    rows = pd.read_csv(csv_path, dtype=str, keep_default_na=False).to_dict("records")
    transcript_totals = Counter(row["Transcript"] for row in rows)
    occurrence = Counter()

    for row in rows:
        sequence = row["AA_seq"].rstrip("*")
        if len(sequence) < min_length:
            continue

        occurrence[row["Transcript"]] += 1
        transcript_id = row["Transcript"]
        if transcript_totals[transcript_id] > 1:
            transcript_id = f"{transcript_id}_dup{occurrence[transcript_id]}"

        unique_id = f"{row['genes']}_{cell_tag}_{transcript_id}"
        header = f">cp|{unique_id}|{row['genes']}_{row['novel_junction_type']}"
        tracking_row = {
            "Header": header,
            "Unique.ID": unique_id,
            "Transcript": row["Transcript"],
            "Gene": row["genes"],
            "Cell.Tag": cell_tag,
            "Junction.Type": row["novel_junction_type"],
            "dPSI": row["novel_junction_dpsi"],
            "NMD": row["NMD"],
            "TSL": row["TSL"],
            "Sequence.Length": len(sequence),
            "Sequence": sequence,
        }
        yield header, sequence, tracking_row


def generate_translated_aa_fasta(args, output_dir: Path):
    """Build a full-length translated-AA FASTA from a splice-junction CSV."""
    if not args.input_csv:
        raise SystemExit("--input_csv is required for library_type=translated_aa")
    if not args.cell_tag:
        raise SystemExit("--cell_tag is required for library_type=translated_aa")

    entries = list(build_translated_aa_entries(args.input_csv, args.cell_tag, args.min_length))

    output_file = output_dir / "translated_aa_library.fasta"
    tracking_file = output_dir / "translated_aa_tracking.tsv"
    with open(output_file, "w") as out, open(tracking_file, "w", newline="") as track:
        tracking_writer = DictWriter(track, fieldnames=TRANSLATED_AA_TRACKING_COLUMNS, delimiter="\t")
        tracking_writer.writeheader()
        for header, sequence, tracking_row in entries:
            out.write(f"{header}\n{sequence}\n")
            tracking_writer.writerow(tracking_row)

    print(f"Wrote {len(entries)} entries from {args.input_csv} to {output_file} and {tracking_file}")


def read_fasta(fasta_file: str):
    """Yield (header, sequence) tuples from a FASTA file; header includes the leading '>'."""
    header = None
    seq_parts = []
    with open(fasta_file, "rt") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line
                seq_parts = []
            else:
                seq_parts.append(line)
        if header is not None:
            yield header, "".join(seq_parts)


def digest_sequence(sequence: str, miss_cleavage: int, min_length: int, max_length: int):
    """Yield (start_pos, peptide) tryptic fragments spanning exactly miss_cleavage missed cleavages."""
    fragments = {}
    last_pos = 0
    for pos, aa in enumerate(sequence):
        if aa in ("K", "R"):
            fragments[last_pos] = sequence[last_pos:pos + 1]
            last_pos = pos + 1
    if last_pos < len(sequence):
        fragments[last_pos] = sequence[last_pos:]

    positions = list(fragments.keys())
    for i in range(len(positions)):
        if i + miss_cleavage < len(positions):
            combo = "".join(fragments[positions[i + offset]] for offset in range(miss_cleavage + 1))
            if min_length <= len(combo) <= max_length:
                yield positions[i], combo


def insert_accession_suffix(header: str, suffix: str) -> str:
    """Append suffix to a pipe-delimited header's accession (its 2nd field, e.g. "cp|ID|desc"
    becomes "cp|ID{suffix}|desc") so the accession alone stays unique; a header with no pipe
    structure gets suffix appended to the whole string instead.
    """
    parts = header.split("|")
    if len(parts) >= 2:
        parts[1] += suffix
        return "|".join(parts)
    return header + suffix


def build_peptide_library(fasta_file: str, miss_cleavage: int, min_length: int, max_length: int):
    """Digest every record in fasta_file into unique tryptic peptides within [min_length, max_length].

    Records that reuse another record's header text are tagged "_dup{N}" (1-based, in
    file order) so their peptides stay traceable to the exact source record.
    """
    records = list(read_fasta(fasta_file))
    header_totals = Counter(header for header, _ in records)

    peptides = {}
    skipped = 0
    occurrence = Counter()
    for header, sequence in records:
        occurrence[header] += 1
        label_header = header if header_totals[header] == 1 else insert_accession_suffix(header, f"_dup{occurrence[header]}")
        if set(sequence) - VALID_AMINO_ACIDS:
            skipped += 1
            continue
        for pos, combo in digest_sequence(sequence, miss_cleavage, min_length, max_length):
            if combo not in peptides:
                start, end = pos + 1, pos + len(combo)
                peptides[combo] = insert_accession_suffix(label_header, f"_{start}-{end}")
    return peptides, skipped


def filter_novel_peptides(peptides: dict, reference_fasta: str):
    """Keep only peptides whose collapsed (I/L) form is not a substring of any reference sequence."""
    content = "\n".join(sequence.replace("I", "L") for _, sequence in read_fasta(reference_fasta))
    return {
        sequence: label
        for sequence, label in peptides.items()
        if sequence.replace("I", "L") not in content
    }


def format_header(label: str) -> str:
    """Apply the cp accession convention: an "sp"/"cp" db tag is normalized to "cp", others are prefixed with "cp|"."""
    body = label[1:]
    if body.startswith("sp|") or body.startswith("cp|"):
        body = "cp" + body[2:]
    else:
        body = "cp|" + body
    return ">" + body


def write_library(peptides: dict, target_path: Path, decoy_path: Path):
    """Write the target FASTA library and its reversed-sequence decoy counterpart."""
    with open(target_path, "w") as target, open(decoy_path, "w") as decoy:
        for sequence, label in peptides.items():
            header = format_header(label)
            target.write(f"{header}\n{sequence}\n")
            decoy.write(f">REVERSE_{header[1:]}\n{sequence[::-1]}\n")


def generate_novel_tryptic(args, output_dir: Path):
    """Digest a FASTA into tryptic peptides and keep only those novel relative to a reference FASTA."""
    if not args.source_fasta:
        raise SystemExit("--source_fasta is required for library_type=novel_tryptic")
    if not args.reference_fasta:
        raise SystemExit("--reference_fasta is required for library_type=novel_tryptic")

    peptides, skipped = build_peptide_library(args.source_fasta, args.miss_cleavage, args.min_length, args.max_length)
    print(f"Digested {len(peptides)} unique peptides ({skipped} source records skipped for invalid characters)")
    novel_peptides = filter_novel_peptides(peptides, args.reference_fasta)
    print(f"{len(novel_peptides)} peptides are novel relative to {args.reference_fasta}")

    target_path = output_dir / "novel_peptides_target.fasta"
    decoy_path = output_dir / "novel_peptides_decoy.fasta"
    write_library(novel_peptides, target_path, decoy_path)
    print(f"Wrote {len(novel_peptides)} target/decoy entries to {target_path} and {decoy_path}")


def parse_args():
    parser = argparse.ArgumentParser(description="Generate a peptide library from exon/splicing data, or a novel tryptic digest from any FASTA")
    parser.add_argument("--library_type", required=True, choices=["orf_fasta", "translated_aa", "novel_tryptic"])
    parser.add_argument("--exon_data_file", default=None, help="Exon-level data file, one row per splicing event (orf_fasta)")
    parser.add_argument("--chromosome_column", default="chromosome", help="Column holding the chromosome (orf_fasta)")
    parser.add_argument("--strand_column", default="strand", help="Column holding the strand (orf_fasta)")
    parser.add_argument("--gene_column", default="gene", help="Column holding the gene ID (orf_fasta)")
    parser.add_argument("--gene_name_column", default="gene_name", help="Column holding the gene name (orf_fasta)")
    parser.add_argument("--id_column", default="ID", help="Column holding the row/event ID (orf_fasta)")
    parser.add_argument("--sequence_inclusion_column", default="sequence_inclusion", help="Column holding the exon-inclusion nucleotide sequence (orf_fasta)")
    parser.add_argument("--sequence_exclusion_column", default="sequence_exclusion", help="Column holding the exon-exclusion nucleotide sequence (orf_fasta)")
    parser.add_argument("--exon_start_column", default="exonStart_0base", help="Column holding the exon start coordinate (orf_fasta)")
    parser.add_argument("--exon_end_column", default="exonEnd", help="Column holding the exon end coordinate (orf_fasta)")
    parser.add_argument("--comparison_column", default="", help="Column holding the comparison name (orf_fasta); blank auto-detects a column ending in .l2fc")
    parser.add_argument("--group_label_column", default="", help="Optional column copied verbatim into the tracking table's Group.Label, e.g. an event-type label (orf_fasta); blank leaves it empty")
    parser.add_argument("--input_csv", default=None, help="Translated splice-junction CSV (translated_aa)")
    parser.add_argument("--cell_tag", default=None, help="Short dataset tag used in FASTA header IDs (orf_fasta, translated_aa)")
    parser.add_argument("--source_fasta", default=None, help="FASTA to digest tryptically (novel_tryptic)")
    parser.add_argument("--reference_fasta", default=None, help="Reference FASTA to check novelty against (novel_tryptic)")
    parser.add_argument("--miss_cleavage", type=int, default=2, help="Number of missed cleavages (novel_tryptic)")
    parser.add_argument("--max_length", type=int, default=30, help="Maximum peptide length (novel_tryptic)")
    parser.add_argument("--min_length", type=int, default=7, help="Minimum peptide/sequence length to keep")
    parser.add_argument("--output_folder", required=True, help="Output directory")
    return parser.parse_args()


def main():
    args = parse_args()
    output_dir = Path(args.output_folder)
    output_dir.mkdir(parents=True, exist_ok=True)

    # @step-if[id=mode]: Library type?
    if args.library_type == "orf_fasta":
        # @step[id=orf,from=mode:orf_fasta]: Generating full-ORF FASTA from exon data
        generate_orf_fasta(args, output_dir)
    elif args.library_type == "translated_aa":
        # @step[id=aa,from=mode:translated_aa]: Generating translated-AA FASTA from a splice-junction CSV
        generate_translated_aa_fasta(args, output_dir)
    else:
        # @step[id=novel,from=mode:novel_tryptic]: Digesting a novel tryptic-peptide target/decoy library
        generate_novel_tryptic(args, output_dir)

    # @step[from=orf+aa+novel]: Library generation complete
    print("Library generation complete.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
