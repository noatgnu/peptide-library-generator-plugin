# Examples

Put sample input files here for manual testing and for CI.

`params.json` holds the `--param`/`--params-json` values `cauldron job run` uses against this
plugin's inputs (see `plugin.yaml`'s `inputs:` section, or run
`cauldron plugin inputs peptide-library-generator` once installed). CI runs this automatically via
`.github/workflows/test-plugin.yml`; it's skipped until this file has real values.

```json
{
  "library_type": "orf_fasta",
  "exon_data_file": "examples/as_results/AS_results_fpkm_JCEC_SE_human.seq.deseq.csv",
  "cell_tag": "demo",
  "min_length": 7
}
```

`as_results/AS_results_fpkm_JCEC_SE_human.seq.deseq.csv` is a small synthetic rMATS-style
skipped-exon results file (2 genes, inclusion/exclusion nucleotide sequences) used to exercise
the Full ORF FASTA library type; its column names happen to match the `*_column` inputs'
defaults, so none need to be set here, but any differently-named exon data file can be used by
overriding those inputs instead of renaming the file. To exercise the other two types manually,
swap in:

- `translated_aa`: `input_csv` pointing at a CSV with `Transcript, AA_seq, genes,
  novel_junction_type, novel_junction_dpsi, NMD, TSL` columns, plus `cell_tag`.
- `novel_tryptic`: `source_fasta` and `reference_fasta` pointing at any two FASTA files, plus
  `miss_cleavage`/`max_length`.
