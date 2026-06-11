# Input contracts

Every GenoBlocks block reads from one of three standardized inputs. Match these
formats and any block in the matching track will accept your data — that is what
makes the blocks interchangeable.

---

## 1. Groupfile (sample metadata) — required by both tracks

A CSV with a header. **The first column is the sample identifier**; every other
column is metadata you can group, color, model, or contrast by.

```csv
sample_name,condition,batch,timepoint,sex
SAMPLE_01,treated,run1,V1,Female
SAMPLE_02,control,run1,V1,Male
```

Rules:
- First column = sample/subject ID. It must match the sample column names in your
  counts file (bulk) or a metadata column in your Seurat object (single-cell).
- Continuous phenotypes (age, a clinical score) go in their own numeric columns —
  the modeling and DESeq2-against-continuous blocks read these directly.
- Missing values are allowed (`NA` / blank); blocks that need a complete column
  drop incomplete rows and tell you which.

A real, filled-in example is `examples/example_pain_metadata.csv` (clinical pain
phenotypes per subject). A blank template is `examples/groupfile_template.csv`.

## 2. Counts file — bulk track

Tab- or comma-delimited, with a header. The first two columns are gene identifiers,
all remaining columns are per-sample counts:

| column | meaning |
|--------|---------|
| `Symbol` | stable gene ID (e.g. ENSG…) — becomes the row name |
| `Gene_Symbol` | human-readable gene symbol (e.g. TP53) |
| `<sample columns>` | one column per sample; header must match the groupfile's first column |

Two flavors are used depending on the block:
- **raw integer counts** → DESeq2 / VST normalization (`02_normalize_cluster`)
- **CPM values**, same layout → QC and z-score clustering (`01_qc`)

Templates: `examples/counts_template_bulk.csv`.

## 3. Single-cell object — single-cell track

A **Seurat `.rds`** object. Blocks expect (and the early blocks create) these
conventions:

- An `RNA` assay; ADT/HTO assays when present are detected automatically.
- Sample identity in a metadata column — default `HTO_maxID` (override with
  `--sample_col` / `--hto_col`).
- Pool / batch in `hash_run` (override with `--pool_col`).
- Group/phenotype label in a metadata column (e.g. `PainGroup`, `Status`) used by
  the comparison blocks.

If your object doesn't yet have these, start at
`singlecell/01_build_object/merge_hashed_replicates.R`, which builds a
compliant merged object (mito/ribo metrics, demux, standardized metadata) from
raw hashed per-pool objects.

---

### Quick compatibility check

| You have… | Start at |
|-----------|----------|
| FASTQ/BAM only | align upstream first (or use the BS-seq / MaxOrigami pipelines), then bring counts here |
| raw bulk counts + groupfile | `bulk/01_qc` → `bulk/02_normalize_cluster` |
| per-pool hashed Seurat objects | `singlecell/01_build_object` |
| a merged Seurat object | `singlecell/02_normalization_variants` or `03_pca_investigation` |
| a labeled lineage subset | `singlecell/04_subclustering` or `05_figures_report` |
