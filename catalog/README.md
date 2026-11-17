# Dataset catalog

`datasets.tsv` is the inventory of the datasets integrated into the resource, one row per dataset.
Only datasets that are processed and in the portal are listed. The portal shows it on the **Datasets** page.

| column | meaning |
|---|---|
| `dataset_id` | stable upper-case ID; also the folder name under `data/<type>/` and `results/<type>/` once ingested |
| `sample_type` | `patient`, `cell_line`, or both (`;`-separated) |
| `population` | `adult`, `pediatric`, `mixed` |
| `modalities` | `;`-separated, from `schema/vocabularies.tsv` (`modality`) |
| `reference` / `source` / `accession` | citation, where it comes from, and its accession, if any |
| `used_by` | portal analysis (or analyses) that use it |
| `status` | `planned` → `downloaded` → `processed` (in `results/`) |

Blank `accession` cells are still to be filled in during the data inventory. Do not guess them.
