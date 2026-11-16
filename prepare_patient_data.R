# Prepare harmonized AML patient-cohort data for the portal.   STATUS: planned (jamboree)
#
#   data/patients/<DATASET_ID>/   raw downloads per cohort (IDs from catalog/datasets.tsv)
#   ->  results/patients/<DATASET_ID>/
#         samples.tsv              one row per sample, fields from schema/sample_annotation.tsv
#                                  (OncoTree / NCIt disease, HGNC mutations, ELN risk, response, survival)
#         datasets.tsv             one row per omics layer (same columns as results/cell_lines/)
#         features/<key>.tsv       feature annotation (label, gene_symbol, entrez_id, ...)
#         values/<key>.tsv.gz      matrix: rows = features, columns = sample_id
#         drug_response/<key>.tsv.gz   ex vivo drug sensitivity, rows = drugs, columns = sample_id
#         manifest.tsv, provenance.tsv
#
# Cohorts (catalog/datasets.tsv, sample_type = patient): CPTAC adult AML, CPTAC-Kids First
# pediatric AML, Beat AML, TCGA LAML, TARGET AML, Wang 2022, Lee 2018, Malani 2022.
# Built-in analysis 1 (drug response biomarkers) uses the cohorts with paired RNA/protein and
# drug sensitivity.
#
# Usage: Rscript prepare_patient_data.R     (from the repo root)

message("prepare_patient_data.R: not implemented yet; see the plan at the top of this file.")
quit(status = 0)

# One section per cohort, each ending in the layout above:
# 1. read raw files from data/patients/<DATASET_ID>/
# 2. map clinical fields to schema/sample_annotation.tsv, values to schema/vocabularies.tsv
# 3. map feature IDs to HGNC symbols (keep the source ID in features/<key>.tsv)
# 4. write results/patients/<DATASET_ID>/ and update the cohort's status in catalog/datasets.tsv
