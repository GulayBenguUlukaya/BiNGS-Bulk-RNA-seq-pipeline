BiNGS Bulk RNA-Seq Pipeline
Welcome to the BiNGS Bulk RNA-seq pipeline — a full workflow for downloading, preprocessing, and analyzing bulk RNA-seq data at the HPC cluster environment.

📦 Repository Structure
bash
Copy
Edit
BiNGS-Bulk-RNA-seq-pipeline/
├── data_rna/
│   ├── raw/
│   ├── preprocessed/
│   └── processed/
├── supporting_files/
│   └── annotation/
├── scripts/
│   ├── preprocessing/
│   ├── clustering/
│   ├── differential_expression/
│   └── functional_analysis/
├── README.md
└── LICENSE
🖥️ How to Set Up on HPC
Clone the repository:

bash
Copy
Edit
cd /your/working/directory/
git clone https://github.com/GulayBenguUlukaya/BiNGS-Bulk-RNA-seq-pipeline.git
Load the R module (R/4.2.0):

bash
Copy
Edit
ml R/4.2.0
.libPaths() inside R scripts should be set to your cluster-specific R library paths.

🔄 Pipeline Steps
1. Data Download and Annotation Setup
Prepare annotation files under supporting_files/annotation/.

Supported organisms:

Homo sapiens

Mus musculus

If multiple organisms exist, you will be asked to specify.

2. Preprocessing
Raw FASTQ data is processed externally (e.g., Salmon quantification).

Processed data are compiled under data_rna/processed/compiled_data/.

3. Clustering and Quality Control
PCA plots colored by condition and batch.

Distance heatmaps across samples.

Batch effect correction (optional).

Figures saved under figures/clustering/.

4. Differential Expression Analysis
Differential expression calculated with DESeq2.

Volcano plots auto-generated.

Results saved under analysis/differential_expression/.

5. Functional Enrichment Analysis
GSEA run against KEGG, Reactome, WikiPathways, and Hallmark gene sets.

Saves:

GSEA result tables (*_gsea_results.csv).

Enrichment curves (per significant pathway).

Top pathways barplots (NES scale; significant pathways asterisked and bolded).

🧹 Good Practices
Always confirm R/4.2.0 modules are loaded.

Ensure all input .fst and .csv files exist in compiled_data/.

Use HPC batch submission when running heavy analyses.

🔖 Metadata

Item	Description
Language	R 4.2.0
Pipeline Owner	Gulay Bengu Ulukaya
License	MIT License
Last Updated	2025
BiNGS Core	Tisch Cancer Institute, Mount Sinai
📫 Questions?
Open an Issue or contact GulayBenguUlukaya directly.

"We don’t just analyze data. We tell the story of biology."

✅ Happy analyzing!
