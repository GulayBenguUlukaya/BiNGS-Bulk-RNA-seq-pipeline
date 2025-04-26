BiNGS Bulk RNA-Seq Pipeline
Welcome to the BiNGS Bulk RNA-seq pipeline — a full workflow for downloading, preprocessing, and analyzing bulk RNA-seq data at the HPC cluster environment. 🚀

📦 Repository Structure
kotlin
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
Clone the repository to your HPC environment:

bash
Copy
Edit
cd /your/working/directory/
git clone https://github.com/GulayBenguUlukaya/BiNGS-Bulk-RNA-seq-pipeline.git
Load R module (ensure R/4.2.0 is available on your HPC):

bash
Copy
Edit
ml R/4.2.0
Set .libPaths() inside the R scripts to point to your HPC-specific R library.

🔄 Pipeline Steps
1. Data Download and Annotation Setup
Prepare organism-specific annotation files under supporting_files/annotation/

Supported organisms:

Homo sapiens

Mus musculus

If multiple organisms exist, the script will ask you which one to use ✅

2. Preprocessing Scripts
Raw data → Preprocessed using tools like STAR, Salmon, etc. (handled outside this repo)

Compiled summary tables will be stored under processed/compiled_data/

3. Clustering and Quality Control (R Script)
Create PCA plots

Batch correction if needed

Distance matrix heatmaps

Save plots under figures/clustering/

4. Differential Expression (R Script)
Uses DESeq2

Results saved to analysis/differential_expression/

Volcano plots automatically generated ✅

5. Functional Enrichment Analysis (R Script)
Runs GSEA against KEGG, Reactome, WikiPathways, and Hallmark gene sets

Saves:

GSEA results tables ✅

GSEA enrichment curves ✅

Top pathways barplots:

NES plotted (−1 = blue ➔ 0 = white ➔ +1 = red)

Statistically significant pathways (padj < 0.01) are asterisked * and bolded ✅

🧹 Good Practices
Always check module versions on the HPC.

Make sure .fst and .csv input files exist in compiled_data/ before running.

Prefer running inside an R batch job script if datasets are large.

🔖 Pipeline Metadata

Item	Description
Language	R 4.2.0
Pipeline Owner	Gulay Bengu Ulukaya
License	MIT License
Last Updated	2025
BiNGS Core	Tisch Cancer Institute, Mount Sinai
📫 Questions?
Feel free to open an Issue or contact GulayBenguUlukaya directly on GitHub!

"We don’t just analyze data. We tell the story of biology."

✅ Happy analyzing! 🎨📊🧬
