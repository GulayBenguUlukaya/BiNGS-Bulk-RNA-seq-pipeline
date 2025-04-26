# BiNGS Bulk RNA-seq Pipeline
A streamlined pipeline for bulk RNA-seq analysis (first release)

Welcome to the BiNGS Bulk RNA-seq pipeline — a full workflow for downloading, preprocessing, and analyzing bulk RNA-seq data at the HPC cluster environment.

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


Pipeline Owner	Gulay Bengu Ulukaya
BiNGS Core	Tisch Cancer Institute, Mount Sinai
📫 Questions?
Open an Issue please and I will get back to you ASAP.


✅ Happy analyzing!


---

### Citation

Citation
To cite material from this pipeline in your publications, please use:

Gulay Bengu Ulukaya. (2025, April 28). BiNGS-Bulk-RNA-seq-pipeline: A streamlined pipeline for bulk RNA-seq analysis (first release). GitHub. https://github.com/GulayBenguUlukaya/BiNGS-Bulk-RNA-seq-pipeline.

A lot of time and effort went into the development of this pipeline. Citations help us understand the needs of the research community, gain recognition for our work, and attract further support for continued development. Thank you for citing this material if it helped you in your data analysis.

---

*These materials have been developed by members of the teaching team at the [The Bioinformatics for Next-Generation Sequencing (BiNGS) Core](https://bings.mssm.edu/). These are open access materials and permitted unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*


