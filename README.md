# BiNGS Bulk RNA-seq Pipeline
A streamlined pipeline for bulk RNA-seq analysis (first release)

BiNGS Core	Tisch Cancer Institute, Mount Sinai

Welcome to the BiNGS Bulk RNA-seq pipeline — a full workflow for downloading, preprocessing, and analyzing bulk RNA-seq data through minerva terminal.

🖥️ How to Set Up on HPC

Clone the repository through terminal:

>> cd /your/working/directory/

>> git clone https://github.com/GulayBenguUlukaya/BiNGS-Bulk-RNA-seq-pipeline.git

🔄 Pipeline Steps

Running .sh scripts in terminal:

>> bash /path/to/script/......sh
  
  
Running .R scripts in terminal:

>> ml R/4.2.0
  
>> R
  
>> source("/path/to/script/......R")
  

🧹 Good Practices

Don't change anything inside the scripts. Answer the questions in terminal and press enter.

🔖 Metadata

When you run the A1.prepare_files.sh script, it will copy over an example metadata and prompt you to edit it. 

Make sure to fill all metadata columns with your dataset info. If your dataset is single-ended, file_path_2 and file_name_2 can be blank.


# 📫 Questions?
Open an **Issue** on this Github repo please and I will get back to you ASAP.




---

### Citation

Citation
To cite material from this pipeline in your publications, please use:

Gulay Bengu Ulukaya. (2025, April 26). BiNGS-Bulk-RNA-seq-pipeline: A streamlined pipeline for bulk RNA-seq analysis (first release). GitHub. https://github.com/GulayBenguUlukaya/BiNGS-Bulk-RNA-seq-pipeline.

A lot of time and effort went into the development of this pipeline. Thank you for citing this material if it helped you in your data analysis.

---

*These materials have been developed by members of the teaching team at the [The Bioinformatics for Next-Generation Sequencing (BiNGS) Core](https://bings.mssm.edu/). These are open access materials and permitted unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*


