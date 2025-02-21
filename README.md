## Airway RNA-seq Dataset
# Dataset Overview

The airway dataset is a RNA-seq experiment dataset that examines the response of human airway smooth muscle cells to glucocorticoid treatment (dexamethasone). The experiment was conducted to understand transcriptional response in asthma treatment. Source
Package: airway
Original Publication: Himes BE, Jiang X, Wagner P, Hu R, Wang Q, Klanderman B, Whitaker RM, Duan Q, Lasky-Su J, Nikolos C, Jester W, Johnson M, Panettieri R Jr, Tantisira KG, Weiss ST, Lu Q. "RNA-Seq Transcriptome Profiling Identifies CRISPLD2 as a Glucocorticoid Responsive Gene that Modulates Cytokine Function in Airway Smooth Muscle Cells." PLoS One. 2014 Jun 13;9(6):e99625. DOI: 10.1371/journal.pone.0099625

# Data Description

Experimental Design: 4 primary human airway smooth muscle cell lines treated with 1 μM dexamethasone or vehicle control (ethanol) Sample Size: 8 samples (4 treated, 4 control) Sequencing: Paired-end RNA-seq Platform: Illumina HiSeq 2000

# Data Structure
The dataset contains: Count Data: Gene-level read counts for each sample Sample Information: Metadata about experimental conditions Gene Annotations: Gene identifiers and additional annotation information
Sample Information Columns

cell: Cell line identifier dex: Treatment status (treated/untreated) 
