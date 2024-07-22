# RNAseq analysis
## Horizontal transfer in CTVT

This pipeline performs the initial RNAseq gene expression element of the horizontal transfer study.

### 
We did PacBio sequencing after we had these RNAseq results. We were able to identify more informative sites
using read variant phasing. Therefore the RNAseq analysis was updated to include the extra variants.
Code to run the final analysis is in the "updated" folder.

#### Requirements
There's no container for this analysis. Instead I set up a Micromamba environment, which
has installed `featureCounts` (from `subread`), `DESeq` and `R` and its required packages.
The `requirements.txt` file contains the list of packages and versions.
