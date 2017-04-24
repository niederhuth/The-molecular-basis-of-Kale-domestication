---
Daily status of the analysis.
---
### To do
1. Add all genes of interest to gof.csv
2. Add KEGG analysis
3. SNP calling and analysis? (Genetic distance between accessions?)
4. How does gene expression PCA compare to genetic PCA?
5. Other QC & global analyses?

### April 24, 2017
1. Update trim_align.sh

### Mar 17, 2017
1. Edits to deg.R adding new analyses
2. Added plotting of genes of interest

### Jan 12, 2017
1. Added deg.R (for calling de genes, etc)
2. Major edits.

#### Dec. 12, 2016
1. Data is phred64
2. Optimized preprocessing and alignment of data.
 * 5-prime and 3-prime trimming increased mapping over 3-trimming alone.
 * Trimming to quality score of 20 performed better than trimming to quality score of 13.
 * High occurance (>90%) of "A" at first base indicated some sort of contamination in read even after trimming. Removing first base from all reads greatly improved mapping further.
 * Keeping reads <30bp (min length 25bp) gives a weird GC content skew. Filter out reads <30bp.
3. Run trim_align.sh on all data.
