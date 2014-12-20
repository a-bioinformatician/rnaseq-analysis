
NGS Plot for plotting density of coverage at features:
https://code.google.com/p/ngsplot/ - 5' UTR and 3' UTR plots
http://crazyhottommy.blogspot.no/2013/04/how-to-make-tss-plot-using-rna-seq-and.html

Analysis methods for comparison of protocols:
http://www.ncbi.nlm.nih.gov/pubmed/20711195

RNA-seq quality assessment (http://www.ncbi.nlm.nih.gov/pubmed/20711195, http://www.biomedcentral.com/1471-2164/15/419):
1. Complexity (read starting alignment locations varied = duplication)
2. Strand specificity
3. Evenness of coverage across genes.
4. Comparison to known transcript structure  - 5' end coverage vs. 3' end coverage and segmentation in coverage across gene.

Metrics:
- coefficient of variation for top 50% expressed genes.
- Percent antisense
- Coverage (%) across gene body (5)
- 5' and 3' end coverage.

- Fraction of bases not covered by reads for each gene in the genome against the fraction of total reads for that gene in the library.
- Number of segments (separated by at least five bases of zero coverage) weighted by average expression of each gene.
