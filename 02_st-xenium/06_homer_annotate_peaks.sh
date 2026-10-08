mfile="/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/01_bulk-mrnaseq/data/motifs/motifs_file.txt"
genefile="/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/data/XeniumPrimeHuman5Kpan_tissue_pathways_metadata.txt"
outdir="/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/results/homer_out"

/Users/mzarodniuk/Documents/Software/HOMER/bin/annotatePeaks.pl tss hg38 \
  -size -400,100 \
  -list "${degfile}" \
  -m "${mfile}" \
  > "${outdir}/DEG_ETS_annotation.txt"

#/Users/mzarodniuk/Documents/Software/HOMER/bin/annotatePeaks.pl tss hg38 \
#  -hist 400 \
#  -ghist \
#  -size -400,100 \
#  -list "${degfile}" \
#  -m "${mfile}" \
#  > "${outdir}/DEG_ETS_annotation_ghist.txt"