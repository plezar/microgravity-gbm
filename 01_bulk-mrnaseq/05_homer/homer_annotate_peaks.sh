mfile="/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/01_bulk-mrnaseq/data/motifs/TEAD_motifs.txt"
degfile="/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/01_bulk-mrnaseq/results/homer_in/U87_DE_DOWN.txt"
outdir="/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/01_bulk-mrnaseq/results/homer_out/"

/Users/mzarodniuk/Documents/Software/HOMER/bin/annotatePeaks.pl tss hg38 \
  -size -400,100 \
  -list "${degfile}" \
  -m "${mfile}" \
  > "${outdir}/DEG_TEAD_annotation.txt"

/Users/mzarodniuk/Documents/Software/HOMER/bin/annotatePeaks.pl tss hg38 \
  -hist 400 \
  -ghist \
  -size -400,100 \
  -list "${degfile}" \
  -m "${mfile}" \
  > "${outdir}/DEG_TEAD_annotation_ghist.txt"



mfile="/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/01_bulk-mrnaseq/data/motifs/motif331.motif.txt"
degfile="/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/01_bulk-mrnaseq/results/homer_in/U87THP1_DE_UP.txt"
outdir="/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/01_bulk-mrnaseq/results/homer_out/"

/Users/mzarodniuk/Documents/Software/HOMER/bin/annotatePeaks.pl tss hg38 \
  -size -400,100 \
  -list "${degfile}" \
  -m "${mfile}" \
  > "${outdir}/DEG_SMAD3_annotation.txt"

/Users/mzarodniuk/Documents/Software/HOMER/bin/annotatePeaks.pl tss hg38 \
  -hist 400 \
  -ghist \
  -size -400,100 \
  -list "${degfile}" \
  -m "${mfile}" \
  > "${outdir}/DEG_SMAD3_annotation_ghist.txt"