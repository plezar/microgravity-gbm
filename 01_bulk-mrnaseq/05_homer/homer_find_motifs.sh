BIN="/Users/mzarodniuk/Documents/Software/HOMER/bin"
IN_DIR="/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/01_bulk-mrnaseq/results/homer_in/"
OUT_DIR="/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/01_bulk-mrnaseq/results/homer_out"

${BIN}/findMotifs.pl "${IN_DIR}/U87_DE_ALL.txt" human "${OUT_DIR}/U87_DE_ALL" -start -400 -end 100 -len 8,10 -p 1
${BIN}/findMotifs.pl "${IN_DIR}/U87_DE_UP.txt" human "${OUT_DIR}/U87_DE_UP" -start -400 -end 100 -len 8,10 -p 1
${BIN}/findMotifs.pl "${IN_DIR}/U87_DE_DOWN.txt" human "${OUT_DIR}/U87_DE_DOWN" -start -400 -end 100 -len 8,10 -p 1


${BIN}/findMotifs.pl "${IN_DIR}/U87THP1_DE_ALL.txt" human "${OUT_DIR}/U87THP1_DE_ALL" -start -400 -end 100 -len 8,10 -p 1
${BIN}/findMotifs.pl "${IN_DIR}/U87THP1_DE_UP.txt" human "${OUT_DIR}/U87THP1_DE_UP" -start -400 -end 100 -len 8,10 -p 1
${BIN}/findMotifs.pl "${IN_DIR}/U87THP1_DE_DOWN.txt" human "${OUT_DIR}/U87THP1_DE_DOWN" -start -400 -end 100 -len 8,10 -p 1
