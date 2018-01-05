# enhancer-like regions
bedtools intersect \
  -a GSE77625/GSE77625_h3k4me1_chow.bed \
  -b GSE77625/GSE77625_h3k27ac_chow.bed \
  > enhancer-like_chow.bed

bedtools intersect \
  -a enhancer-like_chow.bed \
  -b extra/transcripts.bed \
  -v \
  > intergenic_enhancer-like_chow.bed
  
# First sort by chromosome and then start position:
bedtools sort -i intergenic_enhancer-like_chow.bed > intergenic_enhancer-like_chow_sorted.bed
bedtools sort -i extra/transcripts.bed > extra/transcripts_sorted.bed

# Then get the closest genes
bedtools closest \
  -a intergenic_enhancer-like_chow_sorted.bed \
  -b extra/transcripts_sorted.bed \
  -D b \
  > closest_transcripts_to_enhancer_chow.bed

bedtools flank \
  -l 1 \
  -r 0 \
  -s \
  -g extra/mm10.chromsizes \
  -i extra/transcripts.bed \
  > tsses.bed

bedtools intersect \
  -a GSE77625/GSE77625_h3k27ac_hfd.bed \
  -b GSE77625/GSE77625_h3k27ac_chow.bed \
  -v \
  > gained_h3k27ac.bed
  
bedtools intersect \
  -a tsses.bed \
  -b gained_h3k27ac.bed \
  -u \
  > tsses_with_gained_h3k27ac.bed


