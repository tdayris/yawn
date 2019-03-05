#!/bin/awk -f

BEGIN {
  FS="\t"
  ref="../../GRCh38.primary_assembly.genome.chr_only.fa"
}

{
  print "samtools tview -d T -p " $1":"$2-10 " " ??? " --reference " ref
}

END {
  print "DONE!"
}
