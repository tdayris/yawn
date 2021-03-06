#!/bin/awk -f

BEGIN {
  FS="\t"
}

$1 !~ "^#.*" {
  gene=$8
  sub("^.*Gene.refGene=", "", gene)
  sub(";.*$", "", gene)
  sub("\\\\x3b", ";", gene)
  dp=$8
  sub("^.*DP=", "", dp)
  sub(";.*$", "", dp)
  print $1 FS $2 FS $2 + length($5) - 1 FS gene FS dp FS "."
}
