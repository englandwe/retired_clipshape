#!/bin/bash

#little one-liner so the positions of - strand transcripts are reversed - 50 becomes -50 and so on.
awk -F '\t' 'BEGIN {OFS="\t";}{if ($6 == "-") {print -$1,$2,$3,$4,$5,$6,$7,$8,$9} else {print $0}}' "$1" > "$1".reversed

