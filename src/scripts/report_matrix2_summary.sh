#!/bin/sh

zcat $1 | awk -v OFS=',' '{{print $3,$4,$5,$8,$9,$10,$11,$14}}' | sort | uniq -c | awk -v OFS=',' '{{print $2,$1}}' | sed 's/,/\t/g'