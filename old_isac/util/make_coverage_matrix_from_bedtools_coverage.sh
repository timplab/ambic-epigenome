#!/bin/bash
### taken from https://stackoverflow.com/questions/34523920/how-to-merge-specific-columns-from-many-files-in-one-file
input="$@"
for f in $input;do
  dir=$(dirname "$f")
  base=$(basename "$f")
  sample=${base%%.*}
  samples="$samples $sample"
done

# save stuff in temp
header_tmp=$dir/.header.tmp
echo "#chromosome start end$samples" | tr " " "\t" > $header_tmp
rows_tmp=$dir/.rownames.tmp
cut -f1,2,3 ${1} > $rows_tmp
out_tmp=$dir/.matrix.tmp
awk '{a[FNR]=a[FNR]"\t"$4}END{for(i=1;i<=length(a);i++)print a[i]}' $input > $out_tmp
#awk '{a[FNR]=a[FNR]?a[FNR]"\t"$4:$4}END{for(i=1;i<=length(a);i++)print a[i]}' $input > $out_tmp

# output everything
cat $header_tmp
paste -d "" $rows_tmp $out_tmp

# clean up
rm $header_tmp $rows_tmp $out_tmp
