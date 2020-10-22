#!/bin/bash
root=/kyber/Data/Nanopore/projects/ambic/sigma/reads/190225_choSigma_plasmid_target

while IFS=$',' read -r -a line || [[ -n "$line" ]]
do
  base=${line[0]}
  bc=${line[1]}
  for obj in $(find $root -maxdepth 1 -name "$bc*"); do
    newname=${obj/$bc/$base}
    mv $obj $newname
  done
done <<< "$(awk 'NR>1' $root/samples.csv)"
