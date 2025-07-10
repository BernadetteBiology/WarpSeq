#!/bin/bash

# If you download files from the SRA database, they usually end in SRR00000000_1.fastq or SRR00000000_2.fastq. To quickly rename them to SRR00000000_R1.fastq or SRR00000000_R2.fastq use this script.

# Set to 1 to preview file changes (i.e. DRY_RUN=1), 0 to actually rename files (DRY_RUN=0)
DRY_RUN=1

for f in *_1.fastq *_2.fastq; do
  [[ -e "$f" ]] || continue  # skip if no match

  newname=$(echo "$f" | sed -E 's/_([12])\.fastq$/_R\1.fastq/')
  
  if [[ $DRY_RUN -eq 1 ]]; then
    echo "Would rename: $f -> $newname"
  else
    mv "$f" "$newname"
    echo "Renamed: $f -> $newname"
  fi
done
