#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
MLSTDIR="$0"
BLASTDIR="$DIR/../db/blast"
BLASTFILE="$BLASTDIR/mlst.fa"

mkdir -p "$BLASTDIR"
rm -f "$BLASTFILE"

#for N in $(find $MLSTDIR -maxdepth 1 | grep -v '_2$'); do
for N in $(find $MLSTDIR -mindepth 1 -maxdepth 1 -type d); do
  SCHEME=$(basename $N)
  echo "Adding: $SCHEME"
  cat "$MLSTDIR"/$SCHEME/*.tfa \
  	| grep -v 'not a locus'  \
  	| sed -e "s/^>/>$SCHEME./" \
  	>> "$BLASTFILE"
done

makeblastdb -hash_index -in "$BLASTFILE" -dbtype nucl -title "PubMLST" -parse_seqids

echo "Created BLAST database for $BLASTFILE"
