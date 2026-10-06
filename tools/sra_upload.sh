#!/bin/bash
# Run ON HOFFMAN2, on dtn2, inside screen:
#     bash sra_upload.sh uploads/<account_folder> [upload_list.tsv]
#
# Uploads every FASTQ the SRA sheet names (sra/upload_list.tsv: source path,
# name on NCBI, bytes, md5) into <account_folder>/microscale_xqtl/ on NCBI's
# preload FTP. Parent-panel files go up under their A1ref_ names; nothing on
# Hoffman is renamed or moved.
#
# The password is read from ~/.ncbi_netrc (chmod 600), never from the command
# line. Make it once with:
#   read -s -p "NCBI FTP password: " PW; umask 077
#   printf 'machine ftp-private.ncbi.nlm.nih.gov login subftp password %s\n' "$PW" > ~/.ncbi_netrc; unset PW
#
# Safe to re-run: each file resumes where it stopped (-C -), and a file that
# is already complete costs one round trip.
set -uo pipefail
ACCT=${1:?usage: bash sra_upload.sh uploads/<account_folder> [upload_list.tsv]}
LIST=${2:-$(dirname "$0")/upload_list.tsv}
export ACCT SUB=microscale_xqtl HOST=ftp-private.ncbi.nlm.nih.gov NETRC=~/.ncbi_netrc
LOG=~/sra_upload.log
[ -r "$NETRC" ] || { echo "make $NETRC first (see the header of this script)"; exit 1; }
[ -r "$LIST" ]  || { echo "no upload list at $LIST"; exit 1; }

up() {
  local src=$1 tgt=$2 try
  for try in 1 2 3 4 5; do
    if curl -sS --netrc-file "$NETRC" --ftp-create-dirs -C - -T "$src" \
         "ftp://$HOST/$ACCT/$SUB/$tgt"; then
      echo "$(date +%T) ok   $tgt"; return 0
    fi
    sleep $((try * 30))
  done
  echo "$(date +%T) FAIL $tgt"
}
export -f up

N=$(wc -l < "$LIST")
echo "$(date) uploading $N files to $ACCT/$SUB" | tee -a "$LOG"
cut -f1,2 "$LIST" | tr '\t\n' '\0\0' | xargs -0 -n 2 -P 4 bash -c 'up "$0" "$1"' | tee -a "$LOG"

ok=$(grep -c ' ok ' "$LOG"); fail=$(grep -c ' FAIL ' "$LOG")
echo "$(date) done. log: $LOG  (ok lines $ok, FAIL lines $fail; re-run to retry failures)"
echo "files now in the NCBI folder: $(curl -sS --netrc-file "$NETRC" -l "ftp://$HOST/$ACCT/$SUB/" | wc -l) of $N"
