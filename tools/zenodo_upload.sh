#!/bin/bash
# Run ON HOFFMAN2, on the data-transfer node (ssh dtn2.hoffman2.idre.ucla.edu),
# inside screen/tmux: the uploads take hours. Compute nodes have no internet.
#
# Uploads $STAGE (from zenodo_package.sh) into a DRAFT Zenodo record and
# reserves its DOI. It never publishes: the record stays private until someone
# presses Publish in the web UI, after the preprint is posted.
#
#   export ZENODO_TOKEN=...     # or keep it in ~/.zenodo_token, chmod 600;
#                               # zenodo.org -> Applications -> Personal access
#                               # tokens; scopes deposit:write and deposit:actions
#   bash tools/zenodo_upload.sh                  # new draft
#   DEPOSITION_ID=1234567 bash tools/zenodo_upload.sh   # resume into that draft
set -euo pipefail
[ -z "${ZENODO_TOKEN:-}" ] && [ -r ~/.zenodo_token ] && ZENODO_TOKEN=$(<~/.zenodo_token)
: "${ZENODO_TOKEN:?export ZENODO_TOKEN or put it in ~/.zenodo_token (chmod 600)}"
API=${ZENODO_API:-https://zenodo.org/api}   # sandbox: https://sandbox.zenodo.org/api
STAGE=${STAGE:-$SCRATCH/xqtl_zenodo}
README=${README:-$(dirname "$0")/zenodo_README.md}
AUTH=(-H "Authorization: Bearer $ZENODO_TOKEN")
J() { python3 -c "import json,sys; d=json.load(sys.stdin); print($1)"; }

if [ -z "${DEPOSITION_ID:-}" ]; then
  META=$(python3 - "$README" <<'EOF'
import json, sys, html
desc = open(sys.argv[1]).read()
print(json.dumps({"metadata": {
  "upload_type": "dataset",
  "title": "Microscale X-QTL: processed data for runs v7 and v6",
  "description": "<pre>" + html.escape(desc) + "</pre>",
  "license": "cc-by-4.0",
  "access_right": "open",
  "prereserve_doi": True}}))
EOF
)
  R=$(curl -sf "${AUTH[@]}" -H "Content-Type: application/json" -X POST "$API/deposit/depositions" -d "$META")
  DEPOSITION_ID=$(echo "$R" | J 'd["id"]')
  echo "new draft $DEPOSITION_ID, reserved DOI $(echo "$R" | J 'd["metadata"]["prereserve_doi"]["doi"]')"
fi
R=$(curl -sf "${AUTH[@]}" "$API/deposit/depositions/$DEPOSITION_ID")
BUCKET=$(echo "$R" | J 'd["links"]["bucket"]')
HAVE=$(echo "$R" | J '"\n".join(f["filename"] for f in d.get("files", []))')

cp "$README" "$STAGE/README.md"
for f in README.md MD5SUMS.txt INPUTS.tsv $(cd "$STAGE" && ls *.tar.gz); do
  if grep -qx "$f" <<<"$HAVE"; then echo "skip $f (already in draft)"; continue; fi
  echo "$(date +%T) put $f ($(du -h "$STAGE/$f" | cut -f1))"
  curl -sf --retry 3 "${AUTH[@]}" --upload-file "$STAGE/$f" "$BUCKET/$f" > /dev/null
done

## Zenodo's own checksums against ours
curl -sf "${AUTH[@]}" "$API/deposit/depositions/$DEPOSITION_ID/files" \
  | J '"\n".join(f["checksum"]+"  "+f["filename"] for f in d)' | sort -k2 > "$STAGE/zenodo_md5.txt"
if diff <(sort -k2 "$STAGE/MD5SUMS.txt") <(grep -v 'README.md\|MD5SUMS.txt' "$STAGE/zenodo_md5.txt") >/dev/null; then
  echo "checksums match"; else echo "CHECKSUM MISMATCH: compare MD5SUMS.txt and zenodo_md5.txt"; fi
echo "Draft (unpublished): https://zenodo.org/deposit/$DEPOSITION_ID"
echo "Add creators, related identifiers and the preprint link in the web UI. Do not publish yet."
