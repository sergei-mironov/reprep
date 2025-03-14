#!/bin/sh

set -e -x

if ! test -d $PROJECT_ROOT ; then
  echo "PROJECT_ROOT is not set" >&2
  exit 1
fi

FILE=$1
F=$(basename $FILE)
D=$(realpath $(dirname $FILE))
R="$PROJECT_ROOT"
A=$(echo $F | sed 's/\.[^.]*$//').tar.gz
T=$R/_sources/tmp
mkdir -p "$T" || true
cd "$D"

# 1. Use the -p flag with Makefile to print the database of rules and resolutions
make -p $F | expand | grep -v '^[# ]' > $T/alldeps

# 2. Scan the dependencies to get the dependencies of the file $F
cat $T/alldeps | grep "$F:" | awk -F ':' '{print $2}' | {
  for f in $(cat) ; do realpath --relative-to=$R $f ; done ;
} >$T/deps

cd "$R"

# 3. Run gzip to pack them into a .tar.gz archive
rm "$A" || true
tar -czf "$A" -T $T/deps

# 4. Place the archive into the file folder
mv "$A" "$D" || true
test -f "$D/$A"

