#!/bin/sh
# Print lines between start and stop, ignoring empty lines

if test -n "$2" ; then
  B=$1
  E=$2
elif test -n "$1" ; then
  B="{{{ $1"
  E="}}} $1"
else
  echo "Usage: sedlines.sh (PAT|START STOP)" >&2
  exit 1
fi

sed -n "/$B/,/$E/{ /$B/b; /$E/b; p; }" | sed -e :a -e '/^\n*$/{$d;N;ba' -e '}' | sed '/./,$!d'
