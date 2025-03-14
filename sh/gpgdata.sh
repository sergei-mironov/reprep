#!/bin/sh

NAME_EMAIL=$(gpg --list-keys --with-colons | grep "^uid" | awk -F ":" '{print $10}')

# Print name
echo $NAME_EMAIL | sed 's/^\(.*\) <[^>]*>$/\1/'
# Print email
echo $NAME_EMAIL | sed 's/.*<\([^>]*\)>.*/\1/'
