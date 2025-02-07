This repository contains technical surveys on various topics. The notes are formatted as literate
Latex documents and evaluated with the [Litrepl](https://github.com/sergei-mironov/litrepl) tool.
The assistant LLM model employed in these writings was accessed via the
[aicli](https://github.com/sergei-mironov/aicli) interpreter.

<!--
``` sh
(cd tex;
FIRST=y
for f in $(ls -1 *tex | grep -v preamble.tex | grep -v '^_' | grep -v 'template') ; do
  NM=`basename $f .tex`
  TEX=$f
  PDF=$NM.pdf
  if test "$FIRST" = "y" ; then
    FIRST=n
  else
    echo -n " | "
  fi
  echo -n "$NM [[TEX](./tex/$TEX)] [[PDF](./tex/$PDF)]"
done
)
```
-->

<!--result-->
C3-inheritance [[TEX](./tex/C3-inheritance.tex)] [[PDF](./tex/C3-inheritance.pdf)] |
brothers-on-bridge [[TEX](./tex/brothers-on-bridge.tex)] [[PDF](./tex/brothers-on-bridge.pdf)] |
chernoff-bound-vs-qthreshhold [[TEX](./tex/chernoff-bound-vs-qthreshhold.tex)]
[[PDF](./tex/chernoff-bound-vs-qthreshhold.pdf)]
<!--noresult-->
