This repository contains technical surveys on various topics. The notes are formatted as literate
Latex documents and evaluated with the [Litrepl](https://github.com/sergei-mironov/litrepl) tool.
The helper LLM model used in these writings was communicated using the
[aicli](https://github.com/sergei-mironov/aicli) interpreter.

<!--
``` sh
(cd tex;
FIRST=y
for f in $(ls -1 *tex | grep -v preamble.tex) ; do
  NM=`basename $f .tex`
  TEX=$f
  PDF=$NM.pdf
  echo -n "$NM [[TEX](./tex/$TEX)] [[PDF](./tex/$PDF)]"
  if test "$FIRST" = "y" ; then
    echo -n " | "
    FIRST=n
  fi
done
)
```
-->

<!--result-->
C3-inheritance [[TEX](./tex/C3-inheritance.tex)] [[PDF](./tex/C3-inheritance.pdf)] |
brothers-on-bridge [[TEX](./tex/brothers-on-bridge.tex)] [[PDF](./tex/brothers-on-bridge.pdf)]
<!--noresult-->
