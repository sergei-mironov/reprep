This repository contains technical surveys on various topics. The notes are formatted as literate
Latex documents and evaluated with the [Litrepl](https://github.com/sergei-mironov/litrepl) tool.
The helper LLM model used in these writings was communicated using the
[aicli](https://github.com/sergei-mironov/aicli) interpreter.

<!--
``` sh
(cd tex;
for f in $(ls -1 *tex | grep -v preamble.tex) ; do
  NM=`basename $f .tex`
  TEX=$f
  PDF=$NM.pdf
  echo -n "$NM [[TEX](./tex/$TEX)] [[PDF](./tex/$PDF)]"
done
)
```
-->

<!--result-->
brothers-on-bridge [[TEX](./tex/brothers-on-bridge.tex)] [[PDF](./tex/brothers-on-bridge.pdf)]
<!--noresult-->
