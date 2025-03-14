RepRep
------

**RepRep** stands for **Reproducible Reports**. In this repository, we are experimenting with
creating reproducible documents on various technical topics, usually involving computations.

When we refer to a document as reproducible, we imply the following:
- The document contains source files enabling anyone to reconstruct the presentation.
- These sources are authenticated through a cryptographic signature.
- The sources specify build requirements in a universally unique manner.


In this project, we adhere to these principles by: (1) embedding the signed source archive in every
document we create; (2) incorporating [Nix flake](https://wiki.nixos.org/wiki/Flakes) environment
definitions and build instructions within these archives; (3) Leveraging
[Litrepl](https://github.com/sergei-mironov/litrepl) to execute computations and directly validate
results within the document source files.

For documents that incorporate AI assistant dialogues, we employ Litrepl alongside the
[aicli](https://github.com/sergei-mironov/aicli) interpreter as our means of communication.

Documents
---------

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
[[PDF](./tex/chernoff-bound-vs-qthreshhold.pdf)] | maxwell-sim [[TEX](./tex/maxwell-sim.tex)]
[[PDF](./tex/maxwell-sim.pdf)] | mechanical-bloch [[TEX](./tex/mechanical-bloch.tex)]
[[PDF](./tex/mechanical-bloch.pdf)] | pennylane-quantum-lock
[[TEX](./tex/pennylane-quantum-lock.tex)] [[PDF](./tex/pennylane-quantum-lock.pdf)] | stern-gerlach
[[TEX](./tex/stern-gerlach.tex)] [[PDF](./tex/stern-gerlach.pdf)]
<!--noresult-->
