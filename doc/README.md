# Xyce Documentation

## LaTeX

The Xyce documentation is written in
[LaTeX](https://en.wikibooks.org/wiki/LaTeX), so to create PDFs, you will need
a LaTeX distribution (look
[here](https://en.wikibooks.org/wiki/LaTeX/Installation) for more information).
Also, by default the Guides use the Koma-Script "report" class
([scrreprt](https://ctan.org/pkg/scrreprt)), so you will need to install it, if
it is not included in your LaTeX distribution.

## Building the documents

### Required modifications for Open Source users

The official Guides are formatted according to [Sandia](https://www.sandia.gov)
guidelines, using the "SANDreport" class, which is available only within
Sandia. However, users outside of Sandia can still build the guides by making
a couple of small modifications:
* For the Users' Guide, open the `Xyce_UG.tex` file.
* For the Reference Guide, open the `Xyce_RG.tex` file.
* Change the `\documentclass` and `\usepackage` lines at the top of each file
  to appear as follows:
```tex
%\documentclass[11pt,report]{SANDreport}
%\usepackage[sand]{optional}
\documentclass[11pt,letterpaper]{scrreprt}
\usepackage[report]{optional}
```
(note that `%` indicates a comment line in LaTeX).

### Building

In each document directory, there is a `Makefile` for use with [GNU
Make](https://www.gnu.org/software/make/). The file primarily sets some
environment variables, and invokes `latexmk` (which is likely a part of your
LaTeX package). Usage from the command line:
* `make` builds the document
* `make clean` removes auxiliary files created during the build process
* `make realclean` removes all files created during the build process,
  including the final PDF.

## Contributing corrections/improvements

If you wish to offer enhancements/improvements or corrections to the
documentation, please let us know (see the [CONTRIBUTING](../CONTRIBUTING.md)
document for more information). Note that, while the requirements on accepting
documentation changes are looser than accepting code changes, we still may be
limited in what we can incorporate directly.

