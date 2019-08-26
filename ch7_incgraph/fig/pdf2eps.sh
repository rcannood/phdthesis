#!/bin/sh
# $Id: pdf2eps,v 0.01 2005/10/28 00:55:46 Herbert Voss Exp $
# Convert PDF to encapsulated PostScript.
# usage:
# pdf2eps <page number> <pdf file without ext>

pdfcrop "$1.pdf" "$1-temp"
pdftops -eps "$1-temp" "$1.eps"
rm  "$1-temp"
