
# JCAMPDXir.jl

## General description

This package is designed to read and write infrared spectra saved in JCAMP-DX (infrared) format.
It implements the main functionality to read and write JCAMP-DX files according to the `5.01`.   

## About JCAMP-DX file format
For a detailed overview of the JCAMP-DX infrared format, please refer to
[JCAMP-DX for infrared](https://iupac.org/what-we-do/digital-standards/jcamp-dx/)


## Current state of the package
Currently, the package parses JCAMP-DX files written in  `(X++(Y..Y))` and `(XY...XY)` data line formats 
in a single or multiple blocks (each block must be embraced in `##TITLE...##END`).  
Supported data compression methods:
- for reading: no data compression, integer comression, `SQZ`,`PAC`,`DIF` and `DUP` (file can use various combinations of  these compression formats simultaneously)
- for writing: simple integer compression, the package also supports various x- and y- data units conversions

The package was tested for reading all IR-spectra from python package [jcamp](https://github.com/nzhagen/jcamp.git)

## Quick start
```julia
import Pkg 
Pkg.add("https://github.com/Manarom/JCAMPDXir.jl.git")
using JCAMPDXir
(x,y,headers,validation) = read_jdx_file(file_name) 
# to read the file, x - x-values, y - y values, headers - file headers, validation - jcamp specification checks

write_jdx_file(x,y,"MKM","TRANSMITTANCE") 
# to write x - and y - data vectors of the sama size

```
