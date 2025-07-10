# JCAMPDXir
[![Build Status](https://github.com/Manarom/JCAMPDXir.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Manarom/JCAMPDXir.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://manarom.github.io/JCAMPDXir.jl)
## General description

This package is designed to read and write infrared spectra saved in JCAMP-DX (infrared) format.
It implements the main functionality to read and write JCAMP-DX files according to `5.01` JCAMP-DX 
version or lower.   

Full  documentation is available at  [documentation](https://manarom.github.io/JCAMPDXir.jl/)

### About JCAMP-DX file format

For a detailed overview of the JCAMP-DX infrared format, please refer to [JCAMP-DX](https://iupac.org/what-we-do/digital-standards/jcamp-dx/)

## Current state of the package

Currently, the package parses JCAMP-DX files written in  `(X++(Y..Y))` and `(XY...XY)` data line formats 
in a single or multiple blocks (each block must be embraced in `##TITLE...##END`).  
Supported data compression methods:
- for reading: no data compression, integer compression, `SQZ`,`PAC`,`DIF` and `DUP` (file can use various combinations of  these compression formats simultaneously)
- for writing: simple integer compression
For writing, the package also supports various x- and y- units conversions

The package was tested on all IR-spectra examples from python package
[jcamp](https://github.com/nzhagen/jcamp.git)


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
## Contact
To contact me use [GitHub repository](https://github.com/Manarom).