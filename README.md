# JCAMPDXir

[![Build Status](https://github.com/Manarom/JCAMPDXir.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Manarom/JCAMPDXir.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://manarom.github.io/JCAMPDXir.jl)

# General description

This package is designed to read and write infrared spectra saved in JCAMP-DX (infrared) format.
It implements some of the basic functionality to read  and write JCAMP-DX files
according to the `4.24` and also it tested with `5.01`.  

Full  documentation is available at  [documentation](https://manarom.github.io/JCAMPDXir.jl/)

# About JCAMP-DX file format

For a detailed overview of the JCAMP-DX infrared format, please refer to
[JCAMP-DX for infrared 4.24](https://iupac.org/what-we-do/digital-standards/jcamp-dx/)

Other JCAMP-DX parsers:

There is a parser, written in python [jcamp](https://github.com/nzhagen/jcamp.git), 
examples from this repository were used for testing this package

R-languge package [readJDX](https://github.com/bryanhanson/readJDX.git)

## Current state of the package

Currently, the package parses JCAMP-DX files written in  `(X++(Y..Y))` and `(XY...XY)` data line formats 
in a single or multiple blocks (each block must be embraced in `##TITLE...##END`).  
Supported data compression methods:
- for reading: no data comression, integer comression, `SQZ`,`PAC`,`DIF` and `DUP` (file can use various combinations of  these compression formats simultaniously)
- for writing: simple integer comression
For writing, the package also supports various x- and y- units conversions

The package was tested on all IR-spectra examples from python package [jcamp](https://github.com/nzhagen/jcamp.git)


## Quick start

```julia
import Pkg 
Pkg.add("https://github.com/Manarom/JCAMPDXir.jl.git")
using JCAMPDXir
data = read_jdx_file(file_name) # to read the file, data.x - x-values, data.y - y values, data.headers - file headers
# if it is known that all data rows (except the last one) has the same number of numbers, the following version is a liitle faster:
data = read_jdx_file(file_name, fixed_columns_number=true) 
write_jdx_file(x,y,"MKM","TRANSMITTANCE") # to write x - and y - data vectors of the sama size

```

## Contact
To contact me use [GitHub repository](https://github.com/Manarom).