# JCAMPDXir

[![Build Status](https://github.com/Manarom/JCAMPDXir.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Manarom/JCAMPDXir.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://manarom.github.io/JCAMPDXir.jl)

# General description

This package is designed to read and write infrared spectra saved in JCAMP-DX (infrared) format.
The package was written for a specific task and cannot be considered as a complete implementation of all JCAMP-DX format specification. However, it implements some of the basic functionality to read  and write JCAMP-DX files
according to the `JCAMP-DX=4.24`. Full  documentation is available at  [documentation](https://manarom.github.io/JCAMPDXir.jl/)

# About JCAMP-DX file format

For a detailed overview of the JCAMP-DX infrared format, please refer to
[JCAMP-DX for infrared 4.24](https://iupac.org/what-we-do/digital-standards/jcamp-dx/)

There is an opensource parser, written in python [jcamp](https://github.com/nzhagen/jcamp.git), 
examples from this repository were used for testing this package

## Current state of the package

Currently, the package parses JCAMP-DX files written as  '(X++(Y..Y))' format (most common in spectroscopy) 
in a single block (##TITLE...##END), with a fixed number of values per line except for the last line before
the '##END' statement. Data may be written separated by specified delimiter, in PAC format or mixed

Under development:
 - Compound files with multiple data blocks
 - SQZ and DIF data formats


## Quick start

```julia
import Pkg 
Pkg.add("https://github.com/Manarom/JCAMPDXir.jl.git")
using JCAMPDXir
data = read_jdx_file(file_name) # to read the file, data.x - x-values, data.y - y values
write_jdx_file(x,y,"MKM","TRANSMITTANCE") # to write x- and y- data vectors
```

## Contact
To contact me [GitHub repository](https://github.com/Manarom).