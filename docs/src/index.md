
# JCAMPDXir.jl

## General description
This small package is designed to read and write the infrared spectra from/to files written JCAMP-DX format (common file extention is ".jdx").
JCAMP-DX (infrared) file format was developed for the exchange of infrared spectra between different laboratories.
General descripton of JCAMP-DX infrared format can be found in 
[JCAMP-DX for infrared 4.24](https://iupac.org/what-we-do/digital-standards/jcamp-dx/)
In addition to the spectra themselves, the file also stores metadata containing information about the units of measurement
and the conditions under which the spectra were acquired. A detailed specification of the format is provided via the link. 
This package was written for a specific task and cannot be considered as a complete implementation of all JCAMP-DX format specification; 
however, it implements some of the basic functionality to read  and write JCAMP-DX files
according to the `JCAMP-DX=4.24`.

## Contact

To contact me, please do it through the [GitHub repository](https://github.com/Manarom/JCAMPDXir.jl).
