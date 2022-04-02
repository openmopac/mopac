# Input file processing

This directory contains functionality for the processing and parsing of input files. Note that
the input file is copied to a scratch file by `getdat.F90` before any other processing is performed.
The primary parsing of keywords occurs in `wrtkey.F90`, and any new keywords should be introduced
there before adding any other keyword-based control logic associated with that keyword.
