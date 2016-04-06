#!/usr/bin/env python
import os
import sys
import spectrum
import parse_asensetek
import export_csv

if __name__ in "__main__":
    infile = sys.argv[1]
    outfile = sys.argv[2]

    spec = parse_asensetek.parse(infile)
    export_csv.export(spec, outfile)

