#!/bin/bash
~/programs/analysis_test *.scfout
sname=$(grep -oP "System\.Name\s*\K.*" temporal_12345.input)
~/programs/bader *.tden.cube
