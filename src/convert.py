#!/usr/bin/python
# -*- coding: utf8 -*-
# Usage: ./convert.py fileName 0,1,2,...
import sys
import csv
from sas7bdat import SAS7BDAT
cols = sys.argv[2].split(',')
kWriter = csv.writer(sys.stdout)
with SAS7BDAT(sys.argv[1]) as f:
 for row in f:
  kData = []
  for i in cols:
   kData.append(row[int(i)])
  kWriter.writerow(kData)
