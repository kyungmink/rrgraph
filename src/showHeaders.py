#!/usr/bin/python
# -*- coding: utf8 -*-
# Usage: ./showheaders.py fileName
import sys
from sas7bdat import SAS7BDAT
with SAS7BDAT(sys.argv[1]) as f:
    kIter = iter(f)
    print kIter.next()
    print kIter.next()
