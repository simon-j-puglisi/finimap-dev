#!/usr/bin/env python3
"""Copy standard input to standard output, prefixing each line with the seconds elapsed."""
import sys
import time

start = time.time()
for line in sys.stdin:
    if len(line) > 200:
        line = line[:200] + " ...\n"
    sys.stdout.write("%8.1f  %s" % (time.time() - start, line))
    sys.stdout.flush()
