#!/usr/bin/env python3
"""Run a command and record its wall time, CPU time and peak memory.

Usage: timeit.py OUT.txt CMD [ARGS...]
A stand-in for /usr/bin/time on machines that lack it. Exits with the command's exit code.
"""
import resource
import subprocess
import sys
import time

out, cmd = sys.argv[1], sys.argv[2:]
start = time.time()
code = subprocess.call(cmd)
wall = time.time() - start
usage = resource.getrusage(resource.RUSAGE_CHILDREN)
with open(out, "w") as f:
    f.write("wall_s %.1f\n" % wall)
    f.write("user_s %.1f\n" % usage.ru_utime)
    f.write("sys_s %.1f\n" % usage.ru_stime)
    f.write("peak_rss_gb %.2f\n" % (usage.ru_maxrss / 1048576.0))
    f.write("exit_code %d\n" % code)
sys.exit(code)
