import re
import glob
import os

rgx_str = "number of jobs, number of machines, initial seed, upper bound and lower bound :\n[^\d]*(?P<n>\d+)[^\d]*(?P<m>\d+)[^\d]*(?P<seed>\d+)[^\d]*(?P<UB>\d+)[^\d]*(?P<LB>\d+)"
rgx = re.compile(rgx_str)

files = glob.glob("./instances/tai*.txt")

print("instance,n,m,UB,LB")
for file in files:
    filename = os.path.basename(file)
    with open(file, 'r') as f:
        for i, match in enumerate(rgx.finditer(f.read())):
            print(f"{filename}#{i+1},{match.group('n')},{match.group('m')},{match.group('UB')},{match.group('LB')}")