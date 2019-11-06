from bs4 import BeautifulSoup
from urllib.request import urlopen
import re

# Extractor
regex = r"([A-Za-z][A-Za-z0-9]*) +[0-9]+ +([0-9]+) +\([^\)]*\)(?: +([0-9]+))?"
def parse_table(table, result):
    matches = re.finditer(regex, table)
    for match in matches:
        name = match.group(1)
        ub   = match.group(2)
        lb   = ub
        if match.groups == 3:
            lb = match.group(3)
        result.append((name, ub, lb))
    return result

url = "http://anjos.mgi.polymtl.ca/qaplib//inst.html"

b = BeautifulSoup(urlopen(url), features="html.parser")

# Elements 0-2 are the objective function and format example.
tables = b("pre")[3:]
result = []
for table in tables:
    parse_table(table.text, result)


print(f"name\tub\tlb")
for r in result:
    print(f"{r[0]}\t{r[1]}\t{r[1]}")