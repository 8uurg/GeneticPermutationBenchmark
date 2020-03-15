import re

file = "./web20150522.txt"

nmrgx = re.compile(r"^ *(\d+) *x *(\d+)")
tarrgx = re.compile(r"^ta(\d+)-(\d+)")

print("instance,n,m,i,lb,ub")
with open(file, 'r') as f:
    bounds = []
    for line in f:
        if line.startswith('ta'):
            tarry = tarrgx.match(line)
            start = int(tarry.group(1))
            bounds = [(start + i, bound.strip().split('-')) for i, bound in enumerate(line[9:].split('|'))]
        else:
            isnandm = nmrgx.match(line)
            if isnandm is not None:
                n = int(isnandm.group(1))
                m = int(isnandm.group(2))

                filename = f"tai{n}_{m}.txt"
                for num, bound in bounds:
                    if len(bound) > 1:
                        csvbound = f"{bound[0].strip()},{bound[1].strip()}"
                    else:
                        csvbound = f"{bound[0].strip()},{bound[0].strip()}"
                    i = (num - 1)%10+1
                    print(f"{filename}#{i},{n},{m},{i},{csvbound}")
