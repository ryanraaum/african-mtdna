import csv

starts = {}
ends = {}

with open("gomes2015.csv", "rU") as f:
	reader = csv.reader(f)
	reader.next() # skip past header
	for row in reader:
		group = row[3]
		if row[9] == '':
			hvr1 = row[7]
			start, end = hvr1.split('-')
			start = int(start)
			end = int(end)
			if starts.has_key(group):
				starts[group].append(start)
				ends[group].append(end)
			else:
				starts[group] = [start]
				ends[group] = [end]

for group in starts.keys():
	print group, "%d-%d" % (max(starts[group]), min(ends[group]))