def range2region(s):
	if ';' in s:
		parts = s.split(';')
		return range2region(parts[0]) + range2region(parts[1])
	parts = s.split('-')
	parts = [int(x) for x in parts]
	return range(parts[0], parts[1]+1)
