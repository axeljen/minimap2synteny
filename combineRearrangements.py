import argparse

# define a function to check the amount of overlap between two intervals
def overlap(a, b):
	return max(0, min(a[1], b[1]) - max(a[0], b[0]))

# function to parse the rearrangement file
def parse_rearrangements(file):
	rearrangements = []
	with open(file, 'r') as f:
		for i,line in enumerate(f):
			if i == 0:
				# parse to headers/keys
				headers = line.strip().split('\t')
			else:
				# parse to values
				values = line.strip().split('\t')
				# make a dictionary
				rearrangement = dict(zip(headers, values))
				# add to list
				rearrangements.append(rearrangement)
	return rearrangements

# function to sort rearrangements based on ref chrom and type of event
def sort_rearrangements(rearrangements, refgenome):
	# sorted rearrangements
	sorted_rearrangements = []
	# loop over the chromosomes in the reference genome
	for chrom in refgenome:
		# loop over the rearrangements
		for rearrangement in rearrangements:
			# if the rearrangement is on this chromosome
			if rearrangement['refchrom_1'] == chrom:
				# add to the sorted rearrangements
				sorted_rearrangements.append(rearrangement)
	# then sort by refchrom_1, event and last by refstart_1
	sorted_rearrangements = sorted(sorted_rearrangements, key=lambda k: (k['refchrom_1'], k['event'], int(k['refstart_1'])))
	return sorted_rearrangements


# prepare arguments
parser = argparse.ArgumentParser(description='Combine rearrangements from different samples')
parser.add_argument('-i', '--input', help='Input file with rearrangements', required=True, nargs='+')
parser.add_argument('-o', '--output', help='Output file with combined rearrangements', required=True)
parser.add_argument('-s', '--samples', help='Samples labels, if not given the input file names will be used.', nargs='+')
# reference genome file
parser.add_argument('-r', '--refgenome', help='Reference genome file, used to sort rearrangements', required=True)

# parse arguments
args = parser.parse_args()

# if samples are given, make sure the number of samples matches the number of input files
if args.samples:
	labels = args.samples
	if len(labels) != len(args.input):
		raise ValueError('The number of samples does not match the number of input files.')
else:
	# if no samples are given, use the input file names
	labels = [x.split('/')[-1].split('.')[0] for x in args.input]

# make a dictionary with the rearrangements per sample
rearrangements = {}
for i, file in enumerate(args.input):
	rearrangements[labels[i]] = parse_rearrangements(file)

# make a list to store the combined rearrangements
combined_rearrangements = []

# loop over the rearrangements of the first sample
for rearrangement in rearrangements[labels[0]]:
	# add all rearrangements from the first sample to the combined rearrangements
	combined_rearrangements.append({
		'refchrom_1': rearrangement['refchrom_1'],
		'refchrom_2': rearrangement['refchrom_2'],
		'refstart_1': rearrangement['refstart_1'],
		'refstart_2': rearrangement['refstart_2'],
		'refend_1': rearrangement['refend_1'],
		'refend_2': rearrangement['refend_2'],
		'samples': [labels[0]],
		'event': rearrangement['event'],
	})
# loop over the rearrangements of the other samples
for label in labels[1:]:
	for rearrangement in rearrangements[label]:
		# check if this rearrangement is the same as one of the rearrangements in the combined rearrangements
		# if so, add the sample to the samples list
		# if not, add the rearrangement to the combined rearrangements
		for cr in combined_rearrangements:
			# if its a fission, we only care about the first interval
			if rearrangement['event'] == 'fission':
				# if chrom is the same
				if rearrangement['refchrom_1'] == cr['refchrom_1']:
					# if the start and end overlap
					if overlap([int(rearrangement['refstart_1']), int(rearrangement['refend_1'])], [int(cr['refstart_1']), int(cr['refend_1'])]) > 0:
						# add the sample to the samples list
						cr['samples'].append(label)
						# stop looking for matches
						break
			# if its a fusion, inversion or rearrangement we care about both intervals
			else:
				# if chroms are the same
				if rearrangement['refchrom_1'] == cr['refchrom_1'] and rearrangement['refchrom_2'] == cr['refchrom_2']:
					# if the start and end overlap
					if overlap([int(rearrangement['refstart_1']), int(rearrangement['refend_1'])], [int(cr['refstart_1']), int(cr['refend_1'])]) > 0 and overlap([int(rearrangement['refstart_2']), int(rearrangement['refend_2'])], [int(cr['refstart_2']), int(cr['refend_2'])]) > 0:
						# add the sample to the samples list
						cr['samples'].append(label)
						# stop looking for matches
						break
		else:
			# if no match was found, add the rearrangement to the combined rearrangements
			combined_rearrangements.append({
				'refchrom_1': rearrangement['refchrom_1'],
				'refchrom_2': rearrangement['refchrom_2'],
				'refstart_1': rearrangement['refstart_1'],
				'refstart_2': rearrangement['refstart_2'],
				'refend_1': rearrangement['refend_1'],
				'refend_2': rearrangement['refend_2'],
				'samples': [label],
				'event': rearrangement['event'],
			})

# write the combined rearrangements to the output file
with open(args.output, 'w') as f:
	# write the header, first event, then intervals and then a column for each sample
	f.write('event\trefchrom_1\trefstart_1\trefend_1\trefchrom_2\trefstart_2\trefend_2\t' + '\t'.join(labels) + '\n')
	# loop over the rearrangements
	for rearrangement in combined_rearrangements:
		# write the event
		f.write(rearrangement['event'] + '\t')
		# write the first interval
		f.write(rearrangement['refchrom_1'] + '\t' + rearrangement['refstart_1'] + '\t' + rearrangement['refend_1'] + '\t')
		# write second interval, as None if its a fission
		if rearrangement['event'] == 'fission':
			f.write('NA\tNA\tNA\t')
		else:
			f.write(rearrangement['refchrom_2'] + '\t' + rearrangement['refstart_2'] + '\t' + rearrangement['refend_2'] + '\t')
		# write the samples, as 0 or 1 if the sample does or does not have the rearrangement
		for label in labels:
			if label in rearrangement['samples']:
				f.write('1\t')
			else:
				f.write('0\t')
		# write a newline
		f.write('\n')

		






