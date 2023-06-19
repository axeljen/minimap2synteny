""" this script will take alignments against a common reference genome and prepare a 
synteny map in a progressive manner
"""
# import modules
import argparse
import os
import sys
import copy as cp


# define a function that will read the reference genome index file and return a dictionary
def read_ref_index(ref_index):
	# chromosome dictionary
	chrom_dict = {}
	# open the reference genome index file
	with open(ref_index, 'r') as ref:
		# loop over the file
		for line in ref:
			# split the line
			line = line.strip().split()
			# add the chromosome name and length to the dictionary
			chrom_dict[line[0]] = {'length': int(line[1])}
	# return the dictionary
	return chrom_dict

# function to add cumulative starting position of the chromosomes in a genome, adding spacing between chromosomes based on a total length to fit the genome plot
def add_cumulative_starting_position(chrom_dict, plotlen=None, spacing=20000000):
	# start with calculating a spacing between chromosomes
	# total length of the genome
	totlen = sum([i['length'] for i in chrom_dict.values()])
	# number of chromosomes
	numchrom = len(chrom_dict.keys())
	# number of spaces between chromosomes
	numspace = numchrom - 1
	# if plotlen is none, than use spacing of 5000000 and calculate plotlen
	if plotlen == None:
		plotlen = totlen + numspace * spacing
	# calculate the spacing between chromosomes
	spacing = int((plotlen - totlen) / numspace)
	# add the cumulative starting position to the dictionary
	# start with a cumulative starting position of 1
	cumpos = 1
	# cumulative starting position of the first key in the dictionary should be 1
	chrom_dict[list(chrom_dict.keys())[0]]['cumpos'] = 1
	# add the length of this chromosome plus spacing as current cumulative starting position
	cumpos += chrom_dict[list(chrom_dict.keys())[0]]['length'] + spacing
	# loop over the rest of the chromosomes
	for chrom in list(chrom_dict.keys())[1:]:
		# add the cumulative starting position to the dictionary
		chrom_dict[chrom]['cumpos'] = cumpos
		# add the length of this chromosome plus spacing as current cumulative starting position
		cumpos += chrom_dict[chrom]['length'] + spacing
	# return the dictionary
	return chrom_dict

# function to parse the paf alignments to dictionary
def parsePAF(paf):
	# these are all mandatory cols in the paf
	pafcols = ['qseqname','qseqlen','qstart','qend','strand','tseq','tlen','tstart','tend','nbasematch','nbasetot','mapq']
	# and then there are a bunch of potential sam headers, which we'll store as a separate list
	paf_addcols = ['tp','cm','s1','s2','NM','MD','AS','SA','ms','nn','ts','cg','cs','dv','de','rl']
	# list of dictionaries for storing the alignemnts as key, value pairs
	alignments = []
	# loop through all the alignments in the input file and parse to dict
	with open(paf) as f:
		for line in f.readlines():
			alignments.append({})
			# split lines into items
			values = line.strip().split()
			# first 12 cols will always be there
			for i in range(0,len(pafcols)):
				alignments[-1][pafcols[i]] = values[i]
			# rest of the cols, if any, should be split into key,value pairs first
			if len(values) > len(pafcols):
				for v in values[len(pafcols):len(values)]:
					key = v.split(":")[0]
					try:
						vallist = v.split(":")[1:]
						val = ''.join(vallist)
					except:
						#val = str(v.split(":")[1:])
						print("failed to parse {}, {}".format(key, v.split(":")[1:]))
					# check if key is in addcols, otherwise error
					#if not key in paf_addcols:
						#print("Unrecognized tag: {}".format(key))
						# sys.exit(1)
					alignments[-1][key] = val
	return alignments

# define a function to filter paf alignment dictionary based on a minimum alignment length, primary alignment, minimum mapping quality and minimum percent identity
def filterPAF(paf_dict, minlen=100000, primary=True, minmapq=40, minpid=0.9, only_ref_chroms=True, only_query_chroms=True, ref_chroms=None, query_chroms=None):
	# if only ret_chroms is True, then rec_chroms is mandatory, and if only_query_chroms is True, then query_chroms is mandatory, if these are not provided, print an error message and exit
	if only_ref_chroms == True and ref_chroms == None:
		print('ERROR: only_ref_chroms is True but no reference chromosomes are provided')
		sys.exit()
	if only_query_chroms == True and query_chroms == None:
		print('ERROR: only_query_chroms is True but no query chromosomes are provided')
		sys.exit()
	# define a list to store the filtered alignments
	filtered_paf = []
	# make a variable defining if filter is passed
	filter_passed = True
	# loop over the alignments
	for aln in paf_dict:
		# if subject chrom is not in the list of reference chromosomes, skip this alignment
		if only_ref_chroms == True and aln['tseq'] not in ref_chroms:
			continue
		# if query chrom is not in the list of query chromosomes, skip this alignment
		if only_query_chroms == True and aln['qseqname'] not in query_chroms:
			continue
		# if the alignment is not primary, skip this alignment
		if primary == True and aln['tp'] != 'AP':
			continue
		# if the alignment length is shorter than the minimum length, skip this alignment
		if int(aln['qend']) - int(aln['qstart']) < minlen:
			continue
		# if the mapping quality is lower than the minimum mapping quality, skip this alignment
		if int(aln['mapq']) < minmapq:
			continue
		# if the percent identity is lower than the minimum percent identity, skip this alignment
		if int(aln['nbasematch']) / int(aln['nbasetot']) < minpid:
			continue
		# add a key specifying the original reference chromosome, that will be used for plotting later on
		aln['tseq_orig'] = aln['tseq']
		# if the alignment passes all filters, add it to the filtered alignments
		filtered_paf.append(aln)
	# return the filtered alignments
	return filtered_paf

# define a function to return the summed length of all alignments in an alignment dictionary
def sum_aln_len(aln_dict):
	# start with a length of 0
	length = 0
	# loop over the alignments
	for aln in aln_dict:
		# add the length of the alignment to the total length
		length += int(aln['qend']) - int(aln['qstart'])
	# return the total length
	return length

# define a function to find the lowest cumulative starting position of an alignment per chromosome in query genome
def find_lowest_starting_pos(chrom, refgenome, aln_dict):
	# start with an insanely high number that will not be suroassed
	lowest_starting_pos = 100000000000
	# loop over the alignments
	for aln in aln_dict:
		# fetch the raw target starting position
		raw_starting_pos = int(aln['tstart'])
		# translate this to cumulative starting position in the reference genome
		cum_starting_pos = refgenome[aln['tseq']]['cumpos'] + raw_starting_pos
		# if the alignment is on the chromosome of interest
		if aln['qseqname'] == chrom:
			# if the cumulative starting position is lower than the current lowest starting position
			if cum_starting_pos < lowest_starting_pos:
				# update the lowest starting position
				lowest_starting_pos = cum_starting_pos
	# return the lowest starting position
	return lowest_starting_pos

# define a function to sort chromosomes based on lowest mean, median or 10th percentile target position
def sort_chromosomes(aln_dict, sort_by='median'):
	# sort the keys in the dictionary by ascending mean_target_pos
	sorted_chromosomes = sorted(aln_dict.keys(), key=lambda x: aln_dict[x][sort_by])
	# return the sorted chromosomes
	return sorted_chromosomes

# define a function to reorder keys (chromosomes) in a dictionary based on mean_target_pos, median_target_pos or 10th percentile target position
def reorder_chromosomes(aln_dict, sort_by='median'):
	# sort the keys in the dictionary by ascending mean_target_pos
	sorted_chromosomes = sorted(aln_dict.keys(), key=lambda x: aln_dict[x][sort_by])
	# make a new dictionary to store the reordered chromosomes
	reordered_dict = {}
	# loop over the sorted chromosomes
	for i,chrom in enumerate(sorted_chromosomes):
		# add the chromosome to the new dictionary
		reordered_dict[chrom] = aln_dict[chrom]
		# add the new order to the dictionary
		reordered_dict[chrom]['order'] = i + 1
	# return the reordered dictionary
	return reordered_dict

# sort the lowest_starting_pos_dict based on the lowest starting position
# def sort_chromosomes(lowest_starting_pos_dict):
# 	# sort the dictionary based on the lowest starting position
# 	sorted_chromosomes = sorted(lowest_starting_pos_dict.items(), key=lambda x: x[1])
# 	# return the sorted chromosomes
# 	return sorted_chromosomes

# define a function to find the position for a given query position in the reference genome
def find_ref_pos(query_pos, query_chrom, refgenome, aln_dict):
	# define a variable to see if this position is found within an alignment
	found = False
	# loop over the alignments
	for aln in aln_dict:
		# if the alignment is on the query chromosome of interest
		if aln['qseqname'] == query_chrom:
			# if the query position is within the alignment
			if int(aln['qstart']) <= query_pos <= int(aln['qend']):
				# if tseq is None
				if aln['tseq'] == None:
					# return None
					return None, None, None
				# fetch the raw target starting position
				raw_starting_pos = int(aln['tstart'])
				# translate this to cumulative starting position in the reference genome
				cum_starting_pos = refgenome[aln['tseq']]['cumpos'] + raw_starting_pos
				# calculate the position in the reference genome
				ref_pos = cum_starting_pos + (query_pos - int(aln['qstart']))
				ref_chrom = aln['tseq']
				found = True
				# return the position in the reference genome
				return ref_pos, ref_chrom, aln['tseq_orig']
	# if the position is not found within an alignment, return None
	if found == False:
		return None, None, None

# define a function to iterate through windows in the alignment and find the reference position for each query position
def find_ref_pos_window(query_chrom, query_chrom_length, refgenome, aln_dict, window_size=100000):
	# store results in a list of dictionaries, keys being query chrom, query pos, ref chrom and ref cum pos
	results = []
	# loop over the query chromosome
	for i in range(1, query_chrom_length, window_size):
		# find the reference position for this query position
		ref_pos, ref_chrom, original_refchrom = find_ref_pos(i, query_chrom, refgenome, aln_dict)
		# add the results to the list of dictionaries
		results.append({'query_chrom':query_chrom, 'query_pos':i, 'ref_chrom':ref_chrom, 'ref_pos':ref_pos, 'original_refchrom':original_refchrom})
	# return the results
	return results

# define a function to split the alignments by query chromosome
def split_aln_by_query_chrom(aln_dict):
	# make a dictionary to store the alignments, keys being query chroms
	aln_dict_split = {}
	# loop over the alignments
	for aln in aln_dict:
		# if the query chrom is not in the dictionary, add it
		if aln['qseqname'] not in aln_dict_split.keys():
			aln_dict_split[aln['qseqname']] = []
		# add the alignment to the dictionary
		aln_dict_split[aln['qseqname']].append(aln)
	# return the dictionary
	return aln_dict_split

# define a fuction to split the alignments by reference chromosome
def split_aln_by_ref_chrom(aln_dict):
	# make a dictionary to store the alignments, keys being ref chroms
	aln_dict_split = {}
	# loop over the alignments
	for aln in aln_dict:
		# if the ref chrom is not in the dictionary, add it
		if aln['tseq'] not in aln_dict_split.keys():
			aln_dict_split[aln['tseq']] = []
		# add the alignment to the dictionary
		aln_dict_split[aln['tseq']].append(aln)
	# return the dictionary
	return aln_dict_split

# define a function to write the chromdicts, with four columns: chrom, cumpos, length and genome as integer, iterating over a list of chromdicts as input
def write_chromdicts(chromdicts, outfile):
	# open the output file
	with open(outfile, 'w') as out:
		# write the header
		out.write('chrom\tcumpos\tlength\tgenome\n')
		# loop over the chromdicts
		for i,chromdict in enumerate(chromdicts):
			# loop over the chromosomes
			for chrom in chromdict.keys():
				# write the chrom, cumpos, length and genome as integer to the output file
				out.write('{}\t{}\t{}\t{}\n'.format(chrom, chromdict[chrom]['cumpos'], chromdict[chrom]['length'], str(i + 1)))

# define a function to write the syntenymap (results),iterating over a list of resultlists and write the following columns: query_chrom, query_pos, ref_pos, ref_chrom, query_genome as integer (enumerate +1) and target_genome as integer (query_genome - 1)
def write_syntenyMap(resultlists, outfile):
	# open the output file
	with open(outfile, 'w') as out:
		# write the header
		out.write('query_chrom\tquery_pos\tref_pos\tref_chrom\trefgenome_chrom\tquery_genome\ttarget_genome\n')
		# loop over the resultlists
		for i,resultlist in enumerate(resultlists):
			# loop over the results
			for result in resultlist:
				# write the query_chrom, query_pos, ref_pos, query_genome and target_genome to the output file
				out.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(result['query_chrom'], result['query_pos'], result['ref_pos'], result['ref_chrom'], result['original_refchrom'], str(i + 2), str(i + 1)))

# define a function to convert from chromosome and position to cumulative position
def chrompos2cumpos(chrom, pos, chrom_dict):
	# fetch the cumulative starting position for this chromosome
	cumpos = chrom_dict[chrom]['cumpos']
	# add the position to the cumulative starting position
	cumpos += pos
	# return the cumulative position
	return cumpos

# define a function to find the mean and median target positions for a given query chromosomes alignments
def find_mean_median_target_pos(query_chrom, aln_dict, ref_dict):
	# make a list to store the target positions
	target_pos = []
	# loop over the alignments
	for aln in aln_dict:
		# if the query chrom is the query chrom of interest
		if aln['qseqname'] == query_chrom:
			# if tseq is None, skip this alignment
			if aln['tseq'] == None:
				continue
			# convert the target position to cumulative position
			target_pos_cum = chrompos2cumpos(aln['tseq'], int(aln['tstart']), ref_dict)
			# add the target position to the list of target positions
			target_pos.append(int(target_pos_cum))
	# calculate the mean target position
	mean_target_pos = sum(target_pos) / len(target_pos)
	# calculate the median target position
	median_target_pos = sorted(target_pos)[len(target_pos) // 2]
	# calculate the 10th percentile target position
	percentile_10_target_pos = sorted(target_pos)[int(len(target_pos) * 0.1)]
	# return the mean and median target positions and the 10th percentile target position
	return mean_target_pos, median_target_pos, percentile_10_target_pos

# define a function that takes target chrom, start and end positions and returns the query chrom, start and end positions from a list of alignments
def find_query_pos(target_chrom, target_start, target_end, aln_dict_split_by_chrom):
	# loop over the alignments
	for aln in aln_dict_split_by_chrom[target_chrom]:
		# if the start position of the alignment is within the target start and end positions
		if int(aln['tstart']) <= target_start <= int(aln['tend']):
			# check if also the end position of the alignment is within the target start and end positions
			if int(aln['tstart']) <= target_end <= int(aln['tend']):
				# if so, convert the target start and end positions to query start and end positions
				query_start = int(aln['qstart']) + (target_start - int(aln['tstart']))
				query_end = int(aln['qstart']) + (target_end - int(aln['tstart']))
				# return the query chrom, start and end positions
				return aln['qseqname'], query_start, query_end
			# if start but not end position of the alignment is within the target start and end positions
			else:
				# return target end as end of alignment
				query_start = int(aln['qstart']) + (target_start - int(aln['tstart']))
				query_end = int(aln['qend'])
				# return the query chrom, start and end positions
				return aln['qseqname'], query_start, query_end
		# if the end position of the alignment is within the target start and end positions
		elif int(aln['tstart']) <= target_end <= int(aln['tend']):
			# return target start as start of alignment
			query_start = int(aln['qstart'])
			query_end = int(aln['qstart']) + (target_end - int(aln['tstart']))
			# return the query chrom, start and end positions
			return aln['qseqname'], query_start, query_end
	# if no alignment is found, return None
	return None, None, None

# define function to check of majority of alignments are on plus or minus strand
def check_dominating_strand(chrom,aln_dict):
	# define a variable to store the strand
	strand = None
	# define a variable to store the number of plus strand alignments
	plus = 0
	# define a variable to store the number of minus strand alignments
	minus = 0
	# loop over the alignments
	for aln in aln_dict:
		# if the alignment is on the query chromosome of interest
		if aln['qseqname'] == chrom:
			# if the strand is plus
			if aln['strand'] == '+':
				# add 1 to the plus strand variable
				plus += 1
			# if the strand is minus
			elif aln['strand'] == '-':
				# add 1 to the minus strand variable
				minus += 1
	# if the number of plus strand alignments is higher than the number of minus strand alignments
	if plus > minus:
		# set the strand to plus
		strand = '+'
	# if the number of minus strand alignments is higher than the number of plus strand alignments
	elif minus > plus:
		# set the strand to minus
		strand = '-'
	# return the strand
	return strand

# define a function to flip alignment positions from minus to plus strand
def flip_aln_pos(aln_dict, query_chrom, query_chrom_length):
	# loop over the alignments
	for aln in aln_dict:
		# if the query chrom is the query chrom of interest
		if aln['qseqname'] == query_chrom:
			# fetch the target chrom length
			#target_chrom_length = refgenome[aln['tseq']]['length']
			start = int(aln['qstart'])
			end = int(aln['qend'])
			# if the strand is minus
			if aln['strand'] == '-':
				# flip the query start and end positions
				aln['qstart'] = query_chrom_length - end
				aln['qend'] = query_chrom_length - start
			elif aln['strand'] == '+':
				# flip this to minus
				aln['qstart'] = query_chrom_length - end
				aln['qend'] = query_chrom_length - start
	# return the flipped alignments
	return aln_dict

# prepare an argument parser
parser = argparse.ArgumentParser(description='''this script will take alignments against a common reference genome and prepare a
synteny map in a progressive manner''')
# add arguments to parse
# fasta index file of reference genome, using --reference or -r, required
parser.add_argument('--reference', '-r', type=str, required=True, help='the fasta index file of reference genome')
# option to add as many alignment files as wanted, using --alignment or -aln, at least one is required
parser.add_argument('--alignment', '-aln', type=str, required=True, nargs='+', help='alignment in paf format to the reference genome')
# each alignment file should have an index file, using --index or -i, required
parser.add_argument('--index', '-i', type=str, required=True, nargs='+', help='index file of the alignment file')
# output file name, using --output or -o, default is "syntenyMap.txt"
parser.add_argument('--output', '-o', type=str, default='genomes', help='output file prefix, will be used to create two files: <output>.chromdicts.txt and <output>.syntenyMap.txt')
# add argument to supply spacing to plot, defaulting to 20mb
parser.add_argument('--spacing', type=int, default=20000000, help='spacing between reference genome chromosomes in the plot, default is 20mb')
# add argument for minimum alignment length to keep in filter step
parser.add_argument('--minlen', type=int, default=10000, help='minimum alignment length to keep in filter step, default is 10kb')
# minimum mapping quality to keep in filter step
parser.add_argument('--minmapq', type=int, default=40, help='minimum mapping quality to keep in filter step, default is 40')
# minimum percent identity to keep in filter step
parser.add_argument('--minpid', type=float, default=0.9, help='minimum percent identity to keep in filter step, default is 0.9')
# add argument with block size with which to create the syntenymap, default is 100kb
parser.add_argument('--block_size', type=int, default=100000, help='block size with which to create the syntenymap, default is 100kb')
# parse arguments
args = parser.parse_args()

# define a function that tests that all arguments are valid
def test_args(args):
	# check that the reference genome index file exists
	if not os.path.isfile(args.reference):
		# if not, print an error message
		print('ERROR: the reference genome index file does not exist')
		# and exit
		sys.exit()
	# check that the alignment files exist
	for aln in args.alignment:
		if not os.path.isfile(aln):
			# if not, print an error message
			print('ERROR: the alignment file does not exist')
			# and exit
			sys.exit()
	# check that the index files exist
	for index in args.index:
		if not os.path.isfile(index):
			# if not, print an error message
			print('ERROR: the index file does not exist')
			# and exit
			sys.exit()
	# check that the number of alignment files is the same as the number of index files
	if len(args.alignment) != len(args.index):
		# if not, print an error message
		print('ERROR: the number of alignment files is not the same as the number of index files')
		# and exit
		sys.exit()
	# check that the output file name is valid
	if '/' in args.output:
		# if not, print an error message
		print('ERROR: the output file name is not valid')
		# and exit
		sys.exit()

# test the arguments
test_args(args)

# print current date and time to standard error, a few newlines on both sites for readability
sys.stderr.write('\n\n')
sys.stderr.write('Script started at:\n')
os.system('date')
sys.stderr.write('\n\n')

# print commands as executed to standard error
sys.stderr.write('Script executed as:\n')
sys.stderr.write(' '.join(sys.argv) + '\n')


# read the reference genome index file and return a dictionary
# print('reading the reference genome file') to standard error
sys.stderr.write('reading the reference genome file from\n{}\n'.format(args.reference))
refgenome = read_ref_index(args.reference)
# testing
#refgenome = read_ref_index("mmul_chroms.fai")
# add cumulative starting position of the chromosomes in the reference genome, adding spacing between chromosomes based on a total length to fit the genome plot
# print('adding cumulative starting position of the chromosomes in the reference genome') to standard error
refgenome = add_cumulative_starting_position(refgenome, plotlen=None, spacing=args.spacing)
# testing
#refgenome = add_cumulative_starting_position(refgenome, plotlen=None, spacing=20000000)


# make a list to store the alignment dictionaries
aln_dicts = []
# and a list to store the index/genome files
genomes = []
# append parsed index files to the list
for index in args.index:
	genomes.append(read_ref_index(index))

# testing
#for index in ['Cercopithecus_diana_hap1_repeatinput.fa.fai','Allochrocebus_lhoesti_hap1_repeatinput.fa.fai']:
#	genomes.append(read_ref_index(index))
aln_dicts_raw = []
# loop over the alignment files
for i,aln in enumerate(args.alignment):
# testing
#for i,aln in enumerate(['map_Cercopithecus_diana_hap1_asm_ref.paf','map_Allochrocebus_lhoesti_hap1_asm_ref.paf']):
	# print('reading the alignment file') to standard error
	sys.stderr.write('reading the alignment file from\n{}\n'.format(aln))
	# parse the paf alignments to dictionary
	aln_dict = parsePAF(aln)
	# print('filtering the alignments') to standard error
	sys.stderr.write('filtering the alignments\n')
	# filter paf alignment dictionary based on a minimum alignment length, primary alignment, minimum mapping quality and minimum percent identity
	aln_dict = filterPAF(aln_dict, minlen=args.minlen, primary=True, minmapq=args.minmapq, minpid=args.minpid, only_ref_chroms=True, only_query_chroms=False, ref_chroms=refgenome.keys(), query_chroms=genomes[i].keys())
	# testing without args. in argument
	#aln_dict = filterPAF(aln_dict, minlen=10000, primary=True, minmapq=40, minpid=0.9, only_ref_chroms=True, only_query_chroms=False, ref_chroms=refgenome.keys(), query_chroms=genomes[i].keys())
	# add the alignment dictionary to the list of alignment dictionaries
	aln_dicts.append(aln_dict)
	# that will be our main dict that we'll used for all the filtration step, but also append them to a raw dictlist where we'll keep the original alignments
	aln_dicts_raw.append(cp.deepcopy(aln_dict))

# print how long alignments remain after filtering for each alignment
for i,aln_dict in enumerate(aln_dicts):
	alnlen = sum_aln_len(aln_dict)
	sys.stderr.write('after filtering, {} bp alignment remain in alignment file {}\n'.format(alnlen, i + 1))


# the first dict in the dictlist will be directly projected onto the reference genome, so start with processing this one
# print('splitting the alignments by query chromosome') to standard error
sys.stderr.write('splitting the alignments by query chromosome\n')
# split the alignments by query chromosome
aln_dict_split_by_chrom = split_aln_by_query_chrom(aln_dicts[0])
# find the median target position for each query chromosome
# print('finding the median target position for each query chromosome') to standard error
sys.stderr.write('finding the median target position for each query chromosome\n')
# add the median target position to the dictionary
for chrom in aln_dict_split_by_chrom.keys():
# 	# check if the majority of alignments are on the plus or minus strand
# 	strand = check_dominating_strand(chrom, aln_dicts[0])
# 	# if strand is minus, flip the alignments
# 	if strand == '-':
# 		print("flipping chrom {}".format(chrom))
# 		aln_dicts[i] = flip_aln_pos(aln_dicts[0], chrom, genomes[0][chrom]['length'])
	genomes[0][chrom]['median'] = find_mean_median_target_pos(chrom, aln_dict_split_by_chrom[chrom], refgenome)[1]
# sort the chromosomes based on median target position
# print('sorting the chromosomes based on median target position') to standard error
sys.stderr.write('sorting the chromosomes based on median target position\n')
# sort the chromosomes based on median target position
sorted_chromosomes = sort_chromosomes(genomes[0], sort_by='median')
# reorder the chromosomes based on median target position
# print('reordering the chromosomes based on median target position') to standard error
sys.stderr.write('reordering the chromosomes based on median target position\n')
# reorder the chromosomes based on median target position
genomes[0] = reorder_chromosomes(genomes[0], sort_by='median')

# get the plot length by taking the highest cumulative position in the reference genome + the length of the last chromosome
plotlen = refgenome[list(refgenome.keys())[-1]]['cumpos'] + refgenome[list(refgenome.keys())[-1]]['length']
# add the cumulative starting position in relation to reference genome to the alignments in this first alignment dictionary
# print('adding the cumulative starting position in relation to reference genome to the alignments in this first alignment dictionary') to standard error
sys.stderr.write('adding the cumulative starting position in relation to reference genome to the alignments in this first alignment dictionary\n')
# add the cumulative starting position to all the chroms in the first genome
genomes[0] = add_cumulative_starting_position(genomes[0], plotlen=plotlen, spacing=None)

# prep a list of list to store the results
results = [[]]
# get syntenymap for the first query genome
# print('getting syntenymap for the first query genome') to standard error
sys.stderr.write('making the syntenymap for the first query genome\n')
# get syntenymap for the first query genome
for chrom in genomes[0].keys():
	# get the results for this query chromosome
	results[0] += find_ref_pos_window(chrom, genomes[0][chrom]['length'], refgenome, aln_dict_split_by_chrom[chrom], window_size=args.block_size)
	# testing without args. in argument
	#results[0] += find_ref_pos_window(chrom, genomes[0][chrom]['length'], refgenome, aln_dict_split_by_chrom[chrom], window_size=100000)
# loop over the rest of the alignment dictionaries
for i in range(1,len(aln_dicts)):
	# check if we need to flip any chromosomes
	# for chrom in genomes[i]:
	# 	# check if the majority of alignments are on the plus or minus strand
	# 	strand = check_dominating_strand(chrom, aln_dicts[i])
	# 	# if strand is minus, flip the alignments
	# 	if strand == '-':
	# 		print("flipping chrom {}".format(chrom))
	# 		aln_dicts[i] = flip_aln_pos(aln_dicts[i], chrom, genomes[i][chrom]['length'])
	# convert target coordinates to the query genome of the previous alignment dictionary
	# print('converting target coordinates to the query genome of the previous alignment dictionary') to standard error
	sys.stderr.write('converting target coordinates of genome {} to fit genome {}\n'.format(args.alignment[i], args.alignment[i - 1]))
	# test print the same line wihtout args. in argument
	#sys.stderr.write('converting target coordinates of genome {} to fit genome {}\n'.format(genomes[i], genomes[i -1]))
	sys.stderr.write('this may take a while...\n')
	# split the alignments by reference chromosome
	target_dict_split_by_ref_chrom = split_aln_by_ref_chrom(aln_dicts_raw[i - 1])
	for aln in aln_dicts[i]:
		# convert the target start and end positions to query start and end positions
		aln['tseq'], aln['tstart'], aln['tend'] = find_query_pos(aln['tseq'], int(aln['tstart']), int(aln['tend']), target_dict_split_by_ref_chrom)


# sort the remaining genomes based on median target position
# print('sorting the remaining genomes based on median target position') to standard error
sys.stderr.write('progressively sorting the genomes based on median target position\n')
# sort the remaining genomes based on median target position
for i in range(1,len(aln_dicts)):
	# split the alignments by query chromosome
	aln_dict_split_by_chrom = split_aln_by_query_chrom(aln_dicts[i])
	# add the median target position to the dictionary
	for chrom in aln_dict_split_by_chrom.keys():
		#print(chrom, genomes[i - 1])
		genomes[i][chrom]['median'] = find_mean_median_target_pos(chrom, aln_dict_split_by_chrom[chrom], genomes[i - 1])[1]
	# sort the chromosomes based on median target position
	sorted_chromosomes = sort_chromosomes(genomes[i], sort_by='median')
	# reorder the chromosomes based on median target position
	genomes[i] = reorder_chromosomes(genomes[i], sort_by='median')
	# add cumulative starting position to this genome
	genomes[i] = add_cumulative_starting_position(genomes[i], plotlen=plotlen, spacing=None)

# make a syntenymap for the remaining genomes
# print('making a syntenymap for the remaining genomes') to standard error
sys.stderr.write('making a syntenymap for the remaining genomes\n')
# make a syntenymap for the remaining genomes
for i in range(1,len(aln_dicts)):
	# append an empty list to the results
	results.append([])
	# split by chrom again
	aln_dict_split_by_chrom = split_aln_by_query_chrom(aln_dicts[i])
	# loop over chromosomes
	for chrom in genomes[i].keys():
		# get the results for this query chromosome
		results[i] += find_ref_pos_window(chrom, genomes[i][chrom]['length'], genomes[i-1], aln_dict_split_by_chrom[chrom], window_size=args.block_size)
		# testing without args. in argument
		#results[i] += find_ref_pos_window(chrom, genomes[i][chrom]['length'], genomes[i-1], aln_dict_split_by_chrom[chrom], window_size=100000)

# convert the results query pos to cumulative position
for i,resultlist in enumerate(results):
	for result in resultlist:
		result['query_pos'] = chrompos2cumpos(result['query_chrom'], result['query_pos'], genomes[i])

# write the chromdicts to file
# print('writing the chromdicts to file') to standard error
sys.stderr.write('writing the chromdicts to file\n')
# add refgenome to first position in list of genomes
genomes.insert(0, refgenome)
# write the chromdicts to file
write_chromdicts(genomes, args.output + '.chroms.txt')
# test write without args. in argument
#write_chromdicts(genomes, 'genomes.chroms.txt')
# write the syntenymap to file
# print('writing the syntenymap to file') to standard error
sys.stderr.write('writing the syntenymap to file\n')
# write the syntenymap to file
write_syntenyMap(results, args.output + '.syntenyMap.txt')
# test wrute without args. in argument
#write_syntenyMap(results, 'genomes.syntenyMap.txt')

# print that we're done and how long it took
sys.stderr.write('\n\n')
sys.stderr.write('Script finished at:\n')
os.system('date')
# total time
sys.stderr.write('\n\n')
