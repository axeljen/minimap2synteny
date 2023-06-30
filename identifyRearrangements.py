from alnToSyntenyMap import *
import argparse
import datetime
import sys
from tabulate import tabulate
from copy import copy

def flip_negative_alignments(qgenome, alnDict_f):
    for chrom in qgenome.keys():
        strand = check_dominating_strand(chrom, alnDict_f)
        if strand == '-':
            #print("flipping chrom {}".format(chrom))
            alnDict_f = flip_aln_pos(alnDict_f, chrom, qgenome[chrom]['length'])
    return alnDict_f

def create_synteny_map(alnDict_f):
    syntenyMap = aln2syntenyMap(alnDict_f)
    return syntenyMap

def make_synteny_map(qgenome, stepsize, syntenyMap):
    synmap_q = []
    for chrom in qgenome.keys():
        for pos in range(1, qgenome[chrom]['length'], stepsize):
            tchrom, tpos = find_target_pos(chrom, pos, syntenyMap)
            if tpos is not None:
                tpos = int(tpos)
            synmap_q.append({'tchrom': chrom, 'tpos': pos, 'qchrom': tchrom, 'qpos': tpos})
    return synmap_q

def make_synteny_map_target(refgenome, stepsize, syntenyMap):
    synmap_t = []
    for chrom in refgenome.keys():
        for pos in range(1, refgenome[chrom]['length'], stepsize):
            qchrom, qpos = find_target_pos(chrom, pos, syntenyMap)
            if qpos is not None:
                qpos = int(qpos)
            synmap_t.append({'qchrom': qchrom, 'qpos': qpos, 'tchrom': chrom, 'tpos': pos})
    return synmap_t

def remove_none_target_positions(synmap_t):
    synmap_tf = [x for x in synmap_t if x['qpos'] is not None]
    return synmap_tf

# define a function to merge consecutive dicts into blocks if qchrom is the same, and absolute qpos - previous qpos is less tna 1mb away from tpos - previous tpos
def merge_consecutive_dicts(synmap, maxdiff=5000000):
	blocks = []
	# merge consequtive windows into blocks if the target chrom is the same, and the absolute target position - previous target position is less than 1mb away from query position - previous query position
	for i in range(len(synmap)):
		# if the blocks are empty, add this first one
		if len(blocks) == 0:
			blocks.append({'tchrom':synmap[i]['tchrom'],'tstart': synmap[i]['tpos'], 'tend': synmap[i]['tpos'], 'qchrom': synmap[i]['qchrom'], 'qstart': synmap[i]['qpos'], 'qend': synmap[i]['qpos']})
			continue
		# same if the target chrom is not the same as the previous one
		elif synmap[i]['tchrom'] != blocks[-1]['tchrom']:
			blocks.append({'tchrom':synmap[i]['tchrom'],'tstart': synmap[i]['tpos'], 'tend': synmap[i]['tpos'], 'qchrom': synmap[i]['qchrom'], 'qstart': synmap[i]['qpos'], 'qend': synmap[i]['qpos']})
			continue
		# if the target chrom is the same as the previous one and the target chrom is the same too
		if synmap[i]['qchrom'] == blocks[-1]['qchrom']:
			# get the distance between this target pos and the previous target pos
			tdist = abs(synmap[i]['tpos'] - blocks[-1]['tend'])
			# and the distance between this query pos and the previous query pos
			qdist = abs(synmap[i]['qpos'] - blocks[-1]['qend'])
			# if the difference between qdist and tdist is less than , add this position to the block
			if abs(tdist - qdist) < maxdiff:
				blocks[-1]['tend'] = synmap[i]['tpos']
				blocks[-1]['qend'] = synmap[i]['qpos']
			# if the difference between qdist and tdist is more than 1mb, start a new block
			else:
				blocks.append({'tchrom':synmap[i]['tchrom'],'tstart': synmap[i]['tpos'], 'tend': synmap[i]['tpos'], 'qchrom': synmap[i]['qchrom'], 'qstart': synmap[i]['qpos'], 'qend': synmap[i]['qpos']})
		# if qchrom has changed, add a new block
		else:
			blocks.append({'tchrom':synmap[i]['tchrom'],'tstart': synmap[i]['tpos'], 'tend': synmap[i]['tpos'], 'qchrom': synmap[i]['qchrom'], 'qstart': synmap[i]['qpos'], 'qend': synmap[i]['qpos']})
	# return the blocks
	return blocks

# define a function to sort blocks by target position
def sort_blocks_by_target_position(blocks):
	# sort first by target chrom, then by target start
	sorted_blocks = sorted(blocks, key=lambda k: (k['tchrom'], k['tstart']))
	# add keys for the order of the within each target and query chrom
	for chrom in set([x['tchrom'] for x in sorted_blocks]):
		for qchrom in set([x['qchrom'] for x in sorted_blocks if x['tchrom'] == chrom]):
			i = 1
			for block in sorted_blocks:
				if block['tchrom'] == chrom and block['qchrom'] == qchrom:
					block['order'] = i
					i += 1
	# return the sorted blocks
	return sorted_blocks

# function to sort by query position
def sort_blocks_by_query_position(blocks):
	# sort first by query chrom, then by query start
	sorted_blocks = sorted(blocks, key=lambda k: (k['qchrom'], k['qstart']))
	# add keys for the order of the within each target and query chrom
	for chrom in set([x['qchrom'] for x in sorted_blocks]):
		for tchrom in set([x['tchrom'] for x in sorted_blocks if x['qchrom'] == chrom]):
			i = 1
			for block in sorted_blocks:
				if block['qchrom'] == chrom and block['tchrom'] == tchrom:
					block['q-order'] = i
					i += 1
	# return the sorted blocks
	return sorted_blocks

# function to print all dictionary keys and values as a table, from list of dictionaries
def print_dict_as_table(list_of_dicts):
	# print header tabulated
	print(tabulate(list_of_dicts, headers="keys"))

# define a function to identify inversions
def identify_inversions(blocks):
	# if qend - qstart is negative, add a key specifying that this is an inversion
	for block in blocks:
		if block['qend'] - block['qstart'] < 0:
			block['direction'] = 'inversion'
		else:
			block['direction'] = 'forward'
	# keep a separate list with all inversions
	inversions = [x for x in blocks if x['direction'] == 'inversion']
	# return the blocks and inversions
	return blocks, inversions

# function to remove blocks shorter than length_limit
def remove_short_blocks(blocks, length_limit):
	blocks = [x for x in blocks if x['tend'] - x['tstart'] > length_limit and abs(x['qend'] - x['qstart']) > length_limit]
	return blocks


# define a function to add sort position for each tchrom and qchrom, based on starting postions of the respective genomes
def add_sort_position(blocks):
	# first sort the blocks by tchrom and then tstart
	blocks = sorted(blocks, key=lambda k: (k['tchrom'], k['tstart']))
	# add keys for the order of the within each target and query chrom
	for chrom in set([x['tchrom'] for x in blocks]):
		for qchrom in set([x['qchrom'] for x in blocks if x['tchrom'] == chrom]):
			i = 1
			for block in blocks:
				if block['tchrom'] == chrom and block['qchrom'] == qchrom:
					block['order'] = i
					i += 1
	# sort the blocks by tchrom, then qchrom and last qstart
	blocks = sorted(blocks, key=lambda k: (k['tchrom'], k['qchrom'], k['qstart']))
	# add a query order for each query chrom as they're currently sorted in the blocks list
	for chrom in set([x['qchrom'] for x in blocks]):
		for tchrom in set([x['tchrom'] for x in blocks if x['qchrom'] == chrom]):
			i = 1
			for block in blocks:
				if block['qchrom'] == chrom and block['tchrom'] == tchrom:
					block['q-order'] = i
					i += 1
	# sort the blocks by tchrom, then tchrom and last tstart, before returning
	blocks = sorted(blocks, key=lambda k: (k['tchrom'], k['qchrom'], k['qstart']))
	return blocks

# function to swap qstart and qend of inversions
def swap_qstart_qend_inversions(blocks):
	for block in blocks:
		if block['direction'] == 'inversion':
			qstart = block['qstart']
			block['qstart'] = block['qend']
			block['qend'] = qstart
	return blocks


# define a function to check if a target and query chrom is sorted
def check_if_sorted(blocks, tchrom, qchrom):
	# make a list of all blocks with this tchrom and qchrom
	blocks_considered = [x for x in blocks if x['tchrom'] == tchrom and x['qchrom'] == qchrom]
	# check if they are sorted
	if [x['order'] for x in blocks_considered] == sorted([x['order'] for x in blocks_considered]):
		return True
	else:
		return False


# define a function to sort blocks and identify translocations
def identify_translocations(blocks, tchrom, qchrom):
	translocations = []
	# first we check if it is sorted, in which case we just return the blocks and None
	if check_if_sorted(blocks, tchrom, qchrom):
		return blocks, translocations
	# otherwise, make one list of the target order, and one of the query order
	target_order = [x['order'] for x in blocks if x['tchrom'] == tchrom and x['qchrom'] == qchrom]
	query_order = [x['q-order'] for x in blocks if x['tchrom'] == tchrom and x['qchrom'] == qchrom]
	# find the lowest number in the target order that does not have the correct index
	n = len(target_order) + 1
	for y,i in enumerate(target_order):
		#print(target_order[y], y + 1)
		if i < n and target_order[y] != y + 1:
			#print(y,target_order[y], i - 1)
			target_order[y]
			n = i
			index = y
	# now, the index of n should in principle be n - 1, so to check how many steps we need to move, we subtract the index from the target order
	wanted_index = n - 1
	difference = wanted_index - index
	# fetch the block that we want to move
	block_to_move = [x for x in blocks if x['tchrom'] == tchrom and x['qchrom'] == qchrom and x['order'] == n][0]
	# and the blocks we want to swap places with
	blocks_to_swap = [x for x in blocks if x['tchrom'] == tchrom and x['qchrom'] == qchrom and x['order'] > n]
	# find the lowest and highest position in among blocks to swap, and append to translocations
	tchrom = block_to_move['tchrom']
	tstart = min([x['tstart'] for x in blocks_to_swap])
	tend = max([x['tend'] for x in blocks_to_swap])
	# append to translocations
	translocations.append([block_to_move,{'swapped_with': {'tchrom':tchrom, 'tstart':tstart, 'tend':tend}}])
	# add info this being a translocation in the blocks
	block_to_move['translocation'] = True
	# find the index of this block in the full list
	index = [i for i,x in enumerate(blocks) if x['tchrom'] == tchrom and x['qchrom'] == qchrom and x['order'] == n][0]
	# find the index of where to insert it
	insert_index = index + difference
	# remove the block from the list
	blocks.remove(block_to_move)
	# insert it at the correct index
	blocks.insert(insert_index, block_to_move)
	# return the blocks and translocations
	return blocks, translocations


# define a function to identify fission breaking points
def find_fissions(blocks):
	fissions = []
	# make sure that the blocks are sorted by tchrom and the tstart
	blocks = sorted(blocks, key=lambda k: (k['tchrom'], k['tstart']))
	for i in range(1,len(blocks)):
		# if the tchrom is the same, but the qchrom is different, we have a fission
		if blocks[i]['tchrom'] == blocks[i-1]['tchrom'] and blocks[i]['qchrom'] != blocks[i-1]['qchrom']:
			# find the breaking point
			breaking_point = (blocks[i]['previous_tpos'], blocks[i]['tstart'])
			# and their new chroms
			new_chroms = [blocks[i-1]['qchrom'],blocks[i-1]['qstart'],blocks[i-1]['qend'], blocks[i]['qchrom'],blocks[i]['qstart'],blocks[i]['qend']]
			# append to the fissions list
			fissions.append({'tchrom':blocks[i]['tchrom'], 'breaking_point':breaking_point, 'new_chroms':new_chroms})
	return fissions

# define a funciton to find the preceding block for a given target chromosome and position
def find_preceding_block(blocks, tchrom, tpos):
	# make a deepcopy of the blocks
	cblocks = copy.deepcopy(blocks)
	# sort copy by tchrom and tstart
	cblocks = sorted(cblocks, key=lambda k: (k['tchrom'], k['tstart']))
	# find the index of the block with this tchrom and tpos
	index = [i for i,x in enumerate(cblocks) if x['tchrom'] == tchrom and x['tstart'] == tpos][0]
	try:
		# return the block before this one
		# check if the chromosome of the preceeding block is the same
		if cblocks[index - 1]['tchrom'] == tchrom:
			# then the preceeding block pos should be the tend of this block
			tpos = cblocks[index - 1]['tend']
		# otherwise it is at the start of a new chromosome, so we set it to 0
		else:
			tpos = 0
	# if we get an index error, we're at the first block, so we set the tpos to 0
	except IndexError:
		tpos = 0
	# return the tpos
	return tpos

# define a function to find the following block for a given target chromosome and position
def find_following_block(blocks, tchrom, tpos, refgenome):
	# make a deepcopy of the blocks
	cblocks = copy.deepcopy(blocks)
	# sort copy by tchrom and tstart
	cblocks = sorted(cblocks, key=lambda k: (k['tchrom'], k['tstart']))
	# find the index of the block with this tchrom and tpos
	index = [i for i,x in enumerate(cblocks) if x['tchrom'] == tchrom and x['tstart'] == tpos][0]
	try:
		# return the block after this one
		# check if the chromosome of the following block is the same
		if cblocks[index + 1]['tchrom'] == tchrom:
			# then the following block pos should be the tstart of this block
			tpos = cblocks[index + 1]['tstart']
		# otherwise it is at the end of a chromosome, so we set it to the length of the chromosome
		else:
			tpos = refgenome[tchrom]['length']
	# if we get an index error, we're at the last block, so we set the tpos to the length of the chromosome
	except IndexError:
		tpos = refgenome[tchrom]['length']
	# return the tpos
	return tpos
	
# add previous and next target position for the same chromosome to each block
def add_previous_next_target_position(blocks, refgenome):
	for i,block in enumerate(blocks):
		# if it's the first block, we set the previous tpos to 0
		if i == 0:
			block['previous_tpos'] = 0
		# otherwise check if the previous block is on the same chromosome
		elif blocks[i-1]['tchrom'] == block['tchrom']:
			# if it is, we set the previous tpos to the tend of the previous block
			block['previous_tpos'] = blocks[i-1]['tend']
		# otherwise we set it to 0
		else:
			block['previous_tpos'] = 0
		# if it's the last block, we set the next tpos to the length of the chromosome
		if i == len(blocks) - 1:
			block['next_tpos'] = refgenome[block['tchrom']]['length']
		# otherwise check if the next block is on the same chromosome
		elif blocks[i+1]['tchrom'] == block['tchrom']:
			# if it is, we set the next tpos to the tstart of the next block
			block['next_tpos'] = blocks[i+1]['tstart']
		# otherwise we set it to the length of the chromosome
		else:
			block['next_tpos'] = refgenome[block['tchrom']]['length']
	# return the blocks
	return blocks


# define a function to find fusion events
def find_fusions(blocks, refgenome):
	fusions = []
	# sort the blocks by qchrom and then qstart
	blocks = sorted(blocks, key=lambda k: (k['qchrom'], k['qstart']))
	# loop through the blocks
	for i in range(1,len(blocks)):
		# if the qchrom is the same, but the tchrom is different, we have a fusion
		if blocks[i]['qchrom'] == blocks[i-1]['qchrom'] and blocks[i]['tchrom'] != blocks[i-1]['tchrom']:
			# first fusion point is the end of the previous block is direciton is forward, otherwise it is the start of the previous block
			if blocks[i-1]['direction'] == 'forward':
				fusion_point_1 = (blocks[i-1]['tchrom'],blocks[i-1]['tend'],blocks[i-1]['next_tpos'])
			else:
				fusion_point_1 = (blocks[i-1]['tchrom'],blocks[i-1]['previous_tpos'],blocks[i-1]['tstart'])
			# second fusion point is the start of the current block is direciton is forward, otherwise it is the end of the current block
			if blocks[i]['direction'] == 'forward':
				fusion_point_2 = (blocks[i]['tchrom'],blocks[i]['previous_tpos'],blocks[i]['tstart'])
			else:
				fusion_point_2 = (blocks[i]['tchrom'],blocks[i]['tend'],blocks[i]['next_tpos'])
			# find the cumulative position of the fusion points
			cumulative_fusion_point_1 = fusion_point_1[1] + refgenome[fusion_point_1[0]]['cumpos']
			cumulative_fusion_point_2 = fusion_point_2[1] + refgenome[fusion_point_2[0]]['cumpos']
			# if secont point is smaller then first, swap them
			if cumulative_fusion_point_2 < cumulative_fusion_point_1:
				fusion_point_1, fusion_point_2 = fusion_point_2, fusion_point_1
			#print(cumulative_fusion_point_1, cumulative_fusion_point_2)
			#print(fusion_point_1, fusion_point_2)
			#print(blocks[i-1], blocks[i])
			# find the fusion point
			fused_to = (blocks[i]['qchrom'], blocks[i-1]['qend'], blocks[i]['qstart'])
			# and their old chroms
			old_chroms = ("{}:{}".format(blocks[i-1]['tchrom'], blocks[i-1]['tend']), "{}:{}".format(blocks[i]['tchrom'], blocks[i]['tstart']))
			# append to the fusions list
			fusions.append({'qchrom':blocks[i]['qchrom'], 'fusion_point_1':fusion_point_1, 'fusion_point_2': fusion_point_2,
		    'fusion_cords_in_query':fused_to})
	# return the fusions
	return fusions

# define a function to write the genomes to a file, chrom and length and query/ref
def write_genomes(qgenome_sorted, qgenome, refgenome, outfile):
	outfile.write("chrom\tlength\ttype\n")
	for chrom in refgenome.keys():
		outfile.write("{}\t{}\t{}\n".format(chrom, refgenome[chrom]['length'], "ref"))
	for chrom in qgenome_sorted:
		outfile.write("{}\t{}\t{}\n".format(chrom, qgenome[chrom]['length'], "query"))

# define a function to merge consecutive blocks, if their of the same type, on the same chroms and no more than maxdist apart
def merge_consecutive_blocks(blocks, maxdist=5000000):
	blocks_merged = []
	# make sure blocks are sorted by tchrom and tstart
	blocks = sorted(blocks, key=lambda k: (k['tchrom'], k['tstart']))
	# add the first block to the list
	blocks_merged.append(blocks[0])
	# loop through the blocks
	for i in range(1,len(blocks)):
		merged = False
		# if tchrom and qchrom is the same 
		if blocks[i]['tchrom'] == blocks[i-1]['tchrom'] and blocks[i]['qchrom'] == blocks[i-1]['qchrom']:
			# if q-order is one up from the previous block, (or down) from the previous block
			if abs(blocks[i]['q-order'] == blocks_merged[-1]['q-order']) == 1:
				# if direction is the same as the previous block
				if blocks[i]['direction'] == blocks_merged[-1]['direction']:
					# if the maxdist is less than maxdist on both target and query
					if abs(blocks[i]['tstart'] - blocks_merged[-1]['tend']) < maxdist and abs(blocks[i]['qstart'] - blocks_merged[-1]['qend']) < maxdist:
						# merge the blocks
						blocks_merged[-1]['tend'] = blocks[i]['tend']
						blocks_merged[-1]['qend'] = blocks[i]['qend']
						blocks_merged[-1]['q-order'] = blocks[i]['q-order']
						blocks_merged[-1]['order'] = blocks[i]['order']
						merged = True
		# if the blocks were not merged, add the block to the list
		if not merged:
			blocks_merged.append(blocks[i])
	# return the merged blocks
	return blocks_merged

# define a function to parse the different dictionaries into a ready-ro-write table
def parse_blocks_to_table(event, etype):
	if etype == 'fission':
		refchrom_1 = event['tchrom']
		refstart_1 = event['breaking_point'][0]
		refend_1 = event['breaking_point'][1]
		refchrom_2 = "NA"
		refstart_2 = "NA"
		refend_2 = "NA"
		querychrom_1,querystart_1,queryend_1,querychrom_2,querystart_2,queryend_2 = event['new_chroms']
		note = "NA"
	elif etype == 'fusion':
		refchrom_1,refstart_1,refend_1 = event['fusion_point_1']
		refchrom_2,refstart_2,refend_2 = event['fusion_point_2']
		querychrom_1 = event['fusion_cords_in_query'][0]
		querystart_1 = event['fusion_cords_in_query'][1]
		queryend_1 = event['fusion_cords_in_query'][2]
		querychrom_2 = "NA"
		querystart_2 = "NA"
		queryend_2 = "NA"
		note = "NA"
	elif etype == 'translocation':
		t = event[0]
		note = event[1]['swapped_with']['tchrom'] + ":" + str(event[0]['tstart']) + "-" + str(event[0]['tend'])
		refchrom_1,refstart_1,refend_1,refchrom_2,refstart_2,refend_2 = t['tchrom'],t['previous_tpos'],t['tstart'],t['tchrom'],t['tend'],t['next_tpos']
		querychrom_1,querystart_1,queryend_1 = t['qchrom'],t['qstart'],t['qend']
		querychrom_2,querystart_2,queryend_2 = "NA","NA","NA"
	elif etype == 'inversion':
		refchrom_1,refstart_1,refend_1,refchrom_2,refstart_2,refend_2 = event['tchrom'],event['previous_tpos'],event['tstart'],event['tchrom'],event['tend'],event['next_tpos']
		querychrom_1,querystart_1,queryend_1 = event['qchrom'],event['qstart'],event['qend']
		querychrom_2,querystart_2,queryend_2 = "NA","NA","NA"
		note = "NA"
	# add all keys, including etype to a neat dictionary
	parsed_d = {
		'refchrom_1':refchrom_1,
		'refstart_1':refstart_1,
		'refend_1':refend_1,
		'refchrom_2':refchrom_2,
		'refstart_2':refstart_2,
		'refend_2':refend_2,
		'querychrom_1':querychrom_1,
		'querystart_1':querystart_1,
		'queryend_1':queryend_1,
		'querychrom_2':querychrom_2,
		'querystart_2':querystart_2,
		'queryend_2':queryend_2,
		'note':note,
		'event':etype
	}
	# return the dictionary
	return parsed_d









# prepare the argument parser
parser = argparse.ArgumentParser(description='Parse a PAF whole genome alignment file and make a synteny map, and identify large-scale genome rearrangements.')
# alignment file
parser.add_argument('-a', '--aln', type=str, required=True, help='PAF file with whole genome alignment')
# reference genome
parser.add_argument('-r', '--refgenome', type=str, required=True, help='Genome file for the ref genome, with at least two columns: chromosome in first and length in second. Remaining columns are ignored')
# query genome
parser.add_argument('-q', '--querygenome', type=str, required=True, help='Genome file for the query genome, with at least two columns: chromosome in first and length in second. Remaining columns are ignored')
# minimum alignment length
parser.add_argument('-l', '--minlen', type=int, default=10000, help='Minimum alignment length to consider')
# minimum mapping quality
parser.add_argument('-m', '--minmapq', type=int, default=40, help='Minimum mapping quality to consider')
# minimum percent identity
parser.add_argument('-p', '--minpid', type=float, default=0.9, help='Minimum percent identity to consider')
# output prefix
parser.add_argument('-o', '--outprefix', type=str, default="out", help='Prefix for output files')
# min difference query and reference lenghts before breaking a syntenyblock
parser.add_argument('-d', '--mindiff', type=int, default=5000000, help='Minimum difference between query and reference lengths before breaking a syntenyblock')
# min length for a synteny block
parser.add_argument('-b', '--minblock', type=int, default=1000000, help='Minimum length for a synteny block')
# syntenymap step size
parser.add_argument('-s', '--stepsize', type=int, default=100000, help='Step size for the syntenymap')
# parse the arguments
args = parser.parse_args()

# # input alignment
#aln="../syntenyplot/map_Allochrocebus_lhoesti_hap1_asm_ref.paf"
# aln2="../syntenyplot/map_Cercopithecus_diana_hap1_asm_ref.paf"
# # input query genome
#qgenome="../syntenyplot/Allochrocebus_lhoesti_hap1_repeatinput.fa.fai"
# qgenome2="../syntenyplot/Cercopithecus_diana_hap1_repeatinput.fa.fai"
# # reference genome
#rgenome="../syntenyplot/mmul_chroms.fai"
stepsize = args.stepsize
length_limit = args.minlen
aln = args.aln
rgenome = args.refgenome
qgenome = args.querygenome
stepsize = args.stepsize
length_limit = args.minblock

# write to sys.stderr that we've started and how it was run
sys.stderr.write("Started at {} with the following command:\n".format(datetime.datetime.now()))
sys.stderr.write(" ".join(sys.argv) + "\n\n")

# print that we're parsing and filtering the PAF file
sys.stderr.write("Parsing and filtering the PAF file\n")
# parse paf
alnDict = parsePAF(aln)

refgenome = read_ref_index(rgenome)
qgenome = read_ref_index(qgenome)
# add cumulative position to the refgenome
refgenome = add_cumulative_starting_position(refgenome, plotlen=None, spacing=0)

# filter paf
alnDict_f = filterPAF(alnDict, minlen=10000, primary=True, minmapq=40, minpid=0.9, only_ref_chroms=True, only_query_chroms=False, ref_chroms=refgenome.keys(), query_chroms=qgenome.keys())

# print total length of alignmnet after filtring
sys.stderr.write("Total length of alignment after filtering: {}\n".format(sum([int(x['qend']) - int(x['qstart']) for x in alnDict_f])))


# flip negative alignments
alnDict_f = flip_negative_alignments(qgenome, alnDict_f)
# get synteny map
syntenyMap = aln2syntenyMap(alnDict_f)
# print that we're making syntenymaps
sys.stderr.write("Making syntenymaps\n")
# make a syntenymap based on the query genome
#synmap_q = make_synteny_map(qgenome, stepsize, syntenyMap)
# make a syntenymap based on the target genome
synmap_t = make_synteny_map(refgenome, stepsize, syntenyMap)
# remove nones
synmap_tf = remove_none_target_positions(synmap_t)

# merge consecutive dicts into blocks
blocks = merge_consecutive_dicts(synmap_tf, maxdiff=5000000)

# remove blocks shorter than length_limit
blocks = remove_short_blocks(blocks, length_limit)

# add previous and next target position for the same chromosome to each block
blocks = add_previous_next_target_position(blocks, refgenome)

# identify inversions
blocks, inversions = identify_inversions(blocks)

# swap qstart and qend for inversions
blocks = swap_qstart_qend_inversions(blocks)

# add sort order positions of each block
blocks = add_sort_position(blocks)
#print_dict_as_table(blocks)
	
translocations = []
# fetch translocations from the blocks
for chrom in refgenome.keys():
	for qchrom in set([x['qchrom'] for x in blocks if x['tchrom'] == chrom]):
		is_sorted = False
		if check_if_sorted(blocks, chrom, qchrom):
			continue
		else:
			while not is_sorted:
				blocks,t = identify_translocations(blocks, chrom, qchrom)
				translocations.extend(t)
				if check_if_sorted(blocks, chrom, qchrom):
					is_sorted = True

# fetch fissions
fissions = find_fissions(blocks)

# fetch fusions,
fusions = find_fusions(blocks, refgenome)

# merge consecutive blocks
blocks = merge_consecutive_blocks(blocks, maxdist=5000000)

# before writing, we need to sort the query genome by ascending target pos
# get cumulative starting pos by summing the length of all previous chroms
refgenome = add_cumulative_starting_position(refgenome, plotlen=None, spacing=0)

# find the lowest cordinates for each chromosome in the query genome, and order them accordingly in the genome dictionary
for chrom in qgenome.keys():
	# get median tpos for this chromosome
	qgenome[chrom]['median'] = find_mean_median_target_pos(chrom, alnDict_f, refgenome)[1]
# sort the chromsoomes by lowest median position
qgenome_sorted = sorted(qgenome.keys(), key=lambda k: qgenome[k]['median'])



# write that we're writing the files
sys.stderr.write("Writing the output files\n")

# write genomes to a file
with open(args.outprefix + "_genomes.txt", 'w') as outfile:
	write_genomes(qgenome_sorted, qgenome, refgenome, outfile)

# write the target-based syntenymap to a file
with open(args.outprefix + "_syntenymap.txt", 'w') as outfile:
	outfile.write("tchrom\ttpos\tqchrom\tqstart\n")
	for item in synmap_t:
		outfile.write("{}\t{}\t{}\t{}\n".format(item['tchrom'], item['tpos'], item['qchrom'], item['qpos']))


# make a more writeable dict from all different rearrangements
combined_dict = []
# parse all events to the combined dict
for fission in fissions:
	combined_dict.append(parse_blocks_to_table(fission, 'fission'))
for fusion in fusions:
	combined_dict.append(parse_blocks_to_table(fusion, 'fusion'))
for inversion in inversions:
	combined_dict.append(parse_blocks_to_table(inversion, 'inversion'))
for translocation in translocations:
	combined_dict.append(parse_blocks_to_table(translocation, 'translocation'))
# now write this to  the file
with open(args.outprefix + "_rearrangements.txt", 'w') as outfile:
	# start with writing the dict keys as headers
	outfile.write("{}\n".format("\t".join(combined_dict[0].keys())))
	# then write each dict as a line
	for item in combined_dict:
		outfile.write("{}\n".format("\t".join([str(x) for x in item.values()])))
# and then we're done
# print that we're done, and give a summary on how many of each rearrengement type we found
sys.stderr.write("Done at {} with the following summary:\n".format(datetime.datetime.now()))
sys.stderr.write("Found {} fissions\n".format(len(fissions)))
sys.stderr.write("Found {} fusions\n".format(len(fusions)))
sys.stderr.write("Found {} inversions\n".format(len(inversions)))
sys.stderr.write("Found {} translocations\n".format(len(translocations)))
