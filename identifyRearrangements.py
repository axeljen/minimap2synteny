from alnToSyntenyMap import *
import argparse
import datetime
import sys

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
# aln="../syntenyplot/map_Allochrocebus_lhoesti_hap1_asm_ref.paf"
# aln2="../syntenyplot/map_Cercopithecus_diana_hap1_asm_ref.paf"
# # input query genome
# qgenome="../syntenyplot/Allochrocebus_lhoesti_hap1_repeatinput.fa.fai"
# qgenome2="../syntenyplot/Cercopithecus_diana_hap1_repeatinput.fa.fai"
# # reference genome
# rgenome="../syntenyplot/mmul_chroms.fai"

# write to sys.stderr that we've started and how it was run
sys.stderr.write("Started at {} with the following command:\n".format(datetime.datetime.now()))
sys.stderr.write(" ".join(sys.argv) + "\n\n")

# print that we're parsing and filtering the PAF file
sys.stderr.write("Parsing and filtering the PAF file\n")
# parse paf
alnDict = parsePAF(args.aln)

refgenome = read_ref_index(args.refgenome)
qgenome = read_ref_index(args.querygenome)

# filter paf
alnDict_f = filterPAF(alnDict, minlen=10000, primary=True, minmapq=40, minpid=0.9, only_ref_chroms=True, only_query_chroms=False, ref_chroms=refgenome.keys(), query_chroms=qgenome.keys())

# flip negative alignments
for chrom in qgenome.keys():
	strand = check_dominating_strand(chrom, alnDict_f)
	# if strand is minus, flip the alignments
	if strand == '-':
		print("flipping chrom {}".format(chrom))
		alnDict_f = flip_aln_pos(alnDict_f, chrom, qgenome[chrom]['length'])

# get synteny map
syntenyMap = aln2syntenyMap(alnDict_f)

# print that we're making syntenymaps
sys.stderr.write("Making syntenymaps\n")
# make a syntenymap based on the query genome
synmap_q = []

# loop over the chromosomes in the query genome
for chrom in qgenome.keys():
    # find the syntenic position every 100kb
	for pos in range(1, qgenome[chrom]['length'], args.stepsize):
		tchrom,tpos = find_target_pos(chrom, pos, syntenyMap)
		if not tpos == None:
			tpos = int(tpos)
		# append to the synteny map
		synmap_q.append({'qchrom':chrom, 'qpos':pos, 'tchrom':tchrom, 'tpos':tpos})
# make a syntenymap based on the target genome
synmap_t = []
# loop over the chromosomes in the target genome
for chrom in refgenome.keys():
	# find the syntenic position every 100kb
	for pos in range(1, refgenome[chrom]['length'], args.stepsize):
		qchrom,qpos = find_target_pos(chrom, pos, syntenyMap)
		if not qpos == None:
			qpos = int(qpos)
		# append to the synteny map
		synmap_t.append({'qchrom':qchrom, 'qpos':qpos, 'tchrom':chrom, 'tpos':pos})

# length limit to consider a block
length_limit = args.minblock

# remove dict items where the target position is None
synmap_tf = [x for x in synmap_t if x['qpos'] != None]

# merge consexutive dicts into blocks if qchrom is the same, and absolute qpos - previous qpos is less tna 1mb away from tpos - previous tpos
blocks = []
#blocks.append({'tchrom':None,'tstart': 0, 'tend': 0, 'qchrom': None, 'qstart': None, 'qend': None, 'direction': None})

# write that we're calling rearrangements
sys.stderr.write("Calling rearrangements\n")

# merge consequtive windows into blocks if the target chrom is the same, and the absolute target position - previous target position is less than 1mb away from query position - previous query position
for i in range(len(synmap_tf)):
	# if the blocks are empty, add this first one
	if len(blocks) == 0:
		blocks.append({'tchrom':synmap_tf[i]['tchrom'],'tstart': synmap_tf[i]['tpos'], 'tend': synmap_tf[i]['tpos'], 'qchrom': synmap_tf[i]['qchrom'], 'qstart': synmap_tf[i]['qpos'], 'qend': synmap_tf[i]['qpos']})
		continue
	# same if the target chrom is not the same as the previous one
	elif synmap_tf[i]['tchrom'] != blocks[-1]['tchrom']:
		blocks.append({'tchrom':synmap_tf[i]['tchrom'],'tstart': synmap_tf[i]['tpos'], 'tend': synmap_tf[i]['tpos'], 'qchrom': synmap_tf[i]['qchrom'], 'qstart': synmap_tf[i]['qpos'], 'qend': synmap_tf[i]['qpos']})
		continue
	# if the target chrom is the same as the previous one and the target chrom is the same too
	if synmap_tf[i]['qchrom'] == blocks[-1]['qchrom']:
		# get the distance between this target pos and the previous target pos
		tdist = abs(synmap_tf[i]['tpos'] - blocks[-1]['tend'])
		# and the distance between this query pos and the previous query pos
		qdist = abs(synmap_tf[i]['qpos'] - blocks[-1]['qend'])
		# if the difference between qdist and tdist is less than 1mb, add this position to the block
		if abs(tdist - qdist) < 5000000:
			blocks[-1]['tend'] = synmap_tf[i]['tpos']
			blocks[-1]['qend'] = synmap_tf[i]['qpos']
		# if the difference between qdist and tdist is more than 1mb, start a new block
		else:
			blocks.append({'tchrom':synmap_tf[i]['tchrom'],'tstart': synmap_tf[i]['tpos'], 'tend': synmap_tf[i]['tpos'], 'qchrom': synmap_tf[i]['qchrom'], 'qstart': synmap_tf[i]['qpos'], 'qend': synmap_tf[i]['qpos']})
	# if qchrom has changed, add a new block
	else:
		blocks.append({'tchrom':synmap_tf[i]['tchrom'],'tstart': synmap_tf[i]['tpos'], 'tend': synmap_tf[i]['tpos'], 'qchrom': synmap_tf[i]['qchrom'], 'qstart': synmap_tf[i]['qpos'], 'qend': synmap_tf[i]['qpos']})
	#else:
	#	blocks.append({'tchrom':synmap_tf[i]['tchrom'],'tstart': synmap_tf[i]['tpos'], 'tend': synmap_tf[i]['tpos'], 'qchrom': synmap_tf[i]['qchrom'], 'qstart': synmap_tf[i]['qpos'], 'qend': synmap_tf[i]['qpos']})

# if qend - qstart is negative, add a key specifying that this is an inversion
for block in blocks:
	if block['qend'] - block['qstart'] < 0:
		block['direction'] = 'inversion'
	else:
		block['direction'] = 'forward'

# remove blocks shorter than length_limit
blocks = [x for x in blocks if x['tend'] - x['tstart'] > args.minblock]

# keep a separate list with all inversions
inversions = [x for x in blocks if x['direction'] == 'inversion']

# print all the blocks as a table
#for block in blocks:
#	print("%s\t%s\t%s\t%s\t%s\t%s\t%s" % (block['tchrom'], block['tstart'], block['tend'], block['qchrom'], block['qstart'], block['qend'], block['direction']))

# then switch the qstart and qend for all inversions in the original blocks list
for block in blocks:
	if block['direction'] == 'inversion':
		qstart = block['qstart']
		block['qstart'] = block['qend']
		block['qend'] = qstart

# now look for translocations within chromosomes
# if the tchrom is the same, but the qstart is not consecutive with the previous qend, we have a translocation
# start with checking the first item, and we'll always put the translocation on the block with the lowest tstart, for consistency
#if blocks[0]['tchrom'] == blocks[1]['tchrom'] and blocks[0]['qend'] < blocks[1]['qstart']:
#	blocks[0]['order'] = 'translocation'
#else:
#	blocks[0]['order'] = 'linear'
for i in range(0,len(blocks) - 1):
	# if tchrom and qchrom is the same as the next block, but qstart of the next block is smaller than qend of the current, call it a translocation
	if blocks[i]['tchrom'] == blocks[i+1]['tchrom'] and blocks[i]['qchrom'] == blocks[i+1]['qchrom'] and blocks[i]['qend'] > blocks[i+1]['qstart']:
		blocks[i]['order'] = 'translocation'
	else:
		blocks[i]['order'] = 'linear'
# check if the last block does not have an order
if not 'order' in blocks[-1].keys():
	blocks[-1]['order'] = 'linear'


# put translocation in their own list
try:
	translocations = [x for x in blocks if x['order'] == 'translocation']
except:
	print(blocks)
	sys.exit()

# now merge syntenic blocks such that one block for each synteny between tchrom and qchrom is kept
merged_blocks = []

# list with blocks under consideration
blocks_considered = []
# add first block to this
blocks_considered.append(blocks[0])

# loop through the blocks
for i in range(1,len(blocks)):
	# if tchrom and qchrom is the same as the latest block in blocks_considered, add this block to blocks_considered
	if blocks[i]['tchrom'] == blocks_considered[-1]['tchrom'] and blocks[i]['qchrom'] == blocks_considered[-1]['qchrom']:
		blocks_considered.append(blocks[i])
	else:
		# otherwise we try to merge the blocks in blocks_considered, if more than one
		if len(blocks_considered) > 0:
			# find lowest tstart and highest tend
			tstart = min([x['tstart'] for x in blocks_considered])
			tend = max([x['tend'] for x in blocks_considered])
			# find lowest qstart and highest qend
			qstart = min([x['qstart'] for x in blocks_considered])
			qend = max([x['qend'] for x in blocks_considered])
			# now merge these into a new block
			merged_blocks.append({'tchrom':blocks_considered[0]['tchrom'], 'tstart':tstart, 'tend':tend, 'qchrom':blocks_considered[0]['qchrom'], 'qstart':qstart, 'qend':qend})
		else:
			# otherwise, just add the block to merged_blocks
			merged_blocks.append(blocks_considered[0])
		# and reset blocks_considered
		blocks_considered = [blocks[i]]

# add the last block to merged_blocks
if len(blocks_considered) > 0:
	# find lowest tstart and highest tend
	tstart = min([x['tstart'] for x in blocks_considered])
	tend = max([x['tend'] for x in blocks_considered])
	# find lowest qstart and highest qend
	qstart = min([x['qstart'] for x in blocks_considered])
	qend = max([x['qend'] for x in blocks_considered])
	# now merge these into a new block
	merged_blocks.append({'tchrom':blocks_considered[0]['tchrom'], 'tstart':tstart, 'tend':tend, 'qchrom':blocks_considered[0]['qchrom'], 'qstart':qstart, 'qend':qend})


# remove translocations and inversions from the blocks list
blocks_f = [x for x in blocks if x['order'] == 'linear' and x['direction'] == 'forward']

# now we can get fissions and their breaking points
fissions = []
for i in range(1,len(merged_blocks)):
	# if the tchrom is the same, but the qchrom is different, we have a fission
	if merged_blocks[i]['tchrom'] == merged_blocks[i-1]['tchrom'] and merged_blocks[i]['qchrom'] != merged_blocks[i-1]['qchrom']:
		# find the breaking point
		breaking_point = (merged_blocks[i-1]['tend'], merged_blocks[i]['tstart'])
		# and their new chroms
		new_chroms = (merged_blocks[i-1]['qchrom'], merged_blocks[i]['qchrom'])
		# append to the fissions list
		fissions.append({'tchrom':merged_blocks[i]['tchrom'], 'breaking_point':breaking_point, 'new_chroms':new_chroms})

# and fusions
fusions = []
# first, we need to sort the merged blocks by first query chrom and then query start
merged_blocks = sorted(merged_blocks, key=lambda k: (k['qchrom'], k['qstart']))
for i in range(1,len(merged_blocks)):
	# if the tchrom is different, but the qchrom is the same, we have a fusion
	if merged_blocks[i]['tchrom'] != merged_blocks[i-1]['tchrom'] and merged_blocks[i]['qchrom'] == merged_blocks[i-1]['qchrom']:
		# find the fusion point in the query genome
		fusion_point = (merged_blocks[i-1]['qend'], merged_blocks[i]['qstart'])
		# and their old chroms
		old_chroms = ("{}:{}".format(merged_blocks[i-1]['tchrom'], merged_blocks[i-1]['tend']), "{}:{}".format(merged_blocks[i]['tchrom'], merged_blocks[i]['tstart']))
		# append to the fusions list
		fusions.append({'qchrom':merged_blocks[i]['qchrom'], 'fusion_point':fusion_point, 'old_chroms':old_chroms})

# add cumulative start position to each ref chrom
#for chrom in refgenome.keys():
# get cumulative starting pos by summing the length of all previous chroms
refgenome = add_cumulative_starting_position(refgenome, plotlen=None, spacing=0)
# print each chrom and its cumulative starting pos
##	print("{}\t{}".format(chrom, refgenome[chrom]['cumpos']))

# find the lowest cordinates for each chromosome in the query genome, and order them accordingly in the genome dictionary
for chrom in qgenome.keys():
	# get median tpos for this chromosome
	qgenome[chrom]['median'] = find_mean_median_target_pos(chrom, alnDict_f, refgenome)[1]
# sort the chromsoomes by lowest median position
qgenome_sorted = sorted(qgenome.keys(), key=lambda k: qgenome[k]['median'])

# print each chrom and it's median position
#for chrom in qgenome_sorted:
#	print("{}\t{}".format(chrom, qgenome[chrom]['median']))

#print(qgenome_sorted)
# write that we're writing the files
sys.stderr.write("Writing the output files\n")

# write the genomes to a file, chrom and length and query/ref
with open(args.outprefix + "_genomes.txt", 'w') as outfile:
	outfile.write("chrom\tlength\ttype\n")
	for chrom in qgenome_sorted:
		outfile.write("{}\t{}\t{}\n".format(chrom, qgenome[chrom]['length'], "query"))
	for chrom in refgenome.keys():
		outfile.write("{}\t{}\t{}\n".format(chrom, refgenome[chrom]['length'], "ref"))

# write the target-based syntenymap to a file
with open(args.outprefix + "_syntenymap.txt", 'w') as outfile:
	outfile.write("tchrom\ttpos\tqchrom\tqstart\n")
	for item in synmap_t:
		outfile.write("{}\t{}\t{}\t{}\n".format(item['tchrom'], item['tpos'], item['qchrom'], item['qpos']))

# and last write the rearrangements to a file
with open(args.outprefix + "_rearrangements.txt", 'w') as outfile:
	# write the header
	outfile.write("type\trefchrom\trefstart\trefend\tquerychrom\tquerystart\tqueryend\tnote\n")
	# then write the fissions
	for fission in fissions:
		outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("fission", fission['tchrom'], fission['breaking_point'][0], fission['breaking_point'][1], 'Na', 'Na', 'Na', "Split_into:" + fission['new_chroms'][0] + "&" + fission['new_chroms'][1]))
	# then the fusions
	for i,fusion in enumerate(fusions):
		# we write the fusions on two lines each, numbering them to store the fusion event
		fusion_id = "fusion_{}".format(i)
		chrom_startpos = fusion['old_chroms'][0].split(":")[0]
		startpos = fusion['old_chroms'][0].split(":")[1]
		chrom_endpos = fusion['old_chroms'][1].split(":")[0]
		endpos = fusion['old_chroms'][1].split(":")[1]
		outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("fusion", chrom_startpos, startpos, 'Na', fusion['qchrom'], fusion['fusion_point'][0], fusion['fusion_point'][1], fusion_id + "_from:" + fusion['old_chroms'][0] + "&" + fusion['old_chroms'][1]))
		outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("fusion", chrom_endpos, endpos, 'Na', fusion['qchrom'], fusion['fusion_point'][0], fusion['fusion_point'][1], fusion_id + "_from:" + fusion['old_chroms'][0] + "&" + fusion['old_chroms'][1]))
	# then the inversions
	for inversion in inversions:
		outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("inversion", inversion['tchrom'], inversion['tstart'], inversion['tend'], inversion['qchrom'], inversion['qstart'], inversion['qend'], "Na"))
	# then the translocations
	for translocation in translocations:
		outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("translocation", translocation['tchrom'], translocation['tstart'], translocation['tend'], translocation['qchrom'], translocation['qstart'], translocation['qend'], "Na"))

# print that we're done, and give a summary on how many of each rearrengement type we found
sys.stderr.write("Done at {} with the following summary:\n".format(datetime.datetime.now()))
sys.stderr.write("Found {} fissions\n".format(len(fissions)))
sys.stderr.write("Found {} fusions\n".format(len(fusions)))
sys.stderr.write("Found {} inversions\n".format(len(inversions)))
sys.stderr.write("Found {} translocations\n".format(len(translocations)))


# # now let's write the different types of rearrangements to a file
# # start with opening the outfile and write the header
# outfile = open("rearrangements.txt", 'w')
# outfile.write("type\trefchrom\trefstart\trefend\tquerychrom\tquerystart\tqueryend\tnote\n")
# # then write the fissions
# for fission in fissions:
# 	outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("fission", fission['tchrom'], fission['breaking_point'][0], fission['breaking_point'][1], 'Na', 'Na', 'Na', "Split_into:" + fission['new_chroms'][0] + "&" + fission['new_chroms'][1]))
# # then the fusions
# for fusion in fusions:
# 	outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("fusion", 'Na', 'Na', 'Na', fusion['qchrom'], fusion['fusion_point'][0], fusion['fusion_point'][1], "Fused_from:" + fusion['old_chroms'][0] + "&" + fusion['old_chroms'][1]))
# # then the inversions
# for inversion in inversions:
# 	outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("inversion", inversion['tchrom'], inversion['tstart'], inversion['tend'], inversion['qchrom'], inversion['qstart'], inversion['qend'], "Na"))
# # then the translocations
# for translocation in translocations:
# 	outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("translocation", translocation['tchrom'], translocation['tstart'], translocation['tend'], translocation['qchrom'], translocation['qstart'], translocation['qend'], "Na"))

# # close the outfile
# outfile.close()

# # write the unmerged blocks to a file
# outfile = open("unmerged_blocks.txt", 'w')
# outfile.write("refchrom\trefstart\trefend\tquerychrom\tquerystart\tqueryend\tdirection\n")
# for block in blocks:
# 	outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(block['tchrom'], block['tstart'], block['tend'], block['qchrom'], block['qstart'], block['qend'], block['direction']))

# # close the outfile
# outfile.close()
# # sort the list of dictionaries first by tchrom, then by qchrom and last by qstart
# sblocks = sorted(blocks, key=lambda k: (k['tchrom'], k['qchrom'], k['qstart']))

# # print the target and query length for all blocks
# for block in sblocks:
# 	print("%s\t%s\t%s\t%s" % (block['tchrom'], block['tend'] - block['tstart'], block['qchrom'], block['qend'] - block['qstart']))

# # print all items in the blocks as tchrom, tstart, tend, qchrom, qstart, qend
# for block in blocks:
# 	print("%s\t%s\t%s\t%s\t%s\t%s" % (block['tchrom'], block['tstart'], block['tend'], block['qchrom'], block['qstart'], block['qend']))

# # and merge items where tchrom and qchrom are the same, and tstart is greater than previous tend
# blocks_merged = []
# for i in range(len(blocks)):
# 	if i == 0:
# 		blocks_merged.append(blocks[i])
# 		continue
# 	if blocks[i]['tchrom'] == blocks_merged[-1]['tchrom'] and blocks[i]['qchrom'] == blocks_merged[-1]['qchrom'] and blocks[i]['qstart'] > blocks_merged[-1]['qend']:
# 		blocks_merged[-1]['tend'] = blocks[i]['tend']
# 		blocks_merged[-1]['qend'] = blocks[i]['qend']
# 	else:
# 		blocks_merged.append(blocks[i])

# # print as a table
# for b in blocks_merged:
# 	print("{}\t{}\t{}\t{}\t{}\t{}".format(b['tchrom'],b['tstart'],b['tend'],b['qchrom'],b['qstart'],b['qend']))

# # now, let's deal with translocations
# translocations = []
# for i in range(1,len(blocks_merged)):
# 	# if both chromosomes are the same, but qstart is less than previous qend, we have a translocation
# 	if blocks_merged[i]['tchrom'] == blocks_merged[i-1]['tchrom'] and blocks_merged[i]['qchrom'] == blocks_merged[i-1]['qchrom'] and blocks_merged[i]['qstart'] < blocks_merged[i-1]['qend']:
# 		translocations.append((blocks_merged[i-1], blocks_merged[i]))

# # print translocations as a table
# for block in translocations:
# 	print("%s\t%s\t%s\t%s\t%s\t%s" % (block[0]['tchrom'], block[0]['tstart'], block[0]['tend'], block[0]['qchrom'], block[0]['qstart'], block[0]['qend']))
# 	print("%s\t%s\t%s\t%s\t%s\t%s" % (block[1]['tchrom'], block[1]['tstart'], block[1]['tend'], block[1]['qchrom'], block[1]['qstart'], block[1]['qend']))
# 	print("")

# # now, we can identify fissions as blocks where tchrom is the same, but qchrom is different
# fissions = []
# tchrom = blocks_merged[0]['tchrom']
# qchrom = blocks_merged[0]['qchrom']
# tend = blocks_merged[0]['tend']
# tstart = blocks_merged[0]['tstart']
# for i in range(1,len(blocks_merged)):
# 	if blocks_merged[i]['tchrom'] == tchrom and blocks_merged[i]['qchrom'] != qchrom:
# 		# this means we have a fission, with a breaking point between the previous block and this one
# 		breaking_point = (tend, blocks_merged[i]['tstart'])
# 		new_chroms = ()


# # count how many positions are None in the the target positions
# print("Number of positions in the query genome that are not syntenic to the reference genome: %s" % len([x for x in synmap_q if x['tpos'] == None]))