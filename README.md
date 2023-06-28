# minimap2synteny
Create a stacked syntenymap and plot it from one or more whole genome alignment, against the same reference genome.

From one or more paf whole genome alignment files to the same reference genome, the alnToSyntenyMap.py creates syntenymaps in a progressive manner,
such that the first query genomes synteny is displayed against the reference genome, the second against the first query genome etc.

## Usage

Input files needed:
- reference.genomefile.txt, list of chromosomes and their lenghts of the reference genome that was used in the alignments, in two columns. More columns will work too but only the first two are used (i.e. a faidx index file will work).
- alignment_1.paf, at least one alignment file
- query_genome_1.txt, and at least one genome file for the query genome used in the alignment
optionally:
- alignment_2.paf
- query_genome_2.txt
- alignment_3.paf
- query_genome_3.txt etc...

## Step 1
Generate the syntenymaps:

python3 alnToSyntenyMap.py --reference reference.genomefile.txt \
  --alginment alignment_1.paf (alignment_2.paf alignment_3.paf...) \
  --index query_genomefile_1.txt (query_genomefile_2.txt query_genomefile_3.txt) \
  --output output_prefix
  
## Step 2
The above script will create two files: output_prefix.chroms.txt and output_prefix.syntenyMap.txt. To plot this in a ribbon plot, use the rscript plotSyntenyMaps.R:

Rscript plotSyntenyMaps.R output_prefix.chroms.txt output_prefix.syntenyMap.txt (label_refgenome,label_genome1,label_genome2...)

Will output an svg file in the same directory as the syntenymap is stored.


# Identify rearrangements

The script identifyRearrangements.py outputs three tables (largely overlapping with the above scripts) given an alignment, refgenome file and query genome file as input:
- The two genomes sorted by ref coordinates
- syntenymap between the two
- identified rearrangements (translocations, inversions, fissions and fusions)

run as: 
python3 identifyRearrangements.py -a alignment.paf \
  -q query_genomefile.txt \
  -r ref_genomefile.txt \
  -o outprefix

Output cak be plotted with the plotChromRearrangements.R:
Rscript plotChromRearrangements.R outprefix.genomes.txt outprefix.syntenymap.txt outprefix.rearrangements.txt

Note that this script so far only works with pairwise comparisons.