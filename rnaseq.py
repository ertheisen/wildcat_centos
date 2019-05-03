#fastqc_pre.py
from datetime import datetime
import os
import glob
import shutil
import sys
from subprocess import check_output
    
startTime = datetime.now()

# Change directory within container
os.chdir('/data')
print 'Current working directory is:' +  os.getcwd()
os.mkdir('fastqc_pre_out')
os.chdir('/data/fastq')
print 'Current working directory is:' +  os.getcwd()

#Determine if fastq files are present
input_fastq = sorted(glob.glob('*.fastq'))

if len(input_fastq)==0:
	#import input filenames as a list and print out files
	input_files = sorted(glob.glob('*.fastq.gz'))
	if len(input_files)==0:
		print 'No files imported. Exiting...'
		sys.exit()

	print 'Input files:'
	print '\n'.join(input_files)
	print '\n'

	print 'Decompressing .gz files'
	#create fastq filename from .fastq.gz
	temp_list = [f.replace('.gz', '') for f in input_files]
	#run system command to decompress file with zcat
	for i in range(len(input_files)):
		temp_str = 'zcat ' + input_files[i] + ' > ' + temp_list[i]
		check_output(temp_str, shell=True)

	#replace list of input_files with temp_list so future operations are on decompressed files
	input_files = temp_list

else:
	input_files = input_fastq

print 'Decompressed files to be analyzed:'
print '\n'.join(input_files)
print '\n'


for i in range(len(input_files)):
	#create string for system command
	temp_str = 'fastqc --outdir=/data/fastqc_pre_out ' + input_files[i]

	print temp_str

	check_output(temp_str, shell=True)
	print '\n' 

print 'Runtime pre-processing fastq (hh:mm:ss): ' + str(datetime.now() - startTime)

startTime = datetime.now()

# Change directory within container
os.chdir('/data')
print 'Current working directory is:' +  os.getcwd()
print '\n'
os.mkdir('trim_out')
os.mkdir('trim_out/fastqc_post_trim')
os.chdir('/data/fastq')
print 'Current working directory is:' +  os.getcwd()
print '\n'

#Determine if fastq files are present
input_fastq = sorted(glob.glob('*.fastq'))
input_fastq_R1 = sorted(glob.glob('*R1*.fastq'))
input_fastq_R2 = sorted(glob.glob('*R2*.fastq'))

if len(input_fastq)==0:
	#import input filenames as a list and print out files
	input_files_R1 = sorted(glob.glob('*R1*.fastq.gz'))
	if len(input_files)==0:
		print 'No R1 files imported. Exiting...'
		sys.exit()

	print 'Input files R1:'
	print '\n'.join(input_files_R1)
	print '\n'

	print 'Decompressing .gz files'
	#create fastq filename from .fastq.gz
	temp_list_R1 = [f.replace('.gz', '') for f in input_files_R1]
	#run system command to decompress file with zcat
	for i in range(len(input_files_R1)):
		temp_str_R1 = 'zcat ' + input_files_R1[i] + ' > ' + temp_list_R1[i]
		check_output(temp_str_R1, shell=True)

	#replace list of input_files with temp_list so future operations are on decompressed files
	input_files_R1 = temp_list_R1

	#import input filenames as a list and print out files
	input_files_R2 = sorted(glob.glob('*R2*.fastq.gz'))
	if len(input_files)==0:
		print 'No R2 files imported. Exiting...'
		sys.exit()

	print 'Input files R2:'
	print '\n'.join(input_files_R2)
	print '\n'

	print 'Decompressing .gz files'
	#create fastq filename from .fastq.gz
	temp_list_R2 = [f.replace('.gz', '') for f in input_files_R2]
	#run system command to decompress file with zcat
	for i in range(len(input_files_R2)):
		temp_str_R2 = 'zcat ' + input_files_R2[i] + ' > ' + temp_list_R2[i]
		check_output(temp_str_R2, shell=True)

	#replace list of input_files with temp_list so future operations are on decompressed files
	input_files_R2 = temp_list_R2

else:
	input_files_R1 = input_fastq_R1
	input_files_R2 = input_fastq_R2

if len(input_files_R1)!=len(input_files_R2):
	print 'Unpaired files detected. ...Exiting.'
	sys.exit()

print 'Decompressed files to be analyzed:'
print '\n'.join(input_files_R1)
print '\n'.join(input_files_R2)
print '\n'

for i in range(len(input_files_R1)):
#create string for system command
	temp_str = 'trim_galore -o /data/trim_out --paired ' + input_files_R1[i] + ' ' + input_files_R2[i]

	print temp_str

	check_output(temp_str, shell=True)
	print '\n' 

os.chdir('/data/trim_out')
print 'Current working directory is:' +  os.getcwd()
fastqc_input = sorted(glob.glob('*.fq'))

print 'Trimmed files to be analyzed:'
print '\n'.join(fastqc_input)
print '\n'

for i in range(len(fastqc_input)):
#create string for system command
	temp_str = 'fastqc --outdir=fastqc_post_trim ' + fastqc_input[i]

	print temp_str

	check_output(temp_str, shell=True)
	print '\n' 


print 'Runtime trimming (hh:mm:ss): ' + str(datetime.now() - startTime)

startTime = datetime.now()

# Change directory within container
os.chdir('/data/trim_out')
print 'Current working directory is:' +  os.getcwd()

#prepare input files
input_files = sorted(glob.glob('*val*'))
input_files_R1 = []
input_files_R2 = []

if len(input_files) == 0:
	print 'No quality-controlled input files from trim_galore, checking input folder for fastqc output...'
	print '\n'
	input_files = sorted(glob.glob('*.fastq'))
	if len(input_files) == 0:
		print 'No uncompressed fastq files detected, looking for fastq.gz...'
		print '\n'
		input_files = sorted(glob.glob('*.fastq.gz*'))
		if len(input_files) == 0:
			print 'No valid input files detected. Exiting...'	
			sys.exit()
		else:
			print 'Decompressing .gz files'
			#create fastq filename from .fastq.gz
			temp_list = [f.replace('.gz', '') for f in input_files]
			#run system command to decompress file with zcat
			for i in range(len(input_files)):
				temp_str= 'zcat ' + input_files[i] + ' > ' + temp_list[i]
				check_output(temp_str, shell=True)
			input_files = temp_list

for item in input_files:
	if '_R1_' in item and '_R2_' in item:
		print "Input file: " + item + " contains both strings 'R1' and 'R2'. Not including..."
		sys.exit()
	elif '_R1_' in item and '_R2_' not in item:
		input_files_R1.append(item)
	elif '_R1_' not in item and '_R2_' in item:
		input_files_R2.append(item)
	else:
		print "Input file: " + item + "does not contain string 'R1' or 'R2'. Not including..."

if len(input_files_R1) != len(input_files_R2):
	print 'Unequal numbers of files assigned as R1 and R2. Check naming convention. Exiting...'
	sys.exit()

if (len(input_files_R1) + len(input_files_R2)) != len(input_files):
	print 'Not all of input files assigned as R1 or R2. Exiting...'
	sys.exit()



print 'Files assigned as R1:'
print '\n'.join(input_files_R1)
print '\n'

print 'Files assigned as R2:'
print '\n'.join(input_files_R2)
print '\n'

#make sam output names
STAR_prefixes = []

for i in range(len(input_files_R1)):
	STAR_prefix = input_files_R1[i].split('_R1')[0] + '_'
	STAR_prefixes.append(STAR_prefix)

print 'Output STAR prefixes:'
print '\n'.join(STAR_prefixes)
print '\n'

sam_names = []

for i in range(len(input_files_R1)):
	sam_name = input_files_R1[i].split('_R1')[0] + '_Aligned.out.sam'
	sam_names.append(sam_name)

STAR_logs = []

for i in range(len(input_files_R1)):
	STAR_log = input_files_R1[i].split('_R1')[0] + '_Log.out'
	STAR_logs.append(STAR_log)

for i in range(len(input_files_R1)):
	STAR_log = input_files_R1[i].split('_R1')[0] + '_Log.final.out'
	STAR_logs.append(STAR_log)

for i in range(len(input_files_R1)):
	STAR_log = input_files_R1[i].split('_R1')[0] + '_Log.progress.out'
	STAR_logs.append(STAR_log)

gene_counts = []

for i in range(len(input_files_R1)):
	gene_count = input_files_R1[i].split('_R1')[0] + '_ReadsPerGene.out.tab'
	gene_counts.append(gene_count)


print 'Output sam filenames:'
print '\n'.join(sam_names)
print '\n'

#define bowtie2 index for genome
print 'Data will be aligned to reference genome'

genome_index = '/genomes/STAR/'
print 'STAR index for data files found at:'
print genome_index
print '\n'

#run bowtie
print 'Running STAR alignment for reference genome followed by read counts'
print '\n'
for i in range(len(input_files_R1)):
	print 'count = ' +str(i)
	print '\n'
	#create string for system command
	temp_str = 'STAR --runThreadN 16 --genomeDir ' + genome_index + ' --readFilesIn ' + input_files_R1[i] + ' ' + input_files_R2[i] + ' --outFileNamePrefix ' + STAR_prefixes[i] + ' --quantMode GeneCounts --twopassMode Basic --outSAMunmapped Within'

	print temp_str

	check_output(temp_str, shell=True)
	print '\n'


#make new directory for output
os.mkdir('/data/sams')
os.mkdir('/data/sams/logs')
os.mkdir('/data/sams/counts')
print 'Current working directory is:' +  os.getcwd()
print '\n'

#copy files to output folder
output_dir = '/data/sams/'
print 'Moving sam files to output folder'
print '\n'
for i in range(len(sam_names)):
	shutil.move(sam_names[i], output_dir)

output_dir = '/data/sams/logs'
print 'Moving sam files to output folder'
print '\n'
for i in range(len(STAR_logs)):
	shutil.move(STAR_logs[i], output_dir)

output_dir = '/data/sams/counts'
print 'Moving sam files to output folder'
print '\n'
for i in range(len(gene_counts)):
	shutil.move(gene_counts[i], output_dir)


print 'Alignment Runtime (hh:mm:ss): ' + str(datetime.now() - startTime)
print '\n'

###SAM conversion to bam, bedgraph, and BigWig

os.mkdir('/data/bams')
os.chdir('/data/sams')
print 'Current working directory is:' +  os.getcwd()
print '\n'

import pybedtools
from pybedtools import BedTool
import pandas as pd
from pybedtools.helpers import chromsizes
from pybedtools.contrib.bigwig import bedgraph_to_bigwig
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

startTime = datetime.now()


print 'Converting sams to bams, bedgraph, and BigWig without spike-in normalization.'
print '\n'
datafiles = sorted(glob.glob('*.sam'))

print '\n'
print 'Data files loaded:'
print '\n'.join(datafiles)
print '\n'


##convert to bam format
print 'Converting to bam format'
print '\n'

bam_names = []

for i in range(len(sam_names)):
	bam_name = sam_names[i].split('_Aligned')[0] + '.bam'
	bam_names.append(bam_name)

print '\n'
print '\n'.join(bam_names)
print '\n'
print '\n'
print 'SAM to BAM'
print '\n'

bam_string = []

for i in range(len(bam_names)):
        bam_string.append('samtools view -b -S ' + datafiles[i] + ' > /data/bams/' + bam_names[i])

for item in bam_string:
        check_output(item, shell = True)

datafiles = bam_names

##Sort and index
print '\n'
print 'Sorting bams'
print '\n'

sorted_bam_names = []
for i in range(len(bam_names)):
	sorted_bam_name = 'sorted.' + bam_names[i]
	sorted_bam_names.append(sorted_bam_name)

bam_names_sorted = [f.replace('.bam', '') for f in sorted_bam_names]

sort_string = []

for i in range(len(bam_names)):
		sort_string.append('samtools sort /data/bams/' + bam_names[i] + ' /data/bams/' + bam_names_sorted[i])

for item in sort_string:
		check_output(item, shell = True)

print '\n'
print 'Bam files to index:'
print '\n'.join(bam_names_sorted)
print '\n'


print 'Indexing bams'
print '\n'             
os.chdir('/data/bams')
print 'Current working directory is:' +  os.getcwd()

index_string = []

for i in range(len(bam_names_sorted)):
        index_string.append('samtools index ' + sorted_bam_names[i])

for item in index_string:
        check_output(item, shell = True)

sorted_files = sorted(glob.glob('sorted*'))
os.mkdir('/data/bams/sorted_bams')
output_dir = '/data/bams/sorted_bams'
for i in range(len(sorted_files)):
	shutil.move(sorted_files[i], output_dir)

	
##generate bigwig files from bam files
print 'Generating bigwig files'
print '\n'

os.chdir('/data/bams/sorted_bams')

print 'Current working directory is:' +  os.getcwd()
print '\n'

#generate bigwig file names
bw_names = [f.replace('bam', 'bw') for f in datafiles]


bigwig_string = []

for i in range(len(sorted_bam_names)):
		bigwig_string.append('bamCoverage -b '+ sorted_bam_names[i] + ' /data/bams/sorted_bams/' + bw_names[i])

for item in bigwig_string:
		check_output(item, shell = True)


print 'Finished generating bigwig files:'
print '\n'
print 'whole insert bigwig files:' + '\n' + '\n'.join(bw_names)
print '\n'

os.mkdir('/data/bigwigs')

print 'Moving bigwigs to output folder'
print '\n'
output_dir0 = '/data/bigwigs'
for i in range(len(bw_names)):
	shutil.move(bw_names[i], output_dir0)


print 'Finished'
print '\n'
print 'Runtime (hh:mm:ss): ' + str(datetime.now() - startTime)
