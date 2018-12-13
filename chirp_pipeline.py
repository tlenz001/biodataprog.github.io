#!/usr/bin/env python3

import subprocess
import sys
import os
import numpy

# files must be separated into subfolders by read pairs
# control input must end in _ChIRPInput.fastq.gz

def input_dir(): # module to request user input for parent directory of ChiRP files and prevent invalid entries
    i = 0
    indir = ""
    while i < 3: # prevents an infinite loop and forces module exit after 3 invalid attempts
        indir = input("Enter directory of ChiRP files: ") # requests user input
        if not os.path.exists(indir) or os.path.isfile(indir): # checks that directory exists and is not a file
            i += 1
            if i < 3: # printed each time an invalid directory is entered
                print("")
                print("User input is not a valid directory.")
            if i == 1:
                print("2 more chances to enter a valid directory.")
            if i == 2:
                print("1 more chance to enter a valid directory.")
        else: # if a valid directory is entered it is saved as the input of the parent module
            return indir
    if i == 3: # module exits when 3 invalid attempts occur
        print("You have exceeded the maximum allowable attempts. Double-check your directory and try again.")
        print("")
        exit()

directory = input_dir() # user input is saved as the parent directory for ChiRP files used in this program

def input_organism(): # module to request user input for parent directory of fasta file and prevent invalid entries
    i = 0
    org = ""
    while i < 3: # prevents an infinite loop and forces module exit after 3 invalid attempts
        org = input("Enter directory of fasta file: ") # requests user input for fasta directory
        if not os.path.exists(org) or os.path.isfile(org): # checks that directory exists and is not a file
            i += 1
            if i < 3:
                print("")
                print("User input is not a valid directory.") # printed each time an invalid directory is entered
            if i == 1:
                print("2 more chances to enter a valid directory.")
            if i == 2:
                print("1 more chance to enter a valid directory.")
        else:
            return org # if a valid directory is entered it is saved as the input of the parent module
    if i == 3: # module exits when 3 invalid attempts occur
        print("You have exceeded the maximum allowable attempts. Double-check your directory and try again.")
        print("")
        exit()

organism = input_organism()

def bt2_org(): # module saves fasta file for bowtie2
    bt2 = []
    for dir, subdir, file in os.walk(organism):
        for item in file:
            if item.endswith(".fasta"):
                bt2.append(dir + '/' + item)
    bt2_org = ''.join(bt2).strip(".fasta")
    return bt2_org

bt2_index = bt2_org()

def input_chrom(): # module saves .chrom.sizes file for bedtools genomecov
    chrom_file = []
    for dir, subdir, file in os.walk(organism):
        for item in file:
            if item.endswith('.chrom.sizes'):
                chrom_file.append(dir + '/' + item)
    chrom_sizes = ''.join(chrom_file)
    return chrom_sizes

chrom = input_chrom()

def trimmomatic_search(): # module creates path to trimmomatic .jar file for easier loading in java
    trimmomatic_dir_search = "/Users"
    trimmomatic_path = []
    for dir, subdir, file in os.walk(trimmomatic_dir_search):
        for item in file:
            if item == "trimmomatic-0.38.jar":
                trimmomatic_path.append(dir + '/' + item)
    trim = ''.join(trimmomatic_path)
    return trim

trimmomatic = trimmomatic_search()

def picard_search(): # module creates path to picardtools .jar file for easier loading in java
    picard_dir_search = "/Users"
    picard_path = []
    for dir, subdir, file in os.walk(picard_dir_search):
        for item in file:
            if item == "picard.jar":
                picard_path.append(dir + '/' + item)
    picard = ''.join(picard_path)
    return picard

picardtools = picard_search()

for dir, subdir, file in os.walk(directory):
    file_list = []
    for item in file:
        if item.endswith('.fastq.gz'): # checks that the parent directory contains .fastq.gz files
            file_list.append(dir + '/' + item)
    if len(file) < 2: # checks that directory contains pairs of .fastq.gz files
        continue
    file_list.sort()
    subprocess.run(['java', '-Xmx32g', '-jar', '{}'.format(trimmomatic), \
    'PE', '-threads', '12', '{}'.format(file_list[0]), '{}'.format(file_list[1]), \
    '{}'.format(file_list[0]).replace(".fastq.gz", "_paired.fastq"), \
    '{}'.format(file_list[0]).replace(".fastq.gz", "_unpaired.fastq"), \
    '{}'.format(file_list[1]).replace(".fastq.gz", "_paired.fastq"), \
    '{}'.format(file_list[1]).replace(".fastq.gz", "_unpaired.fastq"), \
    'HEADCROP:2', 'CROP:73', 'ILLUMINACLIP:/Users/toddlenz/Desktop/fasta/Illumina_adapters_polys.fasta:2:30:10'])
    # runs a unix subprocess for Trimmomatic software
    # used on paired-end reads
    # trims adapter sequences and bases from beginning and end of reads

print("")
print("==================================================")
print("")

for dir, subdir, file in os.walk(directory):
    paired = []
    for item in file:
        if item.endswith('_paired.fastq'):
            paired.append(dir + '/' + item)
    if len(paired) < 2:
        continue
    paired.sort()
    print(paired)
    subprocess.call(['sickle', 'pe', \
    '-f', '{}'.format(paired[0]), '-r', '{}'.format(paired[1]), '-t', 'sanger', \
    '-o', '{}'.format(paired[0]).replace("_paired.fastq", "_filtered.fastq"), \
    '-p', '{}'.format(paired[1]).replace("_paired.fastq", "_filtered.fastq"), \
    '-s', '{}'.format(paired[0]).replace("_pair1", ""), '-q', '25', '-l', '18'])
    # runs a unix subprocess for sickle software
    # filters paired-end reads by quality score and length

print("Sickle: Completed successfully")
print("==================================================")
print("")

for dir, subdir, file in os.walk(directory):
    filtered = []
    for item in file:
        if item.endswith('_filtered.fastq'):
            filtered.append(dir + '/' + item)
    if len(filtered) < 2:
        continue
    filtered.sort()
    print(filtered)
    subprocess.call(['bowtie2', '-p', '12', '-x', '{}'.format(bt2_index), \
    '-1', '{}'.format(filtered[0]), '-2', '{}'.format(filtered[1]), \
    '-S', '{}'.format(filtered[0]).replace("_pair1", "").replace("_filtered.fastq", "_BT2_aligned.sam")])
    # runs a unix subprocess for bowtie2 software
    # aligns filtered paired reads to the parent genome

print("")
print("Bowtie2: Completed successfully")
print("==================================================")
print("")

for dir, subdir, file in os.walk(directory):
    aligned = []
    for item in file:
        if item.endswith('_BT2_aligned.sam'):
            aligned.append(dir + '/' + item)
    if len(aligned) < 1:
        continue
    print(aligned)
    subprocess.call(['samtools', 'view', '-q', '10', '-f', '0X02', '-F', '0X04', '-@', '12', \
    '-b', '{}'.format(aligned[0]), '-o', '{}'.format(aligned[0]).replace(".sam", ".bam")])
    # runs a unix subprocess for samtools software
    # converts bowtie2 aligned .sam file to a .bam file
    # filters alignments by quality of alignment

for dir, subdir, file in os.walk(directory):
    sorted = []
    for item in file:
        if item.endswith('_BT2_aligned.bam'):
            sorted.append(dir + '/' + item)
    if len(sorted) < 1:
        continue
    print(sorted)
    subprocess.call(['samtools', 'sort', '{}'.format(sorted[0]), \
    '-o', '{}'.format(sorted[0]).replace("_BT2_aligned.bam", "_sorted.bam")])
    # runs a unix subprocess for samtools software
    # sorts reads in .bam file by coordinates in genome

for dir, subdir, file in os.walk(directory):
    flag = []
    for item in file:
        if item.endswith('_sorted.bam'):
            flag.append(dir + '/' + item)
    if len(flag) < 1:
        continue
    print(flag)
    os.system('samtools stats {0} > {1}'.format(flag[0], flag[0].replace(".bam", "_stats.txt")))
    # runs a unix subprocess for samtools software
    # creates a .txt file containing read data for .bam file

print("")
print("Samtools: Completed successfully")
print("==================================================")
print("")

for dir, subdir, file in os.walk(directory):
    dedup = []
    for item in file:
        if item.endswith('_sorted.bam'):
            dedup.append(dir + '/' + item)
    if len(dedup) < 1:
        continue
    print(dedup)
    subprocess.call(['java', '-Xmx32g', '-jar', '{}'.format(picardtools), \
    'MarkDuplicates', 'I={}'.format(dedup[0]), 'O={}'.format(dedup[0]).replace("_sorted.bam", "_dedup.bam"), \
    'M={}'.format(dedup[0]).replace("_sorted.bam", "_dedup_metrics.txt"), 'REMOVE_DUPLICATES=true'])
    # runs a unix subprocess for picardtools software
    # removes duplicate reads from sorted .bam file

for dir, subdir, file in os.walk(directory):
    flag = []
    for item in file:
        if item.endswith('_dedup.bam'):
            flag.append(dir + '/' + item)
    if len(flag) < 1:
        continue
    print(flag)
    os.system('samtools stats {0} > {1}'.format(flag[0], flag[0].replace(".bam", "_stats.txt")))
    # runs a unix subprocess for samtools software
    # creates a .txt file containing read data for deduplicated .bam file

print("")
print("Picardtools: Completed successfully")
print("==================================================")
print("")

for dir, subdir, file in os.walk(directory):
    index = []
    for item in file:
        if item.endswith('_dedup.bam'):
            index.append(dir + '/' + item)
    if len(index) < 1:
        continue
    print(index)
    subprocess.call(['samtools', 'index', '{}'.format(index[0]), \
    '{}'.format(index[0]).replace("_dedup.bam", "_dedup.bam.bai")])
    # runs a unix subprocess for samtools software
    # indexes the deduplicated .bam file

print("")
print("Samtools: Completed successfully")
print("==================================================")
print("")

for dir, subdir, file in os.walk(directory):
    bed = []
    for item in file:
        if item.endswith('_dedup.bam'):
            bed.append(dir + '/' + item)
    if len(bed) < 1:
        continue
    print(bed)
    os.system('Bedtools bamtobed -i {0} > {1}'.format(bed[0], bed[0].replace(".bam", ".bed")))
    # runs a unix subprocess for bedtools software
    # converts deduplicated .bam file to .bed file

for dir, subdir, file in os.walk(directory):
    cov = []
    for item in file:
        if item.endswith('.bed'):
            cov.append(dir + '/' + item)
    if len(cov) < 1:
        continue
    cov.sort()
    print(cov)
    os.system('Bedtools genomecov -i {0} -g {1} -d > {2}'.format(cov[0], chrom, cov[0].replace("_dedup.bed", "_coverage.txt")))
    # runs a unix subprocess for bedtools software
    # creates a .txt file that displays read coverage per nucleotide in the .bed file

print("")
print("Bedtools: Completed successfully")
print("==================================================")
print("")

def input_coverage(): # module normalizes read data in control experiment and saves list of read counts
    input_cov_reads = []
    for dir, subdir, file in os.walk(directory):
        norm_input = []
        nuc_input = []
        for item in file:
            if item.endswith('ChIRPInput_dedup_stats.txt'):
                norm_input.append(dir + '/' + item)
            if item.endswith('ChIRPInput_coverage.txt'):
                nuc_input.append(dir + '/' + item)
        if len(norm_input) < 1:
            continue
        else:
            normalized_input = ''.join(norm_input)
        if len(nuc_input) < 1:
            continue
        else:
            nts_in = ''.join(nuc_input)
        with open(normalized_input, 'r') as fh:
            for line in fh:
                line = line.strip("\n").split("\t")
                if ("reads mapped:") in line:
                    mmr_in = int(line[2])/1000000 # divides total reads mapped by 1000000 to get million mapped reads
            with open(nts_in, 'r') as f:
                for line in f:
                    line = line.strip("\n").split("\t")
                    reads_in = int(line[2])/mmr_in # divides read count at each nuclotide by the million mapped reads
                    input_cov_reads.append(reads_in)
                return input_cov_reads

in_cov = input_coverage()

for dir, subdir, file in os.walk(directory): # module normalizes the read counts per nucleotide in each _coverage.txt file
    norm = []
    nuc = []
    for item in file:
        if item.endswith('_dedup_stats.txt') and not item.endswith('ChIRPInput_dedup_stats.txt'):
            norm.append(dir + '/' + item)
        if item.endswith('_coverage.txt') and not item.endswith('ChIRPInput_coverage.txt'):
            nuc.append(dir + '/' + item)
    if len(norm) < 1:
        continue
    else:
        normalized = ''.join(norm)
    if len(nuc) < 1:
        continue
    else:
        nts = ''.join(nuc)
    norm_reads = (nts).replace("_coverage.txt", "_norm_reads.txt")
    with open(normalized, 'r') as fh:
        for line in fh:
            line = line.strip("\n").split("\t")
            if ("reads mapped:") in line:
                mmr = int(line[2])/1000000
        with open(nts, 'r') as f:
            print(nts)
            genes = []
            nucleotides = []
            cov = []
            nc = open(norm_reads, 'a+') # creates/opens _norm_reads.txt file for editing
            for line in f:
                line = line.strip("\n").split("\t")
                reads = int(line[2])/mmr
                genes.append(str(line[0]))
                nucleotides.append(str(line[1]))
                cov.append(reads)
            normalized_coverage = numpy.array(cov) - numpy.array(in_cov) # subtracts coverage per nucleotide in control from replicates
            norm_cov = numpy.clip(normalized_coverage, 0, None) # sets minimum threshold for read count to 0
            merged = numpy.column_stack((genes, nucleotides, norm_cov)) # formats list of normalized coverage for printing
            for item in merged:
                nc.write("%s\n" % "\t".join(item)) # appends list of normalized reads to .txt file

print("")
print("Normalization: Completed successfully")
print("==================================================")
print("")

for dir, subdir, file in os.walk(directory):
    wigs = []
    for item in file:
        if item.endswith('_norm_reads.txt'):
            wigs.append(dir + '/' + item)
    if len(wigs) < 1:
        continue
    else:
        wig = ''.join(wigs)
    wig_file = (wig).replace("_norm_reads.txt", "_norm_reads.wig")
    with open(wig, 'r') as fh:
        print(wig)
        reads_num = []
        wig_edit = open(wig_file, 'a+')
        reads_num.append("track type=wiggle_0") # sets first line in .wig file for proper formatting
        for line in fh:
            line = line.strip("\n").split("\t")
            if line[1] is "1":
                reads_num.append("fixedStep chrom=%s start=1 step=1" % line[0]) # formats each section header of .wig file
            reads_num.append(line[2])
        wig_edit.write("%s" % "\n".join(reads_num))
        # creates a .wig file from the normalized reads .txt file for viewing in IGV

print("")
print(".wig File Creation: Completed successfully")
print("==================================================")
print("")

for dir, subdir, file in os.walk(directory):
    tdfs = []
    for item in file:
        if item.endswith("_norm_reads.wig"):
            tdfs.append(dir + '/' + item)
    if len(tdfs) < 1:
        continue
    else:
        tdf = ''.join(tdfs)
    print(tdf)
    subprocess.call(['igvtools', 'toTDF', '{}'.format(tdf), '{}'.format(tdf).replace(".wig", ".tdf"), '{}'.format(chrom)])
    # runs a unix subprocess for IGVtools software
    # creates a .tdf file from .wig file for easier viewing in IGV

print("")
print("IGVtools toTDF: Completed successfully")
print("==================================================")
print("")
