#!/usr/bin/env python
import align_util
import os
import samtool_util
import seq_util
import subprocess
import sys

''' call bwa, sampe, groom sam file, use GATK Queue to call variants '''
queue = "pcpgm"
	
''' keywords for configuration file values '''
numseq = "numseq" # number of files to split sam into to sort and merge back
seqleft = "seqleft"
seqright = "seqright"
picard = "picard"
reference = "reference"
prefix = "prefix"
thread = "thread"
thresh = "thresh"
rgid = "RGID"
rgsm = "RGSM"
rglb = "RGLB"
rgpl = "RGPL"
rgpu = "RGPU"
base = "base"

try:
	conf = sys.argv[1]
except IOError as (errno,strerror):
	print "usage: call_variant_series.py config_file"
	
''' load config file into dict '''
opt = dict()
with open(conf) as handle:
	for line in handle.readlines():
		if not line.startswith("//") and line.rstrip() is not "":
			val = line.split("=")
			opt[val[0]] = val[1].rstrip()

# filename intermediates
saileft = opt[base]+".left.sai"
sairight = opt[base]+".right.sai"
samfile = opt[base]+".sam"
bamclean = opt[base]+".clean.bam"
bamfile = opt[base]+".bam"
bamindex = opt[base]+".bai"

''' create objects from utility libraries '''
seq_obj = seq_util.SplitFile(opt[numseq])
bwa_obj = align_util.BWA(opt[prefix])
st_obj = samtool_util.Use()

''' split large fastq file into multiple smaller fastq files '''
left_files = seq_obj.split_fastq(opt[seqleft])
right_files = seq_obj.split_fastq(opt[seqright])

''' do alignments '''
rg = "@RG\tID:"+opt[rgid]+"\tSM:"+opt[rgsm]+"\tLB:"+opt[rglb]+"\tPL:"+opt[rgpl]+"\tPU:"+opt[rgpu]
samfiles = bwa_obj.multifastq_call(left_files,right_files,rg,opt[base])

''' groom sam files () '''
cmd = "groom_sam.py"

bamfiles = list()
jobids = list()
for eachfile in samfiles:
	eachbam = eachfile+".bam"
	bamfiles.append(eachbam)
	
	cmdarg = [cmd,eachfile,eachbam,opt[picard]]
	subprocess.check_call(cmdarg)

''' merge sorted bams back together '''
cmd = "st_merge.py"
cmdarg = [cmd,bamclean]
cmdarg.extend(bamfiles)
subprocess.check_call(cmdarg)

''' remove pcr duplicates from bam file '''
st_obj.rm_dup(bamclean,bamfile)

''' produce a bam index for merged bam '''
st_obj.index(bamfile,bamindex)

''' call GATK with Queue '''
GATKQ="/data/pcpgm/GATK-new-test/standard/tool/Queue/QueueLite.jar"
scalaQ="/PHShome/jje16/svn/variant-reporting/detection/trunk/variantCaller_v2.scala"
cmdarg = ["java","-Xms12g","-Xmx48g","-Djava.io.tmpdir=javatmpdir","-jar",GATKQ,"-S",scalaQ,"-cfg",conf,"-bam",bamfile,"-index",bamindex,"-bsub","-jobQueue",queue,"--disableJobReport","-startFromScratch","-run"]
subprocess.check_call(cmdarg)

"""
''' cleanup intermediate files '''
for seqfile in left_files:
	os.remove(seqfile)
	print seqfile
for seqfile in right_files:
	os.remove(seqfile)
for alignfile in samfiles:
	os.remove(alignfile)
for alignfile in bamfiles:
	os.remove(alignfile)
os.remove(bamclean)
"""
