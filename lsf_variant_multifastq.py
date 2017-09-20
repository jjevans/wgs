#!/usr/bin/env python
import align_util
import lsf_util
import os
import samtool_util
import seq_util
import sys
import time

''' call bwa, sampe, groom sam file, use GATK Queue to call variants '''
queue = "pcpgm"
lsfout = "lsf_out.txt"
lsferr = "lsf_err.txt"
hosts = "-m cmu033 -m cmu034 -m cmu035 -m cmu036 -m cmu037 -m cmu038 -m cmu023 -m cmu024"

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
bamclean = opt[base]+".merge.bam"
bamfile = opt[base]+".bam"
bamindex = opt[base]+".bai"

''' create objects from utility libraries '''
seq_obj = seq_util.SplitFile(opt[numseq])
bwa_obj = align_util.BWA(opt[prefix])

lsf_obj = lsf_util.LSF("-q "+queue+" -o "+lsfout+" -e "+lsferr+" "+hosts)
st_obj = samtool_util.Use()

''' split large fastq file into multiple smaller fastq files '''
print "******* begin splitting left (read 1) fastq: "+time.strftime("%c")+"*******"
left_files = seq_obj.split_fastq(opt[seqleft])
print "******* end splitting left (read 1) fastq: "+time.strftime("%c")+"*******"
print "******* begin splitting right (read 2) fastq: "+time.strftime("%c")+"*******"
right_files = seq_obj.split_fastq(opt[seqright])
print "******* end splitting right (read 2) fastq: "+time.strftime("%c")+"*******"

''' do alignments '''
print "******* begin bwa alignment: "+time.strftime("%c")+"*******"
rg = "@RG\tID:"+opt[rgid]+"\tSM:"+opt[rgsm]+"\tLB:"+opt[rglb]+"\tPL:"+opt[rgpl]+"\tPU:"+opt[rgpu]
samfiles = bwa_obj.multifastq_lsf(lsf_obj,left_files,right_files,rg,opt[base])
print "******* end bwa alignment: "+time.strftime("%c")+"*******"

''' groom sam files () '''
print "******* begin clean, sort of alignment sam: "+time.strftime("%c")+"*******"
cmd = "groom_sam.py"

bamfiles = list()
jobids = list()
for eachfile in samfiles:
	eachbam = eachfile+".bam"
	bamfiles.append(eachbam)
	
	cmdarg = [eachfile,eachbam,opt[picard]]
	jobids.append(lsf_obj.submit(cmd,cmdarg))
lsf_obj.sync(jobids)
print "******* end clean, sort of alignment sams: "+time.strftime("%c")+"*******"

''' merge sorted bams back together '''
print "******* begin merge of sorted bams: "+time.strftime("%c")+"*******"
cmd = "st_merge.py"
cmdarg = [bamclean]
cmdarg.extend(bamfiles)
jobid = lsf_obj.submit(cmd,cmdarg)
lsf_obj.wait(jobid)
print "******* end merge of sorted bams: "+time.strftime("%c")+"*******"

''' remove pcr duplicates from bam file '''
print "******* begin marking duplicates: "+time.strftime("%c")+"*******"
#st_obj.rm_dup(bamclean,bamfile)
cmd = "java"
cmdarg = ["-jar",opt[picard]+"/MarkDuplicates.jar","I="+bamclean,"O="+bamfile,"M="+bamfile+".metric","REMOVE_DUPLICATES=true","AS=true"]
jobid = lsf_obj.submit(cmd,cmdarg)
lsf_obj.wait(jobid)
print "******* end marking duplicates: "+time.strftime("%c")+"*******"

''' produce a bam index for merged bam '''
print "******* begin indexing bam: "+time.strftime("%c")+"*******"
st_obj.index(bamfile,bamindex)
print "******* end indexing bam: "+time.strftime("%c")+"*******"

''' call GATK with Queue '''
print "******* begin GATK: "+time.strftime("%c")+"*******"
GATKQ="/data/pcpgm/GATK-new-test/standard/tool/Queue/QueueLite.jar"
scalaQ="/PHShome/jje16/svn/variant-reporting/detection/trunk/variantCaller_v3.scala"
cmd = "java"
cmdarg = ["-Xms12g","-Xmx48g","-Djava.io.tmpdir=javatmpdir","-jar",GATKQ,"-S",scalaQ,"-cfg",conf,"-bam",bamfile,"-index",bamindex,"-bsub","-jobQueue",queue,"--disableJobReport","-startFromScratch","-run"]
jobid = lsf_obj.submit(cmd,cmdarg)
lsf_obj.wait(jobid)
print "******* end GATK: "+time.strftime("%c")+"*******"

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
