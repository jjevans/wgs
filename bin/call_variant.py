#!/usr/bin/env python
import align_util
import lsf_util
import samtool_util
import subprocess
import sys

''' call bwa, sampe, groom sam file, use GATK Queue to call variants '''

''' lsf and split variables '''
numfile = 10 # number of files to split sam into to sort and merge back
queue = "pcpgm"
#queue = "normal"
lsfout = "lsf_output.txt"

''' keywords for configuration file values '''
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
	print "usage: call_variant.py config_file"
	
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
bammerge = opt[base]+".merge.bam"
bamfile = opt[base]+".bam"
bamindex = opt[base]+".bai"

# initialize lsf
lsf_obj = lsf_util.LSF("-q "+queue+" -o "+lsfout)
jobs = list()
	
''' call bwa, make sam file '''
cmd = "bwa"
cmdarg = ["aln",opt[prefix],opt[seqleft],"-t",opt[thread],"-q",opt[thresh],"-f",saileft]
jobs.append(lsf_obj.submit(cmd,cmdarg))

cmdarg = ["aln",opt[prefix],opt[seqright],"-t",opt[thread],"-q",opt[thresh],"-f",sairight]
jobs.append(lsf_obj.submit(cmd,cmdarg))

lsf_obj.sync(jobs)

# read group line put in in sampe
rg = "@RG\tID:"+opt[rgid]+"\tSM:"+opt[rgsm]+"\tLB:"+opt[rglb]+"\tPL:"+opt[rgpl]+"\tPU:"+opt[rgpu]

cmdarg = ["sampe",opt[prefix],saileft,sairight,opt[seqleft],opt[seqright],"-f",samfile,"-r",rg]
jobid = lsf_obj.submit(cmd,cmdarg)
lsf_obj.wait(jobid)

aln_obj = align_util.BreakFile(samfile)

files = aln_obj.make_files(numfile)

''' groom sam files () '''
cmd = "groom_sam.py"

bamfiles = list()
for eachfile in files:
	eachbam = eachfile+".bam"
	bamfiles.append(eachbam)
	
	cmdarg = [eachfile,eachbam,opt[picard]]
	jobs.append(lsf_obj.submit(cmd,cmdarg))

lsf_obj.sync(jobs)

''' merge sorted bams back together '''
cmd = "st_merge.py"
cmdarg = [bammerge]
cmdarg.extend(bamfiles)
jobid = lsf_obj.submit(cmd,cmdarg)
lsf_obj.wait(jobid)

''' remove duplicates and produce a bam index for merged bam '''
st_obj = samtool_util.Use()
st_obj.rm_dup(bammerge,bamfile)
st_obj.index(bamfile,bamindex)

''' call GATK with Queue '''
cmd = "java"
Queue = "/data/pcpgm/GATK-new-test/standard/tool/Queue/QueueLite.jar"
scala = "/PHShome/jje16/svn/variant-reporting/detection/trunk/variantCaller_v2.scala"
cmdarg = ["-Xms12g","-Xmx48g","-Djava.io.tmpdir=javatmpdir","-jar",Queue,"-S",scala,"-cfg",conf,"-bam",bamfile,"-index",bamindex,"-bsub","-jobQueue",queue,"--disableJobReport","-startFromScratch","-run"]
jobid = lsf_obj.submit(cmd,cmdarg)
lsf_obj.wait(jobid)

