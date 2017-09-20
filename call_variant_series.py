#!/usr/bin/env python
import align_util
import samtool_util
import subprocess
import sys

''' call bwa, sampe, groom sam file, use GATK Queue to call variants '''

''' split variable '''
numfile = 10 # number of files to split sam into to sort and merge back
queue = "pcpgm"

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
bammerge = opt[base]+".merge.bam"
bamfile = opt[base]+".bam"
bamindex = opt[base]+".bai"

''' call bwa, make sam file '''
cmdarg = ["bwa","aln",opt[prefix],opt[seqleft],"-t",opt[thread],"-q",opt[thresh],"-f",saileft]
subprocess.check_call(cmdarg)

cmdarg = ["bwa","aln",opt[prefix],opt[seqright],"-t",opt[thread],"-q",opt[thresh],"-f",sairight]
subprocess.check_call(cmdarg)

# read group line put in in sampe
rg = "@RG\tID:"+opt[rgid]+"\tSM:"+opt[rgsm]+"\tLB:"+opt[rglb]+"\tPL:"+opt[rgpl]+"\tPU:"+opt[rgpu]

cmdarg = ["bwa","sampe",opt[prefix],saileft,sairight,opt[seqleft],opt[seqright],"-f",samfile,"-r",rg]
subprocess.check_call(cmdarg)

aln_obj = align_util.BreakFile(samfile)

files = aln_obj.make_files(numfile)

''' groom sam files () '''
cmd = "groom_sam.py"

bamfiles = list()
for eachfile in files:
	eachbam = eachfile+".bam"
	bamfiles.append(eachbam)
	
	cmdarg = [cmd,eachfile,eachbam,opt[picard]]
	subprocess.check_call(cmdarg)

''' merge sorted bams back together '''
cmd = "st_merge.py"
cmdarg = [cmd,bammerge]
cmdarg.extend(bamfiles)
subprocess.check_call(cmdarg)

''' remove duplicates and produce a bam index for merged bam '''
st_obj = samtool_util.Use()
st_obj.rm_dup(bammerge,bamfile)
st_obj.index(bamfile,bamindex)

''' call GATK with Queue '''
Q = "/data/pcpgm/GATK-new-test/standard/tool/Queue/QueueLite.jar"
scala = "/PHShome/jje16/svn/variant-reporting/detection/trunk/variantCaller_v2.scala"
cmdarg = ["java","-Xms12g","-Xmx48g","-Djava.io.tmpdir=javatmpdir","-jar",Q,"-S",scala,"-cfg",conf,"-bam",bamfile,"-index",bamindex,"-bsub","-jobQueue",queue,"--disableJobReport","-startFromScratch","-run"]
subprocess.check_call(cmdarg)
