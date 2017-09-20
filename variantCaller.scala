//package variantQ

// $Id$

/*
	Run Queue for GATK 2.0.  Modelled after python script for variant calling 
	named align_call_variants.py.

	Jason Evans, 10012012
*/

import scala.io.Source
import java.io.File

import org.broadinstitute.sting.queue.util.Logging
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.extensions.picard._
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine.OUTPUT_MODE

import net.sf.samtools.SAMFileHeader.SortOrder

class variantCaller extends QScript with Logging {
	qscript =>

	// config file locale
	@Input(doc="path to the config file", fullName="config", shortName="cfg", required=true)
	var conffile: File = _
	
	@Input(doc="skip alignment (1), skip doing the metrics (2), 12 skips both", fullName="skip", shortName="skp", required=false)
	var skip: String = _  // fix so can use on command line, BROKEN
	
	////
    /* global vars */
    
	// leave or remove intermediate files on filesystem
	val leaveFile = false

	// variable for string containing java params
	//val javaMem: String = "-Xms12g -Xmx48g"
	val javaMem: String = "-Xmx48g"
	
    // assignment of parameter variables to keys defined in configuration file
    val seqleft = "seqleft"
    val seqright = "seqright"
	val outBase = "base"
    val ref = "reference"
    val knownsnp = "snp"
    val knownindel = "indel"
    val picard = "picard"
    val bwaPrefix = "prefix"
    val bwaThread = "thread"
    val bwaThresh = "thresh"
    val rgid = "RGID"
    val rgsm = "RGSM"
    val rglb = "RGLB"
    val rgpl = "RGPL"
    val rgpu = "RGPU"
    val bait = "bait"
    val roi = "roi"
    val dbsnp = "dbsnp"
    
    // names of the files specified indicating what step to begin with
	val stepName: Seq[String] = Seq("aln","sampe","clean","bam","sort","duplicate","index")
		
    // sort order always coordinate
    val order = SortOrder.coordinate
   	
	//////
    //////
	def script() {
		// Queue process
		
		//logger.debug("skip value: "+skip)
	  
		//////
		
		// covariates for baserecalibrator
		val covs = Seq("ReadGroupCovariate","QualityScoreCovariate","CycleCovariate","ContextCovariate")
	  
		// read and load configuration file containing run parameters
		val configure = new Configure
		val conf: Map[String,String] = configure.loadConf(conffile)

		// sequence files
		val seqs: Seq[String] = Seq(conf(seqleft),conf(seqright)) // fastq files to analyze defined in config
		
		// find step to begin with
		//val startpoint = new startPoint
		//val step: String = startPoint.findStep(conf,stepName)
		
		// INTERMEDIATE FILE EXTENSIONS
    	// intermediate filenames with extensions for alignment, associated classes
    	val sam = conf(outBase) + ".sam"	// bwa makeSam
    	val rg = conf(outBase) + ".rg.sam"  // readGroup
    	val clean = conf(outBase) + ".clean.sam"	//	cleanSam 
    	val bam = conf(outBase) + ".bam"	//	samToBam
    	val sort = conf(outBase) + ".sort.bam"	//sortSam
    	val dup = conf(outBase) + ".dup.bam"	//markDuplicates
    	val dupbai = conf(outBase) + ".bai"	// bamIndex
    		
    	// intermediate filename with extensions for quality recalibration
    	val interval = conf(outBase) + ".bed"	// createTarget
    	val realign = conf(outBase) + ".real"	// indelRealign	
    	val fix = conf(outBase) + ".fix.bam"	// fixMateInfo
    	val fixbai = conf(outBase) + ".fix.bai"	// bamIndex
    	//val rec = conf(outBase) + ".rec"	
    	//val bqsr = conf(outBase) + ".grp"
    	val pre = conf(outBase) + ".pre.bam"
    	val read = conf(outBase) + ".read.bam"
    	val readbai = conf(outBase) + ".pre.bai"
    	val post = conf(outBase) + ".post.bam"
    	
    	// intermediate filenames with extensions for producing qc metrics
    	val indexstat = conf(outBase) + ".is.tbl"
    	val hsmetric = conf(outBase) + ".hs.tbl"
    	val alignmetric = conf(outBase) + ".alignsmry.tbl"
    	val gcbiastbl = conf(outBase) + ".gcbias.tbl"
    	val gcbiaschart = conf(outBase) + ".gcbias.pdf"
    	val gcbiassmry = conf(outBase) + ".gcbias.smry.tbl"
    	val bycycletbl = conf(outBase) + ".bycycle.tbl"
    	val bycyclechart = conf(outBase) + ".bycycle.pdf"
    	val qualdisttbl = conf(outBase) + ".qualdist.tbl"
    	val qualdistchart = conf(outBase) + ".qualdist.pdf"
    	val inserttbl = conf(outBase) + ".insert.tbl"
    	val inserthist = conf(outBase) + ".insert.hist"
    	val complextbl = conf(outBase) + ".complexity.tbl"
    	
    	// intermediate filenames with extensions for calling variants 
    	val indelcall = conf(outBase) + ".indel.vcf"
    	val indelmetric = conf(outBase) + ".indel.tbl"
    	val snpcall = conf(outBase) + ".snp.vcf"
    	val snpmetric = conf(outBase) + ".snp.tbl"
    	val bothcall = conf(outBase) + ".both.vcf"
    	val bothmetric = conf(outBase) + ".both.tbl"
    	val allcall = conf(outBase) + ".all.vcf"
    	val allmetric = conf(outBase) + ".all.tbl"
    	
    	// produce files that will be used in multiple blocks
    	val refFile = new File(conf(ref))
    	val dupFile = new File(dup)
    	val readFile = new File(read)
    			
		
    	////
		// PROCESS
    	//// 
		
		////
		// alignment and post-process, create index
		
		if(!skip.contains("1")){
			val samFile = new File(sam)
	    	//val rgFile = new File(rg)
	    	val cleanFile = new File(clean)
	    	val bamFile = new File(bam)
	    	val sortFile = new File(sort)
	    	val dupbaiFile = new File(dupbai)
    		
	    	// create a sai filename for each inputted file and put in list
	    	var sais: Seq[File] = Nil
    		
	    	var filenum = 0		
	    	for(seq <- seqs) {

	    		val sai = new File(conf(outBase) + "." + filenum + ".sai")
    			
	    		sais = sais :+ sai
 
	    		filenum += 1
    		
	    		add(bwaAlign(seq, sai, conf(bwaPrefix), conf(bwaThread), conf(bwaThresh)))
	    	}
    	
   			// create read group string
   			val rgStr: String = "@RG\tID:"+conf(rgid)+"\tSM:"+conf(rgsm)+"\tLB:"+conf(rglb)+"\tPL:"+conf(rgpl)+"\tPU:"+conf(rgpu)
   			  		
   			// merge and convert to SAM format
   			add(bwaMakeSam(sais, seqs, conf(bwaPrefix), samFile, rgStr),
   			    //readGroup(samFile, rgFile, conf(rgid), conf(rgsm), conf(rglb), conf(rgpl), conf(rgpu)),
    			//cleanSam(rgFile, cleanFile, conf(picard)),
    			cleanSam(samFile, cleanFile, conf(picard)),
    			samToBam(cleanFile, bamFile, conf(picard)),
    			sortSam(bamFile, sortFile, order),
    			markDup(sortFile, dupFile),
    			bamIndex(dupFile, dupbaiFile, conf(picard)))
	
    	}
    		
    	if(!skip.contains("2")){
    	
    		////	
    		// realignment and base recalibrator
    		val intervalFile = new File(interval)
    		val realignFile = new File(realign)
  		  	val fixFile = new File(fix)
    		val fixbaiFile = new File(fixbai)
    		//val recFile = new File(rec)
    		//val bqsrFile = new File(bqsr)
    		val preFile = new File(pre)
    		val readbaiFile = new File(readbai)
    		val postFile = new File(post)
    		val snpFile = new File(conf(knownsnp))
    		val indelFile = new File(conf(knownindel))
    		
    		add(createTarget(dupFile, refFile, indelFile, intervalFile),
    			indelRealign(dupFile, intervalFile, refFile, indelFile, realignFile),
    			fixMateInfo(realignFile, fixFile, order, conf(picard)),
    			bamIndex(fixFile, fixbaiFile, conf(picard)),
    			baseRecal(fixFile, refFile, snpFile, covs, preFile),
    			printRead(fixFile, preFile, refFile, readFile),
    			bamIndex(readFile, readbaiFile, conf(picard)),
    			baseRecal(readFile, refFile, snpFile, covs, postFile))
    		
    	}
    	
    	// command-line option "skip" contains "3" will skip doing the metrics
    	if(!skip.contains("3")){
    		
    		////
    		// perform assessment of qc metrics
    		val isFile = new File(indexstat)
    		val hsMetricFile = new File(hsmetric)
    		val alignMetricFile = new File(alignmetric)
    		val gcBiasTbl = new File(gcbiastbl)
    		val gcBiasPdf = new File(gcbiaschart)
    		val gcBiasSmry = new File(gcbiassmry)
    		val byCycleTbl = new File(bycycletbl)
    		val byCyclePdf = new File(bycyclechart)
    		val qualDistTbl = new File(qualdisttbl)
    		val qualDistPdf = new File(qualdistchart)
    		val insertTbl = new File(inserttbl)

    		
	    	add(indexStat(readFile, isFile, conf(picard)),
    		    hsMetric(readFile, conf(bait), conf(roi), hsMetricFile, conf(picard)),
    		    alignSmryMetric(readFile, refFile, alignMetricFile, conf(picard)),
    		    gcBiasMetric(readFile, refFile, gcBiasTbl, gcBiasPdf, gcBiasSmry, conf(picard)),
    		    qualByCycle(readFile, byCycleTbl, byCyclePdf, conf(picard)),
    		    qualDist(readFile, qualDistTbl, qualDistPdf, conf(picard)))
    		    
    		// if paired end
    		if(seqs.size == 2) {
    			val insertHist = new File(inserthist)
    			val complexTbl = new File(complextbl)
    			
    			add(insertSizeMetric(readFile, insertTbl, insertHist, conf(picard)),
    			    libComplex(readFile, complexTbl, conf(picard)))    
   			}
    	}
   		
   		
   		if(!skip.contains("4")){

	   		////
	   		// call variants
	    	val indelVcf = new File(indelcall)
	    	val indelMetric = new File(indelmetric)
	    	val snpVcf = new File(snpcall)
	    	val snpMetric = new File(snpmetric)
	    	val bothVcf = new File(bothcall)
	    	val bothMetric = new File(bothmetric)
	    	val allVcf = new File(allcall)
	    	val allMetric = new File(allmetric)
    	  
	    	add(callIndel(readFile, refFile, conf(roi), indelVcf, indelMetric), //call indels
	    		callAtSite(readFile, refFile, conf(roi), conf(dbsnp), "SNP", null, snpVcf, snpMetric), // call snps
	    		callAtSite(readFile, refFile, conf(roi), conf(dbsnp), "BOTH", null, bothVcf, bothMetric), // call indels and snps together
	    		callAtSite(readFile, refFile, conf(roi), conf(dbsnp), "BOTH", "EMIT_ALL_CONFIDENT_SITES", allVcf, allMetric)) // call for all bases in roi
		}
	}
	
    ////
    // helper classes
    ////
    
    class Configure {
    	// parse configuration file and put in a map

		def loadConf(conf: String): Map[String,String] = {
			// read in configuration file to a map
	  
			var confMap: Map[String,String] = Map()
			val delim: String = "="
			val comment: String = "//"

			Source.fromFile(conf).getLines().foreach {   
				line => 
					if(!line.startsWith(comment) && !line.isEmpty) {
						val parPair = line.split(delim)
						confMap += parPair(0)->parPair(1)		
					}
			}
	  
			confMap
		}
	}

	/*	
	class startPoint(config: Map){
		// finds the starting point of the tool based on files specified in the config
		
		def findStep(step: Seq) = {
			// go through file strings defined in config file and determine the first file not declared
			
			while(step.tail){
				//if(definedStep(
			}
		}
		
		def definedStep(stepStr: String) = {
			if(config.contains(stepStr)) true
			else false
		}
		
	}
	*/
	
	case class bwaAlign(inFq: File, outSai: File, refPrefix: String, threadCount: String, qualThresh: String) extends CommandLineFunction {
		// call bwa aln
	  
		//include read group info
  
		@Input(doc="fastq file to align")
		val infile = inFq
		@Output(doc="alignment output filename")
		val outfile = outSai
		
		// queue commandLine definition
		def commandLine =	required("bwa") + required("aln") + 
    						required("-f", outfile) +
    						optional("-t", threadCount) + 
    						optional("-q", qualThresh) +
    						required(refPrefix) +
    						required(infile) 
    	
    	this.isIntermediate = qscript.leaveFile
		this.analysisName = outfile + "-bwa_aln"
		this.jobName = outfile + "-bwa_aln"
    
	}
	
	case class bwaMakeSam(inSai: Seq[File], inFq: Seq[File], refPrefix: String, outSam: File, readGroup: String) extends CommandLineFunction {
		// bwa samse/sampe call to merge alignments and format as SAM file
	  
		@Input(doc="alignment files (max 2)")
		val infile = inSai
		@Output(doc="output SAM filename")
		val outfile = outSam
		
		var cmdStr = new String("bwa")

		if(inSai.length == 1) {
			cmdStr += " samse " + "-r '" + readGroup + "' -f " + outfile + " " + refPrefix + " " + infile(0) + " " + inFq(0) 
		}
		else {
			cmdStr += " sampe " + "-r '" + readGroup + "' -f " + outfile + " " + refPrefix + " " + infile(0) + " " + infile(1) + " " + inFq(0) + " " + inFq(1)
		}
		
		def commandLine = cmdStr
		
		this.isIntermediate = qscript.leaveFile
		this.analysisName = outfile + "-bwa_samse_or_sampe"
		this.jobName = outfile + "-bwa_samse_or_sampe"
	}
	
	case class readGroup(inSam: File, outSam: File, identifier: String, sample: String, library: String, platform: String, unit: String) extends AddOrReplaceReadGroups {
		// picard's AddOrReplaceReadGroups
	  
		@Input(doc="SAM file to add or replace read groups")
		val infile = inSam
		@Output(doc="output SAM file")
		val outfile = outSam
	  
		this.input = Seq(infile)
		this.output = outfile
		
		this.RGID = identifier
		this.RGSM = sample
		this.RGLB = library
		this.RGPL = platform
		this.RGPU = unit
		
		//this.sortOrder = order
		
		this.isIntermediate = qscript.leaveFile
		this.analysisName = outfile + "-picard_readgroup"
		this.jobName = outfile + "-picard_readgroup"
	}
	
	case class cleanSam(inSam: File, outSam: File, pathToPicard: String) extends CommandLineFunction {
		// picard's CleanSam.jar
	  
		@Input(doc="SAM file to clean")
		val infile = inSam
		@Output(doc="cleaned SAM filename")
		val outfile = outSam
		
		def commandLine = required("java") + 
						  required(javaMem) +
						  required("-jar", pathToPicard+"/CleanSam.jar") +
						  required("INPUT="+infile) +
						  required("OUTPUT="+outfile)
	
		this.isIntermediate = qscript.leaveFile
		this.analysisName = outfile + "-picard_cleanSam"
		this.jobName = outfile + "-picard_cleanSam"
	}
	
	case class samToBam(inSam: File, outBam: File, pathToPicard: String) extends CommandLineFunction {
		// convert format from sam to bam, picard SamFormatConverter.jar
	  
		@Input(doc="input SAM file")
		val infile = inSam
		@Output(doc="output BAM file")
		val outfile = outBam
		
		def commandLine = required("java") + 
						  required(javaMem) +
						  required("-jar", pathToPicard+"/SamFormatConverter.jar") +
						  required("INPUT="+infile) +
						  required("OUTPUT="+outfile)
		
		this.isIntermediate = qscript.leaveFile
		this.analysisName = outfile + "-picard_samToBam"
		this.jobName = outfile + "-picard_samToBam"  
	}
	
	case class sortSam(inBam: File, outBam: File, order: SortOrder) extends SortSam {
		// picard SortSam
	  
		@Input(doc="BAM file to sort")
		val infile = inBam
		@Output(doc="sorted BAM file as output")
		val outfile = outBam
		
		this.input = Seq(infile)
		this.output = outfile
		this.sortOrder = order
	  
		this.isIntermediate = qscript.leaveFile
		this analysisName = outfile + "-picard_sortSam"
		this.jobName = outfile + "-picard_sortSam"
	}
	
	case class markDup(inBam: File, outBam: File) extends MarkDuplicates {
		// picard MarkDuplicates
	  
		@Input(doc="input BAM file")
		val infile = inBam
		@Output(doc="output BAM file")
		val outfile = outBam
		
		this.input = Seq(infile)
		this.output = outfile
		this.assumeSorted = true
		
		this.isIntermediate = qscript.leaveFile
		this.analysisName = outfile + "-picard_MarkDuplicates"
		this.jobName = outfile + "-picard_MarkDuplicates"
	}
	
	case class bamIndex(inBam: File, outBai: File, pathToPicard: String) extends CommandLineFunction {
		// picard BuildBamIndex
	  
		@Input(doc="BAM file to index")
		val infile = inBam
		@Output(doc="output BAI filename")
		val outfile = outBai
		
		def commandLine = required("java") + 
						  required(javaMem) +
						  required("-jar", pathToPicard+"/BuildBamIndex.jar") +
						  required("INPUT="+infile) + 
						  required("OUTPUT="+outfile)
	  
		this.isIntermediate = qscript.leaveFile
		this.analysisName = outfile + "-picard_bamIndex"
		this.jobName = outfile + "-picard_bamIndex"
	}
	
	case class createTarget(inBam: File, ref: File, knownindel: File, outInterval: File) extends RealignerTargetCreator {
		// GATK RealignerTargetCreator
		
		@Input(doc="BAM file to find intervals with")
		val infile = inBam // need multiple input bams? 
		@Output(doc="interval file to output")
		val outfile = outInterval
		
		this.input_file = Seq(infile)
		this.reference_sequence = ref
		this.out = outfile
		
		// known sites for both snps and indels if defined
		this.known = Seq(knownindel)
		
		this.isIntermediate = qscript.leaveFile
		this.analysisName = outfile + ".target"
		this.jobName = outfile + ".target"
	}
	
	case class indelRealign(inBam: File, inInterval: File, ref: File, knownindel: File, outBam: File) extends IndelRealigner {
		// GATK IndelRealigner
	  
		// include scatterCount = nContigs
	  
		@Input(doc="BAM file to realign")
		val infile = inBam
		@Input(doc="intervals file")
		val interval = inInterval
		@Input(doc="reference sequence fasta file")
		val reffile = ref
		@Output(doc="BAM file of realigned indels")
		val outfile = outBam
		
		this.input_file = Seq(infile)
		this.targetIntervals = inInterval
		this.reference_sequence = ref
		this.out = outfile
		this.known = Seq(knownindel)

		this.isIntermediate = qscript.leaveFile
		this.analysisName = outfile + ".clean"
		this.jobName =  outfile + ".clean"
	}
	
	case class fixMateInfo(inBam: File, outBam: File, order: SortOrder, pathToPicard: String) extends CommandLineFunction {
		// picard's FixMateInformation
	  
		@Input(doc="BAM file to fix")
		val infile = inBam
		@Output(doc="BAM file to output")
		val outfile = outBam
		
		def commandLine = required("java") + 
						  required(javaMem) +
						  required("-jar", pathToPicard+"/FixMateInformation.jar") +
						  required("INPUT="+infile) + 
						  required("OUTPUT="+outfile) + 
						  required("SORT_ORDER="+order)
				
		this.isIntermediate = qscript.leaveFile
		this.analysisName = outfile + "-picard_fixmateinfo"
		this.jobName = outfile + "-picard_fixmateinfo"
	}
	
	case class baseRecal(inBam: File, ref: File, knownsnp: File, covs: Seq[String], outBqsr: File) extends BaseRecalibrator {
		// GATK BaseRecalibrator
		
		// need anything more than infile and cov?, do i define all the outputs if required for next step?
		@Input(doc="BAM file to adjust")
		val infile = inBam
		@Output(doc="output file of recalibrated reads")
		val outfile = outBqsr
		
		// filenames
		this.input_file :+= infile
		this.reference_sequence = ref
		this.knownSites = Seq(knownsnp)
		this.out = outfile
		
		this.covariate ++= covs
    
		this.disable_indel_quals = true
    
		this.isIntermediate = qscript.leaveFile
		this.analysisName = outfile + "-picard_covariates"
		this.jobName = outfile + "-picard_covariates"
		
	}
	
	case class printRead(inBam: File, inBqsr: File, ref: File, outBam: File) extends PrintReads {
		// GATK PrintReads
		
		@Input(doc="input Bam file")
		val infile = inBam 
		@Input(doc="BQSR base qualities file")
		val bqsrfile = inBqsr
		@Output(doc="output Bam file")
		val outfile = outBam
		
		this.input_file :+= infile
		this.reference_sequence = ref
		//this.BQSR = bqsrfile
		this.out = outfile
		
		this.isIntermediate = qscript.leaveFile
		this.analysisName = outfile + "-gatk_recalibration"
    	this.jobName = outfile + "-gatk_recalibration"
	}
	
	case class indexStat(inBam: File, outTbl: File, pathToPicard: String) extends CommandLineFunction {
		// picard BamIndexStats
		// OUTPUT or > ?
	  
		@Input(doc="Bam file to perform stats on")
		val infile = inBam
		@Output(doc="output table of statistics captured from stdout")
		val outfile = outTbl
		
		def commandLine = required("java") + 
						  required(javaMem) +
						  required("-jar", pathToPicard+"/BamIndexStats.jar") +
						  required("INPUT="+infile) //+
						  //required(">") + 
						  //required("OUTPUT="+outfile)
		
		this.isIntermediate = qscript.leaveFile
		this.analysisName = outfile + "-picard_indexstat"
		this.jobName = outfile + "-picard_indexstat"
	}
	
	
	case class hsMetric(inBam: File, baitBed: File, targetBed: File, outTbl: File, pathToPicard: String) extends CommandLineFunction {
		// picard CalculateHsMetrics
	  
		@Input(doc="input BAM files")
		val infile = inBam
		@Input(doc="bait intervals (BED)")
		val baitfile = baitBed
		@Input(doc="target intervals (BED)")
		val targetfile = targetBed
		@Output(doc="output table file")
		val outfile = outTbl
		  
		def commandLine = required("java") + 
						  required(javaMem) +
						  required("-jar", pathToPicard+"/CalculateHsMetrics.jar") +
						  required("BAIT_INTERVALS="+baitfile) +
						  required("TARGET_INTERVALS="+targetfile) +
						  required("INPUT="+infile) +
						  required("OUTPUT="+outfile)
		
		this.isIntermediate = qscript.leaveFile
		this.analysisName = outfile + "-picard_indexstat"
		this.jobName = outfile + "-picard_indexstat"
	}
	
	case class alignSmryMetric(inBam: File, ref: File, outTbl: File, pathToPicard: String) extends CommandLineFunction {
		// picard CollectAlignmentSummaryMetrics
	  
		@Input(doc="Input BAM file")
		val infile = inBam
		@Input(doc="Reference sequence in which alignments were performed")
		val reffile = ref
		@Output(doc="Output file of summary alignment metrics")
		val outfile = outTbl
		  
		def commandLine = required("java") + 
						  required(javaMem) +
						  required("-jar", pathToPicard+"/CollectAlignmentSummaryMetrics.jar") +
						  required("INPUT="+infile) +
						  required("REFERENCE_SEQUENCE="+reffile) +
						  required("OUTPUT="+outfile)
		
		this.isIntermediate = qscript.leaveFile
		this.analysisName = outfile + "-picard_alignsummary"
		this.jobName = outfile + "-picard_alignsummary"
	}
	
	case class gcBiasMetric(inBam: File, ref: File, outTbl: File, outChart: File, outSummary: File, pathToPicard: String) extends CommandLineFunction {
		// picard CollectGcBiasMetrics
	  
		@Input(doc="BAM file for input")
		val infile = inBam
		@Input(doc="Reference sequence used for alignment")
		val reffile = ref
		@Output(doc="Table file of metrics")
		val outfile = outTbl
		@Output(doc="Chart PDF file")
		val outpdf = outChart
		@Output(doc="Output summary file")
		val summaryfile = outSummary
		
		def commandLine = required("java") + 
						  required(javaMem) +
						  required("-jar", pathToPicard+"/CollectGcBiasMetrics.jar") +
						  required("INPUT="+infile) +
						  required("REFERENCE_SEQUENCE="+reffile) + 
						  required("OUTPUT="+outfile) + 
						  required("CHART="+outpdf) +
						  required("SUMMARY_OUTPUT="+summaryfile)
		
		this.isIntermediate = qscript.leaveFile
		this.analysisName = outfile + "-picard_gcbias"
		this.jobName = outfile + "-picard_gcbias"
	}
	
	case class qualByCycle(inBam: File, outTbl: File, outChart: File, pathToPicard: String) extends CommandLineFunction {
		// picard MeanQualityByCycle
	  
		@Input(doc="Input BAM file")
		val infile = inBam
		@Output(doc="Output table file")
		val outfile = outTbl
		@Output(doc="Chart PDF file")
		val outpdf = outChart
		
		def commandLine = required("java") + 
						  required(javaMem) +
						  required("-jar",pathToPicard+"/MeanQualityByCycle.jar") + 
						  required("INPUT="+infile) + 
						  required("OUTPUT="+outfile) +
						  required("CHART="+outpdf)
		
		this.isIntermediate = qscript.leaveFile
		this.analysisName = outfile + "-picard_meanqualbycycle"
		this.jobName = outfile + "-picard_meanqualbycycle"
	}
	
	case class qualDist(inBam: File, outTbl: File, outChart: File, pathToPicard: String) extends CommandLineFunction {
		// picard QualityScoreDistribution
	  
		@Input(doc="Input BAM file")
		val infile = inBam
		@Output(doc="Output table file")
		val outfile = outTbl
		@Output(doc="Chart PDF file")
		val pdffile = outChart
		
		def commandLine = required("java") + 
						  required(javaMem) +
						  required("-jar",pathToPicard+"/QualityScoreDistribution.jar") +
						  required("INPUT="+infile) +
						  required("OUTPUT="+outfile) + 
						  required("CHART="+pdffile)
		
		this.isIntermediate = qscript.leaveFile
		this.analysisName = outfile + "-picard_qualitydist"
		this.jobName = outfile + "-picard_qualitydist"
	}
	
	case class insertSizeMetric(inBam: File, outTbl: File, outHistogram: File, pathToPicard: String) extends CommandLineFunction {
		// picard CollectInsertSizeMetrics
		// only for paired end
	  
		@Input(doc="Input BAM file")
		val infile = inBam
		@Output(doc="Output table file")
		val outfile = outTbl
		@Output(doc="Output histogram file")
		val histfile = outHistogram
		
		def commandLine = required("java") + 
						  required(javaMem) +
						  required("-jar",pathToPicard+"/CollectInsertSizeMetrics.jar") + 
						  required("INPUT="+infile) + 
						  required("OUTPUT="+outfile) +
						  required("HISTOGRAM_FILE="+histfile)
						  
		this.analysisName = outfile + "-picard_insertsize"
		this.jobName = outfile + "-picard_insertsize"
	}
	
	case class libComplex(inBam: File, outTbl: File, pathToPicard: String) extends CommandLineFunction {
		// picard EstimateLibraryComplexity
		// paired end
		// input format?
	  
		@Input(doc="Input BAM file")
		val infile = inBam
		@Output(doc="Output table file")
		val outfile = outTbl
		
		def commandLine = required("java") + 
						  required(javaMem) +
						  required("-jar",pathToPicard+"/EstimateLibraryComplexity.jar") + 
						  required("INPUT="+infile) +
						  required("OUTPUT="+outfile)
						  
		this.isIntermediate = qscript.leaveFile
		this.analysisName = outfile + "-picard_libcomplexity"
		this.jobName = outfile + "-picard_libcomplexity"
	}
	
	case class coverDepth(inBam: File, ref: File, interval: File, outTbl: File) extends DepthOfCoverage {
		// GATK DepthOfCoverage
	  
		val outForm: String = "rtable"
		val numDownsample: Int = 5000
		
		@Input(doc="Input BAM file")
		val infile = inBam
		@Input(doc="Reference sequence used for alignment")
		val reffile = ref
		@Input(doc="Input intervals file")
		val bedfile = interval
		@Output(doc="output table file")
		val outfile = outTbl
		
		this.input_file :+= infile
		this.reference_sequence = ref
		this.intervals = Seq(bedfile)
		this.out = outfile
		
		this.outputFormat = outForm
		this.dcov = numDownsample
		
		this.isIntermediate = qscript.leaveFile
		this.analysisName = outfile + "-gatk_depthofcoverage"
		this.jobName = outfile + "-gatk_depthofcoverage"
	}
	
	case class callIndel(inBam: File, ref: File, intervalBed: File, outVcf: File, outMetric: File) extends UnifiedGenotyper {
		// GATK UnifiedGenotyper
		
		val minFrac: Double = 0.1
		val glmodel: Model = Model.INDEL // genotype likelihood model
			
		@Input(doc="Input BAM file")
		val infile = inBam
		@Input(doc="Reference sequence used for alignment")
		val reffile = ref
		@Input(doc="BED file of intervals")
		val bedfile = intervalBed
		@Output(doc="Output VCF file of called indels")
		val outfile = outVcf
		@Output(doc="Output file of metrics used to call variants")
		val metricfile = outMetric
		
		this.input_file = Seq(infile)
		this.reference_sequence = reffile
		this.out = outfile
		this.metrics = metricfile
		
		if(bedfile.exists) this.intervals = Seq(bedfile)
		
		this.genotype_likelihoods_model = glmodel
		this.min_indel_fraction_per_sample = minIndelFrac
		
		// always provide the output of this class, output vcf
		this.isIntermediate = false // qscript.leaveFile 
		
		this.analysisName = outfile + "-gatk_callindel"
		this.jobName = outfile + "-gatk_callindel"
	}
	
	
	case class callAtSite(inBam: File, ref: File, intervalBed: File, dbSnp: File, glmodel: String, outMode: String, outVcf: File, outMetric: File) extends UnifiedGenotyper {
		// GATK UnifiedGenotyper
		
		val anno: Seq[String] = Seq("AlleleBalance","DepthOfCoverage","SpanningDeletions")
		val callThresh: Int = 30 // standard_min_confidence_threshold_for_calling
		val emitThresh: Int = 30 // standard_min_confidence_threshold_for_emitting
		val downSample: Int = 1000
		val minQual: Int = 10 // min base quality score
		
		@Input(doc="Input BAM file")
		val infile = inBam
		@Input(doc="Reference sequence used for alignments")
		val reffile = ref
		@Input(doc="BED file of intervals")
		val bedfile = intervalBed
		@Input(doc="dbSNP VCF containing rs ids")
		val snpfile = dbSnp
		@Output(doc="Output VCF file of called variants")
		val outfile = outVcf
		@Output(doc="Output file of metrics used to call variants")
		val metricfile = outMetric
		
		this.input_file = Seq(infile)
		this.reference_sequence = reffile
		this.dbsnp = snpfile
		this.out = outfile
		this.metrics = metricfile
		
		if(bedfile.exists) this.intervals = Seq(bedfile)
		
		this.annotation = anno
		this.downsample_to_coverage = downSample
		this.min_base_quality_score = minQual
		this.standard_min_confidence_threshold_for_calling = callThresh
		this.standard_min_confidence_threshold_for_emitting = emitThresh
		
		// better way to do this?
		if(outMode == "EMIT_VARIANTS_ONLY") this.output_mode = OUTPUT_MODE.EMIT_VARIANTS_ONLY
		else if(outMode == "EMIT_ALL_CONFIDENT_SITES") this.output_mode = OUTPUT_MODE.EMIT_ALL_CONFIDENT_SITES
		else if(outMode == "EMIT_ALL_SITES") this.output_mode = OUTPUT_MODE.EMIT_ALL_SITES
		
		// GLM, is there a better way to do this?
		if(glm == "INDEL") this.genotype_likelihoods_model = Model.INDEL
		else if(glmodel == "SNP") this.genotype_likelihoods_model = Model.SNP
		else if(glmodel == "GeneralPloidySNP") this.genotype_likelihoods_model = Model.GeneralPloidySNP
		else if(glmodel == "GeneralPloidyINDEL") this.genotype_likelihoods_model = Model.GeneralPloidyINDEL
		else this.genotype_likelihoods_model = Model.BOTH // defaults to both
		
		//this.isIntermediate = qscript.leaveFile
		this.analysisName = outfile + "-gatk_callatsite"
		this.jobName = outfile + "-gatk_callatsite"
	}
}
