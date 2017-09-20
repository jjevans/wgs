//package variantQ

// $Id$

/*
	Run Queue for GATK 2.0.  Modelled after python script for variant calling 
	named align_call_variants.py.

	Does not run bwa or groom the alignment files.  It just does realignment and 
	variant calling.  Simplified from version 2 in using the default GATK functions 
	instead of defining specific classes for each walker.
	
	Jason Evans, 12032012
*/

import scala.io.Source
import java.io.File

import org.broadinstitute.sting.queue.util.Logging
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine.OUTPUT_MODE

class variantCaller extends QScript with Logging {
	qscript =>

	// config file locale
	@Input(doc="path to the config file", fullName="config", shortName="cfg", required=true)
	var conffile: File = _
	
	// input bam file name
	@Input(doc="path to bam input file", fullName="bamfile", shortName="bam", required=true)
	var bamfile: File = _

	// input bam index file name
	@Input(doc="path to bam index file", fullName="bamindex", shortName="index", required=true)
	var bamindex: File = _
		
	////
    /* global vars */
    
	// leave or remove intermediate files on filesystem
	val leaveFile = false

    // assignment of parameter variables to keys defined in configuration file
	val outBase = "base"
    val ref = "reference"
    val knownsnp = "snp"
    val knownindel = "indel"
    val dbsnp = "dbsnp"
    
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
    		
    	// intermediate filename with extensions for quality recalibration
    	val interval = conf(outBase) + ".bed"	// createTarget
    	val realign = conf(outBase) + ".real.bam"	// indelRealign	
    	val fix = conf(outBase) + ".fix.bam" // fix mate info, samtools
    	val fixbai = conf(outBase) + ".fix.bai"
    	val baserecal = conf(outBase) + ".bqsr.tbl"
    	val read = conf(outBase) + ".read.bam"
    	val readbai = conf(outBase) + ".read.bai"
    	
    	// intermediate filenames with extensions for calling variants 
    	val indelcall = conf(outBase) + ".indel.vcf"
    	val indelmetric = conf(outBase) + ".indel.tbl"
    	val snpcall = conf(outBase) + ".snp.vcf"
    	val snpmetric = conf(outBase) + ".snp.tbl"
    	
    	////
		// PROCESS
    	//// 

    	////	
    	// realignment and base recalibrator
    	val refFile = new File(conf(ref))
    	val intervalFile = new File(interval)
    	val realignFile = new File(realign)
    	val fixFile = new File(fix)
    	val fixbaiFile = new File(fixbai)
  		val bqsrFile = new File(baserecal)
    	val readFile = new File(read)
    	val readbaiFile = new File(readbai)
    	val snpFile = new File(conf(knownsnp))
    	val indelFile = new File(conf(knownindel))
    		

    	// realigner target creator	
    	val rtc = new RealignerTargetCreator	
    	rtc.input_file = Seq(bamfile)
    	rtc.reference_sequence = refFile
    	rtc.out = intervalFile
    	rtc.known = Seq(indelFile) // can also add snp file if desired
    	add(rtc)
    	
    	// indel realignment
    	val real = new IndelRealigner
		real.input_file = Seq(bamfile)
		real.targetIntervals = intervalFile
		real.reference_sequence = refFile
		real.out = realignFile
		real.known = Seq(indelFile) // add snp file if desired
    	add(real)
    	
    	// samtools fixmateinfo, index recalibrated bam
		add(samtoolsFixMate(realignFile, fixFile),
    		samtoolsIndex(fixFile, fixbaiFile))
		
    	// base quality recalibrator
    	val bqsr = new BaseRecalibrator    	
    	bqsr.input_file :+= fixFile
		bqsr.reference_sequence = refFile
		bqsr.knownSites = Seq(snpFile)
		bqsr.out = bqsrFile
		bqsr.covariate ++= covs
		bqsr.disable_indel_quals = true
    	add(bqsr)
    	
    	// print reads to make a recalibrated bam
    	val pr = new PrintReads    	
    	pr.input_file :+= fixFile
		pr.reference_sequence = refFile
		pr.BQSR = bqsrFile
		pr.out = readFile
		add(pr)
		
		// index recalibrated bam
		add(samtoolsIndex(readFile, readbaiFile))
		
	   	////
	   	// call variants
	    val indelVcf = new File(indelcall)
	    val indelMetric = new File(indelmetric)
	    val snpVcf = new File(snpcall)
	    val snpMetric = new File(snpmetric)

		// Unified Genotyper for variant calls
    	val ug0 = new UnifiedGenotyper
    	ug0.input_file = Seq(readFile)
    	ug0.reference_sequence = refFile
    	ug0.min_indel_fraction_per_sample = 0.1
    	ug0.scatterCount = 12
    	
    	// call indels
    	ug0.out = indelVcf
    	ug0.metrics = indelMetric
    	ug0.genotype_likelihoods_model = Model.INDEL
    	// ug.intervals = Seq(conf(roi))
    	add(ug0)
    	
    	// call at sites
    	val ug1 = new UnifiedGenotyper
    	ug1.input_file = Seq(readFile)
    	ug1.reference_sequence = refFile
    	ug1.min_indel_fraction_per_sample = 0.1
    	ug1.scatterCount = 12
    	
    	ug1.dbsnp = conf(dbsnp)
    	ug1.annotation = Seq("AlleleBalance","DepthOfCoverage","SpanningDeletions")
    	ug1.downsample_to_coverage = 1000
    	ug1.min_base_quality_score = 10
    	ug1.standard_min_confidence_threshold_for_calling = 30
    	ug1.standard_min_confidence_threshold_for_emitting = 30
    	
    	// call snps
    	ug1.out = snpVcf 
    	ug1.metrics = snpMetric // outputted metrics
    	ug1.genotype_likelihoods_model = Model.SNP
    	add(ug1)

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
	
	case class samtoolsIndex(inBam: File, outBai: File) extends CommandLineFunction {
		// samtools index
		
		@Input(doc="BAM file to index")
		val infile = inBam
		@Output(doc="Index file to output")
		val outfile = outBai
		
		def commandLine = required("samtools") + 
						  required("index") + 
						  required(infile) + 
						  required(outfile)
						  
		this.isIntermediate = qscript.leaveFile
		this.analysisName = outfile + "-samtools_fixmate"
		this.jobName = outfile + "-samtools_fixmate"
	}
	
	case class samtoolsFixMate(inBam: File, outBam: File) extends CommandLineFunction {
		// samtools fixmate
		
		@Input(doc="BAM file to fix")
		val infile = inBam
		@Output(doc="BAM file to output")
		val outfile = outBam
		
		def commandLine = required("samtools") + 
						  required("fixmate") + 
						  required(infile) + 
						  required(outfile)
						  
		this.isIntermediate = qscript.leaveFile
		this.analysisName = outfile + "-samtools_fixmate"
		this.jobName = outfile + "-samtools_fixmate"
	}
}
