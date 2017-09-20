//package variantQ

// $Id$

/*
	Run Queue for GATK 2.0.  Modelled after python script for variant calling 
	named align_call_variants.py.
	
	Jason Evans, 12032012
*/

import scala.collection.mutable.ListBuffer
import scala.io.Source
import java.io._
	
import org.broadinstitute.sting.queue.util.Logging
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.util.QScriptUtils
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.extensions.picard._
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine.OUTPUT_MODE
import org.broadinstitute.sting.commandline.Hidden

import net.sf.samtools.SAMFileHeader.SortOrder

class variantCaller extends QScript with Logging {
	qscript =>

	// config file locale
	@Input(doc="path to the config file", fullName="config", shortName="cfg", required=true)
	var conffile: File = _
	
	// hidden param for how many intervals (loci) to scatter in unifiedgenotyper
	@Hidden
 	@Input(doc="How many ways to scatter/gather in UnifiedGenotyper", fullName="scattergather", shortName="sg", required=false)
	var ugScatter: Int = 150

	////
    /* global vars */
    
	// leave or remove intermediate files on filesystem
	val rmFile = false

	// heap size for picard
	val javaMem: String = "-Xmx48g"
	
    // assignment of parameter variables to keys defined in configuration file
	val outBase = "base"
    val seqleft = "seqleft"
    val seqright = "seqright"
    val numseq = "numseq"
    val ref = "reference"
    val prefix = "prefix"
    val thread = "thread"
    val thresh = "thresh"
    val rgid = "RGID"
    val rgsm = "RGSM"
    val rglb = "RGLB"
    val rgpl = "RGPL"
    val rgpu = "RGPU"  
    val pic = "picard"  
    val knownsnp = "snp"
    val knownindel = "indel"
    val dbsnp = "dbsnp"
    
    // sort order always coordinate
    val order = SortOrder.coordinate
    
	//////
    //////
	def script() {
		// Queue process
		
		//logger.debug("config value: "+conffile)
	  
		//////
	
		// covariates for baserecalibrator PUT MULTIPLE ENTRIES IN CONFIG WITH A LOOP
		val covs = Seq("ReadGroupCovariate","QualityScoreCovariate","CycleCovariate","ContextCovariate")
	  
		// read and load configuration file containing run parameters
		val configure = new Configure
		val conf: Map[String,String] = configure.loadConf(conffile)
    	
    	// filename of alignment files
    	val mergeBam = conf(outBase) + ".merge.bam"
    	val finalBam = conf(outBase) + ".bam"
    	val finalBai = conf(outBase) + ".bai"
    	
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
    	// intermediate Files
    	val intervalFile = new File(interval) // realigner target creator interval file (BED)
    	val realignFile = new File(realign) // realigned BAM
    	val fixFile = new File(fix) // samtools fixmate BAM
    	val fixbaiFile = new File(fixbai) // fixmate BAM index
  		val bqsrFile = new File(baserecal) // base quality recalibrator table
    	val readFile = new File(read) // print reads BAM
    	val readbaiFile = new File(readbai) // print reads BAM index
    	val snpFile = new File(conf(knownsnp)) // SNP calls VCF
    	val indelFile = new File(conf(knownindel)) // indel calls VCF
    	
    	// sequence fastqs, number of sequences to break fastqs into, and reference sequence file
    	val fastqLeft = new File(conf(seqleft))
    	val fastqRight = new File(conf(seqright))
    	val seqPer = conf(numseq)
    	val refFile = new File(conf(ref)) // reference sequence file    	
    	
    	////
		// PROCESS
    	//// 

		////
		// break left (read 1) fastq into numseq number of sequences per smaller fastqs
		val multiFastq = new breakSeq
		var seqsL: Seq[File] = multiFastq.splitFastq(fastqLeft,seqPer)		
		
		////
		// bwa alignment (left, read 1)
		
		// create a sai filename for each inputted file and put in list
	    var saisL: List[File] = Nil
		var filenum = 0
		
		// left sequences (read 1)
		for(sequence <- seqsL.reverse){
			val sai = new File(conf(outBase) + "." + filenum.toString + ".L.sai")
			
			saisL ::= sai
			
			filenum += 1
			
			add(bwaAlign(sequence, sai, conf(prefix), conf(thread), conf(thresh)))
		}
    	
		////
    	// break right (read 2) fastq into numseq number of sequences per smaller fastqs
    	var seqsR: Seq[File] = multiFastq.splitFastq(fastqRight,seqPer)

		////
		// bwa alignment (right, read 2)
    	// create a sai filename for each inputted file and put in list
    	var saisR: List[File] = Nil
    	filenum = 0
    	
    	// right sequences (read 2)
    	for(sequence <- seqsR.reverse){
			val sai = new File(conf(outBase) + "." + filenum.toString + ".R.sai")
			
			saisR ::= sai
			
			filenum += 1
			
			add(bwaAlign(sequence, sai, conf(prefix), conf(thread), conf(thresh)))		
		}
		
		////
    	// bwa sampe - pair reads
    
    	// create read group string
   		val rgStr: String = "@RG\tID:"+conf(rgid)+"\tSM:"+conf(rgsm)+"\tLB:"+conf(rglb)+"\tPL:"+conf(rgpl)+"\tPU:"+conf(rgpu)
   		
   		var sams: List[String] = Nil
   		for(i <- 0 until saisL.length){
   			val samName = conf(outBase) + "." + i.toString + ".sam"
   			val samFile = new File(samName)
   			
   			add(bwaMakeSam(Seq(saisL(i),saisR(i)), Seq(seqsL(i),seqsR(i)), conf(prefix), samFile, rgStr))
	   		
	   		sams ::= samName
		}    	
    	
    	// clean, sam to bam, sort sam
    	var bams: List[String] = Nil
    	for(sam <- sams){
    		val samFile = new File(sam) // fix as File already made in sampe
    		val cleanFile = new File(sam + ".clean.sam")
    		val cleanBamFile = new File(sam + ".clean.bam")
    		val sortFile = new File(sam + ".sort.bam") 
    	
    		add(cleanSam(samFile,cleanFile,conf(pic)),
    			samtoolsSamToBam(cleanFile,cleanBamFile),
    			sortSam(cleanBamFile,sortFile,order))
    	
    		bams ::= sortFile
    	}
    	
    	// merge bam files to a single bam, picard MarkDuplicates, index bam
    	add(samtoolsMerge(bams,mergeBam),
    		markDup(mergeBam,finalBam),
    		samtoolsIndex(finalBam,finalBai)) // need this or does markDup do it automatically?
    	
		
		// find number of contigs to scatter with realignment, bqsr
		val findScatCnt = new findNumContig
		val nContigs: Int = findScatCnt.getNum(conf(ref)) 
		
    	////	
    	// realignment and base recalibrator
    	////
    	
    	// realigner target creator	
    	val rtc = new RealignerTargetCreator	
    	rtc.input_file = Seq(finalBam)
    	rtc.reference_sequence = refFile
    	rtc.out = intervalFile
    	rtc.known = Seq(indelFile) // can also add snp file if desired
    	rtc.scatterCount = nContigs
    	add(rtc)
    	
    	// indel realignment
    	val real = new IndelRealigner
		real.input_file = Seq(finalBam)
		real.targetIntervals = intervalFile
		real.reference_sequence = refFile
		real.out = realignFile
		real.known = Seq(indelFile) // add snp file if desired
		real.scatterCount = nContigs
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
		bqsr.scatterCount = nContigs
    	add(bqsr)
    	
    	// print reads to make a recalibrated bam
    	val pr = new PrintReads    	
    	pr.input_file :+= fixFile
		pr.reference_sequence = refFile
		pr.BQSR = bqsrFile
		pr.out = readFile
		pr.scatterCount = nContigs
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
    	ug0.scatterCount = ugScatter
    	
    	// call indels
    	ug0.out = indelVcf
    	ug0.metrics = indelMetric
    	ug0.genotype_likelihoods_model = Model.INDEL
    	// ug.intervals = Seq(conf(roi))
    	add(ug0)
    	
    	// call SNPs
    	val ug1 = new UnifiedGenotyper
    	ug1.input_file = Seq(readFile)
    	ug1.reference_sequence = refFile
    	ug1.min_indel_fraction_per_sample = 0.1
    	ug1.scatterCount = ugScatter
    	
    	ug1.dbsnp = conf(dbsnp)
    	ug1.annotation = Seq("AlleleBalance","DepthOfCoverage","SpanningDeletions")
    	ug1.downsample_to_coverage = 1000
    	ug1.min_base_quality_score = 10
    	ug1.standard_min_confidence_threshold_for_calling = 30
    	ug1.standard_min_confidence_threshold_for_emitting = 30
    	
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
	
	class breakSeq {
		// take a large fastq file and break it up into multiple smaller fastq 
		// of inputted number of seqs per file

		def splitFastq(fastq: String, num: String): Seq[File] = {
			var files: List[String] = Nil
    		
    		val numline = 4 // number of lines per sequence
			var linenum = 1
			var numfile = 0
			val entry = new ListBuffer[String]
						
		    for(line <- Source.fromFile(fastq).getLines){
				// getLines is lazy load
				
				entry.append(line)
				
				if(linenum/numline == num.toInt){
					val filename = fastq+"."+numfile.toString+".fastq"
				
					this.writeFastq(entry,filename)
				
					files ::= new File(filename)
				
					numfile = numfile + 1
					linenum = 0
					entry.clear
				}

				linenum = linenum + 1
			}
		
			if(entry.length != 0){
				// print out any remaining lines
				
				val filename = fastq+"."+numfile.toString+".fastq"
				
				this.writeFastq(entry,filename)
				
				files ::= new File(filename)
			
				entry.clear
			}
			
			files
		}
		
		def writeFastq(entry: ListBuffer[String], filename: String){
			// use a buffered writer to print lines in a ListBuffer to a file
			
			val bufferLen: Int = 10000000 // buffer this many chars at a time
			
			val out = new java.io.BufferedWriter(new java.io.FileWriter(filename),bufferLen)
			
			for(line <- entry){
				out.write(line+"\n")
			}
			
			out.close
		}
	}

	class findNumContig {
		// reads the reference file and determines the number of contigs defined in the fasta
		
		def getNum(fasta: String): Int = {
			var num = 0
			
			for(line <- Source.fromFile(fasta).getLines){
				if(isContig(line.take(1))){
					num += 1
				}
			}
		
			num
		}
		
		def isContig(line: String): Boolean = line match {
			// match a fasta >
		    case ">" => true
		    case _ => false
		}

	}
	
	case class bwaAlign(inFq: File, outSai: File, refPrefix: String, threadCount: String, qualThresh: String) extends CommandLineFunction {
		// call bwa aln
	  
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
    	
    	this.isIntermediate = qscript.rmFile
		this.analysisName = outfile + "-bwa_aln"
		this.jobName = outfile + "-bwa_aln"
    
	}
	
	case class bwaMakeSam(inSai: Seq[File], inFq: Seq[File], refPrefix: String, outSam: File, readGroup: String) extends CommandLineFunction {
		// bwa samse/sampe call to merge alignments and format as SAM file
	  
		@Input(doc="alignment files (max 2)")
		val infile = inSai
		@Output(doc="output SAM filename")
		val outfile = outSam
		
		def commandLine = required("bwa") + required("sampe") + 
						  required("-r", readGroup) +
						  required("-f", outfile) +
						  required(refPrefix) +
						  required(infile(0)) +
						  required(infile(1)) +
						  required(inFq(0)) +
						  required(inFq(1))

		this.isIntermediate = qscript.rmFile
		this.analysisName = outfile + "-bwa_sampe"
		this.jobName = outfile + "-bwa_sampe"
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
	
		this.isIntermediate = qscript.rmFile
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
		
		this.isIntermediate = qscript.rmFile
		this.analysisName = outfile + "-picard_samToBam"
		this.jobName = outfile + "-picard_samToBam"  
	}
	
	case class samtoolsSamToBam(inSam: File, outBam: File) extends CommandLineFunction {
		// samtools view sam to bam
		// produces uncompressed bam
		
		@Input(doc="input SAM file")
		val infile = inSam
		@Output(doc="output BAM file")
		val outfile = outBam
	
		def commandLine = required("samtools") + 
						  required("view") + 
						  required("-uS") + 
						  required("-o",outfile) +
						  required(infile)
		
		this.isIntermediate = qscript.rmFile
		this.analysisName = outfile + "-samtools_samToBam"
		this.jobName = outfile + "-samtools_samToBam"
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
	  
		this.isIntermediate = qscript.rmFile
		this analysisName = outfile + "-picard_sortSam"
		this.jobName = outfile + "-picard_sortSam"
	}

	case class samtoolsMerge(inBam: Seq[File], outBam: File) extends CommandLineFunction {
		// merge the multiple bams into a single bam
		
		@Input(doc="BAMs to merge")
		val infile = inBam
		@Output(doc="final merged BAM")
		val outfile = outBam
		
		def commandLine = required("samtools") +
						  required("merge") + 
						  required(outfile) +
						  repeat(infile)
						  
		this.isIntermediate = qscript.rmFile
		this.analysisName = outfile + "-merge_bam"
		this.jobName = outfile + "-merge_bam"	
	}
	
	case class markDup(inBam: File, outBam: File) extends MarkDuplicates {
		// picard MarkDuplicates
	  	val metricFile = new File("MarkDuplicates.metrics")
	  	
		@Input(doc="input BAM file")
		val infile = inBam
		@Output(doc="output BAM file")
		val outfile = outBam
		
		this.input = Seq(infile)
		this.output = outfile
		this.metrics = metricFile // pass in a File or keep simple
		this.assumeSorted = true
		this.REMOVE_DUPLICATES = true
		//this.memoryLimit = 16
		
		this.isIntermediate = qscript.rmFile
		this.analysisName = outfile + "-picard_MarkDuplicates"
		this.jobName = outfile + "-picard_MarkDuplicates"
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
						  
		this.isIntermediate = qscript.rmFile
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
						  
		this.isIntermediate = qscript.rmFile
		this.analysisName = outfile + "-samtools_fixmate"
		this.jobName = outfile + "-samtools_fixmate"
	}
}
