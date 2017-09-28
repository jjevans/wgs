# wgs

Implementation of variant caller for Whole Genome Sequencing.</br></br>

**Highlights:**</br>
1.  Multiple implementations, focus to lsf_variant_multifastq.py</br>
2.  Python library for direct binding to LSF job runner.  No system calls made in process.</br>
3.  Scala implemented Broad Cancer Group tool 'Queue' to manage variant calling.</br>
</br>

**Multiple implementations**</br>
The mature application for alignment and variant calling is in the bin directory called lsf_variant_multifastq.py</br></br>

**usage:** lsf_variant_multifastq.py config</br></br>

**arguments:** pair formatted configuration file; see file variant.conf in the conf directory</br>
</br></br>
**Python for LSF**</br>
</br></br>
See library in lib directory called lsf_util.py</br></br>
![lsf_drmaa_ex](https://user-images.githubusercontent.com/803012/30944202-449715cc-a3c4-11e7-918f-da44b87736fb.png)
</br></br>


**GATK in Scala:**</br>
The application alnVariantCaller_v3_2.scala performs alignment with BWA and GATK variant calling using the 'Queue' framework.
