############
Programs
############

==========
Programs
==========

The list of Python programs in VcfPythonUtils

**vcf_mendel_discrete_filtering.py**

*Filter genotypes according to a Mendelian inheritance pattern*::

	vcf_mendel_discrete_filtering.py --model [dominant|recessive] --ped [file.ped] --filter PASS file.vcf


**variantEvalGenotypeConcordance.py**

*Calculate overall Non Reference Sensitivity (NRS)  and Non Reference Discrepancy (NRD) on a merged evaluation/comparison callset VCF*::

	  variantEvalGenotypeConcordance.py file.vcf

We evaluate a callset by comparing it to another 'gold standard' comparison  callset. The gold comparison calls can be Sanger  or Affymetrix
array derived genotypes. Or if you are comparing variant discovery methods, the comparison calls can be with a different algorithm
from the evaluation calls. Prior to running this program merge the evaluation and comparison VCF files into a single VCF file using
the GATK program `CombineVariants  <http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantutils_CombineVariants.html>`_
variantEvalGenotypeConcordance.py expects the same samples were examined in the evaluation and comparison callsets and in the resulting merged 
VCF file, the odd  sample columns are the evaluation genotype and the even sample columns are the comparison genotype. 

The program proceeeds by building a genotype comparison matrix where the rows are the evaluation calls and columns are the comparison calls.
Iterating through the file it will tally the cell counts in the matrix and then print out the NRS and NRD values to file.variantEval.txt, where
file is the prefix to the input VCF file. It will also write out log files in VCF format of the sites that contribute to NRS, NRD, and concordant
genotypes.


**Per_sample_variantEvalGenotypeConcordance.py**


*Calculate the per-sample Non-reference sensitivity and Non-reference discrepancy in a merge evaluation/comparsion callset*::

	   Per_sample_variantEvalGenotypeConcordance.py file.vcf

This program is very similiar to variantEvalGenotypeConcordance.py. This difference is that it calculates the NRS and NRD 
metrics on a per sample basis and prints ouput and log files for each sample, in addition to the overall NRS and NRD
values. 

**vcf_pysam_allele_piluep.py**

*Given a VCF file and a BAM file containing the sample(s) in the VCF, this program will add additional INFO and FORMAT
tags RA and AA indicating the number of observations of the reference and alternate allele based on the REF and ALT columns
of the input VCF. Results are written to STDOUT with a new VCF*::

       vcf_pysam_allele_pileup.py --bam file.bam [options] file.vcf

This program is meant add information about the total number of reference and alt allele observations to a VCF that doesn't 
include it already. Also, this program relies on the `pysam  <http://www.cgat.org/~andreas/documentation/pysam/contents.html>`_
Python samtools interface. At the time of writing, I was using v0.7.3

**vcf_removeSamples.py**

*Given a VCF file remove sample(s) as provided on the command line*::

       vcf_removeSamples.py [-h] [-vcf VCFILE] sample [sample ...]

Specify the vcf file by the -vcf option and then just list the sample(s) you want to remove. Output is written to STDOUT

**vcf_gt-filter.py**

*Filter records based on genotypes*::

	vcf_gt-filter.py [-h] file.vcf -filter "sample1 0/1" -filter "sample1 1/0" -filter "sample2 0/0" [ -filter FILTER ] [ --no-header ]

Records are filtered based the value given to  -filter. Going dataline by dataline, if all the genotypes pass the filters
for each specified sample, the line is printed. Option to supress the printing of the header lines with --no-header





   