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

*Calculate Non Reference Sensitivity (NRS)  and Non Reference Discrepancy (NRD) on a merged callset VCF*::

	  variantEvalGenotypeConcordance.py file.vcf

We evaluate a callset by comparing it to another 'gold standard' comparison  callset. The gold comparison calls can be Sanger  or Affymetrix
array derived genotypes. Or if you are comparing variant discovery methods, the comparison calls can be with a different algorithm
from the evaluation calls. Prior to running this program merge the evaluation and comparison VCF files into a single VCF file using
the GATK program `CombineVariants  http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantutils_CombineVariants.html`_
variantEvalGenotypeConcordance.py expects the same samples were examined in the evaluation and comparison callsets and in the resulting merged 
VCF file, the odd  sample columns are the evaluation genotype and the even sample columns are the comparison genotype. 

The program proceeeds by building a genotype comparison matrix where the rows are the evaluation calls and columns are the comparison calls.
Iterating through the file it will tally the cell counts in the matrix and then print out the NRS and NRD values to file.variantEval.txt, where
file is the prefix to the input VCF file. It will also write out log files in VCF format of the sites that contribute to NRS, NRD, and concordant
genotypes.
