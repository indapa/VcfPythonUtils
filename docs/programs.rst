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

*Calculate Non Reference Sensitivity  and Non Reference Discrepancy  on a merged callset VCF*::

	  variantEvalGenotypeConcordance.py [options] file.vcf

+------------+------------+-----------+
| Header 1   | Header 2   | Header 3  |
+============+============+===========+
| body row 1 | column 2   | column 3  |
+------------+------------+-----------+
| body row 2 | Cells may span columns.|
+------------+------------+-----------+
| body row 3 | Cells may  | - Cells   |
+------------+ span rows. | - contain |
| body row 4 |            | - blocks. |
+------------+------------+-----------+