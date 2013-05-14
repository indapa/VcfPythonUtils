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

+------------+------------+-----------+-----------+-----------+
|            |     AA	  |    AB     |	    BB    |  no call  |	
+============+============+===========+===========+===========+
|    AA      |     1      |     2     |     3     |    4      |
+------------+------------+-----------+-----------+-----------|
|    AB      |     5      |     6     |     7     |    8      |
+------------+------------+-----------+-----------+-----------+
|    BB      |     9      |     10    |     11    |   12      |
+------------+------------+-----------+----------+------------+
|  no call   |     13     |     14    |     15    |   16      |
+------------+------------+-----------++----------+-----------+
