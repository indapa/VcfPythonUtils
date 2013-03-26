#!/usr/bin/env python
from optparse import OptionParser
from VcfFile import *
from VcfRecord import *
from VcfGenotype import *
import itertools
from common import yield_bedcoordinate
import bx.seq.twobit



def main():

    """ This program adds non-reference positions to a VCF file with variant positions.
        It does this by the following. Given a bed file of non-variant intervals and a 2bit file of the reference genome,
        it retrieves the refernce alllele, and prints out the VCF data line with the ref/ref genotypes. Then 
        it prints a single line from the segregating VCF file, and then start the loop again.
        It assumes that the input vcf is position sorted.
        
        To generate the non-segrgating bed interval file, run the following program from bx-python:
        bed_subtract_basewise.py   reference_genome.bed  segregating.sites.bed
        bed_subtract_basewise.py   ~/software/Pgmsnp/PythonNotebook/simref.1.bed  Simulation1.segregating.bed """

    usage = "usage: %prog [options] file.vcf"
    parser = OptionParser(usage)
    parser.add_option("--bed", type="string", dest="bed", help="bed file with non-variant intervals")
    parser.add_option("--tbf", type="string", dest="tbf", help="2bit file of reference genome", default='/Users/amit/data/MySimulations/Simulation1/Reference/simref.1.2bit')
    
    (options, args)=parser.parse_args()

    try:
        sys.stderr.write("opening twobitfile...\n")
        twobit=bx.seq.twobit.TwoBitFile( open( options.tbf ) )
    except:
        sys.stderr.write("unable to open twobit file!\n")

    segregatingVcf=args[0]
    bedfh=open(options.bed,'r')

    vcfh=open(segregatingVcf,'r')
    vcfobj=VcfFile(segregatingVcf)
    vcfobj.parseMetaAndHeaderLines(vcfh)
    header=vcfobj.returnHeader()
    formatstring="GT"
    print header
    for (chrom,start,end) in yield_bedcoordinate(bedfh):
        start=int(start)
        end=int(end)
        for i in range(start,end):
            begin=i
            end=i+1
            refseq=twobit[chrom][begin:end]
            vrec=VcfRecord(chrom,str(end),'.',refseq,'.','.','.','NS=3')

            vrec.addGenotype( VcfGenotype(formatstring,'0/0') )
            vrec.addGenotype( VcfGenotype(formatstring,'0/0') )
            vrec.addGenotype( VcfGenotype(formatstring,'0/0') )
            print vrec.toStringwithGenotypes()
        vcf_gen=vcfobj.yieldVcfRecordwithGenotypes(vcfh)
        print vcf_gen.next().toStringwithGenotypes()


if __name__ == "__main__":
    main()
