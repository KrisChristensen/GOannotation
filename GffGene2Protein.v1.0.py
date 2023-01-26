##########################################################
### Import Necessary Modules

import argparse                       #provides options at the command line
import sys                       #take command line arguments and uses it in the script
import gzip                       #allows gzipped files to be read
import re                       #allows regular expressions to be used

##########################################################
### Command-line Arguments
parser = argparse.ArgumentParser(description="A script to identify genes and corresponding proteins")
parser.add_argument("-gff", help = "The location of the gff file", default=sys.stdin, required=True)
args = parser.parse_args()

#########################################################
### Open file (object-oriented programming)

class OpenFile():
    def __init__ (self, filename, ty, ):
        """Opens the file (accepts gzipped files)"""
        if re.search(".gz$", filename):
            self.file = gzip.open(filename, 'rb')
        else:
            self.file = open(filename, 'r')
        if ty == "gff":
            sys.stderr.write("Reading gff {} file\n".format(filename))
            self.readLinesGFF(self.file)

    def readLinesGFF(self, f):
        """Reads the lines from gff file and outputs them to global variable"""
        self.filename = f
        self.genePosition = {}
        self.rnaName2GeneName = {}
        self.alreadyUsed = {}
        self.geneCount = 0
        self.linkedProteins = 0
        for line in self.filename:
            try:
                line = line.decode('utf-8')
            except:
                pass        
            line = line.rstrip('\n')
            if not re.search("^#", line):
                self.list_of_columns = line.split("\t")
                ###Finds the link between gene id and protein name (it can be tricky as unique gene names may be the same between species)
                if re.search("^gene", self.list_of_columns[2]) and re.search("protein_coding", self.list_of_columns[8]):
                    self.geneCount += 1
                    self.chrom = self.list_of_columns[0]
                    self.start = self.list_of_columns[3]
                    self.end = self.list_of_columns[4]
                    self.gene_name = "NA"
                    self.gene_id = "NA"
                    self.gene_unique = "NA"
                    ###The GFF file does not contain uniform column 8's########
                    ###########################################################
                    for self.info in self.list_of_columns[8].split(";"):
                        if re.search("^ID=gene", self.info):
                            self.gene_unique = self.info[3:]
                        elif re.search("^Dbxref=GeneID:", self.info):
                            self.gene_id = self.info[14:]
                        elif re.search("^Name=", self.info):
                            self.gene_name = self.info[5:]
                    ###########################################################
                    if self.gene_unique in self.genePosition:              
                        pass
                    else:
                        self.genePosition[self.gene_unique] = "{}\t{}\t{}\t{}\t{}".format(self.gene_id, self.gene_name, self.chrom, self.start, self.end)
                ###The mRNA id is needed to link the protein id to the gene id
                elif re.search("^mRNA", self.list_of_columns[2]):
                    self.rna_id = "NA"
                    self.gene_unique = "NA"
                    ###########################################################
                    for self.info in self.list_of_columns[8].split(";"):
                        if re.search("^ID=rna", self.info):
                            self.rna_id = self.info[3:]
                        elif re.search("^Parent=gene", self.info):
                            self.gene_unique = self.info[7:]
                    ###########################################################
                    if self.rna_id in self.rnaName2GeneName:
                        pass
                    else:
                        self.rnaName2GeneName[self.rna_id] = self.gene_unique
                ###The CDS information contains the protein id
                elif re.search("^CDS", self.list_of_columns[2]):
                    self.parent_rna = "NA"
                    self.gene_id = "NA"
                    self.protein_name = "NA"
                    self.gene_name = "NA"
                    self.product_name = "NA"
                    self.gene_unique = "NA"
                    ###########################################################
                    for self.info in self.list_of_columns[8].split(";"):
                        if re.search("^Parent=rna", self.info):
                            self.parent_rna = self.info[7:]
                        elif re.search("^Dbxref=GeneID:", self.info):
                            self.gene_id = self.info.split(",")[0][14:]
                        elif re.search("^Name=", self.info):
                            self.protein_name = self.info[5:]
                        elif re.search("^gene=", self.info):
                            self.gene_name = self.info[5:]
                        elif re.search("^product=", self.info):
                            self.product_name = self.info[8:]                    
                    ###########################################################
                    if self.parent_rna in self.rnaName2GeneName:
                        self.gene_unique = self.rnaName2GeneName[self.parent_rna]
                    if self.gene_unique == "NA":
                        continue
                    self.gi, self.gn, self.chr, self.str, self.en = self.genePosition[self.gene_unique].split("\t")
                    if (self.parent_rna == "NA" or self.gene_id == "NA" or self.protein_name == "NA" or
                        self.gene_name == "NA" or self.product_name == "NA" or self.gene_unique == "NA"):
                        pass
                    elif "{}\t{}\t{}".format(self.gene_name, self.chr, self.str) not in self.alreadyUsed and self.gene_unique in self.genePosition:
                        self.linkedProteins += 1
                        self.alreadyUsed["{}\t{}\t{}".format(self.gene_name, self.chr, self.str)] = "{}\t{}\t{}\t{}\t{}\t{}\t{}".format(self.gene_unique, self.gene_id, self.gene_name, self.product_name, self.chr, self.str, self.en)
                        print ("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(self.gene_id, self.gene_name, self.protein_name, self.product_name, self.chr, self.str, self.en))
        sys.stderr.write("\tFound {} genes in gff file\n\tFound {} linked proteins\n\n".format(self.geneCount, self.linkedProteins))
        self.filename.close() 

if __name__ == '__main__':
    OpenFile(args.gff, "gff")
