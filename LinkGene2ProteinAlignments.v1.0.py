##########################################################
### Import Necessary Modules

import argparse		               #provides options at the command line
import sys		               #take command line arguments and uses it in the script
import gzip		               #allows gzipped files to be read
import re		               #allows regular expressions to be used

##########################################################
### Command-line Arguments
parser = argparse.ArgumentParser(description="A script to link NCBI gene IDs to Uniprot proteins that have been processed already")
parser.add_argument("-aln", help = "The location of an alignment file with NCBI proteins aligned to uniprot proteins (modified blast outfmt 6)", default=sys.stdin, required=True)
parser.add_argument("-gene", help = "The location of table file with NCBI gene and protein information", default=sys.stdin, required=True)
args = parser.parse_args()


#########################################################
### Open file (object-oriented programming)
class Variables():
    alignments = {}

class OpenFile():
    def __init__ (self, f, typ):
        """Opens a file (gzipped) accepted"""
        if re.search(".gz$", f):
            self.filename = gzip.open(f, 'rb')
        else:
            self.filename = open(f, 'r') 
        if typ == "gene":
            sys.stderr.write("\nOpened gene table file: {}\n".format(f))
            OpenGene(self.filename)
        elif typ == "aln":
            sys.stderr.write("\nOpened alignment file: {}\n".format(f))
            OpenAln(self.filename)


class OpenGene():
    def __init__ (self,f):
        """Reads a gene table file with gene information (unique format)"""
        for self.line in f:
            try:
                line = line.decode('utf-8')
            except:
                pass        
            self.line = self.line.rstrip('\n')
            self.proteinID = self.line.split("\t")[2]
            if self.proteinID in Variables.alignments:
                print("{}\t{}".format(self.line, Variables.alignments[self.proteinID]))
            else:
                print("{}\t{}\t{}\t{}".format(self.line, "NA", "NA", "NA"))
        f.close()


class OpenAln():
    def __init__ (self, aln):
        """Reads an alignment file output by the blast program outfmt 6 (compressed or uncompressed) modified by adding the percent length of the alignment"""          
        for line in aln:
            try:
                line = line.decode('utf-8')
            except:
                pass        
            line = line.rstrip('\n')
            (self.queryID, self.subjectID, self.percentID, self.length, self.mismatch, self.gapopen, self.queryStart, self.queryEnd, 
                self.subjectStart, self.subjectEnd, self.evalue, self.bitscore, self.queryPerLength, self.subjectPerLength) =  line.split("\t")
            Variables.alignments[self.queryID] = "{}\t{}\t{}".format(self.subjectID, self.percentID, self.queryPerLength) 
        aln.close()
        sys.stderr.write("Finished reading alignment file\n\n")


if __name__ == '__main__':
    Variables()
    open_file = OpenFile(args.aln, "aln")
    open_file = OpenFile(args.gene, "gene")
