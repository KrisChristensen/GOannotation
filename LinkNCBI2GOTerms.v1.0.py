##########################################################
### Import Necessary Modules

import argparse		               #provides options at the command line
import sys		               #take command line arguments and uses it in the script
import gzip		               #allows gzipped files to be read
import re		               #allows regular expressions to be used

##########################################################
### Command-line Arguments
parser = argparse.ArgumentParser(description="A script to link NCBI gene IDs to GO terms after")
parser.add_argument("-go", help = "The location of a file linking Uniprot IDs to GO terms", default=sys.stdin, required=True)
parser.add_argument("-ncbi", help = "The location of a file linking NCBI geneIDs to Uniprot IDs", default=sys.stdin, required=True)
parser.add_argument("-sub", help = "The subset GO terms, default: all, option: chr:start-end", default="all")
args = parser.parse_args()


#########################################################
### Open file (object-oriented programming)
class Variables():
    goTerms = {}

class OpenFile():
    def __init__ (self, f, typ):
        """Opens a file (gzipped) accepted"""
        if re.search(".gz$", f):
            self.filename = gzip.open(f, 'rb')
        else:
            self.filename = open(f, 'r') 
        if typ == "ncbi":
            sys.stderr.write("\nOpened ncbi to uniprot file: {}\n".format(f))
            OpenNCBI(self.filename)
        elif typ == "go":
            sys.stderr.write("\nOpened the uniprot to GO file: {}\n".format(f))
            OpenGo(self.filename)


class OpenNCBI():
    def __init__ (self,f):
        """Reads a gene table file with gene information (unique format)"""
        print("GoStat IDs Format Version 1.0")
        for self.line in f:
            try:
                line = line.decode('utf-8')
            except:
                pass        
            self.line = self.line.rstrip('\n')
            self.gi, self.gene, self.prot, self.desc, self.chr, self.start, self.end, self.uniprot, self.pid, self.len = self.line.split("\t")
            self.uniprot = self.uniprot.split("|")[-1]
            if args.sub == "all":
                #sys.stderr.write("\t{}\n".format(self.uniprot))
                if self.uniprot in Variables.goTerms and re.search("^GO", Variables.goTerms[self.uniprot]):
                    print("{}\t{}".format(self.gene, Variables.goTerms[self.uniprot]))            
            else:
                self.tmpChr = args.sub.split(":")[0]
                self.tmpStart = args.sub.split(":")[1].split("-")[0]
                self.tmpEnd = args.sub.split(":")[1].split("-")[1]
                if self.chr == self.tmpChr and int(self.start) >= int(self.tmpStart) and int(self.start) <= int(self.tmpEnd) and int(self.end) >= int(self.tmpStart) and int(self.end) <= int(self.tmpEnd):
                    if self.uniprot in Variables.goTerms and re.search("^GO", Variables.goTerms[self.uniprot]):
                        print("{}\t{}".format(self.gene, Variables.goTerms[self.uniprot]))
        f.close()


class OpenGo():
    def __init__ (self, aln):
        """Reads a uniprot to GO file"""
        self.header = aln.readline()        
        for line in aln:
            try:
                line = line.decode('utf-8')
            except:
                pass
            line = line.rstrip('\n')
            self.id, self.id2, self.uniprot, self.protein, self.gene, self.go, self.synonym = line.split("\t")
            self.go = ",".join(self.go.split("; "))
            Variables.goTerms[self.uniprot] = "{}".format(self.go)
            #sys.stderr.write("\t{}\t{}\n".format(self.uniprot, self.go))
        aln.close()


if __name__ == '__main__':
    Variables()
    open_file = OpenFile(args.go, "go")
    open_file = OpenFile(args.ncbi, "ncbi")
