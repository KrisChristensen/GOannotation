##########################################################
### Import Necessary Modules

import argparse                       #provides options at the command line
import sys                       #take command line arguments and uses it in the script
import gzip                       #allows gzipped files to be read
import re                       #allows regular expressions to be used

##########################################################
### Command-line Arguments
parser = argparse.ArgumentParser(description="Filter an alignment file based on criteria described below.  If a protein has multiple alignments to the same protein, the length of the alignment will be the based on the start and end of all alignments.  The percent identity will be the average.")
parser.add_argument("-aln_file", help = "The location of the alignment file in outfmt=6", default=sys.stdin, required=True)
parser.add_argument("-fastas", help = "The location of the subject fasta file", default=sys.stdin, required=True)
parser.add_argument("-fastaq", help = "The location of the query fasta file", default=sys.stdin, required=True)
parser.add_argument("-min_per", help = "The minimum percent id to keep an alignment (default:0)", default=0)
parser.add_argument("-min_aln_per", help = "The minimum alignment length (percent of each genes' length) to keep an alignment (default:0)", default=0)    
args = parser.parse_args()

#########################################################
### Open file (object-oriented programming)

class OpenFile():
    protein_lengths = {}
    best = {}
    def __init__ (self, filename, fileType):
        """Opens the file and either directs it to a line-reader for alignment files or fasta files"""
        if re.search(".gz$", filename):
            self.filename = gzip.open(filename, 'rb')
        else:
            self.filename = open(filename, 'r')             
        if fileType == "aln":
            sys.stderr.write("\n{}\n\n".format("Opened alignment file"))
            self.readLinesAln()
        elif fileType == "fasta":
            self.readLinesFasta()

    def readLinesAln(self):
        """Reads alignment files, but only those which pass minimum thresholds."""
        self.all_alignments = {}
        self.previous = "NA"
        for line in self.filename:
            try:
                line = line.decode('utf-8')
            except:
                pass         
            if re.search("\S",line):
                line = line.rstrip('\n')
                query,subject,per_id,aln_len,mismatch,gap,qstart,qend,sstart,send,evalue,bit = line.split("\t")
                if query in self.all_alignments and subject in self.all_alignments[query]:
                    self.all_alignments[query][subject] += ",{}".format(line)
                elif query in self.all_alignments and subject not in self.all_alignments[query]:
                    self.processAln()
                    self.all_alignments.clear()
                    self.all_alignments[query] = {}
                    self.all_alignments[query][subject] = "{}".format(line)
                elif query not in self.all_alignments and self.previous != "NA":
                    self.processAln()
                    self.all_alignments.clear()
                    self.all_alignments[query] = {}
                    self.all_alignments[query][subject] = "{}".format(line)
                else:
                    self.all_alignments[query] = {}
                    self.all_alignments[query][subject] = "{}".format(line)
                    self.previous = "Found"
        self.processAln()
        self.filename.close()
        for self.query in OpenFile.best:
            print("{}".format(OpenFile.best[self.query]))

    def readLinesFasta(self):
        """Measures the lengths of the sequences in the fasta """
        self.header = "NA"
        self.seq = ""
        self.number_seqs = 0
        self.total_size = 0
        for line in self.filename:
            try:
                line = line.decode('utf-8')
            except:
                pass         
            line = line.rstrip('\n')
            if re.search("^\>", line):
                if self.header == "NA":
                    self.header = line.split(" ")[0][1:]
                else:
                    if not (int(len(self.header)) > 0 and int(len(self.seq)) > 0):
                        sys.stderr.write("Sequence Found not conforming to fasta: {}\n".format(self.header))
                    OpenFile.protein_lengths[self.header] = int(len(self.seq))
                    self.total_size += int(len(self.seq))
                    self.seq = ""
                    self.header = line.split(" ")[0][1:]
                self.number_seqs += 1
            elif re.search("\w", line):
                self.seq += line
        if not (int(len(self.header)) > 0 and int(len(self.seq)) > 0):
            sys.stderr.write("Sequence Found not conforming to fasta: {}\n".format(self.header))
        OpenFile.protein_lengths[self.header] = int(len(self.seq))
        self.total_size += int(len(self.seq))
        self.seq = ""
        sys.stderr.write("Finished reading protein fasta file: Found {} sequence(s)\n".format(self.number_seqs)) 
        sys.stderr.write("                    Total Amino Acid(s): {}\n\n".format(self.total_size))
        self.filename.close()

    def processAln(self):
        """Processes the alignment(s)."""
        #query,subject,per_id,aln_len,mismatch,gap,qstart,qend,sstart,send,evalue,bit = line.split("\t")
        tquery = "NA"
        tsubjt = "NA"
        tpid = "NA"
        aln_len = "NA"
        tmismatch = "NA"
        tgap = "NA"
        qstar = "NA"
        qen = "NA"
        sstar = "NA"
        sen = "NA"
        tevalue = "NA"
        tbit = "NA"
        for q in self.all_alignments:
            for s in self.all_alignments[q]:
                if int(len(self.all_alignments[q][s].split(","))) > 1:
                    tquery = q
                    tsubjt = s
                    tpida = []
                    qpositions = []
                    spositions = []
                    for l in self.all_alignments[q][s].split(","):
                        qu,su,pe,al,mm,ga,qs,qe,ss,se,ev,bi = l.split("\t")
                        tpida.append(float(pe))
                        qpositions.append(int(qs))
                        qpositions.append(int(qe))
                        spositions.append(int(ss))
                        spositions.append(int(se))
                    tpid = float(sum(tpida))/int(len(tpida))
                    qstar = min(qpositions)
                    qen = max(qpositions)
                    sstar = min(spositions)
                    sen = max(spositions)
                else:
                    tquery,tsubjt,tpid,aln_len,tmismatch,tgap,qstar,qen,sstar,sen,tevalue,tbit = self.all_alignments[q][s].split("\t")
        self.totalQueryLenPercent = 100*float(abs(int(qen)-int(qstar)))/int(OpenFile.protein_lengths[tquery])
        self.totalSubjectLenPercent = 100*float(abs(int(sen)-int(sstar)))/int(OpenFile.protein_lengths[tsubjt])
        if (float(tpid) >= float(args.min_per) and 
            float(self.totalQueryLenPercent) >= float(args.min_aln_per) and 
            float(self.totalSubjectLenPercent) >= float(args.min_aln_per)):
            if tquery in OpenFile.best:
                self.previous = OpenFile.best[tquery].split()
                if float(self.previous[12]) < float(self.totalQueryLenPercent):
                    OpenFile.best[tquery] = ("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(tquery,tsubjt,tpid,aln_len,tmismatch,tgap,qstar,qen,sstar,sen,tevalue,tbit,self.totalQueryLenPercent,self.totalSubjectLenPercent))
            else:
                OpenFile.best[tquery] = ("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(tquery,tsubjt,tpid,aln_len,tmismatch,tgap,qstar,qen,sstar,sen,tevalue,tbit,self.totalQueryLenPercent,self.totalSubjectLenPercent))
       
if __name__ == '__main__':
    open_aln = OpenFile(args.fastas, "fasta")
    open_aln = OpenFile(args.fastaq, "fasta")           
    open_aln = OpenFile(args.aln_file, "aln")
