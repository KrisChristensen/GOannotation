# GOannotation
A generic pipeline to annotate genes with their corresponding gene ontology terms

<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary>Table of Contents</summary>
  <ol>
    <li><a href="#requirements">Requirements</a></li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#license">License</a></li>
  </ol>
</details>

<!-- requirements -->
## Requirements

These scripts have been tested with Python 3.
The script requires the following programs and files.

Programs:<br /><br />
&nbsp;&nbsp;&nbsp;BLAST<br />
&nbsp;&nbsp;&nbsp;Ontologizer<br />
    
Files:<br /><br />
&nbsp;&nbsp;&nbsp;A uniprot protein database (fasta)<br />
&nbsp;&nbsp;&nbsp;A gff file from the organism of interest to assign gene names to protein names<br />
&nbsp;&nbsp;&nbsp;A protein database (fasta) of the organism of interest<br />

<!-- usage -->
## Usage

1) Link protein and gene names for the organism of interest:<br /><br />
&nbsp;&nbsp;&nbsp;python GffGene2Protein.v1.0.py -gff GCF_023373465.1_Oket_V2_genomic.gff.gz > OketV2.Gene2Protein.txt<br /><br />
&nbsp;&nbsp;&nbsp;help (and further explanations): python GffGene2Protein.v1.0.py -h
    
2) Make blast protein database:<br /><br />
&nbsp;&nbsp;&nbsp;makeblastdb -in uniprot_sprot.fasta -dbtype prot<br />
    
3) Align proteins from organism of interest to the uniprot database:<br /><br />
&nbsp;&nbsp;&nbsp;blastp -db uniprot_sprot.fasta -query GCF_023373465.1_Oket_V2_protein.faa -out Oket.vs.uniprot_sprot.aln -max_target_seqs 8 -max_hsps 20 -evalue 0.01 -outfmt 6 -num_threads 12<br />
    
4) Filter the alignments and only keep the best protein:<br /><br />
&nbsp;&nbsp;&nbsp;python FilterAlignmentsBlastFmt6ProteinBest.v1.1.py -aln_file Oket.vs.uniprot_sprot.aln -fastas uniprot_sprot.fasta -fastaq GCF_023373465.1_Oket_V2_protein.faa -min_per 20 -min_aln_per 50 > Oket.vs.uniprot_sprot.lowFiltered.aln<br /><br />
&nbsp;&nbsp;&nbsp;help (and further explanations): python FilterAlignmentsBlastFmt6ProteinBest.v1.1.py -h<br />
    
5) Link proteins and genes:<br /><br />
&nbsp;&nbsp;&nbsp;python LinkGene2ProteinAlignments.v1.0.py -aln Oket.vs.uniprot_sprot.lowFiltered.aln -gene OketV2.Gene2Protein.txt > OketV2.Gene2Protein2Uniprot.lowFilter.txt<br /><br />
&nbsp;&nbsp;&nbsp;help (and further explanations): python LinkGene2ProteinAlignments.v1.0.py -h<br />

6) Use the uniprot mapping software https://www.uniprot.org/id-mapping and downloaded the results (requires certain fields to be compatible with this pipeline: From,Entry,Entry Name,Protein names,Gene Names,Gene Ontology IDs,Gene Names (synonym)).<br /><br />

7) Connect the gene ids to the GO terms for Ontologizer (can also be used to subset regions of the genome)<br /><br />
&nbsp;&nbsp;&nbsp; python LinkNCBI2GOTerms.v1.0.py -go uniprot-compressed_true_download_true_fields_accession_2Cid_2Cprotei-2023.01.24-17.31.15.14.tsv.gz -ncbi OketV2.Gene2Protein2Uniprot.lowFilter.txt > OketV2.Gene2GO.ids<br />
&nbsp;&nbsp;&nbsp;help (and further explanations): python LinkNCBI2GOTerms.v1.0.py -h<br />

8) Download obo file (defines relationships between GO terms): http://purl.obolibrary.org/obo/go.obo<br />

9) Run enrichment analysis:<br /><br />
&nbsp;&nbsp;&nbsp;java -jar Ontologizer.jar -a OketV2.Gene2GO.ids -g go.obo -o folder -s HomeologousRegion.NC_068425.1:26037500-39289791.txt -p OketV2.Gene2GO.ids -m "Benjamini-Hochberg"<br />

<!-- license -->
## License 

Distributed under the MIT License.
