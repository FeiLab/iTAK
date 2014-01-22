iTAK (current version: v1.2 -08/26/11)

Introduction
 
iTAK is a program to identify plant transcription factors (TFs), transcriptional
regulators (TRs) and protein kinases (PKs) from protein or nucleotide sequences 
and then classify individual TFs, TRs and PKs into gene families. Identification
and classification of TFs and TRs are based on the rules (required and forbidden 
pfam protein domains of each gene family) described in Perez-Rodriguez et al 
(2010) (http://nar.oxfordjournals.org/content/38/suppl_1/D822.full). More than 
sixty family of TFs and twenty families of TRs have been characterized in 
plants. A list of these families can be accessed from the PlnTFDB database
(http://plntfdb.bio.uni-potsdam.de/v3.0/).

Plant protein kinases are identified if the sequences have significant hit to 
the protein kinase domain (PF00069) in the Pfam database (http://pfam.sanger.ac.uk/). 
Protein kinase families were adopted from the PlantsP database 
(http://plantsp.genomics.purdue.edu/html/family/class.html) with additions of 
several newly identified families/subfamilies in plants, e.g., WNK like kinase -
with no lysine kinase, Male germ cell-associated kinase (mak). Protein sequences 
of each kinase family were obtained from the PlantsP database and used to build 
Hidden Markov Models (HMMs). The identified plant protein kinases are classified 
into gene families based by comparing their sequences to these HMMs.

System requirement and dependencies

Linux (required)
Perl version 5.10.0 or higher (required). Perl was installed by default on most 
    Linux systems
BioPerl version 1.006 or higher (required). Please check http://www.bioperl.org 
    and wiki/Installing_BioPerl for more details on installation of BioPerl.
HMMER3 (required). Provided in iTAK.
2.0 GB free disk space for installation.

Release notes
iTAK v1.4 - 06/26/11

iTAK v1.3 - 06/26/11

iTAK v1.2 - 06/26/11

iTAK v1.1 - 06/03/11

iTAK v1.0 - 10/10/10

Installation

Installation of iTAK is straightforward. First download the latest version of 
iTAK for your system and uncompress the downloaded file. It will generate a 
folder named "iTAK-1.0.x32" on a 32-bit machine or "iTAK-1.0.x64" on a 64-bit 
machine (we call this folder "iTAK home folder"). iTAK home folder includes 
three subfolders, a "bin" folder containing the HMMER3 executable, a "database" 
folder containing the domain database files and a "doc" folder containing the 
program documentation files. The home folder also contains a perl script, 
iTAK.pl, which is the core script to run the whole iTAK pipeline. Next download 
the formatted domain database files from the download page. Uncompress and move 
the database files to the "database" folder.

Running iTAK

Quick Start

1. Put the protein or nucleotide sequence file in FASTA format in iTAK home 
   folder
2. Go to iTAK home folder and run iTAK with the following command (assuming the 
   input file name is input_seq)
   >perl iTAK.pl -i input_seq
3. The program will generate an output folder named input_seq_output which 
   contains all the output files. See below for the description of the output 
   files.


Parameters
-i [String]  Name of the input sequence file in FASTA format (required)
-s [String]  Type of input sequence. 'p' for protein sequences | 'n' for 
             nucleotide sequences. (default = p)
-m [String]  Type of analysis ('t' for TF identification | 'p' for PK 
             identification | 'b' for both) (default = b)
-a [Integer] number of CPUs used for hmmscan (default = 1)
-o [String]  Name of the output directory (default = "input file name" + 
             "_output")

Output files

Six files will be generated in the output directory.
1. input_seq_tf_seq: transcription factor sequences (FASTA format).
2. input_seq_tf_family: transcription factor classificcation. A tab-delimited 
   txt file containing sequences IDs and their corresponding transcription 
   factor families.
3. input_seq_tf_align: A tab-delimited txt file containing parsed hmmscan result 
   of transcription factors.
4. input_seq_pkseq: protein kinase sequences (FASTA format).
5. input_seq_pkcat: protein kinases classification. A tab-delimited txt file 
   containing sequence IDs and their corresponding protein kinase families.
6. input_seq_pkaln: A tab-delimited txt file containing parsed hmmscan result 
   of protein kinases.


Performance

We run iTAK on a single CPU on a 32-bit laptop with Intel Core2 Duo P9400 @ 
2.40GHz and 4GB memory. It took aout 3.5 hours to identify both transcription 
factors and protein kinases from 33,410 Arabidopsis protein sequences (TAIR9 
release). 

Download 

Current version of iTAK is v1.2. It's available for both 32- and 64-bit linux 
systems. iTAK can be downloaded from the ftp server:
ftp://bioinfo.bti.cornell.edu/pub/program/itak/

Contact
 
For questions and suggestions, please contact us at bioinfo@cornell.edu
