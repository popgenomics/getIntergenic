# getIntergenic  
Produces one fasta file for each noncoding loci (inter genic regions) found in the corresponding gff file.  
Also produces a txt file per fasta file: locus_name, contig_name, first position of the non-coding region (0-based), last position of non-coding region (0-based)  
  
## Compilation  
g++ getIntergenic.cpp -std=c++17 -O3 -o getIntergenic  
  
## Exemple:  
./getIntergenic contig_Hmel201012_ama.fasta Hmel2.gff 10000  
**arg1:** fasta file  
**arg2:** gff file  
**arg2:** minimum size in nucleotides of the intergenic for being printed in a fasta file  
  
The fasta file is an alignement of one contig for different individuals.  
The gff file can be the gff produced for the whole genome.  
  
## Precaution  
The produced **fasta files** are written in 'std::ios::out mode'.  
The produced **txt files** are written in 'std::ios::out mode'.  
Running multiple times getIntergenic with the same fasta input file will overwrite the written fasta output files AND the correspond txt output files  

