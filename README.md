# DNA sequencing by hybridization

Exact and heuristic algorithm to solve the problem of isothermic DNA sequencing by hybridization with negative errors.
Input data form file seq.txt is DNA sequence and first oligonucleotide, from this a biochemical hybridization experiment is made - a spectrum is generated, which consists a set of oligonucleotides, that are a short subsequences of the given DNA fragment. For the experiment the temparature, the percentage of misteakes and repeats can be set. The number of possible repeats here is set to: 0,1,{2,3},{4,5} or more. The aim is to reconstruct the given DNA sequence of a known length on the basis of these generated oligonucleotides. 
