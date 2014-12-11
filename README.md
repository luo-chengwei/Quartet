Quartet
=======

Quartet: Fully automated Embedded Quartet Decomposition Analysis (EQDA) for horizontal gene transfer (HGT) detection

[Dependencies]
1, networkx
2, Biophython
3, PHYLIP
4, NCBI-BLAST standalone 2.2.24+
5, KaKs_Calculator
6, TREE-PUZZLE
7, ClustalW2

[Usage]
python Quartet.py project.config

[config file]
You need to edit a config file for the project. There is a 'template.config' to help you start with.
In general, you need to add:

indir: that is the directory where you put all the .ffn files of the genomes
outdir: all the results will be saved here
clade1, clade2, ...: followed by the .ffn genome files belonging to the clade
blast: that's the directory where BLAST+ executable/binary could be found
clustalW: that's the directory where the binary of clustalW2 could be found
PHYLIP: that's the directory where the binary of PHYLIP package could be found
puzzle: that's the directory where the TREE-PUZZLE binary could be found
KaKsCal: that's the directory where the KaKs_Calculator could be found
Ks: the cutoff for Ks for orthologs to consider (please refer to methods in Luo et al, PNAS, 2011)
nproc: number of processors available for this task


[misc.]

You can also use this tool to construct core genome and the corresponding genome phylogeny
based on concatenated ortholog alignments. You need to un-comment the '#' signs in line 858
and line 873 for that purpose.
