# alignement_transformers
Protein sequence alignement using Transformers

# M2 Bioinformatics - Projet long - Université de paris - Hippolyte MIZINIAK

Link ProTrans (make embedding) : https://github.com/agemagician/ProtTrans

# To use this program, put the data in the same folder (emb : embedding, fasta : fasta files)

# Objectif : Program to perform a pairwise alignment using an embedding produced by a Transformer. This embedding will represent each position of the sequence as a vector of 1024 values describing very finely a position in the context of all the other positions of the sequence from which it comes. 

# Data 

Folder :

- emb :

1a7s.t5emb

1bbr.t5emb

1bina.t5emb

1cg5b.t5emb

1hcga.t5emb

1lh1.t5emb

2fb4l.t5emb

3bjl.t5emb

3hfml.t5emb

- fasta :

1a7s.fasta

1bbr.fasta

1bina.fasta

1hcga.fasta

1lh1.fasta

2fb4l.fasta

3hfm.fasta

# Modules 
- python 3.9
- argparse
- numpy
- matplotlib

# Creating manually conda environment
```
$ conda create env -n alignement_transformers python=3.9
$ conda activate alignement_transformers
$ conda install numpy matplotlib 
```
Command line 
```
python alignement_transformers.py -s1 1bina.fasta -s2 1lh1.fasta -e1 1bina.t5emb -e2 1lh1.t5emb -g -5
```
ARGUMENTS : 

-s1 : fasta file sequence 1

-s2 : fasta file sequence 2

-e1 : embedding sequence 1

-e2 : embdedding sequence 2

-g : gap penality


## Références:
[1] Dzmitry Bahdanau, Kyunghyun Cho, Yoshua Bengio. Neural Machine Translation by Jointly Learning to Align and Translate. arXiv:1409.0473v7 [cs.CL] 2016

[2] Ashish Vaswani, Noam Shazeer, Niki Parmar, Jakob Uszkoreit, Llion Jones, Aidan N. Gomez, Lukasz Kaiser, Illia Polosukhin,. Attention Is All You Need. arXiv:1706.03762v5 [cs.CL] 2017

[3]Ahmed Elnaggar, Michael Heinzinger, Christian Dallago, Ghalia Rehawi, Yu Wang, Llion Jones, Tom Gibbs, Tamas Feher, Christoph Angerer, Martin Steinegger, Debsindhu Bhowmik and Burkhard Rost. ProtTrans: Toward Understanding the Language of Life Through Self-Supervised Learning. IEEE TRANS PATTERN ANALYSIS & MACHINE INTELLIGENCE, VOL. 14, NO. 8,  2021 

[4]  David Mittelman, Ruslan Sadreyev, Nick Grishin. Probabilistic scoring measures for profile–profile comparison yield more accurate short seed alignments. Bioinformatics, Volume 19, Issue 12, 12 August 2003, Pages 1531–1539
