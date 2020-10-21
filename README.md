# whatp53sees
MATLAB code used in "What p53 sees: ATM and ATR activation through crosstalk between DNA damage response pathways."

This repository contains all code used in the above paper that does not use previously published work by other authors. A preprint of this work is available at https://www.biorxiv.org/content/10.1101/2020.08.04.237131v1.article-metrics.


The main function is full_plot. When run as-is, full_plot will generate all figures for the base model; that is, the model with none of the experimental conditions applied.

-----------------

It takes two vectors as optional input. The first vector can contain integers between 1 and 7. These correspond to a list of figures to generate:

1: Fig. 4

2: Fig. 5

3: Fig. 6

4: Fig. S7

5: Fig. S8

6: Fig. S9

7: Fig. S10

The second vector can contain integers between 0 and 7. These correspond to the experimental conditions used in the model:

0: Base model, no changes

1: ATM kinase activity does not activate ATR

2: ATR kinase activity does not activate ATM

3: ATM and ATR do not phosphorylate H2AX

4: Replication stress on MDS activates ATR

5: ssDNA does not break to form DSBs

6: ATM dissociates from end-resected DSBs

7: ATR dissociates from extensively end-resected DSBs

------------------

Here is a short description of the other files:

full_ODEs: ODEs for base model
exp1_ODEs: ODEs for experimental case 1 (ATM kinase activity does not upregulate ATR)
exp4_ODEs: ODEs for experimental case 4 (Multiply damaged sites can activate ATR by stalling replication forks)
exp6_ODEs: ODEs for experimental case 6 (End-resected DSBs cannot activate ATM)
exp7_ODEs: ODEs for experimental case 7 (Extensively end-resected DSBs cannot activate ATR)
