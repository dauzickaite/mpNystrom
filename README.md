# mpNystrom
Code for reproducing results in E. Carson and I. Dauzickaite (2023), Single-pass Nystrom approximation in mixed precision, https://arxiv.org/abs/2205.13355 .
MATLAB R2021a is used.

A modified Matlab function pcg is used to prevent exit due to stagnation. Comment out lines 268 - 270 and 279 - 282 in pcg.m and save as pcg_nostag.m.
Copy the functions in the pcg private directory to the current directory.

Run scriptmpNystrom.m for the synthetic examples, and scriptmpNystromSparseSuit.m for the SuiteSparse problems and the kernel ridge regression. scriptNormWeightPseudoinv.m returns plots for the experiments with the weighted pseudoinverse. 

Note: requires chop (https://github.com/higham/chop) and advanpix toolbox (https://www.advanpix.com/). The kernel ridge regression uses function  libsvm2mat.m from https://notendur.hi.is/tpr/tutorials/svm/hugbunadur.html
