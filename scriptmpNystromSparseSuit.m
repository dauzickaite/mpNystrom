%% choose the problem and parameters
clear
% path to chop
addpath '/Users/id1917/Documents/MATLAB/chop-master' 
% path to advanpix
addpath '/Users/id1917/Documents/MATLAB/AdvanpixMCT-4.8.5.14607' 
% path to SuiteSparse problems
addpath '/Users/id1917/Documents/MATLAB/suite_sparse_problems'

mp.Digits(64);

rnd_runs = 10;
solveOn = true; % switch to solve the linear systems
mu_list = 0.5; % shift for linear systems

% choose u_p
doubleOn = true;
singleOn = true;
halfOn = true;

% problems with u_p = half

mpNystromApproxPlots(1,'bcsstm07',[50:50:200,220,221,350,390],rnd_runs,solveOn,...
    doubleOn,singleOn,halfOn,mu_list)

mpNystromApproxPlots(1,'1138_bus',[33,46,300,400,1000],rnd_runs,solveOn,...
    doubleOn,singleOn,halfOn,mu_list)

% problems without u_p = half
halfOn = false;
mpNystromApproxPlots(1,'nos7',[24,25,51,52,100,300,489,490],rnd_runs,solveOn,...
    doubleOn,singleOn,halfOn,mu_list)

% % krr problem
halfOn = true;
mpNystromApproxPlots(0,'ijcnn1',[30,50,100,500,900],rnd_runs,solveOn,...
    doubleOn,singleOn,halfOn,0.01)


