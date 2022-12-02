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
mu_list = [0,0.1,0.5,1];

% choose u_p
doubleOn = true;
singleOn = true;
halfOn = true;

q43on = false;
q52on = false;

% problems with u_p = half
mpNystromApproxPlots(1,'Journals',10:20:110,rnd_runs,solveOn,...
    doubleOn,singleOn,halfOn,q43on,q52on,mu_list)

mpNystromApproxPlots(1,'bcsstm07',[50:50:200,220,221,350,390],rnd_runs,solveOn,...
    doubleOn,singleOn,halfOn,q43on,q52on,mu_list)

mpNystromApproxPlots(1,'plat362',[50:50:200,258,259,350],rnd_runs,false,...
    doubleOn,singleOn,halfOn,q43on,q52on,mu_list)

mpNystromApproxPlots(1,'494_bus',[9,10,20,21,100,200,300],rnd_runs,solveOn,...
    doubleOn,singleOn,halfOn,q43on,q52on,mu_list)

% problems without u_p = half
halfOn = false;

mpNystromApproxPlots(1,'nos7',[24,25,51,52,100,489,490],rnd_runs,solveOn,...
    doubleOn,singleOn,halfOn,q43on,q52on,mu_list)

mpNystromApproxPlots(1,'bcsstk22',10:20:130,rnd_runs,solveOn,...
    doubleOn,singleOn,halfOn,q43on,q52on,mu_list)

mpNystromApproxPlots(1,'LFAT5',[1,3,4,6,7,10,13],rnd_runs,solveOn,...
    doubleOn,singleOn,halfOn,q43on,q52on,mu_list)

% krr problem
halfOn = true;
mpNystromApproxPlots(0,'ijcnn1',[30,50,100,500],rnd_runs,solveOn,...
    doubleOn,singleOn,halfOn,q43on,q52on,[0,0.01,0.1,1])

% quarter precision problems
q52on = true;
mpNystromApproxPlots(1,'Journals',1:9,rnd_runs,solveOn,...
    0,0,0,0,q52on,mu_list)

mpNystromApproxPlots(1,'bcsstm07',1:9,rnd_runs,solveOn,...
    0,0,0,0,q52on,mu_list)

