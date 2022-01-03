%% choose settings
restoredefaultpath

% add my tools to path
mytools = genpath(['../Matlab_Tools/']);
path(mytools,path);
clear mytools;
addpath('../data/')

clear; cls; clc

%%

name = 'set_fast';

% load baseline settings
[parbase,meth] = baseline_settings();

% modify 
eval(name)

% modify some more if wanted
meth.seed = 1;
% meth.Nportall = 20;

booti = 2;


[samp,parbase,corrmat,meth] = auto_prep(meth,parbase,booti);

% feedback
datetime('now')
disp(['estimating ' meth.name])

tic
[est,objgrid] = one_estimate(meth,parbase,booti);
sec = toc
hrs = sec/60/60    

disp('estimated parameters are')
est.par
err = est.err