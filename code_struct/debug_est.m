%% choose settings


restoredefaultpath

% add my tools to path
mytools = genpath(['../Matlab_Tools/']);
path(mytools,path);
clear mytools;
addpath('../data/')
addpath('../code_exhibits/')

clear; cls; clc


booti = 0
%%

    [parbase,meth] = baseline_settings;

    % modify 
    eval('set_fast')    
    
    %% estimate using p0 grid and quasi newton

    % feedback
    datetime('now')
    disp(['estimating ' meth.name])
 
    tic
    [est,objgrid] = one_estimate(meth,parbase,booti);
    
    toc
    
    %%
    
    est.par