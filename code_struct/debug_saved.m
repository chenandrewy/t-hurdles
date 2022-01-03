%% choose settings
restoredefaultpath

% add my tools to path
mytools = genpath(['../Matlab_Tools/']);
path(mytools,path);
clear mytools;
addpath('../data/')

clear; cls; clc

%%

name = 'set_fast_p0spike2';
outpath = ['../output/' name '/'];

load([outpath 'est.mat'])
load([outpath 'meth.mat'])

%%
Nplot = 11;
parvecname = {'p0'};

booti = 4;

plot_obj(meth,est.par,booti,Nplot,parvecname);