% run andrew's startup
code_dir = pwd;
cd('/cm/chen/')
startup
undock
cd(code_dir)

% add my tools to path
mytools = genpath(['../Matlab_Tools/']);
path(mytools,path);
clear mytools;
addpath('../data')



% turn on all cores 
% this can run into strange problems for some reason, unfortunately.
% feature('numcores')
% if ans > 4
%     parpool('local',feature('numcores'))
% end