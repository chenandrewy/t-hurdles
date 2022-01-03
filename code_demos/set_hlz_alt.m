function [parbase,meth] = set_hlz_alt(booti)
% 2019 03: settings to fit alternative model that matches hlz's moments


% this is the baseline setting that I'll use for bootstrapping
% 2019 03 25 Andrew
% this should be relatively fast, it will probably be less precise than the
% other estimates

% use booti <= 0 to skip bootstrap

% made from
% set_test.m 2019 09 Andrew.  Settings file to share accross main files for
% cleaner code


%% Simulation parameters
clear parbase

% parbase.Nportall = 500; % this is chosen in main code
% meth.Nsim = 48*12*2; % 20
meth.Nsim = 48*12*2; % 20


%% baseline parameters


% 'constant' or 'beta' or 'beta_fixed' or 'constant_fab'
% parbase.corrtype = 'constant_fab'; 
% parbase.rhopar = 0.20;

% parbase.corrtype = 'constant_fab'; 
parbase.corrtype = 'constant'; 
parbase.rhopar = 0.2;


parbase.p0 = 0.444;

parbase.mutype = 'gamma';
parbase.mu1 = 0; % mode
parbase.mu2 = 0.555; % dispersion parameter
parbase.mu3 = 1; % tail parameter

parbase.slope_mu = 0;

parbase.threshtype  = 'step';
parbase.thresh1 = 1.96;  % low / midpoint
parbase.thresh2 = 2.57; % good / slope
parbase.thresh3 = 1/2; % frac missing / nan

parbase.sig = (15/sqrt(12)/sqrt(240)); % p26 + p29
