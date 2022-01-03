% 2021 12: risk premium 

%% estimate slope
load ../data/cz-panel.mat

stats = regstats(czsamp.rpub, czsamp.sigpub);
parbase.slope_mu    = stats.beta(2)/2;

%%
meth.name           = mfilename;
meth.Nportall       = 400;
meth.Nsim           = 50;
meth.thresh3list    = [1/3 0.5 2/3];
meth.boottype       = 'xsection';
meth.seed           = nan;
meth.seed_est       = 1;

