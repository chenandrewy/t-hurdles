% set_fast right now (2021 12) uses Nportall = 400 and Nsim = 50, so let's
% double both.

meth.name           = mfilename;
meth.Nportall       = 800;
meth.Nsim           = 100;
meth.thresh3list    = [1/3 0.5 2/3];
meth.boottype       = 'xsection';
meth.seed           = nan;
meth.seed_est       = 1;

