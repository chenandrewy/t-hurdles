% 2021 12 correlation chosen to fit data

% same as set_fast
meth.name           = mfilename;
meth.Nportall       = 400;
meth.Nsim           = 50;
meth.thresh3list    = [1/3 0.5 2/3];
meth.boottype       = 'xsection';
meth.seed           = nan;
meth.seed_est       = 1;

% difference
parbase.corrtype    = 'constant';
parbase.rho1        = 0.05; % matches median in my sample
parbase.rho2        = nan;
parbase.rho3        = nan;