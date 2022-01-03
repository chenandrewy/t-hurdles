% 2021 12 hump shaped gamma that actually separates nulls from alts that is
% a minimal deviation from HLZ

meth.name           = mfilename;
meth.Nportall       = 400;
meth.Nsim           = 50;
meth.thresh3list    = [1/3 0.5 2/3];
meth.boottype       = 'xsection';
meth.seed           = nan;
meth.seed_est       = 1;

parbase.mu3 = 2;