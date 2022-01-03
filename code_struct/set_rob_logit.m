% 2021 12 logistic selection function


parbase.threshtype  = 'logit';
parbase.thresh1 	= 2.0; % midpoint of logit
parbase.thresh2		= nan; % not used
parbase.thresh3  	= 10; % slope of logit, CZ RAPS has around 10



meth.name           = mfilename;
meth.Nportall       = 400;
meth.Nsim           = 50;
meth.thresh3list    = [3 6 9 12]; % for convenience, allow only a grid
meth.boottype       = 'xsection';
meth.seed           = nan;
meth.seed_est       = 1;