% 2021 12: t alt mu.  better fit for data than normal.  

meth.name           = mfilename;
meth.Nportall       = 400;
meth.Nsim           = 50;
meth.thresh3list    = [1/3 0.5 2/3];
meth.boottype       = 'xsection';
meth.seed           = nan;
meth.seed_est       = 1;

parbase.mutype  	= 't';
parbase.mu3 		= 4; % 4 is estimated by Chen-Zim RAPS
parbase.mu1 		= 0.5;
parbase.mu2 		= 1;


meth.namevec = {'p0','mu1','mu2','musig','sigsig','thresh3'}; 