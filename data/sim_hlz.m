function samp0 = sim_hlz()

%% parameters
rng(0)

Nportall = 1400; % BH is invariant to scale.  Not sure about BY though.
Nsim = 1;

rhopar = 0.2;

p0 = 0.444;

mutype = 'gamma';
mu1 = 0; % mode
mu2 = 0.555; % dispersion parameter
mu3 = 1; % tail parameter

slope_mu = 0;

thresh1 = 1.96;  % low / midpoint
thresh2 = 2.57; % good / slope
thresh3 = 1/2; % frac missing / nan

sigpar = (15/sqrt(12)/sqrt(240)); % p26 + p29


%% generate mus
mu = gamrnd(mu3,mu2,[Nsim Nportall]);

% force nulls to zero
null = rand([Nsim Nportall]) < p0;
mu(null) = 0;

%% generate returns

% Fabian's fast constant rho
sigc = sqrt(rhopar)*sigpar;
sigep = sqrt(sigpar^2-sigc^2);
c = normrnd(0,sigc,[Nsim Nportall]);
ep = normrnd(0,sigep,[Nsim Nportall]);

r = mu + c + ep;

% t-stats
sig = ones(size(mu))*sigpar;
t = r./sig;

%% simulate bias
prob_pub = ones(size(t));
prob_pub(t <= thresh1) = 0;
prob_pub(t > thresh1 & t<= thresh2 ) = thresh3;
pub = rand(size(t)) < prob_pub;

rpub = r(pub);
tpub = t(pub);
sigpub = sig(pub);

%% store
samp0.nobspub = 240*ones(size(tpub));
samp0.rpub = rpub;
samp0.volpub = sigpub*sqrt(240);
samp0.sigpub = sigpub;
samp0.tpub = tpub;
samp0.Npub = sum(pub);
samp0.clearpred = ones(size(tpub));




