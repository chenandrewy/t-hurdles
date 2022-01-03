function [err,sim,simflat,emom] = obj_sim(tarvec,namevec,dat,parbase,corrmat,meth)
% 2019 03: updated to take in transformed parameters to avoid nan

% testing a simulation with correlations 2019 03

% 2021 12 corrmat is now pre-calculated


% fix seed, if asked for
if ~isnan(meth.seed)
    rng(meth.seed,'combRecursive')
end
    
% convert transformed vars to straight vars
parvec = tarv2parv(tarvec,namevec);

% update parambase with parvec and namevec
par = parbase;
for i = 1:length(parvec)
    par.(namevec{i}) = parvec(i);
end

% add-ons
if ~isfield(par,'slope_mu')
    par.slope_mu = 0;
end



if ~isfield(par,'nusig')
    par.nusig = inf;
end
if ~isfield(meth,'sim_agg')
    meth.sim_agg = 'mean';
end


% only find moments if there is data
if isempty(dat)
    findmom = 0;
    [dat.tpub,dat.rpub,dat.sigpub,dat.corrpub] = nanall(1,3);    
else 
    findmom = 1;
end
    


%% setup

% ease of notation
Nsim = meth.Nsim;
Nportall = round(meth.Nportall);


% check for errors
if par.sigsig < 0 || isnan(par.sigsig)
    disp('par.sigsig < 0 or isnan, returning nan')
    [err,sim,simflat] = nanall(1,1);
    return
end

%% simulate all, including unpublished

% get ready for parfor
[mu,sig,r,t,null,pub] = nanall(Nsim,Nportall);
mupub = cell(Nsim,1);
nullpub = cell(Nsim,1);
sigpub = cell(Nsim,1);
rpub = cell(Nsim,1);
tpub = cell(Nsim,1);
corrpub = cell(Nsim,1);
Npub = nanall(Nsim,1);

% initialize simulated moments and ensure they match momfun (this is
% annoyingly complicated, I know)
if  findmom
    [dat.tpub,dat.rpub,dat.sigpub,dat.corrpub] = nanall(1,3);
    
    [tempt,tempr,tempsig,tempcorr] = momfun(dat.tpub,dat.rpub,dat.sigpub,dat.corrpub);

else 
    [tempt,tempr,tempsig,tempcorr] = zerosall(1,1);
    
end

momt = nan(Nsim,length(tempt));
momr = nan(Nsim,length(tempr));
momsig = nan(Nsim,length(tempsig));
momcorr = nan(Nsim,length(tempcorr));

clear temp*;

% sometimes sigsig gets too big and produces inf later. this should head it
% off
if par.sigsig > 7
    err = +inf;
    [sim,simflat] = nanall;
    return
end



%% standard errors and mus

% simulate standard error matrix
%   the eye(Nportall)*eps term ensures psd
if isinf(par.nusig)
    manysig = exp(normrnd(par.musig,par.sigsig,[Nportall Nsim]));  
else
    % added 2019 04 20
    temp = par.sigsig*trnd(par.nusig,[Nportall Nsim]) + par.musig;        
    manysig = exp(temp);
end

% generate true means, fixing modes
% gamrnd uses statistics toolbox
switch par.mutype
    case 'gamma'                   
        manymu = gamrnd(par.mu3,par.mu2,[Nportall Nsim])+ par.mu1;          
    case 't'        
        manymu = par.mu2*trnd(par.mu3,[Nportall Nsim]) + par.mu1;        
end

% add relationship with sig
manymu = manymu + par.slope_mu*manysig;

% force nulls to zero
manynull = rand([Nportall Nsim]) < par.p0;
manymu(manynull) = 0;

%% generate returns and t-stats 
if par.rho1 == 0 && ( strcmp(par.corrtype,'constant_fab') || strcmp(par.corrtype,'constant') )

    % easy version with 0 correlations
    manyr = normrnd(manymu,manysig);

elseif par.sigsig == 0 && strcmp(par.corrtype,'constant_fab')

    %% Fabian's fast constant rho
    % have not full checked this after updating 2019 04
    sigc = sqrt(par.rho1)*exp(par.musig);
    sigep = sqrt(exp(par.musig)^2-sigc^2);

    c = normrnd(0,sigc,[Nportall Nsim]);
    ep = normrnd(0,sigep,[Nportall Nsim]);

    manyr = manymu + c + ep;

else 
    %% do matrix version, no other choice
    manyr = nan(size(manymu));
    manystdnorm = normrnd(0,1,[Nportall Nsim]);
    for simi = 1:Nsim
        % the following 2 lines take all the time
        %   but if we estimate musig and sigsig, there's no way to take this
        %   out    
        tempsig = manysig(:,simi);
        tempmu  = manymu(:,simi);

        Sig = corrmat.*(tempsig*tempsig')+ eye(Nportall)*1e-6;
        [C,err_chol] = chol(Sig);    
    
        if err_chol > 0
            Sig = corrmat.*(tempsig*tempsig') + eye(Nportall)*mean(tempsig.^2)*1e-4;
            [C,err_chol] = chol(Sig);   
            
            if err_chol > 0
                % if still not psd, use nan
                C = nan(size(Sig));
            end
        end

        manyr(:,simi) = tempmu + C'*manystdnorm(:,simi); % use transpose C'!    
    end % for simi

end % if par.rho1 ~= 0

% t-stats
manyt = manyr./manysig;

%% truncation

manytruncp = rand([Nportall Nsim]);

%% find who got published
switch par.threshtype

    case 'step'
        igreat      = manyt>par.thresh2; % keep very significant            
        imarginal   = manyt>par.thresh1 & manyt<par.thresh2;  
        imarginalpub = imarginal & ( manytruncp < (1-par.thresh3) ); % keep 1-par.thresh3 of marginal        
        manypub = igreat | imarginalpub;                       

    case 'logit'
		manypub =1./(1+exp(-par.thresh3.*(manyt-par.thresh1))) > manytruncp;
end    


%% find moments

for simi = 1:Nsim

    %% load stuff
    tempsig = manysig(:,simi);
    tempmu  = manymu(:,simi);
    tempnull    = manynull(:,simi);
    tempr = manyr(:,simi);
    tempt = manyt(:,simi);    
    temppub = manypub(:,simi);
    
    %% examine published
    % make sure these are all rows to match formatting
    % columns are sims, rows are observations
    
    % create pub data
    mupub{simi} = tempmu(temppub)';
    nullpub{simi} = tempnull(temppub)';
    sigpub{simi} = tempsig(temppub)';
    rpub{simi} = tempr(temppub)';
    tpub{simi} = tempt(temppub)';   
 
    if par.rho1 > 0
        corrpub{simi} = unroll_tril(corrmat(temppub,temppub))'; % grab published correlations
    else
        corrpub{simi} = 0;
    end
    
    
    %% find moments if desired
    if findmom
        [momt(simi,:),momr(simi,:),momsig(simi,:),momcorr(simi,:),Npub(simi)] ...
            = momfun(tpub{simi},rpub{simi},sigpub{simi},corrpub{simi});
    end
    
    %% save
    mu(simi,:) = tempmu';
    null(simi,:) = tempnull';
    sig(simi,:) = tempsig';
    r(simi,:) = tempr';
    t(simi,:) = tempt';    
    pub(simi,:) = temppub';
    
end % for simi, all loop

probpub = 100*Npub/meth.Nportall;

% store
sim = packstruct(mu,null,sig,r,t,pub,probpub...
    ,corrmat ...
    ,mupub,nullpub,corrpub ...
    ,rpub,tpub,sigpub ...
    ,momt,momsig,momr,momcorr ...
    ,Npub ...
    );

%% calculate fitting error if you hav what you need

switch meth.sim_agg
    case 'median'        
        aggmomt = median(sim.momt,1);
        aggmomr = median(sim.momr,1);
        aggmomsig = median(sim.momsig,1);
        aggmomcorr = median(sim.momcorr,1);        
    case 'mean'
        aggmomt = mean(sim.momt,1);
        aggmomr = mean(sim.momr,1);
        aggmomsig = mean(sim.momsig,1);
        aggmomcorr = mean(sim.momcorr,1);
end

%% find error 
if findmom
    % if data is available that is
    
    emom.t = meth.wswitch.t*(aggmomt - dat.momt);
    emom.r = meth.wswitch.r*(aggmomr - dat.momr);
    emom.sig = meth.wswitch.sig*(aggmomsig - dat.momsig);
    emom.corr = meth.wswitch.corr*(aggmomcorr - dat.momcorr);
    
    err = sum(meth.wt.*emom.t.^2) ...
        + sum(meth.wr.*emom.r.^2) ...
        + sum(meth.wsig.*emom.sig.^2) ...
        + sum(meth.wcorr.*emom.corr.^2);
    
else 
    
    % otherwise do nan
    err = nan;
    
end % if ~isempty(dat)

if isnan(err)
    err = inf;
end


%% create flat structure for convenience
% flat structure has all fields 1 x something

fn = fieldnames(sim);
for fni = 1:length(fn)    
    temp = sim.(fn{fni});
    
%     if ~strcmp(fn{fni},'corrpub')
    
        if iscell(temp)
            simflat.(fn{fni}) = makerow([temp{:}]);
        else
            simflat.(fn{fni}) = makerow(temp(:));
        end
%     end
    
end
