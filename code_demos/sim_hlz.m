function [sim,simflat] = sim_hlz(par,Nsim)




%% setup

% ease of notation
Nportall = round(par.Nportall);

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


%% correlation matrix

corrmat = par.rhopar*(ones(Nportall)-eye(Nportall)) + eye(Nportall);    



%% standard errors and mus
manymu = gamrnd(par.mu3,par.mu2,[Nportall Nsim]);
% force nulls to zero
manynull = rand([Nportall Nsim]) < par.p0;
manymu(manynull) = 0;

%% generate returns and t-stats 


%% Fabian's fast constant rho
sigc = sqrt(par.rhopar)*par.sig;
sigep = sqrt(par.sig^2-sigc^2);
c = normrnd(0,sigc,[Nportall Nsim]);
ep = normrnd(0,sigep,[Nportall Nsim]);

manyr = manymu + c + ep;

% t-stats
manyt = manyr./par.sig;

%% truncation
manypub = manyt > 1.96;

% manypub = abs(manyt) > 1.96; % testing

%% find moments

% parfor simi = 1:Nsim
for simi = 1:Nsim


    %% load stuff
    tempsig = par.sig;
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
    sigpub{simi} = tempsig;
    rpub{simi} = tempr(temppub)';
    tpub{simi} = tempt(temppub)';   
 

    
    %% save
    mu(simi,:) = tempmu';
    null(simi,:) = tempnull';
    sig(simi,:) = tempsig';
    r(simi,:) = tempr';
    t(simi,:) = tempt';    
    pub(simi,:) = temppub';
    
    Npub(simi) = sum(temppub);
    
end % for simi, all loop

probpub = 100*Npub/par.Nportall;

% store
sim = packstruct(mu,null,sig,r,t,pub,probpub...
    ,corrmat ...
    ,mupub,nullpub,corrpub ...
    ,rpub,tpub,sigpub ...
    ,Npub ...
    );



%% create flat structure for convenience
% flat structure has all fields 1 x something

fn = fieldnames(sim);
for fni = 1:length(fn)    
    temp = sim.(fn{fni});
    
    if iscell(temp)
        simflat.(fn{fni}) = [temp{:}]';
    else
        simflat.(fn{fni}) = temp(:)';
    end

    
end
