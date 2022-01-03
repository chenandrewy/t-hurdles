function [samp,parbase,corrmat,meth] = auto_prep(meth,parbase,booti)

% does automatic prep for objective evaluation / estimation 

%% load empirical data
temp = load('../data/cz-panel.mat');

czmom = temp.czmom;
czsamp = temp.czsamp;

%% bootstrap if desired

if booti == 0
    % point estimate
    samp = czsamp;
    
elseif strcmp(meth.boottype, 'parametric')
    % parameteric bootstrap
    disp(['parametric bootstrapping sample booti = ' int2str(booti)])    
    
    rng(booti*5514, 'combRecursive')
    
    temp = load(['../output/' meth.name '/est']);    
    temppar = temp.est.par;
    
    temp = load(['../output/' meth.name '/rescaled']);    
    tempcorrmat = temp.rescaled.corrmat;
    
    tempmeth = meth;
    tempmeth.Nsim = 1;
    tempmeth.Nportall = temp.rescaled.Nportall;
    tempmeth.seed = randi(1e6);    
    
    [~,~,samp] = obj_sim([],[],[],temppar,tempcorrmat,tempmeth);
    
elseif strcmp(meth.boottype, 'xsection')
    % xsectional boot
    disp(['xsectional bootstrapping sample booti = ' int2str(booti)])        
    
    rng(booti*5514, 'combRecursive')
    
    samp0 = czsamp;
    
    Npub = length(samp0.tpub);
    
    iboot = randi(Npub,[Npub 1]);
    
    tpub = samp0.tpub(iboot);
    rpub = samp0.rpub(iboot);
    sigpub = samp0.sigpub(iboot);
    corrpub = samp0.corrpub;
    
    samp = packstruct(tpub,rpub,sigpub,corrpub);
    
    % corrpub is not used
    
    
end


%% Prepare SMM objective

% find smm targets
% in new code (2019 03), rows are simi, cols are within sim
% so these moments should be row vectors

[samp.momt,samp.momr,samp.momsig,samp.momcorr,samp.Npub] = ...
    momfun(samp.tpub,samp.rpub,samp.sigpub,samp.corrpub);

% --- smm weights
switch meth.wtype

    case 'varinv'                
        
        meth.wt = czmom.se.t.^(-2);
        meth.wr = czmom.se.r.^(-2);
        meth.wsig = czmom.se.sig.^(-2);
        meth.wcorr = czmom.se.corr.^(-2);

end




%% Prepare correlation matrix 
% we're not estimating the correlation matrix 

Nportall = meth.Nportall;
switch parbase.corrtype
    case 'constant' 
        if parbase.rho1 ~= 0
            corrmat = parbase.rho1*(ones(Nportall)-eye(Nportall)) + eye(Nportall);    
        else
            corrmat = 0;
        end
    case 'constant_fab'
        corrmat = parbase.rho1*(ones(Nportall)-eye(Nportall)) + eye(Nportall);    

    case 'saved'
        temp = load('../output/corrsave.mat');
        parbase.rho1 = temp.corrsave.mu;
        parbase.rho2 = temp.corrsave.sigma;
        parbase.rho3 = temp.corrsave.nu;        
        
        corrmat = random_corr(Nportall, parbase.rho1, parbase.rho2, parbase.rho3, 1);
        
    otherwise
        error('bad corrtype')
end

