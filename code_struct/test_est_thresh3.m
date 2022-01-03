% testing estimating thresh3 smoothly
% 2021 12 Andrew

% seems like quasi newton just doesn't work in multiple dimensions.  If
% you're only estimating thresh3, you can get it to work by tweeking
% FiniteDifferenceStepSize.  But once you add one more parameter to
% estiamte, the algorithm likes to return the initial point, or generally
% gets stuck in apparentlylocal minima...
% .... though it depends on the intiial guess.  With a reasonable initial
% guess, it seems to work ok.  


cd ../code_struct/

restoredefaultpath

% add my tools to path
mytools = genpath(['../Matlab_Tools/']);
path(mytools,path);
clear mytools;
addpath('../data/')
addpath('../code_exhibits/')

clear; cls; clc

name = 'set_fast'

[par,meth] = baseline_settings;
eval(name)

meth

booti = 5;

[samp,par,corrmat,meth] = auto_prep(meth,par,booti);

load('../output/set_fast/boot/boot_est0005.mat')

%%

namesmooth = {'p0', 'mu2'};

meth.optnewt = optimoptions('fminunc',...
    'algorithm','quasi-newton',...
    'maxfunctionevaluations',500,...
    'maxiterations',1000*5, ...
    'optimalitytolerance',1e-3, ...
    'steptolerance',1e-3, ...
    'display','iter-detailed', ...
    'FiniteDifferenceStepSize', 1e-2 ...
        );

parguess = est.par;
% parguess.p0 = 0.0

% set up optimization for parameters besides thresh3 and p0 
parvecguess = nan(size(namesmooth));
for pari = 1:length(namesmooth)
    parvecguess(pari) = parguess.(namesmooth{pari});
end

% transform to real line 
tarvec = parv2tarv(parvecguess,namesmooth);

% optimize
minme = @(tarvec)  obj_sim(tarvec,namesmooth,samp,parguess,corrmat,meth);                

[tarvechat,obj,newtflag,newtoutput] = fminunc(minme,tarvec,meth.optnewt);
% [tarvechat,obj] = fminbnd(minme,0,1, optimset('Display','iter'))
% [tarvechat,obj,newtflag,newtoutput] = fminsearch(minme,tarvec);

% transform back
parvechat = tarv2parv(tarvechat,namesmooth); 

namesmooth
parvechat
parvecguess

obj

parhat = parguess;
for pari = 1:length(namesmooth)
    parhat.(namesmooth{pari}) = parvechat(pari);
end

[~,~,simflat] = obj_sim([],[],[],parhat,corrmat,meth);

clf; hold on;
edge = 0:1:15;
histogram(samp.tpub, edge, 'Normalization','probability')
histogram(simflat.tpub, edge, 'Normalization','probability')

  

%% optimize again

namesmooth = {'mu2','musig'};

parguess =  parhat;
objguess = obj;

meth.optnewt = optimoptions('fminunc',...
    'algorithm','quasi-newton',...
    'maxfunctionevaluations',500,...
    'maxiterations',1000*5, ...
    'optimalitytolerance',1e-3, ...
    'steptolerance',1e-3, ...
    'display','iter-detailed', ...
    'FiniteDifferenceStepSize', sqrt(eps) ...
        );


% set up optimization for parameters besides thresh3 and p0 
parvecguess = nan(size(namesmooth));
for pari = 1:length(namesmooth)
    parvecguess(pari) = parguess.(namesmooth{pari});
end

% transform to real line 
tarvec = parv2tarv(parvecguess,namesmooth);

% optimize
minme = @(tarvec)  obj_sim(tarvec,namesmooth,samp,parguess,corrmat,meth);                

[tarvechat,obj,newtflag,newtoutput] = fminunc(minme,tarvec,meth.optnewt);
% [tarvechat,obj] = fminbnd(minme,0,1, optimset('Display','iter'))
% [tarvechat,obj,newtflag,newtoutput] = fminsearch(minme,tarvec);

% transform back
parvechat = tarv2parv(tarvechat,namesmooth); 

namesmooth
parvechat
parvecguess

obj
objguess