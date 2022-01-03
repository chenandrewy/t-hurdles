function se = se_samp(rmat,meth,NbootQ)
% finds standard errors of sample moments by bootstrap

if nargin < 3
    NbootQ = 100;
end


%% bootstrap standarrd errors on sample quantiles

% initialize (this momfun call just grabs dimensions)
momt = nan(NbootQ,length(tempt));
momr = nan(NbootQ,length(tempr));
momsig = nan(NbootQ,length(tempsig));
momcorr = nan(NbootQ,length(tempcorr));

disp('bootstrapping se of quantiles')
tic
for booti = 1:NbootQ
    
    if mod(booti,20)==0
        disp(['booti = ' int2str(booti)])
    end

% shuffle data for bootstrap if desired
% seems like I need to resample both time and signals
% https://conference.iza.org/conference_files/pada2009/hounkannounon_b5157.pdf
% there is surely room for improvement in this bootstrap, and in the (hacked) use of
% the random seed, but given the theoretical results it seems like
% improving the bootstrap will have little effect (also given that the
% average correlations are around zero anyway) Andrew 2021 12.  

    r = samp0.r;

    % draw sample    
    rng(booti*5512,'CombRecursive');
    
    
    % shuffle (not sure about this right now)
    Tboot = size(r,1);
    Nboot = size(r,2);
    
    imonth = randi(Tboot, [Tboot 1]);
    isignal = randi(Nboot, [Nboot 1]);    
    
    tempr = r(imonth,:);
    tempr = tempr(:,isignal);

    r = tempr;
    
    % summary stats 
    rpub = nanmean(r);
    nobspub = sum(~isnan(r));
    sigpub = nanstd(r)./sqrt(nobspub);
    tpub = rpub./sigpub;    
    
    corrpubmat = corr(r, 'rows', 'pairwise');
    corrpub = unroll_tril(corrpubmat);    
    
    
    [momt(booti,:),momr(booti,:),momsig(booti,:),momcorr(booti,:)] ...
        = momfun(tpub,rpub,sigpub,corrpub,meth);


end
toc



%% finalize
se.t = std(momt);
se.r = std(momr);
se.sig = std(momsig);
se.corr = std(momcorr);


