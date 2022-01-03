% estimates correlation parameters.  reused in point estimates and
% bootstraps later.  make sure this is right before running the bootstraps!


restoredefaultpath

% add my tools to path
mytools = genpath(['../Matlab_Tools/']);
path(mytools,path);
clear mytools;

clear

%% ==== USER CHECK FOR NU ==== 
% note: hand-selected corr parameters are here.  nu from here will be saved
% for both point estimates and bootstraps later
 
% import data
temp = load('../data/cz-panel.mat');
r = table2array(temp.insamp(:,2:end)); 

corremp = corr(r, 'rows', 'pairwise');
% corremp(isnan(corremp)) = nan;
corremplong = unroll_tril(corremp);
corremplong = corremplong(~isnan(corremplong));

% hand-selected model
Nportall = 1000;

muguess = 0.1;
sigmaguess = 0.5;
nu = 0.01;
corrmod = random_corr(Nportall,muguess,sigmaguess,nu,1);
corrmodlong = unroll_tril(corrmod);

% compare
edge = -1:0.05:1;
clf; hold on;
histogram(corremplong,edge,'normalization','probability')
histogram(corrmodlong,edge,'normalization','probability')
% histogram(unroll_tril(corrhatnew),edge,'normalization','probability')
legend('data','model')


%% ==== ESTIMATE MU AND SIGMA ====
% 2021 12 Andrew: the quasi-newton often got stuck near the initial guess.
% It seems fminbnd + grid is necessary
% takes about 10 minutes

tic

seed = 1;
mulist = 0:0.02:0.2;

% optimize for each mui in list
clear sigmahatlist errlist
for mui = 1:length(mulist)
    disp(['solving for mu = ' num2str(mulist(mui))])

    minme = @(sigma) obj_corr(mulist(mui),sigma,nu,Nportall,corremp,seed);
    options = optimset('Display','iter');
    [sigmahatlist(mui),errlist(mui)] = fminbnd(minme,0.01,2.0, options);
    
end

% choose best from list
[errhat,muihat] = min(errlist);
mu = mulist(muihat);
sigma = sigmahatlist(muihat);

toc


%%  ==== SAVE TO DISK

clear corrsave
corrfit = random_corr(Nportall,mu,sigma,nu,seed);
corrsave.mat = corrfit;
corrsave.mu = mu;
corrsave.sigma = sigma;
corrsave.nu = nu;

corrsave.mulist = mulist;
corrsave.sigmahatlist = sigmahatlist;
corrsave.errlist = errlist;

corrsave
disp('saving corrsave')
save('../output/corrsave.mat','corrsave')


%% check


corrfit = random_corr(Nportall,mu,sigma,nu,seed);
corrfitlong = unroll_tril(corrfit);


edge = -1:0.05:1;
clf; hold on;
subplot(1,2,1);
plot(mulist, errlist, 'x')
subplot(1,2,2); hold on;

histogram(corremplong,edge,'normalization','probability')
histogram(unroll_tril(corrfit),edge,'normalization','probability')
legend('Data','Model')


[
    mean(corremplong) mean(corrfitlong)
    median(corremplong) median(corrfitlong)    
    var(corremplong) var(corrfitlong)        
    ]





%% === FUNCTIONS 



function err = obj_corr(mu,sigma,nu,Nportall,corremp,seed)
    corrmod = random_corr(Nportall,mu,sigma,nu,seed);
    
    corremplong = unroll_tril(corremp);
    corremplong = corremplong(~isnan(corremplong));
    
    qlist = 0.1:0.1:0.9;
    err = mean(( ...
        quantile(corremplong, qlist) ...
        - quantile(unroll_tril(corrmod), qlist) ...
        ).^2);    

%     err = mean( ...
%         (mean(corremplong) - mean(unroll_tril(corrmod))).^2 ...
%         + (std(corremplong) - std(unroll_tril(corrmod))).^2 ...
%     );

end
