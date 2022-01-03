% use this to run fast settings that take less than, say, 12 hours.
% 2021 12 - Andrew
% right now all main results use this.

% whole thing should take about 12 hours

%% choose settings
restoredefaultpath

cd ../code_struct/

% add my tools to path
mytools = genpath(['../Matlab_Tools/']);
path(mytools,path);
clear mytools;
addpath('../data')

clear; cls; clc;

% name = 'set_fast'; % main results setting, run 2021-12-15 ish
% name = 'set_rob_rho_const'; % robustness: constant corr fit to data
% name = 'set_rob_rho_hlz'; % robustness: hlz's constant corr
% name = 'set_rob_prem'; % robustness: risk premium
% name = 'set_rob_logit'; % robustness: logit selection function
% name = 'set_rob_gamma2'; % robustness: hump shaped gamma alt
% name = 'set_rob_normal'; % robustness: normal alt
name = 'set_rob_t04'; % robustness: t-dist following CZ RAPS

startlist = {  '1', '333', '666'};
endlist   = {'332', '665', '1000'};
pathproj = '/cm/chen/t-hurdles/2021-07-PostOpenAP/';

%% check settings
[parbase,meth] = baseline_settings;
eval(name)

meth

if exist('parbase','var')
    parbase
end


disp('please check the settings')
disp('does it look ok?')
temp = input('type 1 to continue: ');

if temp ~= 1
    error('aborting')    
end

%% check correlation stuff
if ~strcmp(parbase.corrtype, 'saved')  

    corrinfo = dir('../output/corrsave.mat');
    load('../output/corrsave.mat')

    % load empirical data
    temp = load('../data/cz-panel.mat');
    r = table2array(temp.insamp(:,2:end)); 
    corremp = corr(r, 'rows', 'pairwise');

    Nportall = size(corrsave.mat,1);
    corrfit = random_corr(Nportall,corrsave.mu,corrsave.sigma,corrsave.nu,1);
    corrfitlong = unroll_tril(corrfit);


    clf; hold on;
    subplot(1,2,1);
    plot(corrsave.mulist, corrsave.errlist, 'x')
    setplot('','mu','error')

    subplot(1,2,2); hold on;
    edge = -1:0.05:1;
    histogram(unroll_tril(corremp),edge,'normalization','probability')
    histogram(unroll_tril(corrfit),edge,'normalization','probability')
    legend('Data','Model')

    disp('please check the corrsave figure')
    disp('does it look ok?')
    temp = input('type 1 to continue: ');

    if temp ~= 1
        error('aborting')    
    end

end

%% point estimate 
booti = 0;

% estimate
[est, objgrid] = one_estimate(meth,parbase,booti);

%% calculate rescaled model 

[~,~,corrmat,meth] = auto_prep(meth,est.par,booti);

% simulate to find scaling
[~,~,simflat] = obj_sim([],[],[],est.par,corrmat,meth);    
rescaled.pubrate = sum(~isnan(simflat.rpub))/length(simflat.mu);
rescaled.Nportall = round(meth.Nportall/rescaled.pubrate);    
if strcmp(parbase.corrtype, 'saved')
    rescaled.corrmat = random_corr(rescaled.Nportall,est.par.rho1,est.par.rho2,est.par.rho3,1);    
else
    rescaled.corrmat = parbase.rho1*(ones(rescaled.Nportall)-eye(rescaled.Nportall)) + eye(rescaled.Nportall);    
end

%% save

dirname = ['../output/' meth.name '/'];
mkdir(dirname)
save([dirname 'meth'],'meth')
save([dirname 'est'],'est')
save([dirname 'rescaled'],'rescaled')
copyfile([meth.name '.m'] , dirname)    


%% bootstrap fast for testing



for jobi = 1:length(startlist)
    bootistart = startlist{jobi};
    bootiend   = endlist{jobi};    

    % create list of system commands to run
    logname = [pathproj 'output/' name '/hpc' bootistart '-' bootiend];
    runme = {
        ['cd ' pathproj 'code_struct/; ']
        ['mkdir ' pathproj 'output/' name '; ']
        ['sbatch ']
        [' --job-name=boot' bootistart '-' bootiend ]
        [' --error=' logname '.err']
        [' --output=' logname '.out ']
        [' submit-to-hpc.sh ' name ' ' bootistart ' ' bootiend]
        };

    % submit jobs
    system(strjoin(runme(:)))

end



%% ==== CHECK BOOTSTRAP PROGRESS ====
edit(['../output/' name '/hpc1-332.out'])