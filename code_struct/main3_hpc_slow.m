% run the slow version of bootstrap to check
% 2021 12 - andrew
% right now the idea is to use main3_hpc_fast to do the main results in
% like 12 hours, and then do the super slow thing to say, hey, it dones't
% matter if you let the computer think for 2 weeks on it.

% note: parametric bootstrap won't work here, since it all is run at the
% same time.

%% choose settings
restoredefaultpath

% add my tools to path
mytools = genpath(['../Matlab_Tools/']);
path(mytools,path);
clear mytools;

dock

clear; cls; clc

% name = 'set_slow';
% startlist = {  '0', '333', '666'};
% endlist   = {'332', '665', '1000'};

% for adding more smoothness to main results

name = 'set_fast';
startlist = {'1000', '1333', '1666'};
endlist   = {'1332', '1665', '2000'};


pathproj = '/cm/chen/t-hurdles/2021-07-PostOpenAP/';


%% check settings
eval(name)

meth


disp('please check the settings')
disp('does it look ok?')
temp = input('type 1 to continue: ');

if temp ~= 1
    error('aborting')    
end

clear meth* 


%% check data processing

temp = load('../data/cz-panel.mat');
r = table2array(temp.insamp(:,2:end)); 
t = nanmean(r)./nanstd(r).*sqrt(sum(~isnan(r)));

clf;
histogram(t)
setplot('','t-stat','n')
min_t = min(t)
length_t = length(t)


disp('please check the data processing')
disp('does it look ok?')
temp = input('type 1 to continue: ');

if temp ~= 1
    error('aborting')    
end


%% check correlation stuff
corrinfo = dir('../output/corrsave.mat');
load('../output/corrsave.mat')

% load empirical data
temp = load('../data/cz-panel.mat');
r = table2array(temp.insamp(:,2:end)); 
corremp = corr(r, 'rows', 'pairwise');

corrfit = corrsave.mat;
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

corrsave

disp('please check the corrsave figure and the corrsave struct')
disp('does it look ok?')
temp = input('type 1 to continue: ');

if temp ~= 1
    error('aborting')    
end

%% bootstrap baseline!!!!

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


