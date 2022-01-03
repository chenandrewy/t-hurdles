% 2019 04 Andrew

% takes about 20 minutes to load bootstraps
cd '../code_exhibits/'

restoredefaultpath

% add my tools to path
mytools = genpath(['../Matlab_Tools/']);
structfun = genpath(['../code_struct/']);
datafun   = genpath(['../data/']);
path(mytools,path);
path(structfun,path);
path(datafun,path);

clear


name_est = 'set_fast_ok2';
name_boot = name_est;
nbootmax = 2000;
rescale = 1;

%% load data

load('../data/cz-panel.mat')


%% load selected estimate

% these mats are loaded as the filenames (est.mat => est)
load(['../output/' name_est '/est.mat'])
load(['../output/' name_est '/meth.mat'])
est.par.err = est.err; % for symmetry

% rescale (if desired)
load(['../output/' name_est '/rescaled.mat'])
meth2 = meth;
if rescale
    meth2.Nportall = rescaled.Nportall;
    meth2.corrmat  = rescaled.corrmat; 
end

%% load sample, additional stats
[samp,~,corrmat,~] = auto_prep(meth,est.par,0);

% compute additional stats
[sim,htest,shrink,simflat,emom] = addstats(samp,est.par,meth2.corrmat,meth2);


%% load bootstrap
matlist = dir(['../output/' name_boot '/boot/*.mat']);
nboot = min(length(matlist), nbootmax);

bootpar = []; 
for i = 1:nboot
    
    temp = load([matlist(i).folder '/' matlist(i).name]);
    temppar = temp.est.par;   
    temppar.err = temp.est.err;    
    
    bootpar  = [bootpar; struct2table(temppar)];    
    
end
 

%% ==== DISTRIBUTION OF BOOTSTRAPPED PARAM (FOR CHECKING) ====

i = bootpar.thresh3 == 0.5;

clf;
subplot(1,2,1)

tempvar = bootpar.p0;
edge = 0:0.05:1;

    
    center = edge(1:end-1)+mean(diff(edge))/2;
    thresh3list = unique(bootpar.thresh3);
    clear n
    for omegai = 1:length(thresh3list)
        id = bootpar.thresh3 == thresh3list(omegai); 
        n(omegai,:) = histcounts(tempvar(id),edge);
        f = n/length(bootpar.p0);
    end
    
        % draw bars
    h = bar(center,f','stacked');
    setplot('Panel A','\pi_F','Frequency')
    
    
    ylim([0 0.30])
    
        hleg = legend([h] ...
        ,{
        '$\hat{\omega}$ = 1/3'
        '$\hat{\omega}$ = 1/2'
        '$\hat{\omega}$ = 2/3'
        }...
        ,'interpreter','latex');

subplot(1,2,2)
tempvar = bootpar.mu2;
edge = 0.2:0.025:0.5;

    center = edge(1:end-1)+mean(diff(edge))/2;
    thresh3list = unique(bootpar.thresh3);
    clear n
    for omegai = 1:length(thresh3list)
        id = bootpar.thresh3 == thresh3list(omegai); 
        n(omegai,:) = histcounts(tempvar(id),edge);
        f = n/length(bootpar.p0);
    end
    
        % draw bars
    h = bar(center,f','stacked');
    setplot('Panel A','\lambda','Frequency')
    
    ylim([0 0.30])
    
        hleg = legend([h] ...
        ,{
        '$\hat{\omega}$ = 1/3'
        '$\hat{\omega}$ = 1/2'
        '$\hat{\omega}$ = 2/3'
        }...
        ,'interpreter','latex');

quantile(bootpar.p0, [0.5 0.25 0.5 0.75 0.95])



%% ==== TABLE: data moments ====
% assumes meth.wtype = 'varinv'

num = [
    10:10:90
    czmom.est.corr
    -1*czmom.se.corr
    czmom.est.t(2:end)
    -1*czmom.se.t(2:end)
    czmom.est.r
    -1*czmom.se.r;
    czmom.est.sig
    -1*czmom.se.sig;    
    nan(1,9)
    czmom.est.t(1) nan(1,8)
    -1*czmom.se.t(1) nan(1,8)    
    ];

left = {
    'percentile'
    'corr'
    '(se)'
    't'
    '(se)'
    'r'
    '(se)'
    'sig'
    '(se)'    
    ''
    't pct marg'
    '(se)'
};

tab_data = [left num2cell(num)];
tab_data = array2table(tab_data)


writetable(tab_data, '../output_exhibits/tables from matlab.xlsx', 'sheet', 'data')

%% ==== TABLE: PARAMETER ESTIMATES ===
% just write as data, we'll use excel to make a nice table
plist = [5 25 50 75 95];


% process point estimate
point = struct2cell(est.par);
parname = fieldnames(est.par);

% table(parname, point)

% process bootstrap
temp = table2cell(bootpar);
for i = 1:size(temp,1)
    for j = 1:size(temp,2)
        if isstr(temp{i,j})
            temp{i,j} = nan;
        end
    end
end
bootnum = cell2mat(temp);
bootpct = ptile(bootnum, plist);

tab_par = cell2table([parname point num2cell(bootpct')]);


for i = 1:length(plist)
    pliststr{i} = ['pct_' int2str(plist(i))];
end

tab_par.Properties.VariableNames = [
    {'parname' 'point'} pliststr
    ];

tab_par

writetable(tab_par, '../output_exhibits/tables from matlab.xlsx', 'sheet', 'par')

%% === FIGURES: CORRELATION FIT ====



edge = -1:0.05:1;
clf; hold on;
histogram(samp.corrpub,edge,'normalization','probability')
histogram(unroll_tril(meth2.corrmat),edge,'normalization','probability')
legend('Data','Model (incl unpublished)')
setplot('', 'Pairwise Correlation','frequency')

set(gcf,'position', [ 784.3333  372.6667  557.5000  378.3333])

export_fig('../output_exhibits/fit_corr_all.pdf');

%% ==== FIGURE: FIT UNI ====
scale = 1.15;

clf; hold on;
subplot(1,3,1); hold on;
    edges = 0:1:20;
    hdatt = histogram(samp.tpub,edges,'normalization','probability');
    histogram(simflat.tpub,edges,'normalization','probability')
    setplot('Panel A','t-stat','Frequency',scale); addspace
    set(gca,'linewidth', scale);
    legend('Data','Model')
    xlim([-0 18])
    ylim([0 0.30])


subplot(1,3,2); hold on;
    edges = 0:0.2:3;
    histogram(samp.rpub,edges,'normalization','probability')
    histogram(simflat.rpub,edges,'normalization','probability')
    setplot('Panel B','Sample Mean Return','Frequency',scale); addspace
    set(gca,'linewidth', scale);


subplot(1,3,3); hold on;
    edges = 0:0.05:0.6;
    histogram(samp.sigpub,edges,'normalization','probability')
    histogram(simflat.sigpub,edges,'normalization','probability')
    setplot('Panel C','Standard Error','Frequency',scale); addspace
    set(gca,'linewidth', scale);

set(gcf,'position',[ 1094         595           849         308 ])
    
export_fig('../output_exhibits/fit_uni.pdf');


%% ==== FIT: OTHER ====
clf

subplot(1,2,1); hold on;
    edges = -1:0.1:1;
    histogram(samp.corrpub,edges,'normalization','probability')
    histogram(simflat.corrpub,edges,'normalization','probability')
    setplot('Panel D','Pairwise Correlation','Frequency',scale); addspace
    set(gca,'linewidth', scale);
    % legend('Data','Model')
    


subplot(1,2,2); hold on;
    itemp = 1:min(156,length(simflat.sigpub));
    hdat = plot(samp.sigpub,samp.rpub,'o');
    hsim = plot(simflat.sigpub(itemp),simflat.rpub(itemp),'x');
    setplot('Panel E','Standard Error','Sample Mean Return',scale); addspace
    set(gca,'linewidth', scale);

    legend([hdat hsim], {'Data','Model'})

    xlim([0 0.8])


set(gcf,'position',[   1396         534         735         315 ])
    
    % export    
export_fig('../output_exhibits/fit_other.pdf');


