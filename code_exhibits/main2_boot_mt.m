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

%% calculate bootstrap stats
tic

[t_fdr05, t_fdr01, pubfdr, s_mean, s_smooth, drate05, drate01] = nanall([nboot 1]);
parfor i = 1:nboot
    
    if mod(i,10) == 0
        disp(['processing booti = ' int2str(i)])
    end    

    temppar = table2struct(bootpar(i,:));    
    % compute additional stats 
    [~,temphtest,tempshrink,~] = addstats(samp,temppar,meth2.corrmat,meth2);    
    
    % save parameters
    t_fdr05(i) = temphtest.t_fdr05;
    t_fdr01(i) = temphtest.t_fdr01;
    pubfdr(i) = temphtest.pubfdr;
    s_mean(i) = tempshrink.samp.mean;
    s_smooth(i) = tempshrink.samp.smooth;
    drate05(i) = temphtest.drate05;
    drate01(i) = temphtest.drate01;
    
end
bootstat = table(t_fdr05,  t_fdr01, pubfdr, s_mean, s_smooth, drate05, drate01);

disp('done loading bootstrap')
toc



%% save temporarily
save('../output_exhibits/temp_main_exhibits.mat')

%% load
load('../output_exhibits/temp_main_exhibits.mat')

%% ==== FIGURE: PROP 1 INTUITION ====

clf; hold on;

edge = 0:0.05:0.8;
center = edge(1:end-1)+mean(diff(edge))/2;    

subplot(1,2,1); hold on;

    n = histcounts(bootpar.p0, edge);
    f1 = n/length(bootpar.p0);

    n = histcounts(bootstat.drate05, edge);
    f2 = n/length(bootpar.p0);

    hb = bar(center, [f1' f2'], 1.5);

    setplot('Panel A', 'Estimate', 'Frequency')
    hleg = legend( {'$\pi_F$','$Pr(|t|>1.96)$'} ...
        ,'interpreter','latex');

subplot(1,2,2);
    edge = -0.5:0.05:0.7;
    center = edge(1:end-1)+mean(diff(edge))/2;    

    % plot
    n = histcounts(bootpar.p0-bootstat.drate05, edge);
    f1 = n/length(bootpar.p0);
    hb = bar(center, f1', 0.8);
    
    % annotate
    ylim([0 0.16])
    hv = vline(0,'-');  hv.Color = [1 1 1]*0.4;
    text( 0.02,0.15,'raise hurdle ->','color','k','fontsize',10)
    text(-0.43,0.15,'<- lower hurdle','color','k','fontsize',10)

    setplot('Panel B','$\hat{\pi}_F - \widehat{Pr}(|t|>1.96)$','Frequency')

    ha = gca;
    ha.XLabel.Interpreter = 'latex';

set(gcf, 'position', [ 791.8333  303.5000  651.6667  436.6667])

export_fig('../output_exhibits/prop-intuition.pdf')


%% ==== FIGURE: BOOTSTRAPPED T-HURDLES ====
% = user =
% choose how to break up the chart
grey = [1 1 1]*0.5;
breakvar = bootpar.p0;
tempbreak = [0.1 0.4];

% histogram breakpoints
edge = [0:0.25:3.5];

breakid = ones(size(bootpar.p0));
for i = 1:length(tempbreak)
    breakid(breakvar > tempbreak(i)) = i+1;
end


clf       
    
subplot(2,1,1)
    % settings
    tempvar = bootstat.t_fdr05;
    
    
    % prep bars
    center = edge(1:end-1)+mean(diff(edge))/2;    
    clear n
    for breaki = 1:3
        n(breaki,:) = histcounts(tempvar(breakid==breaki),edge);
        f = n/length(bootpar.p0);
    end
    
    % draw bars
    h = bar(center,f','stacked');
    setplot('Panel A','t-hurdle for FDR = 5%','Frequency')
    
    ylim([0 0.35])
    
    % traditional t-hurdle
    hv = vline(-1*norminv(0.05/2),'k--');
    text(2,0.3,'classical t-hurdle','color','k','fontsize',10)
    
    
    hleg = legend([h] ...
        ,{
        ['$\hat{\pi}_F <$  ' num2str(tempbreak(1)) ]
        ['$\hat{\pi}_F \in ($' num2str(tempbreak(1)) ',' num2str(tempbreak(2)) ']']
        ['$\hat{\pi}_F >$  ' num2str(tempbreak(2)) ]
        }...
        ,'interpreter','latex');
    hleg.FontSize = 12;
    hleg.Position = [ 0.2657    0.7652    0.1984    0.1055];
    
    
    addspace(0.05)
    
subplot(2,1,2)
    tempvar = bootstat.t_fdr01;    
    
    center = edge(1:end-1)+mean(diff(edge))/2;
    clear n
    for breaki = 1:3
        n(breaki,:) = histcounts(tempvar(breakid==breaki),edge);
        f = n/length(bootpar.p0);
    end
    
    bar(center,f','stacked')
    setplot('Panel B','t-hurdle for FDR = 1%','Frequency')
    
    ylim([0 0.35])


    % classical t-hurdle
    hv = vline(-1*norminv(0.01/2),'k--');
    text(2.6 ,0.3,'classical t-hurdle','color','k','fontsize',10)
    
    
set(gcf, 'position',     [876.8333  200.1667  613.3333  497.5000])
       


export_fig('../output_exhibits/boot-t-hurdles.pdf')

%% ==== FIGURE: ID'D STUFF ====

% = user =
% choose how to break up the chart
grey = [1 1 1]*0.5;
breakvar = bootpar.p0;
tempbreak = [0.10 0.40];

breakid = ones(size(bootpar.p0));
for i = 1:length(tempbreak)
    breakid(breakvar > tempbreak(i)) = i+1;
end

clf       
    
subplot(2,1,1)
    % settings
    tempvar = bootstat.pubfdr;
    edge = 0:2.5:25;
        
    % prep bars
    center = edge(1:end-1)+mean(diff(edge))/2;   
    clear n
    for breaki = 1:3
        n(breaki,:) = histcounts(tempvar(breakid==breaki),edge);
        f = n/length(bootpar.p0);
    end
    
    % draw bars
    h = bar(center,f','stacked');
    setplot('Panel A','','Frequency')
    xlabel('FDR Among Published Predictors (%)')
    
    ylim([0 0.35])          
    
    
    addspace(0.05)
    
subplot(2,1,2)

    % settings
    tempvar = bootstat.s_smooth;    
    edge = 10:2.5/2:30;
    
    center = edge(1:end-1)+mean(diff(edge))/2;
    clear n
    for breaki = 1:3
        n(breaki,:) = histcounts(tempvar(breakid==breaki),edge);
        f = n/length(bootpar.p0);
    end
    
    bar(center,f','stacked')
    setplot('Panel B','','Frequency')
    xlabel('Smooth Shrinkage for Published Predictors (%)')
    
    ylim([0 0.35])

    
    hleg = legend([h] ...
        ,{
        ['$\hat{\pi}_F <$  ' num2str(tempbreak(1)) ]
        ['$\hat{\pi}_F \in ($' num2str(tempbreak(1)) ',' num2str(tempbreak(2)) ']']
        ['$\hat{\pi}_F >$  ' num2str(tempbreak(2)) ]
        }...
        ,'interpreter','latex');
    hleg.FontSize = 12;
    
    
set(gcf, 'position',     [876.8333  200.1667  613.3333  497.5000])
       

hv = vline(26,'-');  hv.Color = [1 1 1]*0.4;
text(26.2,0.3,{'McLean-Pontiff','Upper Bound'},'color','k','fontsize',10)

export_fig('../output_exhibits/boot-pubstats.pdf')



%% bootstrap bias adjustments quantiles for reference
plist = [5 25 50 75 95];


% c = boot.rhopar;
% sdrho =      2*c./( 2*c.*sqrt(2*c+1) )   


tempboot = boot; 
datboot = [
    tempboot.t_fdr05
    tempboot.t_fdr01
    tempboot.pubfdr
    tempboot.smooth_s    
    tempboot.median_s    
    tempboot.mean_s  
    ];
numboot = ptile(datboot,plist,2)

sdboot = std(datboot,[],2)



open numboot
sum(tempboot.t_fdr05 > 1.96)/length(tempboot.t_fdr05)


