% 2019 03 Andrew
% checks results from bootstrap
% if run on hpc, edit ../hpc/check_boot.sh to select the output name


clear

%% user


outname = 'deleteme';
Nbootmax = 100;
Nshrinkmax = 100;

%% setup

% load poitn estimate stuff
folder = ['../output/' outname '/'];
load([folder 'samp.mat']);
load([folder 'meth.mat']);
load([folder 'est.mat']);
load([folder 'corrmat.mat']);

namevec = meth.namevec; % parameters that got estimated
par1 = est.par;


% find boot files
bootpath = ['../output/' outname '/boot/'];
filelist = dir([bootpath 'boot_est*.mat']);
filelist = {filelist.name}';


% limit the number of boots
Nboot= min(Nbootmax,length(filelist))




%%
% list stats to save
testname = {'t_fdr01' 't_fdr05' 't_fwer01' 't_fwer05' 'pubfdr'};
shrinkname = {'mean_s' 'median_s' 'smooth_s' 'mean_muhat' 'median_muhat'};
allname = [namevec testname];


    % for debugging
%     filelist = filelist(1:10);




    if 0
        disp('TESTING !!!!!!!!!!!!!!!!!!!!!')
        disp('TESTING !!!!!!!!!!!!!!!!!!!!!')    
        disp('TESTING !!!!!!!!!!!!!!!!!!!!!')

        meth.Nsim = 1;
        

    end    
    
%% feedback
disp(['checking bootstrap with meth'])


meth

Nboot
    

%% Simulate bootstrapped estimations


clear boot

clear par htest shrinkpub shrinksim
% load bootstrap results and store in boot struct
disp('importing bootstrap and calculating additional stats')
tic
% parfor esti = 1:Nboot
for esti = 1:Nboot
    if mod(esti,10) == 0
        esti
    end
    temp = load([bootpath filelist{esti}]);    
    par(esti) = temp.est.par;
    
%     % for t257, alter truncation to match HLZ baseline
%     if strcmp(name,'t257_500_96_300')
%         temp.est.par.thresh1 = 1.96;
%         temp.est.par.thresh3 = 0.50;
%     end


    
    % calculate multiple testing stats
    [err,sim,simflat] = obj_sim([],[],[],temp.est.par,corrmat,meth);   
    htest(esti) = test_hypotheses(sim,meth.mugood);
    shrinkpub(esti) = shrink_samp(temp.est.par,samp.rpub,samp.sigpub,Nshrinkmax);
    
    % aggregate simulated stuff
    shrinksim(esti).sim_mean_mupub      = mean(simflat.mupub);
    shrinksim(esti).sim_median_mupub    = median(simflat.mupub);
    shrinksim(esti).sim_mean_rpub      = mean(simflat.rpub);
    shrinksim(esti).sim_median_rpub    = median(simflat.rpub);    
    shrinksim(esti).sim_smooth_s        = 100*(1-mean(simflat.mupub)/mean(simflat.rpub));
    shrinksim(esti).sim_smooth_s_mix    = 100*(1-mean(simflat.mupub)/mean(samp.rpub));    
    
    
end

%% covert to boot struct
clear boot
for i = 1:length(namevec)
    boot.(namevec{i}) = [par(:).(namevec{i})];    
end

for i = 1:length(shrinkname)
    boot.(shrinkname{i}) = [shrinkpub(:).(shrinkname{i})];    
end

for i = 1:length(testname)
    boot.(testname{i}) = [htest(:).(testname{i})];    
end

shrinksimname = fieldnames(shrinksim(1));
for i = 1:length(shrinksimname)
    boot.(shrinksimname{i}) = [shrinksim(:).(shrinksimname{i})];    
end



sec = toc
minutes = sec/60


%% calculate stuff vs t-stat
plist = [5 25 50 75 95];
t = htest(1).t;

[fdrmat,mumat,rmat] = nanall(Nboot,length(t));
for esti = 1:Nboot
    fdrmat(esti,:) = htest(esti).fdr;
    mumat(esti,:) = htest(esti).mu_mean;
    rmat(esti,:) = htest(esti).r_mean;
end
smat = 100*(1-mumat./rmat);


tfun.t = t;
tfun.fdr = ptile(fdrmat,plist,1);
tfun.s = ptile(smat,plist,1);
tfun.mu = ptile(mumat,plist,1);
tfun.r = ptile(rmat,plist,1);
tfun.plist  = plist;



%% add thresh3 if it's not there

if ~isfield(boot, 'thresh3')
    boot.thresh3 = repmat(par1.thresh3, size(boot.pubfdr));
end


%% figure: boot: stat hist by omega
figure(1)
plist = [5 95];

clf

    
subplot(2,2,1)
%
    tempvar = boot.t_fdr05;
    edge = [0:0.25:4];
    
    center = edge(1:end-1)+mean(diff(edge))/2;
    thresh3list = unique(boot.thresh3);
    clear n
    for omegai = 1:length(thresh3list)
        id = boot.thresh3 == thresh3list(omegai); 
        n(omegai,:) = histcounts(tempvar(id),edge);
        f = n/length(boot.p0);
    end
    
    % draw bars
    h = bar(center,f','stacked');
    setplot('','t-hurdle for FDR = 5%','Frequency')
    
    ylim([0 0.35])
    
    % 5, 95% CI
    ci = ptile(tempvar,plist,2);       
    hv = vline(ci,'r--');
    text(ci(1)+0.03,0.3,'5%','color','red','fontsize',10)
    text(ci(2)-0.5,0.3,'95%','color','red','fontsize',10)
    
    
    hleg = legend([h] ...
        ,{
        '$\hat{\omega}$ = 1/3'
        '$\hat{\omega}$ = 1/2'
        '$\hat{\omega}$ = 2/3'
        }...
        ,'position',[0.3374 0.7508 0.1169 0.0868] ...
        ,'interpreter','latex');
    
subplot(2,2,2)
    tempvar = boot.t_fdr01;
    edge = [0:0.25:4];
    
    center = edge(1:end-1)+mean(diff(edge))/2;
    thresh3list = unique(boot.thresh3);
    clear n
    for omegai = 1:length(thresh3list)
        id = boot.thresh3 == thresh3list(omegai); 
        n(omegai,:) = histcounts(tempvar(id),edge);
        f = n/length(boot.p0);
    end
    
    bar(center,f','stacked')
    setplot('','t-hurdle for FDR = 1%','Frequency')
    
    ylim([0 0.35])

    % 5, 95% CI
    ci = ptile(tempvar,plist,2);       
    hv = vline(ci,'r--');
    text(ci(1)+0.03,0.3,'5%','color','red','fontsize',10)
    text(ci(2)-0.6,0.3,'95%','color','red','fontsize',10)
    
    
    

subplot(2,2,3)
    tempvar = boot.mean_s;
    edge = [5:2.5:40];
    
    center = edge(1:end-1)+mean(diff(edge))/2;
    thresh3list = unique(boot.thresh3);
    clear n
    for omegai = 1:length(thresh3list)
        id = boot.thresh3 == thresh3list(omegai); 
        n(omegai,:) = histcounts(tempvar(id),edge);
        f = n/length(boot.p0);
    end
    
    bar(center,f','stacked')
    setplot('','Mean Shrinkage for Published (%)','Frequency')
    

    ylim([0 0.30])

    % 5, 95% CI
    ci = ptile(tempvar,plist,2);       
    hv = vline(ci,'r--');
    text(ci(1)+1  ,0.25 ,'5%','color','red','fontsize',10)
    text(ci(2)-5  ,0.25 ,'95%','color','red','fontsize',10)

    
subplot(2,2,4)
    tempvar = boot.pubfdr;
    edge = [0:2.5:40];
    
    center = edge(1:end-1)+mean(diff(edge))/2;
    thresh3list = unique(boot.thresh3);
    clear n
    for omegai = 1:length(thresh3list)
        id = boot.thresh3 == thresh3list(omegai); 
        n(omegai,:) = histcounts(tempvar(id),edge);
        f = n/length(boot.p0);
    end
    
    bar(center,f','stacked')

    setplot('','FDR for Published (%)','Frequency')
    
    ylim([0 0.30])

    % 5, 95% CI
    ci = ptile(tempvar,plist,2);       
    hv = vline(ci,'r--');
    text(ci(1)+3  ,0.25 ,'5%','color','red','fontsize',10)
    text(ci(2)-5  ,0.25 ,'95%','color','red','fontsize',10)
    



%% figure: boot: stat hist by omega: 6 panel
grey = [1 1 1]*0.5;
figure(2)

clf

   
    
subplot(3,2,1)
%
    tempvar = boot.t_fdr05;
    edge = [0:0.25:4];
    
    center = edge(1:end-1)+mean(diff(edge))/2;
    thresh3list = unique(boot.thresh3);
    clear n
    for omegai = 1:length(thresh3list)
        id = boot.thresh3 == thresh3list(omegai); 
        n(omegai,:) = histcounts(tempvar(id),edge);
        f = n/length(boot.p0);
    end
    
    % draw bars
    h = bar(center,f','stacked');
    setplot('Panel A','t-hurdle for FDR = 5%','Frequency')
    
    ylim([0 0.35])
    
    % traditional t-hurdle
    hv = vline(-1*norminv(0.05/2),'k--');
    text(2,0.3,'traditional t-hurdle','color','k','fontsize',10)
    
    
    hleg = legend([h] ...
        ,{
        '$\hat{\omega}$ = 1/3'
        '$\hat{\omega}$ = 1/2'
        '$\hat{\omega}$ = 2/3'
        }...
        ,'interpreter','latex');
    hleg.Position = [  0.1451    0.8283    0.1158    0.0860];
    
subplot(3,2,2)
    tempvar = boot.t_fdr01;
    edge = [0:0.25:4];
    
    center = edge(1:end-1)+mean(diff(edge))/2;
    thresh3list = unique(boot.thresh3);
    clear n
    for omegai = 1:length(thresh3list)
        id = boot.thresh3 == thresh3list(omegai); 
        n(omegai,:) = histcounts(tempvar(id),edge);
        f = n/length(boot.p0);
    end
    
    bar(center,f','stacked')
    setplot('Panel B','t-hurdle for FDR = 1%','Frequency')
    
    ylim([0 0.35])


    % traditional t-hurdle
    hv = vline(-1*norminv(0.01/2),'k--');
    text(0.7 ,0.3,'traditional t-hurdle','color','k','fontsize',10)
    
    
    
subplot(3,2,3)
    tempvar = boot.pubfdr;
    edge = [0:2.5:40];
    
    center = edge(1:end-1)+mean(diff(edge))/2;
    thresh3list = unique(boot.thresh3);
    clear n
    for omegai = 1:length(thresh3list)
        id = boot.thresh3 == thresh3list(omegai); 
        n(omegai,:) = histcounts(tempvar(id),edge);
        f = n/length(boot.p0);
    end
    
    bar(center,f','stacked')

    setplot('Panel C','FDR for Published (%)','Frequency')
    
    ylim([0 0.30])

    % median
    hv = vline(median(tempvar),'r--');
    text(9, 0.2,'median','color','red','fontsize',10)
    
    

subplot(3,2,4)
    tempvar = boot.smooth_s;
    edge = [5:2.5:40];
    
    center = edge(1:end-1)+mean(diff(edge))/2;
    thresh3list = unique(boot.thresh3);
    clear n
    for omegai = 1:length(thresh3list)
        id = boot.thresh3 == thresh3list(omegai); 
        n(omegai,:) = histcounts(tempvar(id),edge);
        f = n/length(boot.p0);
    end
    
    bar(center,f','stacked')
    setplot('Panel D','Smooth Shrinkage for Published (%)','Frequency')
    

    ylim([0 0.30])

    % median
    hv = vline(median(tempvar),'r--');
    text(20, 0.2,'median','color','red','fontsize',10)
    
    


subplot(3,2,5)
    tempvar = boot.median_s;
    edge = [5:2.5:40];
    
    center = edge(1:end-1)+mean(diff(edge))/2;
    thresh3list = unique(boot.thresh3);
    clear n
    for omegai = 1:length(thresh3list)
        id = boot.thresh3 == thresh3list(omegai); 
        n(omegai,:) = histcounts(tempvar(id),edge);
        f = n/length(boot.p0);
    end
    
    bar(center,f','stacked')
    setplot('Panel E','Median Shrinkage for Published (%)','Frequency')
    

    ylim([0 0.30])

    % median
    hv = vline(median(tempvar),'r--');
    text(20, 0.2,'median','color','red','fontsize',10)
    
    
subplot(3,2,6)
    tempvar = boot.mean_s;
    edge = [5:2.5:40];
    
    center = edge(1:end-1)+mean(diff(edge))/2;
    thresh3list = unique(boot.thresh3);
    clear n
    for omegai = 1:length(thresh3list)
        id = boot.thresh3 == thresh3list(omegai); 
        n(omegai,:) = histcounts(tempvar(id),edge);
        f = n/length(boot.p0);
    end
    
    bar(center,f','stacked')
    setplot('Panel F','Mean Shrinkage for Published (%)','Frequency')
    

    ylim([0 0.30])

    % median
    hv = vline(median(tempvar),'r--');
    text(20, 0.2,'median','color','red','fontsize',10)
        



    

%% some summ stats
plist = [5 50 95];
pname = {'p05', 'p50', 'p95'};

% htest stats
clear num
for i = 1:length(testname)
    num(i,:) = ptile(boot.(testname{i}),plist);
end

tab_htest = array2table(num,'rownames',testname,'variablenames',pname)


% shrink stats
clear num
for i = 1:length(shrinkname)
    num(i,:) = ptile(boot.(shrinkname{i}),plist);
end

tab_shrink = array2table(num,'rownames',shrinkname,'variablenames',pname)



%% export
    
slurm_ppn = str2num(getenv('SLURM_CPUS_ON_NODE'));

if ~isempty(slurm_ppn)

	disp(sprintf('Running on SLURM: auto exporting'));
    
    
    folder = ['../output/boot_' name '/results/']

    undock
    figure(1)
    set(gcf,'position', [    2132         106           670         512 ])
    export_fig([folder 'boot_4pan_' meth.name '.pdf'])
    
    
    figure(2)
    set(gcf,'position', [  1512         151         818         878 ])
    export_fig([folder 'boot_6pan_' meth.name '.pdf'])
    

    
    
    
    writetable(tab_htest,[folder 'boot_' meth.name '.xlsx'] ...
        ,'sheet','sheet1','WriteRowNames',true,'range','A2')  
    writetable(tab_shrink,[folder 'boot_' meth.name '.xlsx'] ...
        ,'sheet','sheet1','WriteRowNames',true,'range','A10')     
    
    save([folder 'check_boot_' meth.name '.mat'])
    
end

