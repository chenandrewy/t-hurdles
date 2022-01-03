%% choose settings
restoredefaultpath

% add my tools to path
mytools = genpath(['../Matlab_Tools/']);
path(mytools,path);
clear mytools;
addpath('../data');
addpath('../code_exhibits');

clear; cls; clc

%% setup

name = 'set_fast';

% Load Baseline Settings
[par,meth] = baseline_settings;
eval(name)

% modify, if desitred
meth.seed = 'shuffle';
meth.Nportall = 200;
par.corrtype = 'saved';
% par.rho1 = 0;

% plot settings
parvecname = {'p0'};
Nplot = 11;
booti = 3;
Nsimlist = 100 + [0 100 200];


[samp,par,corrmat,meth] = auto_prep(meth,par,booti);
[~,htest,shrink,~,~] = addstats(samp,par,corrmat,meth);


%% calculate objective
clear alldata leglab
for Nsimi = 1:length(Nsimlist)
    Nsim = Nsimlist(Nsimi)
    
    tempmeth = meth;
    temppar = par;
    tempmeth.Nsim = Nsim;
    
    alldata{Nsimi} = plot_obj(tempmeth,temppar,booti,Nplot,parvecname);

    leglab{Nsimi} = ['Nsim = ' int2str(Nsimlist(Nsimi))];
end

%% plot together

clf; hold on;

clist = linspace(0,0.9,length(Nsimlist));

for pari = 1:length(alldata{1})
    figure; hold on;
    for Nsimi = 1:length(Nsimlist)    
        plotdat = alldata{Nsimi};    

        
%         subplot(1,3,pari); 
%         hold on;
        plot(plotdat{pari}.val,plotdat{pari}.obj,'x-')
        vline(par.(parvecname{pari}))
        xlabel(parvecname{pari}); ylabel('GMM obj')
        setplot        
        legend(leglab, 'location','best')
    end    
end




ylim([35 45])

%% ==== END NSIM CHECK ====
return

%% check Nportall
meth.Nsim = 50;
parlist = [200:100:1000];

%% calculate objective
clear alldata leglab
parfor i = 1:length(parlist)
    Nsim = parlist(i)
    
    tempmeth = meth;
    temppar = par;
    tempmeth.Nportall = parlist(i);
    
    alldata{i} = plot_obj(tempmeth,temppar,booti,Nplot,parvecname);    
end





%% plot together
cls

clf; hold on;

clear leglab

for i = 1:length(parlist)
    leglab{i} = ['Nportall = ' int2str(parlist(i))];
end

for pari = 1:length(alldata{1})
    figure; hold on;
    for Nsimi = 1:length(parlist)    
        plotdat = alldata{Nsimi};    

        
%         subplot(1,3,pari); 
%         hold on;
        plot(plotdat{pari}.val,plotdat{pari}.obj,'x-')
        vline(par.(parvecname{pari}))
        xlabel(parvecname{pari}); ylabel('GMM obj')
        setplot        
        legend(leglab, 'location','best')
    end    

end



ylim([35 45])