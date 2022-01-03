% 2019 03 Andrew
% first file of new 2019 WFA draft

grey = [1 1 1]*0.5;


clear
% simple example using hlz

% restoredefaultpath

blueauto = [     0    0.4470    0.7410];
redauto  = [  0.8500    0.3250    0.0980];
green = [0 0.6 0];

%% User

rng(2)


[parbase,meth] = set_hlz_alt(-1);


%% simulate Np0 versions


meth.mugood = 'strict';

% p0 and mu are from Table 5, page 30.  thresh is found in table footnote
p0list      = [
    0.444
    0
    0.444
    ];

% this one fits better
mulist = {
    {'gamma', [0 0.555 1.00]}
    {'gamma', [0 0.555 1.00]} % old: 0.5050
    {'gamma', [0 0.250 1.00]}    
};


threshlist = {
    {'step', [1.96 2.57 0]}
    {'step', [1.96 2.57 0]}
    {'step', [1.96 2.57 0]}
};
    



Nalllist = [
    1378
    900
    1378
];    


disp('calculating model sim and stats')

all = {};
clear err pub mom par sim simflat
tic
for pari = 1:3
    
    % load parameters    
    temppar = parbase;
    
    temppar.p0      = p0list(pari);
    temppar.mutype = mulist{pari}{1};
    temppar.mu1 = mulist{pari}{2}(1);
    temppar.mu2 = mulist{pari}{2}(2);
    temppar.mu3 = mulist{pari}{2}(3);    
    
    temppar.threshtype  = threshlist{pari}{1};
    temppar.thresh1 = threshlist{pari}{2}(1);
    temppar.thresh2 = threshlist{pari}{2}(2);
    temppar.thresh3 = threshlist{pari}{2}(3);    
    
    temppar.Nportall = Nalllist(pari);
    
    % simulate
    [sim{pari},simflat{pari}] = sim_hlz(temppar,meth.Nsim);
    
    
    % store par
    par{pari} = temppar;   
    

end
toc



%% additional calcs

% adjust sample following HLZ
clear simadj
for pari = 1:3
    
    for simi = 1:meth.Nsim
    
        % adjust
        idmarginal = find(sim{pari}.tpub{simi} < 2.57);
        idgood = find(sim{pari}.tpub{simi}  > 2.57);

%             idmarginal = find(abs(sim{pari}.tpub{simi}) < 2.57); % testing
%             idgood = find(abs(sim{pari}.tpub{simi})  > 2.57);
        
        Nmarginal = length(idmarginal);
        idkeep = [idmarginal(1:floor(Nmarginal/2)) idgood];

        % shuffle so we dont' have to do so many integrations
        idkeep = idkeep(randperm(length(idkeep)));


        simadj{pari} = sim{pari};

        simadj{pari}.rpub{simi} = sim{pari}.rpub{simi}(idkeep);
        simadj{pari}.sigpub{simi} = sim{pari}.sigpub{simi};       
        simadj{pari}.mupub{simi} = sim{pari}.mupub{simi}(idkeep);           
        simadj{pari}.nullpub{simi} = sim{pari}.nullpub{simi}(idkeep);       
        
    end
    
    % flatten
    fn = fieldnames(simadj{pari});
    for fni = 1:length(fn)    
        temp = simadj{pari}.(fn{fni});
        if iscell(temp)
            simadjflat{pari}.(fn{fni}) = [temp{:}]';
        else
            simadjflat{pari}.(fn{fni}) = temp(:)';
        end
    end    
    
end

%%

Nshrinkmax = 5000;

% hypothesis testing
clear htest shrink
tic
for pari = 1:3
    htest{pari} = test_hypotheses(sim{pari},meth.mugood);
end
toc

%% table: parameters

par{1}

tab_par = [
    par{1}.Nportall par{1}.p0 par{1}.mu2 
    par{2}.Nportall par{2}.p0 par{2}.mu2
    ]

open tab_par


%% table: moments

datamom = [353 2.39 3.16 6.34]; % HLZ p 28 middle


for pari = 1:2
    for simi = 1:meth.Nsim
        sim{pari}.momt(simi,:) = quantile(sim{pari}.tpub{pari},[0.2 0.5 0.9]);
    end
end

simmom = [
    median(sim{1}.Npub) median(sim{1}.momt)
    median(sim{2}.Npub) median(sim{2}.momt)    
    ];

momerr = simmom-datamom;
   
% HLZ p 28 bottom par
weight = [1 10000*ones(1,3)];     
err = sum(weight.*(momerr.^2),2);

tab_mom = [[datamom; simmom] nan(3,1) [nan; err] ]


%%
return

%% stats for reference

tab_stat = [
    htest{1}.t_fdr05 htest{1}.t_fdr01 htest{1}.pubfdr 
    htest{2}.t_fdr05 htest{2}.t_fdr01 htest{2}.pubfdr 
]

tab_mu = [
    mean(simflat{1}.mupub)  median(simflat{1}.mupub)    
    mean(simflat{2}.mupub)  median(simflat{2}.mupub)
    ]

%% plot 2 panel paper fig



clf

% all t-stats    
subplot(2,1,1); hold on;

edge = [-3:0.25:10]; % t
center = mean([edge(1:end-1); edge(2:end)],1);
    
    c1 = histcounts([simflat{1}.t],edge);
    h1 = bar(center,c1/meth.Nsim); h1.FaceAlpha = 0.7; h1.FaceColor = blueauto;
    
    c2 = histcounts([simflat{2}.t],edge);
    h2 = bar(center,c2/meth.Nsim); h2.FaceAlpha = 0.7; h2.FaceColor = redauto;
    
    setplot('Panel A: All Factors'' t-stats','t-stat $t_i$','# of Factors')
    xlhand = get(gca,'title'); set(xlhand, 'interpreter', 'latex');
    xlhand = get(gca,'xlabel'); set(xlhand, 'interpreter', 'latex');
    
    legend({...
        ['HLZ baseline'] ...
        ,['Alternative'] ...
        })
    
    
    ylim([0 120])    
    hv = vline(1.96,'r-'); hv.LineWidth = 4;    
    text(2.3,110,'Observed ->','fontsize',11,'color','r')    
    text(-0.8,110,'<- Unobserved','fontsize',11,'color','r')
    
    
    
    

subplot(2,1,2); hold on;

    muedge = [0:1/20:3.5];
    histogram(simadjflat{1}.mupub,muedge,'Normalization','probability');
    histogram(simadjflat{2}.mupub,muedge,'Normalization','probability');
    setplot('Panel B: Observed Factors'' Expected Returns','Expected Return $\mu_i$ (\%, monthly)','Frequency')
    xlhand = get(gca,'title'); set(xlhand, 'interpreter', 'latex');
    xlhand = get(gca,'xlabel'); set(xlhand, 'interpreter', 'latex');
    
%     xlim([-0.1 2.5])
    
    legend({...
        ['HLZ baseline'] ...
        ,['Alternative'] ...
        })    
    

set(gcf,'position', [     1139         232         717         549  ])

if 0 
    % export   
    export_fig('../output_demos/hlz_demo_2pan.png','-m8')
end







%% plot 1 panel slides 1



clf; hold on;

% all t-stats    


edge = [-3:0.25:10]; % t
center = mean([edge(1:end-1); edge(2:end)],1);
    
    c1 = histcounts([simflat{1}.t],edge);
    h1 = bar(center,c1/meth.Nsim); h1.FaceAlpha = 0.7; h1.FaceColor = blueauto;
    
%     c2 = histcounts([simflat{2}.t],edge);
%     h2 = bar(center,c2/meth.Nsim); h2.FaceAlpha = 0.7; h2.FaceColor = redauto;
    
    
    setplot('','t-stat','# of Factors')

    
    hleg = legend({...
         ['Set 1: HLZ Estimate: p_0 = 44%, \lambda = 0.56%'] ...
        ,['Set 2: Alternative    : p_0 =    0%, \lambda = 0.51%'] ....
        },'fontsize',11);
    hleg.Position = [ 0.5615    0.5564    0.3114    0.1006];
    
    
    ylim([0 120])    
%     hv = vline(1.96,'r-'); hv.LineWidth = 4;    
%     text(2.3,110,'Reported ->','fontsize',11,'color','r')    
%     text(-0.8,110,'<- Unreported','fontsize',11,'color','r')
    
    
set(gcf,'position', [       1143         438         713         343
  ])    
    

if 0 
    % export
    export_fig('../output_demos/hlz_t_1.pdf')
end




%% plot 1 panel slides 2



clf; hold on;

% all t-stats    


edge = [-3:0.25:10]; % t
center = mean([edge(1:end-1); edge(2:end)],1);
    
    c1 = histcounts([simflat{1}.t],edge);
    h1 = bar(center,c1/meth.Nsim); h1.FaceAlpha = 0.7; h1.FaceColor = blueauto;
    
    c2 = histcounts([simflat{2}.t],edge);
    h2 = bar(center,c2/meth.Nsim); h2.FaceAlpha = 0.7; h2.FaceColor = redauto;
    
    
    setplot('','t-stat','# of Factors')

    
    hleg = legend({...
         ['Set 1: HLZ Estimate: p_0 = 44%, \lambda = 0.56%'] ...
        ,['Set 2: Alternative    : p_0 =    0%, \lambda = 0.51%'] ....
        },'fontsize',11);
    hleg.Position = [ 0.5615    0.5564    0.3114    0.1006];
    
    
    ylim([0 120])    
%     hv = vline(1.96,'r-'); hv.LineWidth = 4;    
%     text(2.3,110,'Reported ->','fontsize',11,'color','r')    
%     text(-0.8,110,'<- Unreported','fontsize',11,'color','r')
    
    
set(gcf,'position', [       1143         438         713         343
  ])    
    

if 0 
    % export
    export_fig('../output_demos/hlz_t_2.pdf')
end




%% plot 1 panel slides 3



clf; hold on;

% all t-stats    


edge = [-3:0.25:10]; % t
center = mean([edge(1:end-1); edge(2:end)],1);
    
    c1 = histcounts([simflat{1}.t],edge);
    h1 = bar(center,c1/meth.Nsim); h1.FaceAlpha = 0.7; h1.FaceColor = blueauto;
    
    c2 = histcounts([simflat{2}.t],edge);
    h2 = bar(center,c2/meth.Nsim); h2.FaceAlpha = 0.7; h2.FaceColor = redauto;
    
    
    setplot('','t-stat','# of Factors')

    
    hleg = legend({...
         ['Set 1: HLZ Estimate: p_0 = 44%, \lambda = 0.56%'] ...
        ,['Set 2: Alternative    : p_0 =    0%, \lambda = 0.51%'] ....
        },'fontsize',11);
    hleg.Position = [ 0.5615    0.5564    0.3114    0.1006];
    
    
    ylim([0 120])    
    hv = vline(1.96,'r-'); hv.LineWidth = 4;    
%     text(2.3,110,'Reported and Observed->','fontsize',11,'color','r')    
%     text(-0.8,110,'<- Unreported','fontsize',11,'color','r')
    
    
set(gcf,'position', [       1143         438         713         343
  ])    
    

if 0 
    % export
    export_fig('../output_demos/hlz_t_3.pdf')
end




%% plot 1 panel slides 3



clf; hold on;

% all t-stats    


edge = [-3:0.25:10]; % t
center = mean([edge(1:end-1); edge(2:end)],1);
    
    c1 = histcounts([simflat{1}.t],edge);
    h1 = bar(center,c1/meth.Nsim); h1.FaceAlpha = 0.7; h1.FaceColor = blueauto;
    
    c2 = histcounts([simflat{2}.t],edge);
    h2 = bar(center,c2/meth.Nsim); h2.FaceAlpha = 0.7; h2.FaceColor = redauto;
    
    
    setplot('','t-stat','# of Factors')

    
    hleg = legend({...
         ['Set 1: HLZ Estimate: p_0 = 44%, \lambda = 0.56%'] ...
        ,['Set 2: Alternative    : p_0 =    0%, \lambda = 0.51%'] ....
        },'fontsize',11);
    hleg.Position = [ 0.5615    0.5564    0.3114    0.1006];
    
    
    ylim([0 120])    
    hv = vline(1.96,'r-'); hv.LineWidth = 4;    
    text(2.3,110,'Reported and Observed->','fontsize',11,'color','r')    
    text(-0.8,110,'<- Unreported','fontsize',11,'color','r')
    
    
set(gcf,'position', [       1143         438         713         343
  ])    
    

if 0 
    % export
    export_fig('../output_demos/hlz_t_3b.pdf')
end


%% plot 1 panel slides 4



clf; hold on;

% all t-stats    


edge = [-3:0.25:10]; % t
center = mean([edge(1:end-1); edge(2:end)],1);
i = center > 1.96;
    
    c1 = histcounts([simflat{1}.t],edge);
    h1 = bar(center(i),c1(i)/meth.Nsim); h1.FaceAlpha = 0.7; h1.FaceColor = blueauto;
    
    c2 = histcounts([simflat{2}.t],edge);
    h2 = bar(center(i),c2(i)/meth.Nsim); h2.FaceAlpha = 0.7; h2.FaceColor = redauto;
    
    setplot('','t-stat','# of Factors')

    
    hleg = legend({...
         ['Set 1: HLZ Estimate: p_0 = 44%, \lambda = 0.56%'] ...
        ,['Set 2: Alternative    : p_0 =    0%, \lambda = 0.51%'] ....
        },'fontsize',11);
    hleg.Position = [ 0.5615    0.5564    0.3114    0.1006];
    
        
    ylim([0 120]); xlim([-4 10]);
    hv = vline(1.96,'r-'); hv.LineWidth = 4;    
    text(2.3,110,'Reported and Observed->','fontsize',11,'color','r')    
    text(-0.8,110,'<- Unreported','fontsize',11,'color','r')
    
    
set(gcf,'position', [       1143         438         713         343
  ])    
    

if 0 
    % export
    export_fig('../output_demos/hlz_t_4.pdf')
end



%% plot 1 panel slides 5



clf; hold on;

% all t-stats    


edge = [-3:0.25:10]; % t
center = mean([edge(1:end-1); edge(2:end)],1);
    
    c1 = histcounts([simflat{1}.t],edge);
    h1 = bar(center,c1/meth.Nsim); h1.FaceAlpha = 0.7; h1.FaceColor = blueauto;
    
    c2 = histcounts([simflat{2}.t],edge);
    h2 = bar(center,c2/meth.Nsim); h2.FaceAlpha = 0.7; h2.FaceColor = redauto;
    
    c3 = histcounts([simflat{3}.t],edge);
    h3 = bar(center,c3/meth.Nsim); h3.FaceAlpha = 0.7; h3.FaceColor = green;
    
    setplot('','t-stat','# of Factors')

    
    hleg = legend({...
         ['Set 1: HLZ Estimate: p_0 = 44%, \lambda = 0.56%'] ...
        ,['Set 2: Alternative    : p_0 =    0%, \lambda = 0.51%'] ....
        ,['Set 3: Small \lambda            : p_0 = 44%, \lambda = 0.25%'] ....
        },'fontsize',11);
    hleg.Position = [ 0.5615    0.5564    0.3114    0.1006];
    
    
    ylim([0 120])    
    hv = vline(1.96,'r-'); hv.LineWidth = 4;    
    text(2.3,110,'Reported and Observed->','fontsize',11,'color','r')    
    text(-0.8,110,'<- Unreported','fontsize',11,'color','r')
    
    
set(gcf,'position', [       1143         438         713         343
  ])    
    

if 0 
    % export
    export_fig('../output_demos/hlz_t_5.png','-m8')
end



%% plot slides consistent panel


clf; hold on;
   

    muedge = [0:1/20:3.5];
    histogram(simadjflat{1}.mupub,muedge,'Normalization','probability');
    histogram(simadjflat{2}.mupub,muedge,'Normalization','probability');
    setplot('Published Factors Only','Expected Return $\mu_i$ (\%, monthly)','Frequency')
    xlhand = get(gca,'title'); set(xlhand, 'interpreter', 'latex');
    xlhand = get(gca,'xlabel'); set(xlhand, 'interpreter', 'latex');
    
%     xlim([-0.1 2.5])
    
    hleg = legend({...
         ['Set 1: HLZ Estimate: p_0 = 44%, \lambda = 0.56%'] ...
        ,['Set 2: Alternative    : p_0 =    0%, \lambda = 0.51%'] ....
        },'fontsize',11);  
    
set(gcf,'position', [       1143         438         713         343
  ])    
    
if 0 
    % export   
    export_fig('../output_demos/hlz_mu.pdf')
end




%% plot 1 panel slides intro



clf; hold on;

% all t-stats    


edge = [-3:0.25:10]; % t
center = mean([edge(1:end-1); edge(2:end)],1);

    
    c1 = histcounts([simflat{1}.t],edge);
    h1 = bar(center,c1/meth.Nsim); h1.FaceAlpha = 0.7;
    
    c2 = histcounts([simflat{2}.t],edge);
    h2 = bar(center,c2/meth.Nsim); h2.FaceAlpha = 0.7;
    
    setplot('','t-stat','# of Factors')

    
    hleg = legend({...
         ['Model 1: Many Insignificant Predictors'] ...
        ,['Model 2: Few  Insignificant Predictors'] ....
        },'fontsize',11);
    hleg.Position = [ 0.5615    0.5564    0.3114    0.1006];
    
    
    ylim([0 120]); xlim([-3.1750   10.1750]);
    hv = vline(1.96,'r-'); hv.LineWidth = 4;    
    text(2.3,110,'Published and Observed ->','fontsize',11,'color','r')    
    text(-0.8,110,'<- Unreported','fontsize',11,'color','r')
    
    
    
    
set(gcf,'position', [       1143         438         713         343
  ])    
    

if 0 
    % export
    export_fig('../output_demos/hlz_t_intro_1.png','-m8')
end


%% plot 1 panel slides intro 2



clf; hold on;

% all t-stats    


edge = [-3:0.25:10]; % t
center = mean([edge(1:end-1); edge(2:end)],1);
i = center > 1.96;

    c1 = histcounts([simflat{1}.t],edge);
    h1 = bar(center(i),c1(i)/meth.Nsim); h1.FaceAlpha = 0.7;
    
    c2 = histcounts([simflat{2}.t],edge);
    h2 = bar(center(i),c2(i)/meth.Nsim); h2.FaceAlpha = 0.7;
    
    setplot('','t-stat','# of Factors')

    
    hleg = legend({...
         ['Model 1: Many Insignificant Predictors'] ...
        ,['Model 2: Few  Insignificant Predictors'] ....
        },'fontsize',11);
    hleg.Position = [ 0.5615    0.5564    0.3114    0.1006];
    
    
    ylim([0 120]); xlim([-3.1750   10.1750]);
    hv = vline(1.96,'r-'); hv.LineWidth = 4;    
    text(2.3,110,'Published and Observed ->','fontsize',11,'color','r')    
    text(-0.8,110,'<- Unreported','fontsize',11,'color','r')
    
    
set(gcf,'position', [       1143         438         713         343
  ])    
    

if 0 
    % export
    export_fig('../output_demos/hlz_t_intro_2.png','-m8')
end



%% test


clf; hold on;

edge = [0:1/20:2.0];
center = mean([edge(1:end-1); edge(2:end)],1);

    c1 = histcounts([simadjflat{1}.mu],edge);
    h1 = bar(center,c1/meth.Nsim); h1.FaceAlpha = 0.7;
    
    c2 = histcounts([simadjflat{2}.mu],edge);
    h2 = bar(center,c2/meth.Nsim); h2.FaceAlpha = 0.7;
    
    

    
    hleg = legend({...
         ['Set 1: HLZ Estimate: p_0 = 44%, \lambda = 0.56%'] ...
        ,['Set 2: Alternative    : p_0 =    0%, \lambda = 0.51%'] ....
        },'fontsize',11);    

%% t-hurdle vs shrinkage slide



clf

% all t-stats    
subplot(1,2,1); hold on;

edge = [0:1/20:2.0];
center = mean([edge(1:end-1); edge(2:end)],1);

    c1 = histcounts([simadjflat{1}.mu],edge);
    h1 = bar(center,c1/meth.Nsim); h1.FaceAlpha = 0.7;
    
    c2 = histcounts([simadjflat{2}.mu],edge);
    h2 = bar(center,c2/meth.Nsim); h2.FaceAlpha = 0.7;
    
    setplot('All Factors','Expected Return (%, Monthly)','# of Factors')
    xlhand = get(gca,'title'); set(xlhand, 'interpreter', 'latex');
    xlhand = get(gca,'xlabel'); set(xlhand, 'interpreter', 'latex');
    
    legend({...
        ['HLZ baseline'] ...
        ,['Alternative'] ...
        })
    
    
%     ylim([0 120])    
%     hv = vline(1.96,'r-'); hv.LineWidth = 4;    
%     text(2.3,110,'Observed ->','fontsize',11,'color','r')    
%     text(-0.8,110,'<- Unobserved','fontsize',11,'color','r')
    
    
    
    

subplot(1,2,2); hold on;

    muedge = [0:1/20:3.0];
    histogram(simadjflat{1}.mupub,muedge,'Normalization','probability');
    histogram(simadjflat{2}.mupub,muedge,'Normalization','probability');
    setplot('Published Only','Expected Return $\mu_i$ (\%, monthly)','Frequency')
    xlhand = get(gca,'title'); set(xlhand, 'interpreter', 'latex');
    xlhand = get(gca,'xlabel'); set(xlhand, 'interpreter', 'latex');
    
%     xlim([-0.1 2.5])
    
    legend({...
        ['HLZ baseline'] ...
        ,['Alternative'] ...
        })    
    

set(gcf,'position', [     1139         232         717         549  ])

if 0 
    % export   
    export_fig('../output_demos/hlz_demo_2pan.png','-m8')
end




%% testing: BH and BY adjustments

t = sim{1}.tpub{1};

t = normrnd(5,1,[1 400]);
N = length(t);
p = normcdf(-t);

ps = sort(p);
is = (1:N);
c = sum(1./is);

c = 1;

histogram(ps)

istar = find(ps <= (is./(N*c)*0.05), 1, 'last');
pstar = ps(istar);

tstar = -norminv(pstar)

%% BH BY matrix form

simi = 1;

t = sim{simi}.t;
t(~sim{simi}.pub) = nan;
Npub = sim{simi}.Npub;
Nall = size(t,2);


p = normcdf(-t);
ps = sort(p,2);
is = ones(meth.Nsim,1)* (1:size(t,2));

Nmat = repmat(Npub,[1 Nall]);

temp = sort(sim{simi}.pub,2,'descend');
c = sum(1./is.*temp,2);

fdr_ub = 0.05;

reject = ps <= is./(Nmat.*c)*fdr_ub;
pstar = max(ps.*reject,[],2);

tstar = -norminv(pstar);

clf
plot(tstar)

%% check some simple stats from hlz data



for simi = 1
    mean(simadj{1}.nullpub{1})
    mean(simadj{1}.tpub{1})
end

%% ==== 2021 11 REVISIT ====

% HLZ demo

clf; hold on;

% all t-stats    


edge = [-3:0.25:10]; % t
center = mean([edge(1:end-1); edge(2:end)],1);

    
    c1 = histcounts([simflat{1}.t],edge);
    h1 = bar(center,c1/meth.Nsim); h1.FaceAlpha = 0.7;
    
    c2 = histcounts([simflat{2}.t],edge);
    h2 = bar(center,c2/meth.Nsim); h2.FaceAlpha = 0.7;
    
    setplot('','|t-stat|','# of Factors')

    
    hleg = legend({...
         ['Model 1: HLZ''s Baseline (\pi_F = 0.44)'] ...
        ,['Model 2: HLZ''s Baseline w/ \pi_F = 0'] ....
        },'fontsize',11);
    hleg.Position = [ 0.5615    0.5564    0.3114    0.1006];
    
    
    ylim([0 120]); xlim([-3.1750   10.1750]);
    hv = vline(1.96,'r-'); hv.LineWidth = 4;    
    text(2.3,110,'Published and Observed ->','fontsize',11,'color','r')    
    text(-0.8,110,'<- Unreported','fontsize',11,'color','r')
    
    
    
    
set(gcf,'position', [       1143         438         713         343
  ])    
    

if 0 
    % export
    export_fig('../output_demos/hlz-demo-simple.pdf')
end
