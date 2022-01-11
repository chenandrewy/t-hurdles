% 2021 12: new demo for intuition 
% moving from R code to matlab for matching formatting

nsim = 200;
nfac = 1e3;

red = [ 0.8500    0.3250    0.0980];
blue = [   0    0.4470    0.7410];
gray = [1 1 1]*0.3;

rng(1)

% parameters hand-selected so that pr_tgt2 = 0.5 and pr_tgt3_tgt2 = 0.5
% too lazy too code up numerical solver
tfalse = abs(normrnd(0,1, [nfac nsim]));
ttrue  = gamrnd(2,0.8, [nfac nsim])+1.65; 
tall = [tfalse ttrue];


clf; 

subplot(1,2,1)

    % prep
    edge = 0:0.5:5;
    center = edge(1:end-1) + mean(diff(edge))/2;
    n = histcounts(tfalse,edge);
    nf = n/nsim;

    % plot
    hb = bar(center,nf); hb.FaceColor = red;
    
    hv = vline(1.96,'-');  
    hv.LineWidth = 3;
    hv.Color = gray;
    text(2.1, 300, 't-hurdle = 1.96','fontsize',11)    
    
    % annotate
    x = [0.34 0.30]; 
    y = [0.29 0.20];
    annotation('textarrow',x,y,'String','FDR = 100%')    
    
    % format
    legend('False Factor','fontsize',11)
    setplot('Panel A','|t-statistic|','Number of Factors')       

    addspace(0.05)
    
subplot(1,2,2)
    
    % prep
    edge = 0:0.5:10;
    center = edge(1:end-1) + mean(diff(edge))/2;
    n = histcounts(tfalse,edge); nf = n/nsim;
    n = histcounts(ttrue,edge); nt = n/nsim;    

    % plot
    hb = bar(center,[nf; nt]', 'stacked');
    hb(1).FaceColor = red; hb(2).FaceColor = blue;
        
    hv = vline(1.96,'-');  
    hv.LineWidth = 3;
    hv.Color = gray;
    text(2.1, 300, 't-hurdle = 1.96','fontsize',11)
        
    x = [0.78 0.67]; 
    y = [0.33 0.19];
    annotation('textarrow',x,y,'String','FDR = 5%')
    
    % format
    legend('False Factor', 'True Factor','fontsize',11)
    setplot('Panel B','|t-statistic|','Number of Factors')    
    
    addspace(0.05)
    
    
set(gcf, 'position', [807.6667  501.8333  630.8333  255.0000])    
    

pr_tgt2 = sum(tall(:)>1.96)/(2*nsim*nfac)

pr_tgt3_tgt2 = sum(tall(:)>3) / sum(tall(:)>1.96)

export_fig('../output_exhibits/intro_fig.pdf')