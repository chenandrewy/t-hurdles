% just desperately trying to get lln to work
% Andrew 2021 12

Nport = 500;
Nsimmax = 5000;
scale = 0.5;



clf;

for i = 1:20
    
    mu = exprnd(scale, Nport, Nsimmax);

    medmu = median(mu,1);
    meanmedmu = cumsum(medmu)./(1:Nsimmax);

    meanmu = mean(mu,1);
    meanmeanmu = cumsum(meanmu)./(1:Nsimmax);


    subplot(2,1,1); hold on;
    plot(meanmedmu)
    medtheory = log(2) * scale;
    hline(medtheory);
    ylim(medtheory + [-0.01 0.01])


    subplot(2,1,2); hold on;
    plot(meanmeanmu)
    medtheory = scale;
    hline(medtheory);
    ylim(medtheory + [-0.01 0.01])

end

%%

