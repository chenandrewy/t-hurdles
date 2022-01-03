function htest = test_hypotheses(sim,mugood,sided)
% 2019 03 Andrwe
% finds hypo testing stats for cell-struct versions of simulated data

% this one is more careful to find stats in every sim and then takes
% medians

% it also doesn't accept cells, at least not right now

if nargin <= 2 % 2-sided by default! updated 2021 12 (prev 1-sided)
    sided = 2;
end

% setup
htest.mugood = mugood;
htest.t = 0:0.01:14;


switch sided
    case 1
        t2 = sim.t;
        h01 = -norminv(0.01);
        h05 = -norminv(0.05);
    otherwise
        t2 = abs(sim.t);
        h01 = -norminv(0.01/2);
        h05 = -norminv(0.05/2);        
end

%% fdrs among all

% find false discovery proportion
if ~ischar(mugood)
    false = sim.mu <= mugood;
elseif strcmp(mugood,'strict')
    false = sim.null;
end

[fdp,fwe,mu,r,drate] = nanall(size(false,1),length(htest.t));
for i = 1:length(htest.t)   
    
    discovery = t2 > htest.t(i);            
                
    false_and_disc = sum( false & discovery, 2 );
    Ndiscovery = sum(discovery,2);
    
    
    % following BH 1995: FDP = 0 if no discoveries
    if Ndiscovery > 0 
        fdp(:,i) = 100*false_and_disc./Ndiscovery;          
    else
        fdp(:,i) = 0;        
    end
    
    
    % family wise error
    fwe(:,i) = any(false_and_disc,2);    
    mu(:,i) = sum( discovery.* sim.mu, 2 )./Ndiscovery;    
    r(:,i) = sum( discovery.* sim.r, 2 )./Ndiscovery;        
    
end

% fdr is E[fdp], fwer is E[fwe] = prob(fdp>0)
fdr = mean(fdp,1); % already in pct
fwer = mean(fwe,1)*100; % need to make pct

mu_mean = mean(mu,1);
r_mean = mean(r,1);

% thurdles
t_fdr01 = htest.t(find(fdr<01,1));
t_fdr05 = htest.t(find(fdr<05,1));
t_fwer01 = htest.t(find(fwer<01,1));
t_fwer05 = htest.t(find(fwer<05,1));

% discovery rate
drate01 = mean(mean(t2 > h01,2));
drate05 = mean(mean(t2 > h05,2));



%% fdr among published


pubfdp = nanall(size(sim.mu,1),1);
for simi = 1:length(sim.tpub)
    if ~ischar(mugood)
        Nfalse = sum( sim.mupub{simi} <= mugood );
    elseif strcmp(mugood,'strict')
        Nfalse = sum( sim.nullpub{simi} );
    end
        
    pubfdp(simi) = 100*Nfalse/length(sim.mupub{simi});
end

pubfdr = mean(pubfdp,1);



htest = add2struct(htest,fdp,fwe,fdr,fwer...
    ,mu_mean,r_mean ...
    ,t_fdr05,t_fdr01,t_fwer05,t_fwer01 ...
    ,pubfdp,pubfdr ...
    ,drate01,drate05);

