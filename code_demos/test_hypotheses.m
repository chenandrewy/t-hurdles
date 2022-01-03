function htest = test_hypotheses(sim,mugood)
% 2019 03 Andrwe
% finds hypo testing stats for cell-struct versions of simulated data

% this one is more careful to find stats in every sim and then takes
% medians

% it also doesn't accept cells, at least not right now


% setup

htest.mugood = mugood;
htest.t = 0:0.01:10;

%% fdrs among all

% find false discovery proportion
if ~ischar(mugood)
    false = sim.mu <= mugood;
elseif strcmp(mugood,'strict')
    false = sim.null;
end

[fdp,fwe] = nanall(size(false,1),length(htest.t));
for i = 1:length(htest.t)
    
    discovery = sim.t > htest.t(i);
        
%         discovery = abs(sim.t) > htest.t(i); % testing
    
    false_and_disc = sum( false & discovery, 2 );
    Ndiscovery = sum(discovery,2);
    fdp(:,i) = 100*false_and_disc./Ndiscovery;  
    
    % family wise error
    fwe(:,i) = any(false_and_disc,2);

    
end

% fdr is E[fdp], fwer is E[fwe] = prob(fdp>0)
fdr = mean(fdp,1); % already in pct
fwer = mean(fwe,1)*100; % need to make pct

% thurdles
t_fdr01 = htest.t(find(fdr<01,1));
t_fdr05 = htest.t(find(fdr<05,1));
t_fwer01 = htest.t(find(fwer<01,1));
t_fwer05 = htest.t(find(fwer<05,1));

t_fdr01_2s = htest.t(find(fdr<01/2,1));
t_fdr05_2s = htest.t(find(fdr<05/2,1));


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
    ,t_fdr05,t_fdr01,t_fwer05,t_fwer01 ...
    ,t_fdr05_2s,t_fdr01_2s ...
    ,pubfdp,pubfdr);

