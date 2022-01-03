function plotdat = plot_obj(meth,par,booti,Nplot,parvecname)
% Nplot = 21 works
% parvecname =  {'p0','mu2','musig','sigsig','thresh3'} works


[samp,par,corrmat,meth] = auto_prep(meth,par,booti);

% set up limits of parameter space
clear parlim
for pari = 1:length(parvecname)
    switch parvecname{pari}
        case 'p0'
            parlim(pari,:) = [0 0.8];
        case 'thresh3'
            if strcmp(par.threshtype,'step')
                parlim(pari,:) = [1/4 3/4];            
            else
                parlim(pari,:) = meth.thresh3list([1 end]);
            end
        otherwise
            parlim(pari,1) = par.(parvecname{pari})*0.7;
            parlim(pari,2) = par.(parvecname{pari})*1.3;
    end
end


% plot
Nsubrow = ceil(length(parvecname)/3);
for pari = 1:length(parvecname)
    temppar = par;
    
    cparlist = linspace(parlim(pari,1),parlim(pari,2),Nplot-1);
    cparlist = sort([cparlist temppar.(parvecname{pari})]); % add optimium to list

    % evaluate objective
    objlist = nan(Nplot,1);
    for vali = 1:Nplot                

        % transform to real line 
        tarvec = parv2tarv(cparlist(vali),{parvecname{pari}});
        
        % evalute objective!
        objlist(vali) = obj_sim(tarvec,{parvecname{pari}},samp,temppar,corrmat,meth);
        
    end

    if length(parvecname)>1
        Nsub = 3;
    else 
        Nsub = 1;
    end
     subplot(Nsubrow,Nsub,pari)
    
    plot(cparlist,objlist,'x-')
    xlabel(parvecname{pari}); ylabel('GMM obj');
    
    vline(temppar.(parvecname{pari}))
    
    % save for output
    plotdat{pari}.val = cparlist;
    plotdat{pari}.obj = objlist;
    plotdat{pari}.name = parvecname{pari};    
    
end
