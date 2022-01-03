function many_estimates(set_name,bootistart,bootiend)


%% estimate and bootstrap
for booti = bootistart:bootiend
    
    if booti == 0
        disp('point estimating')
    else    
        disp(['bootstrap number ' int2str(booti)])
    end
    
    %% settings
    
    % load baseline
    [parbase,meth] = baseline_settings;

    % modify 
    eval(set_name)    
    
    %% estimate using p0 grid and quasi newton

    % feedback
    datetime('now')
    disp(['estimating ' meth.name])
 
    tic
    [est,objgrid] = one_estimate(meth,parbase,booti);
    sec = toc
    hrs = sec/60/60    

    disp('estimated parameters are')
    est.par
    err = est.err
    

    %% save
    
    % make directories
    dirname = ['../output/' meth.name '/'];
    mkdir(dirname)


    disp('saving to ')
    dirname

    if booti == 0

        %% calculate rescaled model 

        % simulate to find scaling
        [samp,~,corrmat,meth] = auto_prep(meth,parbase,booti);
        [~,~,simflat] = obj_sim([],[],samp,est.par,corrmat,meth);    
        rescaled.pubrate = mean(simflat.Npub/meth.Nportall);
        rescaled.Nportall = round(meth.Nportall/rescaled.pubrate);    
        rescaled.corrmat = random_corr(rescaled.Nportall,est.par.rho1,est.par.rho2,est.par.rho3,1);    

        %% save
        save([dirname 'meth'],'meth')
        save([dirname 'est'],'est')
        save([dirname 'corrmat'],'corrmat')    
        save([dirname 'rescaled'],'rescaled')
        copyfile([meth.name '.m'] , dirname)    

    else

        dirname = ['../output/' meth.name '/boot/'];
        mkdir(dirname)

        filename = [dirname 'boot_est' sprintf('%04d',booti) '.mat'];
        save(filename,'est')

    end
        
end


