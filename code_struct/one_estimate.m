function [est,objgrid] = one_estimate(meth,parbase,booti)
% 2018 09 Andrew
% estimates model by smm 
% optimizes parameterts in namevec for each value of p0list, thresh3list, then chooses
% the best among that optimization

rng(meth.seed_est, 'combRecursive') 
if isnan(meth.seed)
    % set obj_seed randomly, but reprodroducibly
    meth.seed = randi(50000)*booti; 
end


% prepare sample, corrmat, etc
[samp,parbase,corrmat,meth] = auto_prep(meth,parbase,booti);


%% create grids for p0 and thresh3 if desired
if any(strcmp(meth.namevec,'p0'))
    p0list = meth.p0list;
else
    p0list = parbase.p0;
end

if any(strcmp(meth.namevec,'thresh3'))
    thresh3list = meth.thresh3list;
else
    thresh3list = parbase.thresh3;
end

%% create list of names for smooth opt

namesmooth = meth.namevec;
namesmooth(strcmp(namesmooth,'p0')) = [];
namesmooth(strcmp(namesmooth,'thresh3')) = [];
    
    
%% set up rolled out double loop
% rolling out helps with parfor, I hope

p0ilist = (1:length(p0list))';
thresh3ilist = (1:length(thresh3list))';
thresh3p0ilist = allcomb(thresh3ilist,p0ilist);

%% estimate using p0 grid and quasi newton

% rolled out double loop
objlist = nan(size(thresh3p0ilist,1),1);
parvechatlist = nan(length(thresh3p0ilist),length(namesmooth));


parfor thresh3p0i = 1:size(thresh3p0ilist,1)    
    
    % load up parameters (thresh3 and p0)
    thresh3i    = thresh3p0ilist(thresh3p0i,1);
    p0i         = thresh3p0ilist(thresh3p0i,2);       

    temppar = parbase;
    temppar.thresh3 = thresh3list(thresh3i);
    temppar.p0 = p0list(p0i);


    % set up optimization for parameters besides thresh3 and p0 
    parvec = nan(size(namesmooth));
    for pari = 1:length(namesmooth)
        parvec(pari) = temppar.(namesmooth{pari});
    end
    
    % transform to real line 
    tarvec = parv2tarv(parvec,namesmooth);
    
    % optimize    
    minme = @(tarvec)  obj_sim(tarvec,namesmooth,samp,temppar,corrmat,meth);                
    if isfinite(minme(tarvec))
        % initial guess is OK, so optimize        
        [tarvechat,obj,newtflag,newtoutput] = fminunc(minme,tarvec,meth.optnewt);
    else
        % initial guess is bad, probably p0 is too high
        tarvechat = tarvec;
        obj = inf;
    end    
    
    % transform back
    parvechat = tarv2parv(tarvechat,namesmooth); 

    % save
    objlist(thresh3p0i) = obj;
    parvechatlist(thresh3p0i,:) = parvechat;
    

    if strcmp(meth.feed,'full')
        disp(['p0 = ' num2str(temppar.p0)  ' thresh3 = ' num2str(temppar.thresh3) ' obj = ' num2str(obj) ])
    end
    
    
end

%% optimize over grids and save

% optimize over grid
[err,thresh3p0ihat] = min(objlist);

thresh3i    = thresh3p0ilist(thresh3p0ihat,1);
p0i         = thresh3p0ilist(thresh3p0ihat,2);

% store in struct
par = parbase;
par.thresh3 = thresh3list(thresh3i);
par.p0 = p0list(p0i);

for pari = 1:length(namesmooth)
    par.(namesmooth{pari}) = parvechatlist(thresh3p0ihat,pari);
end
    


%% roll up for easy interpretation
for thresh3p0i = 1:length(objlist)
    thresh3i    = thresh3p0ilist(thresh3p0i,1);
    p0i         = thresh3p0ilist(thresh3p0i,2);    
    objgrid(thresh3i,p0i) = objlist(thresh3p0i);
end



est = packstruct(par,err);

%% debug
% figure(99)
% clf; hold on;
% plot(p0list,objgrid,'x-')
% ylim([0 400])
% setplot('obj at different thresh3','p0','err')


