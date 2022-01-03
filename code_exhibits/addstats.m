function [sim,htest,shrink,simflat,emom] = addstats(samp,par,corrmat,meth)
% find additiional stats.  Just because sometimes you don't want to
% estimate and everything too

if ~isfield(par,'slope_mu')
    par.slope_mu = 0;
end


[~,sim,simflat,emom] = obj_sim([],[],samp,par,corrmat,meth);


% htest
htest = test_hypotheses(sim,meth.mugood);


% shrinkage
shrinksamp = shrink_samp(par,samp.rpub,samp.sigpub,meth.Nshrinkmax);
% shrinksim  = shrink_samp(par,[sim.rpub{:}],[sim.sigpub{:}],meth.Nshrinkmax);


clear shrink
shrink.samp.s    = shrinksamp.s;
shrink.samp.muhat    = shrinksamp.muhat;
shrink.samp.mean = mean(shrinksamp.s);
shrink.samp.median = median(shrinksamp.s);
shrink.samp.smooth = 100*(1-mean(shrinksamp.muhat)./mean(samp.rpub));

% shrink.sim.s    = shrinksim.s;
% shrink.sim.muhat    = shrinksim.muhat;
% shrink.sim.mean = mean(shrinksim.s);
% shrink.sim.median = median(shrinksim.s);

% prepare sim's cell arrays the easy way
% temprpub = [sim.rpub{:}];
% tempNmax = min(length(temprpub),meth.Nshrinkmax);
% tempmupub = [sim.mupub{:}];

% shrink.sim.smooth = 100*(1-mean(shrinksim.muhat)./mean(temprpub(1:tempNmax)));
% shrink.sim.smoothez = 100*(1-mean(tempmupub)/mean(temprpub)); % note this is always fast
