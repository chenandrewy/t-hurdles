function [parbase,meth] = baseline_settings()
% testing random seeding

% baseline settings
% other settings will only modify these

%% ==== USER ENTRY ====

%% seeeding

% controls obj_sim seed.  use nan to not reseed in obj
meth.seed = nan; 

% controls one_estimate seed.
meth.seed_est = 1; 


%% SMM target choices

% SMM
% momtype quant or hist
meth.momtype = 'Qplus';
meth.Qbreak = 0.1:0.1:0.9; 
    
% wtype = 'eye' 'eyepct', 'varinv', 
meth.wtype = 'varinv';

% mean median
meth.sim_agg = 'mean';

% choose to turn off classes of moments
meth.wswitch.t      = 1;
meth.wswitch.r      = 1;
meth.wswitch.sig    = 1; % need this!
meth.wswitch.corr   = 0; % estimating rho doesn't work well


%% optimization 

% 'fminunc', 'qnewton_chen'
meth.optcode = 'fminunc';

switch meth.optcode
    case 'fminunc'
        % set algo and tolerances
        % is 1e-3 enough? eyeballing the baseline point estiamte it looks
        % good enough        
        meth.optnewt = optimoptions('fminunc',...
            'algorithm','quasi-newton',...
            'maxfunctionevaluations',500,...
            'maxiterations',1000*5, ...
            'optimalitytolerance',1e-3, ...
            'steptolerance',1e-3, ...
            'display','none' ...
                );

        % uncomment for full feedback
            % meth.optnewt = optimoptions(meth.optnewt,'display','iter-detailed');
    
    case 'qnewton_chen'
        optset('qnewton','ShowIters',1)
        
end % switch meth.optcode
    
% meth.optgss = optimset('display','iter-detailed');            
meth.optgss = optimset('display','none');            
    
% other
% 'full','none'
meth.feed = 'full';  
meth.Nshrinkmax = 500;
meth.mugood = 'strict'; % percent or strict



%% other feedback
meth.newfig = 1;
% if meth.newfig; cls; end

% export fig automatically if on slurm
slurm_ppn = str2num(getenv('SLURM_CPUS_ON_NODE'));

if ~isempty(slurm_ppn)
    meth.exportme = 1;    
else
    meth.exportme = 0;
end

%% Which parameters to estimate and gridding


meth.p0list      = 0:0.05:1.0;
meth.thresh3list = [1/3 1/2 2/3];
meth.namevec = {'p0','mu2','musig','sigsig','thresh3'}; 


%% Simulation parameters

% these are the goal
meth.Nportall = 1000; 
meth.Nsim = 192; 

meth.boottype = 'parametric';

%% baseline parameters

% 'constant' 'constant_fab' or 'saved'
% saved loads results from main1_fitcorr.m
parbase.corrtype = 'saved'; 

parbase.p0 = 0.444;

parbase.mutype = 'gamma';
parbase.mu1 = 0; % mode
parbase.mu2 = 0.3183; % dispersion parameter
% parbase.mu2 = 1.555; % dispersion parameter
parbase.mu3 = 1; % tail parameter

parbase.slope_mu = 0;

parbase.threshtype  = 'step';
parbase.thresh1 = 1.50;  % low / midpoint
parbase.thresh2 = 2.57; % good / slope
parbase.thresh3 = 1/2; % frac missing / nan

parbase.musig = -1.7075;
parbase.sigsig = 0.4896;


