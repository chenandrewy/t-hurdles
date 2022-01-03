function shrink = shrink_samp(par,r,sig,Nshrinkmax)
% 2018 09 Andrew: shrink sample returns
    % alternatively you can use something like:
%     kernel_dirac = @(x) normpdf(rcurr,x,sigcurr).*par.p0.*dirac(x);
%     denom = double(int(kernel_dirac(x),-inf,inf)) + integral(kernel_exp,-inf,inf);


if ~isfield(par,'slope_mu')
    par.slope_mu = 0;
end

% precalculate location parameter b
%   updated with slope_mu 2019 03
switch par.mutype
    case 'gamma'           
        b =  - (max(par.mu3,1)-1)*par.mu2 + par.mu1 + par.slope_mu*sig;       
    case 't'        
        b = par.mu1 + par.slope_mu*sig;
    case 'gamma_noncen'           
        b =  par.mu1 + par.slope_mu*sig;               
end


%%

Nshrink = min(length(r),Nshrinkmax);

[s,muhat] = nanall([Nshrink 1]);
for porti = 1:Nshrink

    % turn off unnecessary warning about kernel_pos being overwritten
    warning('off','MATLAB:mir_warning_maybe_uninitialized_temporary');

    rcurr = r(porti); 
    sigcurr = sig(porti);

    % define kernel pos (part of integral kernal that does not correspond to
    % dirac term).
    switch par.mutype
        case 'gamma'                       
            kernel_pos = @(mu) normpdf(rcurr,mu,sigcurr).*(1-par.p0) ...
                .*gampdf(mu-b(porti),par.mu3,par.mu2);              
        case 't'        
            kernel_pos = @(mu) normpdf(rcurr,mu,sigcurr).*(1-par.p0) ...
                .*(1/par.mu2).*tpdf((mu-b(porti))/par.mu2,par.mu3);                              
        case 'gamma_noncen'                       
            kernel_pos = @(mu) normpdf(rcurr,mu,sigcurr).*(1-par.p0) ...
                .*gampdf(mu-b(porti),par.mu3,par.mu2);                 
    end    


    % denominator: normpdf comes from dirac term
    denom = normpdf(rcurr,0,sigcurr).*par.p0 + integral(kernel_pos,-inf,inf);
    % numerator note the dirac term doesn't matter here, since mu appears
    % in the kernel, and the dirac term = 0 unless mu = 0
    kernel_pos_mu = @(mu) mu.*kernel_pos(mu);
    num = integral(kernel_pos_mu,-inf,inf);

    muhat(porti) = num/denom;
    s(porti) = 100*(1-muhat(porti)/rcurr);

end

shrink.muhat = muhat;
shrink.s = s;



%% find stats

shrink.mean_s = mean(shrink.s);
shrink.median_s = median(shrink.s);
shrink.smooth_s = 100*(1-mean(shrink.muhat)./mean(r));
shrink.mean_muhat = mean(shrink.muhat);
shrink.median_muhat = median(shrink.muhat);

