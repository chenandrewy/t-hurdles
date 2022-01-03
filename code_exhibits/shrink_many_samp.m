function shrinkpar = shrink_many_samp(par,simpart,Nshrinkmax,pubonly)
% Andrew 2019 03 
% applies shrink_samp to a struct vector or a cell array of struct vectors

if nargin <= 3
    pubonly = 1;
end

% Nshrinkmax is used to keep the code from taking too long

if iscell(simpart)
    for pari = 1:length(simpart)
        
        switch pubonly
            case 1
                r   = [simpart{pari}.rpub];
                sig = [simpart{pari}.sigpub];
            case 0
                r   = [simpart{pari}.r];
                sig = [simpart{pari}.sig];
                
        end
        shrinkpar{pari} = shrink_samp(par{pari},r,sig,Nshrinkmax);
        
    end % for pari
    
else % iscell
    

    switch pubonly
        case 1
            r   = [simpart.rpub];
            sig = [simpart.sigpub];
        case 0
            r   = [simpart.r];
            sig = [simpart.sig];

    end
    shrinkpar = shrink_samp(par,r,sig,Nshrinkmax);    
    
end % iscell

