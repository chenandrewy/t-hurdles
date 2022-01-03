function tarv = parv2tarv(parv,namevec)
% converts transformaed parameters to straight paramaters
% 2019 03 Andrew

% by default don't do anything
tarv = parv;

% transform if necessary
for i = 1:length(parv)
    
    switch namevec{i}
        case 'sigsig'
            tarv(i) = log(parv(i)); % smoosh [0 inf] onto [-inf inf]
        case 'mu1'
            tarv(i) = log(parv(i)); % smoosh [0 inf] onto [-inf inf]            
    end % switch
    
end % for i