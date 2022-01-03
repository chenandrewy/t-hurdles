function parv = tarv2parv(tarv,namevec)
% converts transformaed parameters to straight paramaters
% 2019 03 Andrew

% didn't do nuthin by deafult
parv = tarv;

for i = 1:length(tarv)
    
    switch namevec{i}
        case 'sigsig'
            parv(i) = exp(tarv(i)); % make sigsig positive          
        case 'mu1'
            parv(i) = exp(tarv(i)); % mean alt return should always be pos
    end % switch
    
end % for i