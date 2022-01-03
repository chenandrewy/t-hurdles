function parvec = par2parvec(par,namevec)

parvec = nan(size(namevec));
for i = 1:length(namevec)
    parvec(i) = par.(namevec{i});    
end