function par = parvec2par(parvec,namevec,par)

for i = 1:length(parvec)
    par.(namevec{i}) = parvec(i);
end
