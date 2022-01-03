function [momt,momr,momsig,momcorr,Npub] = momfun(tpub,rpub,sigpub,corrpub)
% find moments 2019 03
% corrpub can be []

% 2021 12: added abs transform 
% philosophy is to do abs() as late as possible for generality

pbreak = 10:10:90;

tpub = abs(tpub);
rpub = abs(rpub);

Npub = length(tpub); 


pmarginal = sum(tpub >= 1.50  & tpub <= 2.57)/Npub*100;
momt = ptile(tpub,pbreak);

% combine
momt = [pmarginal momt];        

momsig = ptile(sigpub,pbreak);
momr = ptile(rpub,pbreak);     
momcorr = ptile(corrpub,pbreak);

        
