function corrmat = random_corr(d,mu,sigma,nu,seed)
% 2019 03 Andrew
% copied from
% https://stats.stackexchange.com/questions/124538/how-to-generate-a-large-full-rank-random-correlation-matrix-with-some-strong-cor
% and based on Lewandowski, Kurowicka, and Joe 2009

% generates random correlation matrix (d x d) based on vines and extended onion method
% partial correlations are beta distributed on -1 1 with parameters alpha =
% beta = 1

% Actually, it seems this is vased on the C-vine (page 9, Section 2). The Onion method (3.2) is apparently much faster (Table 5)

% 2021 12 Andrew: updated to use truncated location scale t distribution
% instead of beta, since this fits the new data better (fat tails in the
% correlations)



rng(seed,'CombRecursive');

% pdist = makedist('tLocationScale','mu',0.05,'sigma',0.3,'nu',0.001);
pdist = makedist('tLocationScale','mu',mu,'sigma',sigma,'nu',nu);
pdist = truncate(pdist,-1,1);
P = random(pdist, [d d]);

%     % old code
%     P = betarnd(betaparam,betaparam, [d d]); %// sample partial correlations from beta
%     P = (P-0.5)*2;     %// linearly shifting to [-1, 1]

corrmat = eye(d);
for k = 1:d-1
    for i = k+1:d        
        p = P(k,i);
        for l = (k-1):-1:1 %// converting partial correlation to raw correlation
            p = p * sqrt((1-P(l,i)^2)*(1-P(l,k)^2)) + P(l,i)*P(l,k);
        end
        corrmat(k,i) = p;
        corrmat(i,k) = p;
    end
end


% 2021 12: since t-dist is no longer permutation invariant (see Joe),
% there's no point to this and this can confuse people.  Not that there was
% really a point to permutation invariance anyway

    % permuting the variables to make the distribution permutation-invariant
    % permutation = randperm(d);
    % corrmat = corrmat(permutation, permutation);
