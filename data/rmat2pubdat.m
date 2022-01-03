function [tpub,rpub,sigpub,corrpub,corrpubmat,nobspub] = rmat2pubdat(rmat)
% 2021 12 andrew: cleaning up copy pasta


nobspub = sum(~isnan(rmat));
rpub = nanmean(rmat);
sigpub = nanstd(rmat)./sqrt(nobspub);
tpub = rpub./sigpub;
corrpubmat = corr(rmat, 'rows', 'pairwise');
corrpub = unroll_tril(corrpubmat);