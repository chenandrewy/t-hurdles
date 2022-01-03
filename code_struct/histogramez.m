function h = histogramez(x,edge)

if nargin <= 1
    h = histogram(x,'normalization','probability');
else
    h = histogram(x,edge,'normalization','probability');
end

