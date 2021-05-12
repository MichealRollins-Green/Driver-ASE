function fdr = fdr_mut(mut)

[ss,tt] = size(mut.data);
s = sum(mut.data,2);
fdr = ones(ss,1);

for n=1:ss
    fdr(n) = sum(s >= sum(mut.data(n,:))) / length(s);
end
