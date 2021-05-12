function ss=gt_plot_tumor_only(t,min_cov,pbinom_adj)

% min cov and min pbinom_adj
% pbinom_adj: -log(0.01)/log(10)=2 and -log(0.05)/log(10)=1.3
mc=(t.a+t.b)>=min_cov & abs(t.p_binom)>pbinom_adj;

% distribution
aa=sum(mc,2);
[~,ll]=size(mc);

ss=zeros(size(ll));

for n=1:ll
    ss(n)=sum(aa>n);
end

% calculate the total of non-NAN value, i.e., 
% only consider sample with het geno;
mhet=(t.a+t.b)>=min_cov & isfinite(t.p_binom);
% distribution

