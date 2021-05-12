function [r,pp]=my_corr(ase,mut)

%if calculating ASE with formula: ASE=|log2(T.a/T.b)|, this is necessary;
%In the manuscript, we used the formula ASE=|2*(T.a/(T.a+T.b))-1|, which is
%always less or equal to 1, thus the following code will not affect ase
%data.;
ase.data=abs(ase.data);
ase.data(ase.data>10)=10;

[ss,~]=size(ase.data);
r=nan(ss,1);
pp=nan(ss,1);

for n=1:ss
    [r(n),pp(n)]=my_nancorr(ase.data(n,:)',mut.data(n,:)');
end

%disp(sum(r<0));
%pp(r<0)=1;

pp(isnan(pp))=1;
pp=-log10(pp);
pp(~isfinite(pp))=1;
