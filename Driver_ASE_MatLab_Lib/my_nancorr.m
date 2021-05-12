function [c,p]=my_nancorr(a,b,varargin)
%Replace default pearson correlation p value with ranksum p value;
NotMissing=isfinite(a)&isfinite(b);
aa=a(NotMissing);
bb=b(NotMissing);

if isempty(aa)||isempty(bb)
    c=nan;
    p=nan;
elseif sum(bb(~isnan(aa)))< 3 % not really enough mutations and ase to make a call
    c=nan;
    p=nan;
else
    %if ~isempty(varargin)
     if ~isempty(varargin)
        %Just get correlation R2 but not pearson correlation p;
        %[c,~]=corr(aa,bb,varargin{1},varargin{2});
        [c,p]=corr(aa,bb);
        %use ranksum to replace pearson correlation P value;
        %[p,~]=ranksum(aa(bb==1),aa(bb==0));
        %[~,p]=corr(aa,bb);
    else
        [c,p]=corr(aa,bb);
    end
end


