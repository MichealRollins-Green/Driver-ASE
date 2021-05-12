function p=TA(aa,bb)

% p=TA(a_matrix,b_matrix) i.e. (allele1,allele2)
% m is a precomputated matrix (1000*1000) for ....;
% In the dialog, values are all 0.5;
% other values are equal for row(a+1) and col(b+1);
load m;
[ss,tt]=size(aa);

p=nan(ss,tt);

for ii=1:ss
    for jj=1:tt
        a=aa(ii,jj);
        b=bb(ii,jj);
        
        if isfinite(a)&&isfinite(b)
            % downsample
            if a>1000||b>1000
                rr=max([a;b])/1000;
                a=round(a/rr);
                b=round(b/rr);
            end
            p(ii,jj)=m.data(a+1,b+1);
        end
    end
end

