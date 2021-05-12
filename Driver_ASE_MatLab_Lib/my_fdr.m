function fdr=my_fdr(pp,p)
% p should be matrix of (iter,genes_p_or_r_info);
% pp is the real data but p is the permutated data;
% Make sure p and pp are log10P values;

%[ss,X]=size(p);
p_s=p;

%Only used in old code;
%for n=1:ss
%    p_s(n,:)=sort(p(n,:),'descend');
%end

fdr=ones(size(pp));

%new code based on tomas's code;
tt0=p_s(:);
tt=tt0(~isnan(tt0));%exclude these nan values;

%Need to sort pp descendly for new code!;
[pps,p_idx]=sort(pp,'descend');

%Just get the largest FDR in each row for calculating FDR;
for n=1:length(pp)
    %fdr(n)=sum(p_s(:,1)>pp(n))/ss;
    
    %Compare real P with permuted P values genome-wide;
    %fdr(n)=sum(p_s(:)>pp(n))/(X*ss);
    
    %updated based on Tomas's code;
    %Align new fdr with the order of original P values;
    fdr(p_idx(n))=sum((tt>=pps(n))/length(tt))/(n/length(fdr));
end


