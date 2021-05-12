function [rr,ppr]=my_corr_shuff(ase,mut,iter,min_samples_powered,min_cov,varargin)

if ~isempty(varargin)
    coding_only=varargin{1};
else
    coding_only=0;
end

%r and p will be matrix of (iter,genes_p_or_r_info);
%updated based on tomas's new code;
[~,m]=subset_ase_mut_by_total_available_ase_mut(ase,mut,coding_only,min_samples_powered,min_cov);
[sc,ss]=size(m.data);
rr=zeros(sc,iter);
%ppr=zeros(sc,iter);
ppr=nan(sc,iter);%will exclude these nan values when calculating FDR;
hh=0;

for ii=1:iter  
    [~,Ir]=sort(rand(ss,1));%sort column-wide, instead of genome-wide;
    mutr=downSizeTo(mut,Ir,ss);
    [ar,mr]=subset_ase_mut_by_total_available_ase_mut(ase,mutr,coding_only,min_samples_powered,min_cov);
     mr.data(mr.data>0)=1;
     [rr_temp,pp_temp]=my_corr(ar,mr);

    if length(pp_temp)> sc
      ppr(:,ii)=pp_temp(1:sc);
      rr(:,ii)=rr_temp(1:sc);
    else
      ppr(1:length(pp_temp),ii)=pp_temp;
      rr(1:length(pp_temp),ii)=rr_temp;
    end
    
    if hh==10
    fprintf('.');
    hh=0;
    end
    hh=hh+1;    
    
end
fprintf('');%print a blank line to separe the above output with new result!;
