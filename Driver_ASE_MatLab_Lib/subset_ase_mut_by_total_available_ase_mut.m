function [a,m]=subset_ase_mut_by_total_available_ase_mut(ase,mut,coding,min_samples_powered,at_x_reads,varargin)

% [a,m]=subset_ase_mut(ase,mut,coding,min_samples_powered,at_x_reads,look)

if ~isempty(varargin)
  look=varargin{1};  
end

% OLD CODE
% % powered ase and mutation
% x=ase.cov;
% s=sum(x>=at_x_reads,2);
% off=find(s>=min_samples_powered);
% 
% % min muts
% x=mut.data;
% s=sum(x>0 & isfinite(ase.data),2);
% off2=find(s>=min_mut_ase);
% off=intersect(off,off2);

% coding
% if coding==1
%     [X,I]=vlookup_list(mut.gene,look.tx,look.tx);
%     x=find(look.coding(I));
%     off=intersect(off,x);
% end



% NEW CODE

% ASE powered
x=ase.cov;
off_ase=x>=at_x_reads;

% mut powered
off_mut=logical(mut.data);

% min samples
off=sum(off_ase&off_mut,2)>=min_samples_powered;

% coding
if (coding==1)
    [~,I]=vlookup_list(mut.gene,look.tx,look.tx);
    x=find(look.coding(I));
    off=intersect(off,x);
end


a=ase;
m=mut;

a.data=a.data(off,:);
a.cov=a.cov(off,:);
a.p=a.p(off,:);
a.gene=a.gene(off);

m.data=m.data(off,:);
m.gene=m.gene(off);

% also dump measurements that don't meet threshold
a.data(a.cov<at_x_reads)=nan;


