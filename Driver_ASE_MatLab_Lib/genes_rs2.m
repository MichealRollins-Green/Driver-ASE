function [mut,ase]=genes_rs2(T,M,v2g,vars)

% [genes,r]=genes_rs(T,S,V,var2gene,look,vars)
% vars are variants to consider for analysis can come from V
% mut is the mutation sum matrix for each gene and ase is the downsized ase
% matrix

% 1 - collapse vars to genes (may be more than one gene/var)
% downsize M to be more manageable
m=M;
[m.rowlabels,oo]=intersect(M.rowlabels,vars);
m.data=M.data(oo,:);
m.syn=m.syn(oo);

% which genes have any vars at all?
% collapse v2g to genes with relevant vars
[v2g.var,aa]=intersect(v2g.var,vars);
if (sum(aa)==0)
    warning('No intersect between mutations and v2g mutations')
    mut=[];
    ase=[];
else
v2g.data=v2g.data(aa,:);
off=sum(v2g.data)>0;
v2g.gene=v2g.gene(off);
v2g.data=v2g.data(:,off);

% collapse v2g to vars thar are disruptive
% disrptv=s.rowlabels(s.syn>1);
% [v2g.var,aa]=intersect(v2g.var,disrptv);
% v2g.data=v2g.data(aa,:);
% off=sum(v2g.data)>0;
% v2g.gene=v2g.gene(off);
% v2g.data=v2g.data(:,off);

%
[~,tt]=size(M.data);
mut.collabels=M.tcga;
mut.data=zeros(length(v2g.gene),tt);

% go through each gene in gg, and pull sum of mutations from S
mut.gene=v2g.gene;

% update m and v2g to be aligned on vars
[zz,aa,bb]=intersect(v2g.var,m.rowlabels);
m.rowlabels=zz;
v2g.var=zz;
m.data=m.data(bb,:);
v2g.data=v2g.data(aa,:);

% collapse mutation data - only consider non-syn mutations
for n=1:length(v2g.gene)
    my_dd=full(m.data(logical(v2g.data(:,n)),:));
    mut.data(n,:)=sum(my_dd,1);
end

% mut has gene and mutations, make ase, reorder both to be consistent
[cols,aa,bb]=intersect(mut.collabels,T.collabels);
mut.collabels=cols;
mut.data=mut.data(:,aa);
ase.collabels=cols;
ase.data=T.lr(:,bb);
ase.cov=T.cov(:,bb);
ase.p=T.p_true_ase(:,bb);
ase.gene=T.tx;

% same genes
[~,aa,bb]=intersect(mut.gene,ase.gene);
mut.gene=mut.gene(aa);
mut.data=mut.data(aa,:);
ase.gene=ase.gene(bb);
ase.data=ase.data(bb,:);
ase.cov=ase.cov(bb,:);
ase.p=ase.p(bb,:);
end


