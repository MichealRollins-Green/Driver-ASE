function ase_mut_assoc=extract_ase_mut(look,n,fdr_iter)

% clear directories
mkdir('ase_mut_genes/');
warning('off','all');
fprintf(['working on ' num2str(n) '\n']);
eval(['load mut_ase_auto/mut_ase_' num2str(n) ';']);

[a,m]=subset_ase_mut(ase,mut,0,10,50,0,look);
m.data(m.data>0)=1;
ase_mut_assoc.ase_data=a.data;
ase_mut_assoc.mut_data=m.data;
ase_mut_assoc.collabels=a.collabels;
ase_mut_assoc.tx=a.gene;
ase_mut_assoc.gene=vlookup_list(a.gene,look.tx,look.gene);
ase_mut_assoc.ase_trueASEp=a.p;

%my_corr and my_corr_shuff will generate -log10(p value);
%The above two functions will change ase data, which is log2(Allele_A_counts/Allele_B_counts);
%a.data=abs(a.data);
%a.data(a.data>10)=10;

[X,pp]=my_corr(a,m);
[X,p_shuff]=my_corr_shuff(a,m,fdr_iter);

ase_mut_assoc.assoc_fdr=my_fdr(pp,p_shuff);
ase_mut_assoc.assoc_p=pp;
ase_mut_assoc.fm=fdr_mut(m);

