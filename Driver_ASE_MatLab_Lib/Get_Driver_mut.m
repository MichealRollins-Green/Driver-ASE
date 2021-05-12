function [out,Mut]=Get_Driver_mut(iter,look,mut_ase_look)
Mut.data=[];
Mut.rowlabels=[];
Mut.collabels=[];
Mut.anno_type=[];
Mut.fdr=[];
Mut.fm=[];
Mut.p=[];
for n=1:18
  if (exist(['mut_ase_auto/mut_ase_' num2str(n) '.mat'],'file')==2)
    fprintf(['working on ' num2str(n) '\n']);
    eval(['load mut_ase_auto/mut_ase_' num2str(n) ';']);
    [a,m]=subset_ase_mut(ase,mut,0,10,50,0,look);
    m.data(m.data>0)=1;
    [~,pp]=my_corr(a,m);
    [~,p_shuff]=my_corr_shuff(a,m,iter);
    fdr=my_fdr(pp,p_shuff);
    fm=fdr_mut(m);
    o=m;
    o.fdr=fdr;
    o.fm=fm;
    o.pp=pp;
    %Change cutoff if only want to keep fdr<=0.8;
    fdr_Point8idx=find(fdr<=1);
    %o.anno_type=mut_ase_look.name(n);
    sn=length(fdr_Point8idx);
    o.anno_type=repmat(mut_ase_look.type(n),sn,1);
    Mut.fdr=[Mut.fdr;o.fdr(fdr_Point8idx)];
    Mut.p=[Mut.p;o.pp(fdr_Point8idx)];
    Mut.fm=[Mut.fm;o.fm(fdr_Point8idx)];
    
    
    %Get mut data;
    Mut.data=[Mut.data;m.data(fdr_Point8idx,:)];
    gene=vlookup_list(m.gene(fdr_Point8idx),look.tx,look.gene);
    Mut.rowlabels=[Mut.rowlabels;gene];
    Mut.anno_type=[Mut.anno_type;o.anno_type];
    Mut.collabels=m.collabels;
    %Output all info for later usage;
    out(n)=o;
  else
    fprintf(['Not exist for ' num2str(n) '\n']);
  end
end
 Mut.gene=Mut.rowlabels;
 Mut.rowlabels=strcat(Mut.rowlabels,'-',Mut.anno_type);