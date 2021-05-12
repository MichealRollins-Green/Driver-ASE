function [out,Mut,features]=Get_Driver_mut_regex4feature(iter,look,mut_ase_look,regex4feature)
Mut.data=[];
Mut.rowlabels=[];
Mut.tx=[];
Mut.collabels=[];
Mut.anno_type=[];
Mut.fdr=[];
Mut.fm=[];
Mut.p=[];

mut_ase_look_tmp=FieldExtractor(mut_ase_look,'type',0,0,regex4feature);
if isempty(mut_ase_look_tmp.type)
    disp('')
    error(['No match for ' regex4feature ,'\n']);
end
features=mut_ase_look_tmp;
for n=1:length(mut_ase_look_tmp.type)
  mut_ase_n=mut_ase_look_tmp.num(n);
  if (exist(['mut_ase_auto/mut_ase_' num2str(mut_ase_n) '.mat'],'file')==2)
    disp('')
    fprintf(['\nworking on ' num2str(mut_ase_n) ' for ' char(mut_ase_look_tmp.type(n)) '\n']);
    eval(['load mut_ase_auto/mut_ase_' num2str(mut_ase_n) ';']);
    %[a,m]=subset_ase_mut(ase,mut,0,10,50,0,look);
    %set the min_mut >3;
    %Set ase.cov>=30; updated on Dec-2-2018;
    %[a,m]=subset_ase_mut_by_total_available_ase_mut(ase,mut,0,10,30,3,look);
    
    %updated on Dec-24-2019;
    min_cov=30;%Keep this parameter consistent with that used in the function my_corr and my_corr_shuff;
    min_samples=3;
    coding_only=0;
    [a,m]=subset_ase_mut_by_total_available_ase_mut(ase,mut,coding_only,min_samples,min_cov);
    if (size(a.data,1)>0 && size(m.data,1)>0)
    m.data(m.data>0)=1;
    [~,pp]=my_corr(a,m);
    %Make sure to use subseted ase and mut matrics, instread of ase and
    %mut, as the final distribution of p_shuff would be different;
    [~,p_shuff]=my_corr_shuff(a,m,iter,min_samples,min_cov);
    fdr=my_fdr(pp,p_shuff);
    f_idx=fdr>1;
    fdr(f_idx)=1;
    %use matlab mafdr to calcuate fdr;
    %fdr=mafdr(10.^-pp);
    fm=fdr_mut(m);
    o=m;
    o.ase_data=a.data;
    o.fdr=fdr;
    o.fm=fm;
    o.pp=pp;
    %Change cutoff if only want to keep fdr<=0.8;
    fdr_Point8idx=find(fdr<=1);
    %o.anno_type=mut_ase_look_tmp.name(n);
    sn=length(fdr_Point8idx);
    o.anno_type=repmat(mut_ase_look_tmp.type(n),sn,1);
    Mut.fdr=[Mut.fdr;o.fdr(fdr_Point8idx)];
    Mut.p=[Mut.p;o.pp(fdr_Point8idx)];
    Mut.fm=[Mut.fm;o.fm(fdr_Point8idx)];
    
    
    %Get mut data;
    Mut.data=[Mut.data;m.data(fdr_Point8idx,:)];
    gene=vlookup_list(m.gene(fdr_Point8idx),look.tx,look.gene);
    tx=m.gene(fdr_Point8idx);
    Mut.tx=[Mut.tx;tx];
    Mut.rowlabels=[Mut.rowlabels;gene];
    Mut.anno_type=[Mut.anno_type;o.anno_type];
    Mut.collabels=m.collabels;
    %Output all info for later usage;
    out(n)=o;
    else
     disp('')
     disp('No mut and ase after QC with subset_ase_mut');
    end
  else
    disp('')
    fprintf(['Not exist for ' num2str(mut_ase_n) '\n']);
  end
end
 Mut.gene=Mut.rowlabels;
 Mut.rowlabels=strcat(Mut.rowlabels,'-',Mut.anno_type);