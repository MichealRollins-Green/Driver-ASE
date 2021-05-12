function S=load_mutations(mutations)
if (exist(mutations,'dir')==0)
    error('Cannot find the directory')
end
S.rowlabels=textread(strcat(mutations,'/rowlabels.txt'),'%s','delimiter','\t');
S.collabels=textread(strcat(mutations,'/collabels.txt'),'%s','delimiter','\t');
S.data=load(strcat(mutations,'/matrix.tab'));
S.data(isnan(S.data))=0;
S.data=sparse(S.data);
