function S=load_annotations(annotations)
npars=nargin;

if (exist(annotations,'dir')==0)
    disp('Cannot find the directory, searching current directory')
    S.data=load('matrix.tab');
    S.rowlabels=textread('rowlabels.txt','%s','delimiter','\t');
    [~,S.aa,S.syn]=textread('syn.tab','%s%s%n','delimiter','\t');
    S.collabels=textread('collabels.txt','%s','delimiter','\t');
    
else
   S.data=load(strcat(annotations,'/matrix.tab'));
   S.rowlabels=textread(strcat(annotations,'/rowlabels.txt'),'%s','delimiter','\t');
   [~,S.aa,S.syn]=textread(strcat(annotations,'/syn.tab'),'%s%s%n','delimiter','\t');
   S.collabels=textread(strcat(annotations,'/collabels.txt'),'%s','delimiter','\t');
end;

S.data(isnan(S.data))=0;
S.data=sparse(S.data);


