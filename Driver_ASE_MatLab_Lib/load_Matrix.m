function S=load_Matrix

if (exist('annotations','dir')==0)
    disp('Cannot find the dir in current path: annotations, and will search matrix data in current dir')
     if (exist('matrix.tab','file')==2)
         S.data=load('matrix.tab');
         S.rowlabels=textread('rowlabels.txt','%s','delimiter','\t');
         S.collabels=textread('collabels.txt','%s','delimiter','\t');
         S.data(isnan(S.data))=0;
         S.data=sparse(S.data);

     elseif (exist('SparseMatrix.tab','file')==2)
         S.data=load('SparseMatrix.tab');
         S.data=spconvert(S.data);
         S.rowlabels=textread('rowlabels.txt','%s','delimiter','\t');
         S.collabels=textread('collabels.txt','%s','delimiter','\t');  
     else
         warning('No matrix.tab or SparseMatrix.tab');
     end
else
    if (exist('annotations/matrix.tab','file')==2)
        S.data=load('annotations/matrix.tab');
        S.rowlabels=textread('annotations/rowlabels.txt','%s','delimiter','\t');
        S.collabels=textread('annotations/collabels.txt','%s','delimiter','\t');
        S.data(isnan(S.data))=0;
        S.data=sparse(S.data);

    elseif (exist('annotations/SparseMatrix.tab','file')==2)
        S.data=load('annotations/SparseMatrix.tab');
        S.data=spconvert(S.data);
        S.rowlabels=textread('annotations/rowlabels.txt','%s','delimiter','\t');
        S.collabels=textread('annotations/collabels.txt','%s','delimiter','\t');  
    else
         warning('No annotations/matrix.tab or annotations/SparseMatrix.tab');
     end
end



