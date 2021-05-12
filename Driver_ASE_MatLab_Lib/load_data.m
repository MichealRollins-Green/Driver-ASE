function T=load_data(matrix)
if (exist(matrix,'dir')==0)
    error('Cannot find the directory')
end
T.a=load(strcat(matrix,'/tumor.a'));
T.b=load(strcat(matrix,'/tumor.b'));
% T.agree=load(strcat(matrix,'/tumor.agree'));
% T.num_snps=load(strcat(matrix,'/tumor.num_snps'));
T.p_binom=load(strcat(matrix,'/tumor.p'));
% %T.problem_cnv=load('matrix/problem.cnv');
% %T.problem_snp=load('matrix/problem.snp');
% 
% %T.problem_cnv(isnan(T.problem_cnv))=0;
% %T.problem_cnv=sparse(T.problem_cnv);
% %T.problem_snp(isnan(T.problem_snp))=0;
% %T.problem_snp=sparse(T.problem_snp);

% N.a=load(strcat(matrix,'/normal.a'));
% N.b=load(strcat(matrix,'/normal.b'));
% N.agree=load(strcat(matrix,'/normal.agree'));
% N.num_snps=load(strcat(matrix,'/normal.num_snps'));
% N.p_binom=load(strcat(matrix,'/normal.p'));

T.gene=textread('matrix/rowlabels.txt','%s','delimiter','\t');
% N.gene=T.gene;

% % T.tx=textread('matrix/tx','%s','delimiter','\t');
% % N.tx=T.tx;
% % T.gene=textread('matrix/genes','%s','delimiter','\t');
% % N.gene=T.gene;

T.collabels=textread(strcat(matrix,'/tumor.collabels'),'%s','delimiter','\t');
% N.collabels=textread(strcat(matrix,'/normal.collabels'),'%s','delimiter','\t');
