function [All_Assoc,top_mut_dsd]=M1_Import_ASE_and_Mutation_data(cancer_type,regex4feature)

directory = fileparts(which(mfilename)); %get directory path of exexuting function
cd ../MatLab_Variables;
load 'directory_variables'

disp(strcat({'Loaded variables'}, matlab_vars, somatic_calls, matrix, mutations, annotations, mut_ase_auto, hits, driver_beds));

cwd=strcat('../../../Driver_ASE_MatLab/MatLab_Analysis','/',cancer_type);
cd(cwd);
%get fullpath;
cwd=pwd;

disp('Changed to directory: ');
disp(pwd);

T = load_data(matrix);

disp(strcat('Loaded array T:'));
disp(T);

[look.tx,look.gene,look.coding] = textread(strcat(somatic_calls,'/tx2gene2coding.lookup.txt'),'%s%s%n','delimiter','\t');

disp('Loaded variable look:');
disp(look);

T.tx = T.gene;
disp('Copied T.gene to new T.tx');

T.gene = vlookup_list(T.tx,look.tx,look.gene);
fprintf('\nRan function vlookup_list and assigned results to T.gene\n');

T.collabels = regexprep(T.collabels,'(TCGA-[^-]+-[^-]+).*$','$1');

%load mutation data;
A = load_annotations(annotations);
M = load_mutations(mutations);
M.syn = A.syn;
M.collabels=regexprep(M.collabels,'.*:(TCGA.{8}).*','$1');

cd(somatic_calls);
disp(strcat({'Changed to directory: '},pwd)); 
[var2gene.gene,var2gene.var] = textread('var2gene.tab','%s%s','delimiter','\t');

[look_id2tcga.id,look_id2tcga.tcga,look_id2tcga.uutc] = textread('look.tab','%s%s%s','delimiter','\t:');
M.tcga = vlookup_list(M.collabels,look_id2tcga.uutc,look_id2tcga.tcga);
fprintf('\nRan vlookup_list to get tcga IDs for mutations\n');

%Only focus on samples available of mutations;
kidx=ismember(T.collabels,M.tcga);
T=FieldExtractorByIdx(T,1,kidx);
%remove duplicate samples;
[~,unqidx]=unique(T.collabels);
%for demo, only 30 unique samples were selected for ASE-Mut association;
T=FieldExtractorBySubIdx(T,1,unqidx);

%QC for tumor ASE data;
% gt_plot_tumor_only is a function to check relationship between the number of harvested genes and the number of tumor samples harboring ASE gene among all tumor samples
% The threshold of total reads, 20, 50, 75, 100 and 200, and p_binom>2, are applied for each gene in tumor samples
% The total number of genes passed the above filters are assigned as y-axis values
% At different tumor sample size, the above gene number is calculated

ss = gt_plot_tumor_only(T,20,2);
hold off
plot(ss,'b','linewidth',2);
hold on
ss=gt_plot_tumor_only(T,30,2);
plot(ss,'g','linewidth',2);
hold on
ss = gt_plot_tumor_only(T,50,2);
plot(ss,'k','linewidth',2);
ss = gt_plot_tumor_only(T,75,2);
plot(ss,'r','linewidth',2);
ss = gt_plot_tumor_only(T,100,2);
plot(ss,'c','linewidth',2);
ss = gt_plot_tumor_only(T,200,2);
plot(ss,'m','linewidth',2);
grid on
legend('20','30','50','75','100','200')
xlabel('Number of samples in which gene is detected');
xlim([0,length(T.collabels)]);
ylabel('Num of Genes with ASE can be measured in tumors');


T.lr = log2(T.a./T.b);
fprintf('\nRan function log2 to divide T.a by T.b and assigned values to T.lr\n');
%TA is a function to adjust read counts of ref|alt to 1000 and calcuate ASE p ???
%T.p_true_ase = TA(T.a,T.b);
T.p_true_ase=T.p_binom;
%disp('Ran function TA and assigned results to T.p_true_ase');
T.cov = T.a + T.b;

[mut_ase_look.num,mut_ase_look.name,mut_ase_look.type] = textread(strcat(somatic_calls,'/mut_ase_look.txt'),'%n%s%s','delimiter','\t');
% create var2gene matrix-style lookup
% read function var2gene_matrix in the attachment
% var2gene_matrix creates v2g.data based on unique vars and genes
% this part is very slow;
%v2g = var2gene_matrix(var2gene.var,var2gene.gene);

v2g=load_Matrix;
v2g.var=v2g.rowlabels;
v2g.gene=v2g.collabels;
v2g=rmfield(v2g,{'rowlabels' 'collabels'});

disp('Ran var2gene_matrix and assigned results to variable v2g');

%Go to the output directory
cd(strcat(directory,'/',cancer_type));
disp(strcat({'Changed to directory: '},pwd));

if(exist(mut_ase_auto,'dir') == 0)
    mkdir (mut_ase_auto)
end
%make sure to look ase and mut data according to mut_ase_look;

run_mut_ase4diff_features_with_Rgx(mut_ase_auto,T,M,v2g,A,mut_ase_look,regex4feature);
fprintf('\nRan function run_mut_ase4diff_features\n');

%Get all significant genes with fdr <=0.2;
%Need to change the above cutoff if want to get all significant genes;

%Note: Get_Driver_mut will generated -log10P for association!
[All_Assoc,top_mut_dsd] = Get_Driver_mut_regex4feature(1000,look,mut_ase_look,regex4feature);
fprintf('\nRan function Get_Driver_mut_regex4feature\n');

% %use fdr and raw assoc p to filter the above dsd;
% -log10(0.05)=1.3 for the top_mut_dsd.p;
% idx = top_mut_dsd.fdr <= 0.2 & top_mut_dsd.p >= 1.3;
% top_mut_dsd.rowlabels(idx);

d = top_mut_dsd;
d.rowlabels = strcat(top_mut_dsd.gene,'-',top_mut_dsd.anno_type);
if(exist(hits,'dir') == 0)
    mkdir(hits);
end
print_excel_sheet(d,strcat(hits,'/mut_all.tab'))
d.data=d.fm;
d.collabels = {'FDR-mut'};
print_excel_sheet(d,strcat(hits,'/fm_all.tab'))
d.data = d.fdr;
d.collabels = {'FDR'};
print_excel_sheet(d,strcat(hits,'/fdr_all.tab'))
d.data=d.p;
d.collabels = {'-log10(p)'};
print_excel_sheet(d,strcat(hits,'/assoc_P_all.tab'));

%Make UCSC beds for top drivers;
d = ExtrAssocFromTopMutDSD(top_mut_dsd,'anno_type','.*','fdr',0.2);
d1 = ExtrAssocFromTopMutDSD(d,'anno_type','.*','p',0.05);
d = d1;
%update the top_mut_dsd;
top_mut_dsd=d;

if(exist(driver_beds,'dir') == 0)
    mkdir(driver_beds);
end
cd (driver_beds);
disp(strcat({'Changed to directory '},pwd));

for i = 1:length(d.anno_type)
   reg_idx = RegexInCell(mut_ase_look.type,d.anno_type(i));
   anno = mut_ase_look.name(reg_idx);
   gene = d.gene(i);
   %for multiple match!
   for ii = 1:length(anno)
     %ip=find(~cellfun(@isempty,regexpi(V.collabels,anno)));
     [ip,~] = RegexInCell(A.collabels,anno(ii));
     
     % Function print_vars(look,V,var2gene,gene_id,V_off)
     %print_vars(look,V,var2gene,gene,ip)
       
     %Generate bed files for deciphering underlying mechanism
     output = sprintf(['%s.%s.out.bed'],char(gene),char(d.anno_type(i)));
     print_vars4bedgraph(look,A,var2gene,gene,ip,output)
   end
end

tidx = ismember(strcat(top_mut_dsd.gene,':',top_mut_dsd.anno_type),strcat(d.gene,':',d.anno_type));
top = FieldExtractorBySubIdx(top_mut_dsd,0,tidx);
OutPutDriverParametersIntoTxt(top,'driver_top_fdr0.2.txt');


% save(strcat('T_',cancer_type),'T','-v7.3');
% disp(strcat({'Saved T as T_'},cancer_type,{' in directory '},pwd));
% save(strcat('v2g_',cancer_type), 'v2g');
% disp(strcat({'Saved v2g as v2g_'},cancer_type,{' in directory '},pwd));
% save look look
% disp(strcat({'Saved look in directory '},pwd));
% save(strcat('var2gene_',cancer_type), 'var2gene');
% disp(strcat({'Saved var2gene as var2gene_'},cancer_type,{' in directory '},pwd));
% save(strcat('M_',cancer_type), 'M');
% disp(strcat({'Saved M as M_'},cancer_type,{' in directory '},pwd));
% save(strcat('A_',cancer_type), 'A', '-v7.3');
% disp(strcat({'Saved A as A_'},cancer_type,{' in directory '},pwd));
% save(strcat('look_id2tcga_',cancer_type), 'look_id2tcga');
% disp(strcat({'Saved look_id2tcga as id2tcga_'},cancer_type,{' in directory '},pwd));
% save(strcat('mut_ase_look_',cancer_type), 'mut_ase_look');
% disp(strcat({'Saved mut_ase_look as mut_ase_look_'},cancer_type,{' in directory '},pwd));
% save(strcat('Mut_o_',cancer_type), 'All_Assoc','top_mut_dsd');
% disp(strcat({'Saved All_Assoc and top_mut_dsd as Mut_o_'},cancer_type,{' in directory '},pwd));

disp('Your ASE-Mut association results are saved in: ');
disp(pwd);
cd(cwd);
disp(strcat({'Changed to directory '},pwd));

end