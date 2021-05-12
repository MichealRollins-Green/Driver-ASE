function [ase_gene,mut_gene,anova_p,ase_mean]=extract_ase_mut4gene(genesymbol,look,ase_mut_assoc)

g_idx=find(strcmp(ase_mut_assoc.gene,genesymbol));
%Get subset genes
ase_gene=abs(ase_mut_assoc.ase_data(g_idx,:));
ase_gene(ase_gene>10)=10;
mut_gene=ase_mut_assoc.mut_data(g_idx,:);

%get mean of ase by mut groups;
ase_mean= group_mean(ase_gene,mut_gene);
[X,gidx]=sort(ase_mean(:,1));

%make boxplot
figure('Color',[1 1 1]);
boxplot(ase_gene,mut_gene,'Symbol','o','OutlierSize',2,'colors',[0 0 0]);
box off
xlabel('mutation occurrence (0=No 1=Yes)');
ylabel([genesymbol,' ASE value -log2(A/B)']);

hold on
plot(ase_mean(gidx,2),'LineStyle','none','Marker','d','MarkerSize',3,'MarkerFaceColor','k','Color',[0 0 0]);
%get anova p value
anova_p=anova1(ase_gene,mut_gene,'off');
