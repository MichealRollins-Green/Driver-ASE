function run_mut_ase4diff_features(T,M,v2g,A,mut_ase_look,iter)

for n = 1:iter
    feature = mut_ase_look.type(n);
    fidx = find(strcmp(mut_ase_look.name(n),A.collabels));
    disp('')
    fprintf(['working on ' num2str(n) ' for regulatory feature ' char(feature) '\n']);
    [mut,ase] = genes_rs2(T,M,v2g,A.rowlabels(logical(A.data(:,fidx))));
    eval(['save mut_ase_auto/mut_ase_' num2str(n) ' mut ase;']);
    disp('')
end
