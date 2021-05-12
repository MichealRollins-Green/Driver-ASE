function run_mut_ase4diff_features_with_Rgx(mut_ase_auto,T,M,v2g,A,mut_ase_look,feature_rgx)
tidx=RegexInCell(mut_ase_look.type,feature_rgx);
for i=1:length(tidx)
    n=tidx(i);
    feature=mut_ase_look.type(n);
    IsMatchedFeature=~isempty(regexp(cellstr(feature),feature_rgx,'ONCE'));
    if (IsMatchedFeature)
    bed=mut_ase_look.name(n);
    beds=A.collabels;
    fidx= strcmpi(bed,beds);
    disp('')
    fprintf(['working on ' num2str(n) ' for regulatory feature ' char(feature) '\n']);
    [mut,ase]=genes_rs2(T,M,v2g,A.rowlabels(logical(A.data(:,fidx))));
    if (size(ase,1)>0 && size(mut,1))
      eval(['save ' strcat(mut_ase_auto,'/mut_ase_',num2str(n)) ' mut ase;']);
    end
    end
    disp('')
end
