function print_vars4bedgraph(look,A,var2gene,gene_id,A_off,outbed)

% tx
[~,aa]=intersect(look.gene,gene_id);
tx=look.tx(aa);

% all vars in gene
[~,off]=ismember(var2gene.gene,tx);
vars=var2gene.var(logical(off));

% vars in gene AND type of var
[~,v_off]=ismember(vars,A.rowlabels);
comm=intersect(v_off,find(A.data(:,A_off)));
if length(comm)<1
    fprintf(['\nWrong with this driver\n' char(gene_id)]);
    return;
end
vars_filt=A.rowlabels(comm);

% print
tn=length(vars_filt);
pos=nan(tn,1);
for n=1:tn;
    tmp=strsplit(vars_filt{n},'-');
    ps=str2double(char(tmp(2)));
    pos(n)=ps;
end;
pos=sort(pos);
CHR=regexprep(vars_filt{1},'-.*','');
fileID = fopen(outbed,'w');
fprintf(['Please get the bed file in the file: ' outbed '\n']);
%fprintf(['track name="driver mutations" description="' char(gene_id) ' (' tx{1} ') for ' A.collabels{A_off} '" itemRgb="On"\n']);
fprintf(fileID,['track name="driver mutations" description="' char(gene_id) ' (' tx{1} ')' '" itemRgb="On"\n']);

fprintf(fileID,['browser position ' CHR ':' num2str(pos(1)-1000) '-' num2str(pos(tn)+1000) '\n']);
fprintf(fileID,'browser hide all\nbrowser dense refGene\nbrowser full altGraph\n');
for n=1:tn;
    varInfo=strsplit(vars_filt{n},'-');
    chr=varInfo{1};
    st=str2double(char(varInfo(2)))-1;
    ed=varInfo{2};
    %fprintf(fileID,[chr '\t' num2str(st) '\t' ed '\t' vars_filt{n} '-' char(gene_id)  '\n']);
    %Make the output simpler;
    fprintf(fileID,[chr '\t' num2str(st) '\t' ed '\t' vars_filt{n} '\n']);
end
fclose(fileID);