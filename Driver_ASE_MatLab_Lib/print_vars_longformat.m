function print_vars_longformat(look,V,var2gene,gene_id,V_off,output_fid)

% tx
[~,aa]=intersect(look.gene,gene_id);
tx=look.tx(aa);

% all vars in gene
[~,off]=ismember(var2gene.gene,tx);
vars=var2gene.var(logical(off));

% vars in gene AND type of var
[~,v_off]=ismember(vars,V.rowlabels);
vars_filt=V.rowlabels(intersect(v_off,find(V.data(:,V_off))));

% print
for n=1:length(vars_filt);
    fprintf(output_fid,[char(gene_id) '\t' tx{1} '\t' V.collabels{V_off} '\t']);
    fprintf(output_fid,[vars_filt{n} '\n']);
end

