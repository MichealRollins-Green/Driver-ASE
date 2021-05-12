function v2g=var2gene_matrix(var,gene)

% step 1 - collapse vars
[v2g.var,bb,cc]=unique(var);
ll=max(cc);
dd=sparse(ll,length(var));

fprintf([num2str(length(bb)/10000) '\n']);
hh=0;
for n=1:ll
    dd(n,n==cc)=1;
    if hh==10000
        fprintf('.');
        hh=0;
    end
    hh=hh+1;
end

fprintf('.');

% step 2 - collapse genes
[v2g.gene,bb,cc]=unique(gene);

v2g.data=sparse(length(v2g.var),length(v2g.gene));

fprintf([num2str(length(bb)/1000) '\n']);
hh=0;

for n=1:length(bb)
    
    xx=dd(:,n==cc);
    xx=sum(xx,2);
    xx=full(xx);
    v2g.data(:,n)=xx;
    if hh==1000
        fprintf('.');
        hh=0;
    end
    hh=hh+1;
end
