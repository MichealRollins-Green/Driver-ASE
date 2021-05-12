function [target_out,off]=vlookup_list(reference,target,target_values)

% will return off, a position in reference where target is found, for each
% e.g. usage
% ucsc_expression=vlookup(ucsc_genes,refseq_genes,refseq_expression)
% target(off)=reference

off=nan(length(reference),1);
target_out=cell(length(reference),1);

for n=1:length(reference)
    tt=find(strcmp(reference(n),target));
    if isempty(tt)
        continue
    else
        off(n)=tt(1);
    end
end

xx=isfinite(off);
target_out(xx)=target_values(off(xx));

for n=1:length(target_out)
    if isempty(target_out{n})
        target_out{n}='';
    end
end





