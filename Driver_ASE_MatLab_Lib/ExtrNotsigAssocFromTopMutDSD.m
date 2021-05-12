function d=ExtrNotsigAssocFromTopMutDSD(dsd,Var,PerlRegexp,varargin)
%PerlRegexp will be applied on Var;
if (~ismember(Var,fieldnames(dsd)))
    error('No field for the dsd');
end
d=dsd;
flt_idx=find(~cellfun(@isempty,regexp(getfield(dsd,Var),PerlRegexp)));
if (nargin>3)
 if (strcmp('fdr',varargin{1}))
    fdr_idx=find(d.fdr>varargin{2});
    [off_idx,~]=intersect(fdr_idx,flt_idx);

 elseif (strcmp('p',varargin{1}))
    p_idx=find(d.p <= -log10(varargin{2}));
    [off_idx,~]=intersect(p_idx,flt_idx);
 end
else
    off_idx=flt_idx;
end


%off_idx=~cellfun(@isempty,regexp(d.anno_type,'(chip-seq|UTR|1kb.1kb|miRNA|CpG|Intron|DNAse_HS)'));
%filter with sample size,make sure there are more than 3 indivs containing
%mutations for ASE analysis;
off=find(sum(d.data,2)>=3);
off_idx=intersect(off,off_idx);
d.data=d.data(off_idx,:);
d.rowlabels=d.rowlabels(off_idx);
d.gene=d.gene(off_idx);
d.fdr=d.fdr(off_idx);
d.anno_type=d.anno_type(off_idx);
d.fm=d.fm(off_idx);
d.p=d.p(off_idx);
