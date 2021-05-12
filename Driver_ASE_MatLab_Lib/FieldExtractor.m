function d=FieldExtractor(dsd,Var,Colwide,FlterType,Flter)
%Usage:
% new_d=FieldExtractor(d,'data',1,1,'data(:,1)<0.05');
% new_d=FieldExtractor(C4,'rowlabels',1,0,'rs11');


%Extract all fields of d based on the filter applied to Var;
%Colwide means applying the filter columwide, otherwise in Rowwide;
%FlterType: 0=>character filter by using perl regexpression Flter for Var
%           1=>numeric filter, i.e, <=0.05 or others, for numeric var;
%PerlRegexp will be applied on Var;

%For debugging;
disp('.........................................................................')
disp('It is going to filter the bellowing Cell Data with FieldExtractor');
dsd
disp('.........................................................................')

tag=0;
if (~ismember(Var,fieldnames(dsd)))
    disp('No field for the dsd');
    tag=1;
end

d=dsd;

if (FlterType==0)
   %Flter=strcat(Flter,'\>');%make exact match!
   flt_idx=~cellfun(@isempty,regexpi(getfield(dsd,Var),Flter));
   if (sum(flt_idx)==0)
       disp(['Flter is ' Flter]);
       disp('No match for the it');
       tag=1;
   end;
 end;   
if (FlterType==1)
    flt_idx=eval(['d.' Flter]);
 end;
 
if tag==1
    d=[];
else
if (Colwide==1)
    d=FieldExtractorByIdx(d,1,flt_idx);
else
    d=FieldExtractorByIdx(d,0,flt_idx);    
end
end

