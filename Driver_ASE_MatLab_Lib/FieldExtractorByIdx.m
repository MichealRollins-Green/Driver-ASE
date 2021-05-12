function d=FieldExtractorByIdx(dsd,Colwide,idx_in)
%Usage:
%Note: if extracting colwide, make sure to these vector in struct as
%1*n format,except for collabels, which will be handled by the function;

% new_d=FieldExtractorByIdx(d,1,idx_in);
% new_d=FieldExtractorByIdx(C4,0,idx_in);
%Extract all fields of d based on the filter applied to all Vars with MATCHED dimansion;
%Make sure the size(idx_in,1) is matchable with potential Vars in dsd;
%Colwide means applying the filter columwide, otherwise in Rowwide;
disp('')
disp('....................................Start running.....................................')
disp('It is going to filter the bellowing Cell Data with FieldExtractorByIdx');
disp(dsd)


d=dsd;
%check the idx_in type: double (0 or 1) or index number (1,2,3,..,n).
tmp_sum=sum(idx_in==1)+sum(idx_in==0);
if (tmp_sum==length(idx_in))
    flt_idx=logical(idx_in);%make sure to change vector from double to logical;
elseif (any(idx_in==0))
    idx_in(1:10)
    error('You provided index vector is not correct')
else
    flt_idx=idx_in;
end
        
fldnames=fieldnames(d);
tot=size(fldnames,1);
fsize=size(flt_idx,1);
if (Colwide==1)
    for fn=1:tot
      fl=char(fldnames(fn));
       if (strcmp(fl,'collabels'))
              %Commonly, d.collables will be n*1 vector
              %Maybe only d.data is in the n*x matrix format;
              [cs1,cs2]=size(eval(['d.' fl]));
              if (cs1>cs2)
                  warning('You dsd collabels is in the n*1 vector format, we will filter it anyhow!');
                  cs=cs1;
              else
                  cs=cs2;
              end
       else
              cs=size(eval(['d.' fl]),2); 
       end
      
      if (cs > 1 && cs==fsize && isfield(d,fl) && ~strcmp(fl,'rowlabels'))
          disp(['Going to filter ' upper(fl) ' with idx_in in column-wide'])
          if (strcmp(fl,'collabels'))
              %Commonly, d.collables will be n*1 vector
              %Maybe only d.data is in the n*x matrix format;
              %eval(['d.' fl '=d.' fl '(flt_idx,:);']);
                  d.collabels=d.collabels(flt_idx);
          else
              eval(['d.' fl '=d.' fl '(:,flt_idx);']);
          end
      else
          if (cs==1 && size(eval(['d.' fl]),2)==1)
             disp(['Not going to filter ' upper(fl) ' with idx_in in column-wide, as the dimansion of col and row are both single'])
          else
          disp(['Not going to filter ' upper(fl) ' with idx_in in column-wide, as the dimansion is not the same between flt_idx and target cell'])
          end
      end
      
    end
else
    for fn=1:tot
      fl=char(fldnames(fn));
      cs=size(eval(['d.' fl]),1);%make sure to use size(x,1) here;
      %Not extract data from collabels;
      if (cs > 1 && cs==fsize && isfield(d,fl) && ~strcmp(fl,'collabels'))
          disp(['Going to filter ' upper(fl) ' with idx_in in row-wide'])
           eval(['d.' fl '=d.' fl '(flt_idx,:);'])
      else
          if (cs==1 && size(eval(['d.' fl]),1)==1)
             disp(['Not going to filter ' upper(fl) ' with idx_in in row-wide, as the dimansion of col and row are both single'])
          else
             disp(['Not going to filter ' upper(fl) ' with idx_in in row-wide, as the dimansion is not the same between flt_idx and target cell'])
          end
      end
      
    end
    
end;
disp('....................................Complete.....................................')
disp('')

