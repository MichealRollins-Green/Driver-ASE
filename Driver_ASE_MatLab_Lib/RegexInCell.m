function [idx,matched]=RegexInCell(input_cell,regexp)
%regexp=strcat(regexp,'\>');%make exact match!
idx=find(~cellfun(@isempty,regexpi(input_cell,regexp)));
matched=input_cell(idx);

