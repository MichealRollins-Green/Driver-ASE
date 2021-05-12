function print_excel_sheet(all_data,filename)

% print_excel_sheet(all_data,filename)
% data format: all_data
%         data: [1 2]
%   collabels: {'zc'  'xxc'}
%    rowlabels: {'r'}

pid=fopen(filename,'w');

fprintf(pid,['Rowlabels']);

for n=1:length(all_data.collabels)
    fprintf(pid,['\t' char(all_data.collabels(n)) ]);
end
fprintf(pid,'\n');

for n=1:length(all_data.rowlabels)
    fprintf(pid,[char(all_data.rowlabels(n)) ]);
    
    for j=1:length(all_data.collabels)
        
        fprintf(pid,['\t' num2str(all_data.data(n,j)) ]);
        
    end
    fprintf(pid,'\n');
end

fclose(pid);
