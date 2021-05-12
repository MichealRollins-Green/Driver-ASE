function a=downSizeTo(A,numbered_offsets,original_list_len)

% a=downSizeTo(A,numbered_offsets,original_list_len)

ff=fieldnames(A);
a=A;

for n=1:numel(ff)
    
    tA=A.(ff{n});
    [ss,tt]=size(tA);
    
    % downsize y
    if ss==original_list_len
        a.(ff{n})=tA(numbered_offsets,:);
    elseif tt==original_list_len
        a.(ff{n})=tA(:,numbered_offsets);
    end
end
