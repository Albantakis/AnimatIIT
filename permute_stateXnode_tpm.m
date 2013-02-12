function permuted_tpm = permute_stateXnode_tpm(tpm)

permuted_tpm = zeros(size(tpm));

for i = 1:size(tpm,1)
    
    permuted_row_index = trans10(flipud(trans2(i-1,size(tpm,2))));
    permuted_tpm(permuted_row_index,:) = tpm(i,:);
    
end
