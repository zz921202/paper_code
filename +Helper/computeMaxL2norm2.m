function max_norm_squared = computeMaxL2norm2(cell_array)
    max_norm = 0;
    k = length(cell_array);
    for ind = 1:k
        max_norm = max(norm(cell_array(ind)), max_norm);
    end
    max_norm_squared = max_norm * max_norm;
end