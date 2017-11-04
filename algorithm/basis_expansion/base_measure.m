function v = base_measure(grid, s, base_measure_type, base_measure_args)
    if(strcmp(base_measure_type,'gaussian'))
        v = -base_measure_args*grid.^2/2;
    elseif(strcmp(base_measure_type,'laplace'))
        v = -base_measure_args*abs(grid);
    elseif(strcmp(base_measure_type,'uniform'))
        v = zeros(size(grid));
    end
end
