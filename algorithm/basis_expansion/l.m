% negative node conditional log likelihood function
% psi - basis expansion matrices
% s : node at which to compute the conditional log likelihood
% params : parameters for node s.
function v = l(psi, s, params, args)
    [~,n,~] = size(psi);
    tmp = sum(permute(repmat(params, 1,1,n), [1,3,2]).*psi, 3);
    fs = tmp(s,:)';
    fns = sum(tmp,1)';
    fns = 1 + fns - fs;
    v = -mean(fs.*fns);
    
    % log normalization constant
    % TODO: need a better integration technique.
    fs_tmp = args.grid_basis_matrix*params(s,:)';
    y = exp(bsxfun(@plus, fs_tmp*fns', args.grid_base_measure'));
    v = v + mean(log(trapz(args.grid, y)));
    
    if isnan(v) ||isinf(v)
        %warning('nan or inf detected in loss function');
    end
end