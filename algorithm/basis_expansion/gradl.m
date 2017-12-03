% gradient of negative node conditional log likelihood
% psi - basis expansion matrices
% s : node at which to compute the gradient
% opt_s : flag denoting if the gradient should be computed for variables corresponding to node s.
% params : parameters for node s.
% returns d x 1 matrix.

function v = gradl(psi, s, opt_s, params, args)
    % for numerical integration
    [~,n,d] = size(psi);
    tmp = sum(permute(repmat(params, 1,1,n), [1,3,2]).*psi, 3);
    fs = tmp(s,:)';
    fns = sum(tmp,1)';
    fns = 1 + fns - fs;
    
    if opt_s == 1
        v = -fns'*squeeze(psi(s,:,:))/n;
        
        % derivative of log normalization constant
        psi_s = args.grid_basis_matrix;
        y = bsxfun(@plus, (psi_s*params(s,:)')*fns', args.grid_base_measure');
        maxy = max(y, [], 1);
        y = exp(bsxfun(@minus, y, maxy));
        dens = trapz(args.grid, y)';
        
        nums = zeros(n, d);
        for i = 1 : d
            tmp = bsxfun(@times, y, psi_s(:,i));
            nums(:,i) = trapz(args.grid, tmp)'.*fns./dens;
        end
        
        v = v + mean(nums,1);
    else
        % derivative of log normalization constant
        fs_tmp = args.grid_basis_matrix*params(s,:)';
        y = bsxfun(@plus, fs_tmp*fns', args.grid_base_measure');
        maxy = max(y, [], 1);
        y = exp(bsxfun(@minus, y, maxy));
        
        dens = trapz(args.grid, y)';
        tmps = trapz(args.grid, bsxfun(@times, y, fs_tmp))'./dens;
        v = zeros(size(params));
        for t = 1 : size(params,1)
            if t == s
                continue;
            end
            v(t,:) = (tmps-fs)'*squeeze(psi(t,:,:))/n;
        end
        v(s,:) = [];
    end
    
    
    if any(isnan(v(:))) || any(isinf(v(:)))
        warning('nan or inf detected in gradient');
    end
end