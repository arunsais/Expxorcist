% Alternating optimization of the regularized node conditional log likelihood

function params = alt_min(psi, s, params, lambda, args)
alt_min_iters = 25;

for i = 1 : alt_min_iters
    %% optimization of node specific params
    % create anonymous functions
    f = @(Z) l_h(psi, s, 1, params, Z, args);
    gradf = @(Z) gradl_h(psi, s, 1, params, Z, args);
    prox = @(Z, eta) prox_gl(Z, lambda * eta);
    
    % acc_prox_gd
    params(s,:) = accproxgd(f, gradf, prox, params(s,:), 40, 1);
    %fprintf('%d, %f\n', s, f(params(s,:)) + lambda * norm(params(s,:)));
    
    %% optimization of interaction params
    % create anonymous functions
    f = @(Z) l_h(psi, s, 0, params, Z, args);
    gradf = @(Z) gradl_h(psi, s, 0, params, Z, args);
    prox = @(Z, eta) prox_gl(Z, lambda * eta);
    
    % acc_prox_gd
    params([1:s-1,s+1:end],:) = accproxgd(f, gradf, prox, params([1:s-1,s+1:end],:), 40, 1);
end
end

function v = l_h(psi, s, opt_s, params, Z, args)
if opt_s == 1
    params(s,:) = Z;
else
    params([1:s-1,s+1:end],:) = Z;
end
v = l(psi, s, params, args);
end

function v = gradl_h(psi, s, opt_s, params, Z, args)
if opt_s == 1
    params(s,:) = Z;
else
    params([1:s-1,s+1:end],:) = Z;
end
v = gradl(psi, s, opt_s, params, args);
end