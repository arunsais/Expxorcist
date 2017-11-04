% proximal operator for group lasso

function v = prox_gl(x, lambda)
norms = sum(x.^2, 2).^0.5;
shrinkage = 1 - lambda./norms;
shrinkage(norms <= lambda) = 0;
v = sparse(1:numel(shrinkage), 1:numel(shrinkage), shrinkage)*x;
end
