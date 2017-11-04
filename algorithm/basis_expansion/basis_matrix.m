function v = basis_matrix(X, d, xmins, xmaxs)
    [n, p] = size(X);
    v = zeros(p,n,d);
    for i = 1 : p
        if numel(xmins) == 1
            xmin = xmins;
            xmax = xmaxs;
        else
            xmin = xmins(i);
            xmax = xmaxs(i);
        end
        v(i,:,:) = basis_matrix_h(X(:,i), d, xmin, xmax);
    end
    
    if p == 1
        v = squeeze(v(1,:,:));
    end
end

% do basis expansion. 
% what basis to use depends on the problem.
% using cosine basis : http://www.stat.cmu.edu/~larry/=sml/functionspaces.pdf
function v = basis_matrix_h(X,d, xmin, xmax)
    idx = 1:d;
    v = sqrt(2/(xmax-xmin))*cos(pi*(X-xmin)/(xmax-xmin)*idx);
end
