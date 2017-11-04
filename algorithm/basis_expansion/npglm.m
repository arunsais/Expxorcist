%% INPUTS
% X: data matrix - n by p matrix; n - number of samples, p - number of nodes in the graph.
% Xtest: test data matrix - ntest by p matrix
% lambda: regularization parameter for group sparse regularized node conditional MLEs.
% d: truncation parameter for basis expansion.
% base_measure_type: specifies which base measure to use. look at base_measure.m file for a list of options
% base_measure_args: a float for specifying the base measure.

%% OUTPUTS
% params: returns the parameters learned by the algorithm
% adj_or: adjacency matrix of the graph learned. uses ‘or’ rule to stitch the estimates from different nodes
% adj_and: adjacency matrix of the graph learned. uses ‘and’ rule to stitch the estimates from different nodes
% train_nllk: p x 1 matrix - node conditional likelihoods of the training data evaluated at all the nodes.
% test_nllk : p x 1 matrix - node conditional likelihoods of the test data evaluated at all the nodes.

% IMPORTANT INFO
% *** current implementation assumes p > 1 ***
% *** the code infers the range of the data from X***
% *** the code uses the same range for all the p variables***
% *** currently the code only works for cosine basis functions ***



function [params, adj_or, adj_and, train_nllk, test_nllk] = npglm(X, Xtest, lambda, d, base_measure_type, base_measure_args)
    %% Problem setup
    [n, p] = size(X);
    
    % Infer the range of the data
    xmin = min(X(:)); xmax = max(X(:)); 


    % grid for numerical integration	
    % pre-evaluate the basis functions at the grid points, for faster numerical integration
    grid = xmin:(xmax-xmin)/100:xmax;
    grid_basis_matrix = basis_matrix(grid', d, xmin, xmax);
    args.grid = grid; args.grid_basis_matrix = grid_basis_matrix;
    
    % params - stores the parameters estimated at each node
    % params is a p x p x d matrix
    % params(i,:,:) - stores parameter estimates obtained by solving node
    % conditional maximum likelihood estimation at node i.
    params = zeros(p,p,d);
    
    % pre-compute the basis matrices at the training data points
    % psi is a p x n x d matrix. with psi(i,:,:) corresponds to node i.
    psi = basis_matrix(X, d, xmin, xmax);
    
    %% compute node conditional ML estimates at each  node
    fprintf('started estimation. n = %d, p = %d, lambda = %f\n', n, p, lambda);
    for s = 1 : p
        %% Perform Alternating Minimization
        args.grid_base_measure = base_measure(args.grid, s, base_measure_type, base_measure_args);
        params(s,:,:) = alt_min(psi, s, squeeze(params(s,:,:)), lambda, args);
        
        % TODO : use these params to initialize optimization for node (s+1).
        fprintf('finished estimation of node %d\n', s);
    end
    
    % compute stats
    train_nllk = zeros(p,1);
    test_nllk = zeros(p,1);
    
    for s = 1 : p
        args.grid_base_measure = base_measure(args.grid, s, base_measure_type, base_measure_args);
        train_nllk(s) = l(psi, s, squeeze(params(s,:,:)), args);
    end
    
    psi = basis_matrix(Xtest, d, xmin, xmax);
    for s = 1 : p
        args.grid_base_measure = base_measure(args.grid, s, base_measure_type, base_measure_args);
        test_nllk(s) = l(psi, s, squeeze(params(s,:,:)), args);
    end
    
    [adj_or, adj_and] = getEdges(params);
end

function [adj_or, adj_and] = getEdges(params)
    [p, ~, ~] = size(params);
    adj_or = zeros(p,p);
    adj_and = zeros(p,p);
    for i = 1 : p
        for j = i+1 : p
            e1 = norm(squeeze(params(i,j,:)));
            e2 = norm(squeeze(params(j,i,:)));
            if isnan(e1) || isnan(e2)
                adj_or(i,j) = 1;
                adj_or(j,i) = 1;
                continue;
            end
            
            if max(e1, e2) < 1e-7
                continue;
            end
            
            adj_or(i,j) = 1;
            adj_or(j,i) = 1;
        end
    end
    
    for i = 1 : p
        for j = i+1 : p
            e1 = norm(squeeze(params(i,j,:)));
            e2 = norm(squeeze(params(j,i,:)));
            if isnan(e1) || isnan(e2)
                adj_and(i,j) = 1;
                adj_and(j,i) = 1;
                continue;
            end
            
            if min(e1, e2) < 1e-7
                continue;
            end
            
            adj_and(i,j) = 1;
            adj_and(j,i) = 1;
        end
    end
end
