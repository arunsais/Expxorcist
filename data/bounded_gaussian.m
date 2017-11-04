% method to sample data from a gaussian graph.
function [] = bounded_gaussian(n, p, exNum, graph)
rng(p);
% generate covariance matrix
if(strcmp(graph, 'chain'))
    sigma = eye(p);
    rho = -0.8;
    for i = 1 : p
        for j = 1 : p
            sigma(i,j) = rho^abs(i - j);
        end
    end
    omega = inv(sigma);
    omega(abs(omega) < 1e-5) = 0;
elseif(strcmp(graph, 'grid'))
    omega = eye(p,p);
    for i = 1 : p
        if rem(i, 10) ~= 1
            omega(i, i - 1) = 0.2;
        end
        if rem(i, 10) ~= 0
            omega(i, i + 1) = 0.2;
        end
        if i > 10
            omega(i, i - 10) = 0.2;
        end
        if i <= p - 10
            omega(i, i + 10) = 0.2;
        end
    end
    fprintf('determinant: %f\n',det(omega));
    sigma = inv(omega);
end

adj = omega; adj(omega ~= 0) = 1;
adj(1:p+1:p*p) = 0;

d = diag(omega).^-0.5;
omega = diag(d)*omega*diag(d);
sigma = inv(omega);
mu = ones(p,1);

rng(n);
xTrain = cell(1, exNum);
xTest = cell(1, exNum);
for i = 1 : exNum
    % generate latent variables
    xTrain{i} = mvnrnd(mu, sigma, n);
    xTest{i} = mvnrnd(mu, sigma, n);
    disp(max(abs(xTrain{i}(:))));
    assert(max(abs(xTrain{i}(:))) < 15);
end

% save to a file
TrueModel = [ 'gaussian_' graph];
if ~exist(TrueModel, 'dir')
    mkdir(TrueModel);
end
fileName = [TrueModel '/data_' num2str(n) '_' num2str(p) '.mat'];
save(fileName, 'xTrain', 'xTest', 'omega', 'sigma', 'adj','mu');

end