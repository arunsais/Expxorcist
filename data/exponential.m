% instead sample from a multivariate exponential distribution
function exponential(n, p, exNum, graph)
rng(p);
% generate interaction terms matrix
if(strcmp(graph, 'chain'))
    omega = zeros(p,p);
    rho = 0.1;
    for i = 1 : p-1
        omega(i,i+1) = rho;
        omega(i+1,i) = rho;
    end
elseif(strcmp(graph, 'grid'))
    omega = zeros(p,p);
    for i = 1 : p
        if rem(i, 10) ~= 1
            omega(i, i - 1) = 0.2;
        end
        if rem(i, 10) ~= 0
            omega(i, i + 1) = 0.2;
        end
        if i > 10
            omega(i, i - 10) = -0.2;
        end
        if i <= p - 10
            omega(i, i + 10) = -0.2;
        end
    end
end

adj = omega; adj(omega ~= 0) = 1;
mu = 2*ones(p,1);

rng(n);
xTrain = cell(1, exNum);
xTest = cell(1, exNum);
for i = 1 : exNum
    % perform Gibbs sampling. 
    % use rejection sampling to sample from each node conditional
    % distribution
    xTrain{i} = gibbsSampler(mu, omega, n);
    xTest{i} = gibbsSampler(mu, omega, n);
end

% save to a file
TrueModel = [ 'exponential_' graph];
if ~exist(TrueModel, 'dir')
    mkdir(TrueModel);
end
fileName = [TrueModel '/data_' num2str(n) '_' num2str(p) '.mat'];
save(fileName, 'xTrain', 'xTest', 'omega', 'adj','mu');

end

function X = gibbsSampler(mu, omega, n)
p = numel(mu);
X = ones(n, p);
x = ones(1, p);
burnin =  500; thinning = 1;
for i = 1 : n + burnin
    for k = 1 : thinning
        for j = 1 : p
            % sample from node conditional distribution
            fns = mu(j)-dot(omega(j,:), abs(x));
            x(j) = exp_sample(fns);
        end
    end
    
    if i > burnin
        thinning = 10;
        X(i - burnin, :) = x;
    end
end
end

function v = exp_sample(fns)
    xmin = -15; xmax = 15;
    if fns >= 0
        sig = rand;
        if(sig < 0.5)
            sig = -1;
        else
            sig = 1;
        end
        
        while true
            mag = exprnd(fns^-1, 1,1);
            if mag < xmax
                break;
            end
        end
        v = sig*mag;
    else
        assert(false);
    end
end