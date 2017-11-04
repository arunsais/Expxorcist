% sample data from non gaussian distributions
function bounded_nongaussian(n, p, exNum, graph)
rng(p);
% generate interaction terms matrix
if(strcmp(graph, 'chain'))
    omega = zeros(p,p);
    rho = 1;
    for i = 1 : p-1
        omega(i,i+1) = rho;
        omega(i+1,i) = rho;
    end
elseif(strcmp(graph, 'grid'))
    omega = zeros(p, p);
    rho = 1;
    for i = 1 : p
        if rem(i, 10) ~= 1
            omega(i, i - 1) = rho;
        end
        if rem(i, 10) ~= 0
            omega(i, i + 1) = rho;
        end
        if i > 10
            omega(i, i - 10) = rho;
        end
        if i <= p - 10
            omega(i, i + 10) = rho;
        end
    end
end

adj = omega; adj(omega ~= 0) = 1;
mu = 1*ones(p,1);

c = 0;
fs_type = ones(p,1); 

rng(n);
xTrain = cell(1, exNum);
xTest = cell(1, exNum);
for i = 1 : exNum
    % perform Gibbs sampling. 
    % use rejection sampling to sample from each node conditional
    % distribution
    xTrain{i} = gibbsSampler(mu, omega, n, c, fs_type);
    xTest{i} = gibbsSampler(mu, omega, n, c, fs_type);
    fprintf('done %d\n', i);
end

% save to a file
TrueModel = [ 'exp_' graph];
if ~exist(TrueModel, 'dir')
    mkdir(TrueModel);
end
fileName = [TrueModel '/data_' num2str(n) '_' num2str(p) '.mat'];
save(fileName, 'xTrain', 'xTest', 'omega', 'adj','mu', 'c', 'fs_type');

end

function X = gibbsSampler(mu, omega, n, c, fs_type)
p = numel(mu);
X = ones(n, p);
x = ones(1, p);
burnin =  500; thinning = 1;
for i = 1 : n + burnin
    for k = 1 : thinning
        for j = 1 : p
            % sample from node conditional distribution
            fns = mu(j)+dot(omega(j,:), eval_fs(x, fs_type));
            x(j) = sample(c, fns, fs_type(j));
        end
    end
    
    if i > burnin
        thinning = 10;
        X(i - burnin, :) = x;
    end
end
end

function v = sample(c, fns, fs_type)
    xmin = -1; xmax = 1;
    if fs_type == 1
        proposal_max = exp(abs(fns));
        while true
            v = (xmax-xmin)*rand + xmin;
            v2 = rand*proposal_max;
            pv = exp(sin(4*pi*v)*fns + c * abs(v));
            if pv >= v2
                break;
            end
        end
    elseif fs_type == 2
        proposal_max = exp(1.1*abs(fns));
        
        while true
            v = (xmax-xmin)*rand + xmin;
            v2 = rand*proposal_max;
            pv = exp((2*(exp(-20*(v-0.5).^2) + exp(-20*(v+0.5).^2)) - 1)*fns + c * abs(v));
            if pv >= v2
                break;
            end
        end
    elseif fs_type == 3
        if fns > 0
            proposal_max = exp(fns*log(2));
        else
            proposal_max = exp(abs(fns)*log(10));
        end
        while true
            v = (xmax-xmin)*rand + xmin;
            v2 = rand*proposal_max;
            pv = exp(log(1.1+v)*fns + c * abs(v));
            if pv >= v2
                break;
            end
        end
    end
end

function v = eval_fs(x, fs_type)
fs1 = sin(4*pi*x);
fs2 = 2*(exp(-20*(x-0.5).^2) + exp(-20*(x+0.5).^2)) - 1;

v = fs1;
v(fs_type == 2) = fs2(fs_type == 2); 

end