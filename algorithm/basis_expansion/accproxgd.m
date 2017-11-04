% FISTA - accelerated proximal gradient descent
% x - starting point
% T - number of iterations
% eta - step size
function x = accproxgd(f, gradf, prox, x, T, eta)   
    gammaO = 1;
    x0 = x;
    [eta, ~] = btls(f, gradf, prox, x, eta, f(x));
    for i = 1 : T
        t = eta/sqrt(i);
        xt = prox(x - t * gradf(x), t);
        gamma = 0.5 + 0.5 * sqrt(1 + 4 * gammaO^2);
        alpha = (1 - gammaO)/gamma;
        x = (1 - alpha) * xt + alpha * x0;

        % save variables for next iteration
        gammaO = gamma;
        x0 = xt;
    end
end
