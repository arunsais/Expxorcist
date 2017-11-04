% back tracking line search
% computes a step size that guarantees decrease in function value
% details of the algorithm can be found at 
%   http://people.eecs.berkeley.edu/~elghaoui/Teaching/EE227A/lecture18.pdf (slide 29)  

% f - function
% gradf - gradient of the function
% prox - proximal operator
% x - current point
% t - step size
% f_val - function evaluated at x
function [t, x_new, f_new_val] = btls(f, gradf, prox, x, t, f_val)
    beta = 0.7;
    if t == []
        t = 1;
    end
    
    gradfx = gradf(x);
    %fprintf('%f,', norm(gradfx));
    gtx = (x - prox(x - t * gradfx, t))/t;
    f1 = f(x - t * gtx);
    f2 = f_val - t * dot(gradfx(:), gtx(:)) + t/2 * norm(gtx(:), 2)^2;
    while f1 > f2
        t = beta * t;
        gtx = (x - prox(x - t * gradfx, t))/t;
        f1 = f(x - t * gtx);
        f2 = f_val - t * dot(gradfx(:), gtx(:)) + t/2 * norm(gtx(:), 2)^2;
    end
    x_new = x - t * gtx;
    f_new_val = f1;
end