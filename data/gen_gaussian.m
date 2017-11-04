% bounded_gaussian(50,100,10,'chain');
% bounded_gaussian(100,100,10,'chain');
% bounded_gaussian(200,100,10,'chain');
% bounded_gaussian(500,100,10,'chain');
% 
% bounded_gaussian(50,100,10,'grid');
% bounded_gaussian(100,100,10,'grid');
% bounded_gaussian(200,100,10,'grid');
% bounded_gaussian(500,100,10,'grid');

p = [20, 30, 50];
for i = 1 : numel(p)
    n = p(i)*[1, 2, 4, 8, 12, 16, 20];
    for j = 1 : numel(n)
        bounded_gaussian(n(j), p(i), 20,'chain');
        bounded_gaussian(n(j), p(i), 20,'grid');
    end
end

% bounded_gaussian(50,10,10,'chain');
% bounded_gaussian(100,10,10,'chain');
% bounded_gaussian(200,10,10,'chain');
% bounded_gaussian(500,10,10,'chain');
% 
% bounded_gaussian(50,10,10,'grid');
% bounded_gaussian(100,10,10,'grid');
% bounded_gaussian(200,10,10,'grid');
% bounded_gaussian(500,10,10,'grid');