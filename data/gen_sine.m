
p = [50];
for i = 1 : numel(p)
    n = p(i)*[1,2,4,10];
    parfor j = 1 : numel(n)
        bounded_nongaussian(n(j), p(i), 10,'chain');
        bounded_nongaussian(n(j), p(i), 10,'grid');
    end
end
