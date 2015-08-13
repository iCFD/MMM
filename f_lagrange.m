function output = f_primary(x, xi, f)

ns=numel(xi);
p=ones(1,ns);

for j = 1:ns
    k = j;
    for j1 = 1:ns
        if j1 ~= k
            p(j) = p(j) * (x-xi(j1)) / (xi(k)-xi(j1));
        end 
    end % j1
end % j

output=dot(p,f);