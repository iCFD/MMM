function output = bdc_f(f,nx)

output = f;
output(1) = output (nx-1);
output(nx) = output (2);

%output = [f(nx),f,f(1)];