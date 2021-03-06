function output = fx_primary(x, xi, f)

output=(2.*x-xi(2)-xi(3))/(xi(1)-xi(2))/(xi(1)-xi(3))*f(1)...
      +(2.*x-xi(1)-xi(3))/(xi(2)-xi(1))/(xi(2)-xi(3))*f(2)...
      +(2.*x-xi(1)-xi(2))/(xi(3)-xi(1))/(xi(3)-xi(2))*f(3);