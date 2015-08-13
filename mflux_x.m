function fm_x = mflux_x(u,v,xi,dx,nx,ischeme)

% Initialize Arrays
f_l =zeros(1,nx);
f_r =zeros(1,nx);
fx_l=zeros(1,nx);
fx_r=zeros(1,nx);

% Compute Boundary values and derivatives
for e=1:nx % for every element
    fe = v(:,e).*u(:,e); % compute f = v.*u;
    f_l(e) = f_lagrange(-1,xi,fe);
    f_r(e) = f_lagrange(+1,xi,fe);
    fx_l(e)=fx_lagrange(-1,xi,fe);
    fx_r(e)=fx_lagrange(+1,xi,fe);
end 

% impose BCs
f_l =bdc_f( f_l,nx);
f_r =bdc_f( f_r,nx);
fx_l=bdc_f(fx_l,nx);
fx_r=bdc_f(fx_r,nx);

% Riemann and DRiemann
for i = 2:nx-1
    if v(1,i) >= 0
        fb(i)=f_r(i);   fb_x(i)=fx_r(i);
    else
        fb(i)=f_l(i+1); fb_x(i)=fx_l(i+1);
    end 
end 

% impose BCs
fb=bdc_f(fb,nx); fb_x=bdc_f(fb_x,nx);

% Derivative of the Hermite Interpolation function 
% evaluated at solutions points \xi=[-1,0,1]
for i = 1:nx-1
    switch ischeme
        case 1 % MCV3
            fm_x(1,i) = 2.0/dx * fb_x(i);
            fm_x(2,i) = 2.0/dx * (3.0*(fb(i+1)-fb(i)) - fb_x(i) - fb_x(i+1))/4.0;
            fm_x(3,i) = 2.0/dx * fb_x(i+1);
        case 2
            fm_x(:,i) = zeros(3,1);
        case 3
            fm_x(:,i) = zeros(3,1);
    end
end

% impose BCs
fm_x = bdc(fm_x,nx);