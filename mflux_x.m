function fm_x = mflux_x(u,v,xi,dx,nx,ischeme)

% Initialize Arrays
fb = zeros(1,nx);
fb_x=zeros(1,nx);
f_l =zeros(1,nx);
f_r =zeros(1,nx);
fx_l=zeros(1,nx);
fx_r=zeros(1,nx);

% compute flux values
f = v.*u;

% Compute Boundary values and derivatives
for e=1:nx % for every element
    f_l(e) = f_lagrange(-1,xi,f(:,e));
    f_r(e) = f_lagrange(+1,xi,f(:,e));
    fx_l(e)=fx_lagrange(-1,xi,f(:,e));
    fx_r(e)=fx_lagrange(+1,xi,f(:,e));
end 

% impose BCs
f_l =bdc_f( f_l,nx);
f_r =bdc_f( f_r,nx);
fx_l=bdc_f(fx_l,nx);
fx_r=bdc_f(fx_r,nx);

% Riemann and DRiemann
for i = 2:nx-1
    if v(1,i) >= 0
        fb(i)=f_r(i-1); fb_x(i)=fx_r(i-1);
    else
        fb(i-1)=f_l(i-1); fb_x(i-1)=fx_l(i-1);
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
            fm_x(2,i) = 2.0/dx * (3.0*(fb(i+1)-fb(i))-fb_x(i)-fb_x(i+1))/4.0;
            fm_x(3,i) = 2.0/dx * fb_x(i+1);
        case 2
            fm_x(1,i) = 2.0/dx * (2.0*(f(1,i)+f(2,i))-0.5*(7.0*fb(i)+fb(i+1)));
            fm_x(2,i) = 2.0/dx * (f(3,i)-f(1,i))/2;
            fm_x(3,i) = 2.0/dx * (-2.0*(f(2,i)+f(3,i))+0.5*(7.0*fb(i+1)+fb(i)));
        case 3
            s3=sqrt(3);
            fm_x(1,i) = 2.0/dx * (s3*(3/4*f(1,i)+5/6*f(2,i)-1/12*f(3,i))+3/4*(1.5-s3)*fb(i+1)-3/4*(1.5+s3)*fb(i));
            fm_x(2,i) = 2.0/dx * (f(3,i)-f(1,i))*s3/3;
            fm_x(3,i) = 2.0/dx * (s3*(1/12*f(1,i)-5/6*f(2,i)-3/4*f(3,i))+3/4*(1.5+s3)*fb(i+1)-3/4*(1.5-s3)*fb(i));
    end
end

% impose BCs
fm_x = bdc(fm_x,nx);