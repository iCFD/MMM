function dfdxi = mmc_FRv2(a,u,L,nE,scheme)
% A modified multi-moment reconstruction using Rational Polynomials

% Build discrete flux
f = a.*u;

% Interpolate flux and dflux values at the boundaries of Ij
fL=L.l*f; dfL=L.dl*f; fR=L.r*f; dfR=L.dr*f;

% Initialize indexes
ip=[2:nE,1];    % i+1
im=[nE,1:nE-1]; % i-1

% Riemann data at the boundaries
if a(1)>0
    fb=fR(im); dfb=dfR(im);
else
    fb=fL; dfb=dfL;
end

% Compute coefs
s=fb(ip)-fb;    %beta=(s-dfb)/(dfb(ip)-s);
%a=fb;   
b=dfb;  c=3*s-2*dfb-dfb(ip);    d=dfb(ip)+dfb-2*s;

% Compute solution array for df/dxi
dfdxi=zeros(3,nE);
switch scheme
    case 1 % MCV3 - original
        dfdxi(1,:)=b;
        dfdxi(2,:)=b+c+3*d/4;
        dfdxi(3,:)=b+2*c+3*d;
    case 2 % MCV3 - rational modified
        dfdxi(1,:)=b;
        dfdxi(2,:)=b+c+3*d/4;
        dfdxi(3,:)=b+2*c+3*d;
end