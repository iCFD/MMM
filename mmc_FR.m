function dfdxi = mmc_FR(a,u,L,nE,scheme)
% Multi-moment Flux Reconstruction 
% MCV3, MCV3_UPCC & MCV3_CPCC

% Build discrete flux
f = a.*u;

% Interpolate flux and dflux values at the boundaries of Ij
fL=L.l*f; dfL=L.dl*f; fR=L.r*f; dfR=L.dr*f;

% Initialize indexes
ip=[2:nE,1];    % i+1
im=[nE,1:nE-1]; % i-1

% Build Numerical fluxes across faces
%f_nface = [0,fR]; % -side
%f_pface = [fL,0]; % +side

% Apply Periodic BCs
%f_nface(1) = f_nface(end); % left BD
%f_pface(end) = f_pface(1); % right BD

% LF numerical flux
%alpha=abs(a);
%nflux=0.5*(flux(u_nface)+flux(u_pface)-alpha*(u_pface-u_nface));
%nfluxL=nflux(1:end-1); nfluxR = nflux(2:end);

% Riemann data at the boundaries
if a(1)>0
    fb=fR(im); dfb=dfR(im);
else
    fb=fL; dfb=dfL;
end

% Compute solution array for df/dxi
dfdxi=zeros(3,nE);

switch scheme
    case 1 % MCV3
        dfdxi(1,:)=dfb;
        dfdxi(2,:)=(3*fb(ip)-3*fb-dfb-dfb(ip))/4;
        dfdxi(3,:)=dfb(ip);
    case 2 % MCV3_UPCC
        dfdxi(1,:)= 2.0*(f(1,:)+f(2,:))-0.5*(7.0*fb+fb(ip));
        dfdxi(2,:)=(f(3,:)-f(1,:))/2;
        dfdxi(3,:)=-2.0*(f(2,:)+f(3,:))+0.5*(7.0*fb(ip)+fb);
    case 3 % MCV3_CPCC
        s3=sqrt(3);
        dfdxi(1,:)=s3*(3/4*f(1,:)+5/6*f(2,:)-1/12*f(3,:))+3/4*(1.5-s3)*fb(ip)-3/4*(1.5+s3)*fb;
        dfdxi(2,:)=(f(3,:)-f(1,:))*s3/3;
        dfdxi(3,:)=s3*(1/12*f(1,:)-5/6*f(2,:)-3/4*f(3,:))+3/4*(1.5+s3)*fb(ip)-3/4*(1.5-s3)*fb;
end