function dfdxi = mmc_FR(a,u,L,nE,scheme)
% Multi-moment Flux Reconstruction 
% MCV3, MCV3_UPCC & MCV3_CPCC

% Solution array for df/dxi
dfdxi=zeros(3,nE);

% Build discrete flux
f = a.*u;

% Interpolate flux and dflux values at the boundaries of Ij
fL=L.l*f; dfL=L.dl*f; fR=L.r*f; dfR=L.dr*f;

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
    fb=fR([nE,1:nE-1]); dfb=dfR([nE,1:nE-1]);
else
    fb=fL; dfb=dfL;
end

% Compute df/dxi
switch scheme
    case 1 % MCV3
        dfdxi(1,:)=dfb;
        dfdxi(2,:)=(3*fb([2:nE,1])-3*fb-dfb-dfb([2:nE,1]))/4;
        dfdxi(3,:)=dfb([2:nE,1]);
    case 2 % MCV3_UPCC
        %dfdxi(1,:)=(2*(f(i,1)+f(i,2))-0.5*(7*fb(i)+fb(i+1)));
        %dfdxi(2,:)=(f(i,3)-f(i,1))/2;
        %dfdxi(3,:)=(-2*(f(i,2)+f(i,3))+0.5*(7*fb(i+1)+fb(i)));
    case 3 % MCV3_CPCC
        %s3=sqrt(3);
        %dfdxi(1,:)=(s3*(3/4*f(i,1)+5/6*f(i,2)-1/12*f(i,3))+3/4*(1.5-s3)*fb(i+1)-3/4*(1.5+s3)*fb(i));
        %dfdxi(2,:)=(f(i,3)-f(i,1))*s3/3;
        %dfdxi(3,:)=(s3*(1/12*f(i,1)-5/6*f(i,2)-3/4*f(i,3))+3/4*(1.5+s3)*fb(i+1)-3/4*(1.5-s3)*fb(i));
end