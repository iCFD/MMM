function fm_x = mm_FR(a,u,L,dx,nE,K,scheme)
% Multi-moment Flux Reconstruction 
% MCV3, MCV3_UPCC & MCV3_CPCC
 
% Build discrete flux
f = a.*u;

% Interpolate u and flux values at the boundaries of Ij
switch scheme
    case {1,2}
        u_lbd = u(1,:);
        u_rbd = u(end,:);
        f_lbd = f(1,:);
        f_rbd = f(end,:);
    case {3}
        u_lbd = L.lcoef*u;
        u_rbd = L.rcoef*u;
        f_lbd = L.lcoef*f;
        f_rbd = L.rcoef*f;
end
% Build Numerical fluxes acroos faces
u_pface = [u_lbd,0]; % + side
u_nface = [0,u_rbd]; % - side

% Apply Periodic BCs
u_nface(1) = u_nface(end); % left BD
u_pface(end) = u_pface(1); % right BD

% LF numerical flux
alpha = max(max(abs(dflux(u))));
nflux = 0.5*(flux(u_nface)+flux(u_pface)-alpha*(u_pface-u_nface));
nfluxL = nflux(1:end-1); nfluxR = nflux(2:end);

%call Riemann(u,v,fb,fb_x,xi,dx,nE,K)

for i=1:nE
    switch scheme
        case 1 % MCV3
            fm_x(1,i)=fb_x(i);
            fm_x(2,i)=(3*(fb(i+1)-fb(i))-fb_x(i)-fb_x(i+1))/4;
            fm_x(3,i)=fb_x(i+1);
        case 2 % MCV3_UPCC
            fm_x(1,i)=(2*(f(i,1)+f(i,2))-0.5*(7*fb(i)+fb(i+1)));
            fm_x(2,i)=(f(i,3)-f(i,1))/2;
            fm_x(3,i)=(-2*(f(i,2)+f(i,3))+0.5*(7*fb(i+1)+fb(i)));
        case 3 % MCV3_CPCC
            s3=sqrt(3);
            fm_x(1,i)=(s3*(3/4*f(i,1)+5/6*f(i,2)-1/12*f(i,3))+3/4*(1.5-s3)*fb(i+1)-3/4*(1.5+s3)*fb(i));
            fm_x(2,i)=(f(i,3)-f(i,1))*s3/3;
            fm_x(3,i)=(s3*(1/12*f(i,1)-5/6*f(i,2)-3/4*f(i,3))+3/4*(1.5+s3)*fb(i+1)-3/4*(1.5-s3)*fb(i));
    end
end

fm_x = 2/dx*fm_x;