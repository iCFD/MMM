%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Solving 1-D wave equation with 
%        3-point multi-moment constrained collocation schemes
%            namely MCV3, MCV3_UPCC and MCV3_CPCC schemes
%
%               du/dt + df/dx = S(x),  for x \in [a,b]
%               where f = v(x,t)*f: linear/quasilinear
%
%              coded by Manuel Diaz, NTU, 2015.07.19
%                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Refs: 
% [1] A flux reconstruction approach to high-order schemes including
%     Discontinuous Galerkin methods. H.T. Huynh, AIAA 2007.
% [2] Feng Xiao, Satoshi Ii, Chungang Chen, Xingliang Li, A note on the
%     general multi-moment constrained flux reconstruction formulation for
%     high order schemes, Applied Mathematical Modelling, Volume 37, 
%     Issue 7, 2013, ISSN 0307-904X, DOI 10.1016/j.apm.2012.10.050.    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes:  Implementation following code by Dr. Feng Xiao.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; %close all; clc;

%% Parameters
fluxfun = 'linear'; 
    cfl = 0.40;	% CFL condition.
   tEnd = 2.00;	% final time.
      K = 3;	% degree of accuaracy (default value).
     nE = 200;	% number of elements.
 scheme = 3;	% (1)MCV3, (2)MCV3_UPCC and (3)MCV3_CPCC.

% Build Solutions Points
switch scheme
    case {1,2}  % 'MCV3' or 'MCV3-UPCC'
        xi = 2.*((1:K)'-1)/(K-1)-1;      % Uniformly Distributed
    case {3}	% 'MCV3-CPCC'
        xi = -cos(((1:K)'-0.5)/K*pi);    % Chebyshev Nodes
    otherwise 
        error('not a listed scheme');
end

% Build Mesh
a=-1; b=1; dx=(b-a)/nE; xc=(a-dx/2):dx:(b+dx/2); 
nE = nE+2; x=ones(3,1)*xc+(dx/2)*xi*ones(1,nE);

% Define velocity fields functions
switch fluxfun
    case 'linear'
        advect = @(x) 1*ones(size(x));
    case 'quasilinear' 
        advect = @(x) 1.0+0.5*sin(2*pi*x);
end

% Build discrete velocity field
v=advect(x);

% Build IC
ICcase = 1;	% (1)Testing, (2)Gaussian.
switch ICcase
    case 1	% Testing IC
        u0 = TestingIC(x); % Jiang and Shu IC
    case 2	% Guassian IC
        u0 = IC(x,1); % Gaussian wave: u_0(x) = exp(-20*x.^2)
end

% set BCs 
u0(:,1)=0; u0(:,nE)=0;

% Exact solution
ue=u0; 

% Set plot range
d=0.2; plotrange = [a,b,min(u0(:))-d,max(u0(:))+2*d];

%% Solver Loop

% set initial time step
dt0=cfl*dx/max(v(:));

% Set initial time & load IC
t=0; u=u0; it=0; dt=dt0;

% index
i=2:nE-1;

while t < tEnd
%for nn = 1:20
    
    % Correction for final time step
    if t+dt>tEnd, dt=tEnd-t; end
    
    % Initialization
    uo = u;
    
    % step 1    
    fm_x=mflux_x(u,v,xi,dx,nE,scheme); f1=fm_x;
    for k=1:K; u(k,i)=uo(k,i)-dt*f1(k,i); end 
    u = bdc(u,nE);
    
    % step 2
    fm_x=mflux_x(u,v,xi,dx,nE,scheme); f2=fm_x;
    for k=1:K; u(k,i)=uo(k,i)-dt*(f1(k,i)+f2(k,i))/4; end
    u = bdc(u,nE);

    % step 3
    fm_x=mflux_x(u,v,xi,dx,nE, scheme); f3=fm_x;
    for k=1:K; u(k,i)=uo(k,i)-dt*(f1(k,i)+f2(k,i)+4*f3(k,i))/6; end 
    u = bdc(u,nE);
    
    % Compute cell averages (for output)
    u_bar=(u(1,:)+4*u(2,:)+u(3,:))/6;        
    
    % increment time and counter
    t=t+dt; it=it+1;
    
    % Plot u
    if rem(it,40)==0;
        figure(1); plot(x,u0,'-k',x,u,'-',xc,u_bar,'sr'); 
        axis(plotrange); grid on; daspect([1.5,2,1]); drawnow;
    end
end

%% Post process
% Save results to file
%mkdir('AdvectionRK33'); save('AdvectionRK33/Plot.mat','x','xc','u','ue','u_bar','dx','dt0','it','tEnd','P','cfl','nE');

% Plot solution
h=plot(x(:),ue(:),'-k',x(:),u(:),'-+r',xc(:),u_bar(:),'sb'); axis(plotrange);
legend(h,'Exact','MCV3','Cell Averages'); legend boxoff; grid on; daspect([1.5,2,1]);
title('MMC-FR','interpreter','latex','FontSize',18);
xlabel('$\it{x}$','interpreter','latex','FontSize',14);
ylabel({'$\it{u(x)}$'},'interpreter','latex','FontSize',14);