function output = bdc(u,nx)

output=u; output(:,1)=output(:,nx-1); output(:,nx)=output(:,2);
