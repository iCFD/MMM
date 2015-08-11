function mm_FR(u,v,fm_x,xi,dx,nx,ns,scheme)
%     MCV3 & MCV3_UPCC 
%      IMPLICIT REAL*8(A-H,O-Z)
%      real*8 u(0:nx+1,ns),v(0:nx+1,ns),fm_x(0:nx+1,ns)
%      real*8 f(0:nx+1,ns),fb(0:nx+1),fb_x(0:nx+1)
%      real*8 xi(ns),dx
%      integer nx,ns
 
      for i=0:nx+1
       for k=1:ns
       f(i,k)=v(i,k)*u(i,k)
       end
      end
 
      call Riemann(u,v,fb,fb_x,xi,dx,nx,ns)

      for i=1,nx
          switch scheme
              case 1 % MCV3
                  fm_x(i,1)=2/dx*fb_x(i);
                  fm_x(i,2)=2/dx*(3*(fb(i+1)-fb(i))-fb_x(i)-fb_x(i+1))/4;
                  fm_x(i,3)=2/dx*fb_x(i+1);
                  
              case 2 % MCV3_UPCC
                  fm_x(i,1)=2/dx*(2*(f(i,1)+f(i,2))-0.5*(7*fb(i)+fb(i+1)));
                  fm_x(i,2)=1/dx*(f(i,3)-f(i,1));
                  fm_x(i,3)=2/dx*(-2*(f(i,2)+f(i,3))+0.5*(7*fb(i+1)+fb(i)));
                  
              case 3 % MCV3_CPCC
                  s3=sqrt(3);
                  fm_x(i,1)=2/dx*(s3*(3/4*f(i,1)+5/6*f(i,2)-1/12*f(i,3))+3/4*(1.5-s3)*fb(i+1)-3/4*(1.5+s3)*fb(i));
                  fm_x(i,2)=2/dx*(f(i,3)-f(i,1))*s3/3;
                  fm_x(i,3)=2/dx*(s3*(1/12*f(i,1)-5/6*f(i,2)-3/4*f(i,3))+3/4*(1.5+s3)*fb(i+1)-3/4*(1.5-s3)*fb(i));
          end
      end

      call bdc(fm_x,nx,ns)
      
      end

      subroutine bdc(u,nx,ns)
      IMPLICIT REAL*8(A-H,O-Z)
      real*8 u(0:nx+1,ns)
         do k=1,ns
            u(0,k) = u(nx,k)
            u(nx+1,k) = u(1,k)
         enddo

      return
      end