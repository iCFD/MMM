subroutine Riemann(u,v,fb,fb_x,xi,dx,nx,ns)

      do i=1,nx
       do k=1,ns
       fe(k)=v(i,k)*u(i,k)
       enddo
       f_l(i)=f_primary(-1,xi,ns,fe)
       f_r(i)=f_primary( 1,xi,ns,fe)
       fx_l(i)=fx_primary(-1,xi,ns,fe)
       fx_r(i)=fx_primary( 1,xi,ns,fe)
      enddo

      call bdc_f(f_l,nx)
      call bdc_f(f_r,nx)
      call bdc_f(fx_l,nx)
      call bdc_f(fx_r,nx)

      do i=1,nx
       if(v(i,1).ge.0.0) then
       fb(i)=f_r(i-1)
       fb_x(i)=fx_r(i-1)
       else
       fb(i)=f_l(i)
       fb_x(i)=fx_l(i)
       endif
      enddo
      call bdc_f(fb,nx)
      call bdc_f(fb_x,nx)

      return
      end

      subroutine bdc_f(f,nx)
      IMPLICIT REAL*8(A-H,O-Z)
      real*8 f(0:nx+1)
       f(0) = f(nx)
       f(nx+1) = f(1)
       return
       end
       
%       real*8 function f_gauss_c(x,f1,f2,f3)
%       IMPLICIT REAL*8(A-H,O-Z)     
%       cc=dsqrt(3.0d0)/2.
%       f_guass_c=2./3.*x*(x-cc)*f1-4./3.*(x+cc)*(x-cc)*f2 + 2./3.*(x+cc)*x*f3
%       return
%       end
% 
%       real*8 function fx_gauss_c(x,f1,f2,f3)
%       IMPLICIT REAL*8(A-H,O-Z)     
%       cc=dsqrt(3.0d0)/3.
%       fx_guass_c=(4./3.*x-cc)*f1-8./3.*x*f2+(4./3.*x+cc)*f3
%       return
%       end

      real*8 function f_primary(x,xi,ns,f)
      IMPLICIT REAL*8(A-H,O-Z)
      real*8 xi(ns),f(ns),p(ns)

      do j=1,ns 
      p(j)=1.0
      k=j
      do j1=1,ns 
      if (j1.ne.k) then
      p(j)=p(j)*(x-xi(j1))/(xi(k)-xi(j1))
      end if
      end do
      end do

      f_primary=0.0
      do j=1,ns 
       f_primary=f_primary+p(j)*f(j)
      end do
      return
      end

      real*8 function fx_primary(x,xi,ns,f)
      IMPLICIT REAL*8(A-H,O-Z)
      real*8 xi(ns),f(ns),p(ns),q(ns)

      do j=1,ns 
      p(j)=0.0
      q(j)=1.0
      k=j
      do j1=1,ns 
      if (j1.ne.k) then
      p(j)=p(j)+(x-xi(j1))
      q(j)=q(j)/(xi(k)-xi(j1))
      end if
      end do
      end do

      fx_primary=0.0
      do j=1,ns 
      fx_primary=fx_primary+p(j)/q(j)*f(j)
      end do

      fx_primary=(2.*x-xi(2)-xi(3))/(xi(1)-xi(2))/(xi(1)-xi(3))*f(1)...
                +(2.*x-xi(1)-xi(3))/(xi(2)-xi(1))/(xi(2)-xi(3))*f(2)...
                +(2.*x-xi(1)-xi(2))/(xi(3)-xi(1))/(xi(3)-xi(2))*f(3);
      return
      end