! TVD split MHD code
! copyright (C) 2001,2003, Ue-Li Pen
! written November 2001 by Ue-Li Pen, pen@cita.utoronto.ca
!This program is free software; you can redistribute it and/or
!modify it under the terms of the GNU General Public License
!as published by the Free Software Foundation; either version 2
!of the License, or (at your option) any later version.
!
!This program is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License
!along with this program; if not, write to the Free Software
!Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
!
! debugged to pass alven test March 29, 2002
!
! release 0.1 May 2003
! release 0.2 Sept 2009 : fixed b2 interpolation in cfl
!
! General Notes:
! restrictions: requires nx=nz or nz=1
! see http://arxiv.org/astroph/abs/astro-ph/0305088
! or http://www.cita.utoronto.ca/~pen/MHD
! for questions contact pen@cita.utoronto.ca
!
implicit none
integer, parameter :: n=100,nx=n,ny=1,nz=1
real, dimension(5,nx,ny,nz) :: u
real, dimension(3,nx,ny,nz) :: b
real t,dt,tf
integer iter,i
! B field is stored on the left side of each cell


call init(u,b,nx,ny,nz)
tf=nx*4
t=0
iter=0
do
   iter=iter+1
   if (t>=tf) exit
   dt=0.9*cfl(u,b)
   dt=min(dt,(tf-t)/2)
   t=t+2*dt
   call fluidx(u,b,nx,ny,nz,dt)
   call advectbyzx(u,b,nx,ny,nz,dt)
! the y sweep
   call transpose12(u,b,u,b,nx,ny,nz)
   call fluidx(u,b,ny,nx,nz,dt)
   call advectbyzx(u,b,ny,nx,nz,dt)
! z sweep
   call transpose13(u,b,u,b,ny,nx,nz)
   call fluidx(u,b,nz,nx,ny,dt)
   call advectbyzx(u,b,nz,nx,ny,dt)
   call advectbyzx(u,b,nz,nx,ny,dt)
   call fluidx(u,b,nz,nx,ny,dt)
   
! back
   call transpose13(u,b,u,b,nz,nx,ny)
   call advectbyzx(u,b,ny,nx,nz,dt)
   call fluidx(u,b,ny,nx,nz,dt)
! x again
   call transpose12(u,b,u,b,ny,nx,nz)
   call advectbyzx(u,b,nx,ny,nz,dt)
   call fluidx(u,b,nx,ny,nz,dt)
   if (mod(iter,10) .eq. 1) write(*,*) 't=',t,iter,u(2,nx/4,1,1)
end do
call output

contains
  subroutine output
    integer i,j
    open(10,file='u.dat')
    open(20,file='e.dat')
    do j=1,1
    do i=1,nx
       write(10,'(i5,7(1x,e9.3))') i,u(1,i,j,1)-1,u(2:4,i,j,1)/u(1,i,j,1),b(:,i,j,1)-(/1,0,0/)
    end do
    do i=1,nx
       write(20,'(i5,7(1x,e30.20))') i,u(5,i,j,1)-sum(u(2:4,i,j,1)**2)/u(1,i,j,1)/2-sum(b(:,i,j,1)**2)/2
    end do
    end do
    close(20)
    close(10)
  end subroutine output


  function cfl(u,b)
    real, dimension(:,:,:,:)::u,b
    real cfl
! locals
    real v,ps,p,c,gamma,bx,by,bz,b2
    integer i,j,k,ip,jp,kp

    gamma=5./3.
    c=0
    !$omp parallel do private(j,i,kp,jp,ip,bx,by,bz,v,ps,p,b2) reduction(max:c)
    do k=1,nz
       do j=1,ny
          do i=1,nx
             kp=mod(k,nz)+1
             jp=mod(j,ny)+1
             ip=mod(i,nx)+1
             bx=(b(1,i,j,k)+b(1,ip,j,k))/2
             by=(b(2,i,j,k)+b(2,i,jp,k))/2
             bz=(b(3,i,j,k)+b(3,i,j,kp))/2
             v=maxval(abs(u(2:4,i,j,k)/u(1,i,j,k)))
             b2=bx*bx+by*by+bz*bz 
             ps=(u(5,i,j,k)-sum(u(2:4,i,j,k)**2,1)/u(1,i,j,k)/2)*(gamma-1)+(2-gamma)*b2/2
             p=ps-b2/2
!   c=max(c,v+sqrt(abs(  (sum(b(:,i,j,k)**2)*(gamma-1)+gamma*p)/u(1,i,j,k))))
             c=max(c,v+sqrt(abs(  (b2*2+gamma*p)/u(1,i,j,k))))
          end do
       end do
    end do
    cfl=1/c
  end function cfl


end

subroutine init(u,b,nx,ny,nz)
  implicit none
  integer nx,ny,nz
  real u(5,nx,ny,nz),b(3,nx,ny,nz)
  real p0
  integer i
  real, parameter :: epsilon=.1
  
  p0=3./5.
  u=0
  b=0
  b(1,:,:,:)=1
  u(1,:,:,:)=1
  u(5,:,:,:)=1.5*p0
  u(2,:,1,1)=0.0001*sin( (/ (2*3.14159*i/nx,i=1,nx) /) )
  u(1,:,1,1)=u(1,:,1,1)+u(2,:,1,1)
  u(5,:,:,:)=u(5,:,:,:)+sum(b**2,1)/2+u(2,:,:,:)**2/u(1,:,:,:)/2
  u(5,:,1,1)=u(5,:,1,1)+p0*u(2,:,1,1)*5./3./(2./3.)
!    return
! circularly polarized alven wave:
! background : rho=B_x=1 all others zero
  u=0
  b=0
  b(1,:,:,:)=1
  u(1,:,:,:)=1
  u(5,:,:,:)=0.1  ! to keep things stable
  do i=1,nx
     u(3,i,:,:)=epsilon*sin( 2*3.14159*i/nx )
     u(4,i,:,:)=epsilon*cos( 2*3.14159*i/nx )
  enddo
  b(2,:,:,:)=-u(3,:,:,:)
  b(3,:,:,:)=-u(4,:,:,:)
  b(2,:,:,:)=(b(2,:,:,:)+cshift(b(2,:,:,:),-1,2))/2
  b(3,:,:,:)=(b(3,:,:,:)+cshift(b(3,:,:,:),-1,3))/2
  
  u(5,:,:,:)=u(5,:,:,:)+sum(b**2,1)/2+sum(u(2:4,:,:,:)**2,1)/u(1,:,:,:)/2
!  return
  
! alven wave:
! background : rho=B_x=1 all others zero
! \dot{v_y}+\div_x(-B_y)=0
! \dot{B_y}+\div_x (       -   v_y)     = 0
! let v_y=\epsilon sin(2 \pi (x-t)/L)
! then B_y=-v_y
  u=0
  b=0
  b(1,:,:,:)=1
  u(1,:,:,:)=1
  u(5,:,:,:)=.001  ! to keep things stable
  do i=1,nx
     u(3,i,:,:)=0.1*sin( 2*3.14159*i/nx )
  enddo
  b(2,:,:,:)=-u(3,:,:,:)
  b(2,:,:,:)=(b(2,:,:,:)+cshift(b(2,:,:,:),-1,1))/2
  u(5,:,:,:)=u(5,:,:,:)+sum(b**2,1)/2+u(3,:,:,:)**2/u(1,:,:,:)/2
 ! return
  ! magnetic
  u=0
  b=0
  b(2,:,:,:)=1
  u(1,:,:,:)=1
  u(5,:,:,:)=.000000001  ! to keep things stable
  do i=1,nx
     u(2,i,:,:)=0.001*sin( 2*3.14159*i/nx )
  enddo
  b(2,:,:,:)=b(2,:,:,:)+u(2,:,:,:)
  u(1,:,:,:)=b(2,:,:,:)
  u(5,:,:,:)=u(5,:,:,:)+sum(b**2,1)/2+sum(u(2:4,:,:,:)**2,1)/u(1,:,:,:)/2
  write(*,*) 'magnetic'
!  return
  ! magnetosonic
  p0=3./5./2.
  u=0
  b=0
  b(2,:,:,:)=1./sqrt(2.)
  u(1,:,:,:)=1
  u(5,:,:,:)=p0*1.5
  do i=1,nx
     u(2,i,:,:)=0.0001*sin( 2*3.14159*i/nx )
  enddo
  b(2,:,:,:)=b(2,:,:,:)+u(2,:,:,:)/sqrt(2.)
  u(1,:,:,:)=u(1,:,:,:)+u(2,:,:,:)
  u(5,:,:,:)=u(5,:,:,:)+sum(b**2,1)/2+sum(u(2:4,:,:,:)**2,1)/u(1,:,:,:)/2
  u(5,:,:,:)=u(5,:,:,:)+p0*u(2,:,:,:)*5./3./(2./3.)
  write(*,*) 'magnetosonic'
end subroutine init

subroutine fluidx(u,b,nx,ny,nz,dt)
  implicit none
  integer nx,ny,nz
  real u(5,nx,ny,nz),b(3,nx,ny,nz),dt
  real, dimension(3,nx) :: b3x
  real, dimension(5,nx) :: u1x
  integer j,k,jp,kp

  !$omp parallel do private(j,u1x,b3x,jp,kp)
  do k=1,nz
     do j=1,ny
        b3x=b(:,:,j,k)/2
        jp=mod(j,ny)+1
        kp=mod(k,nz)+1
        b3x(1,:)=b3x(1,:)+cshift(b3x(1,:),1)
        b3x(2,:)=b3x(2,:)+b(2,:,jp,k)/2
        b3x(3,:)=b3x(3,:)+b(3,:,j,kp)/2
        call tvd1(u(1,1,j,k),b3x,nx,dt)
     end do
  end do
end subroutine fluidx

subroutine advectbyzx(u,b,nx,ny,nz,dt)
  implicit none
  integer nx,ny,nz
  real u(5,nx,ny,nz),b(3,nx,ny,nz),dt
  real, dimension(nx) :: fluxbx,b1x,vx
  integer j,k,jm,km

  !$omp parallel do private(j,jm,vx,b1x,fluxbx)
  do k=1,nz
     do j=1,ny
        jm=mod(j+ny-2,ny)+1
        vx=(u(2,:,jm,k)+u(2,:,j,k))/(u(1,:,jm,k)+u(1,:,j,k))
        vx=(cshift(vx,-1)+cshift(vx,1)+2*vx)/4
        b1x=b(2,:,j,k)
        call tvdb(fluxbx,b1x,vx,nx,dt)
        b(2,:,j,k)=b1x
        b(1,:,j,k)=b(1,:,j,k)-cshift(fluxbx,-1)
        b(1,:,jm,k)=b(1,:,jm,k)+cshift(fluxbx,-1)
     end do
  end do
  do k=1,nz
     !$omp parallel do private(km,vx,b1x,fluxbx)
     do j=1,ny
        km=mod(k+nz-2,nz)+1
        vx=(u(2,:,j,km)+u(2,:,j,k))/(u(1,:,j,km)+u(1,:,j,k))
        vx=(cshift(vx,-1)+cshift(vx,1)+2*vx)/4
        b1x=b(3,:,j,k)
        call tvdb(fluxbx,b1x,vx,nx,dt)
        b(3,:,j,k)=b1x
        b(1,:,j,k)=b(1,:,j,k)-cshift(fluxbx,-1)
        b(1,:,j,km)=b(1,:,j,km)+cshift(fluxbx,-1)
     end do
  end do
end subroutine advectbyzx

recursive subroutine tvdb(flux,b,vg,n,dt)
  integer i,n,ip,ipp,im
  real dt
  real, dimension(n) :: flux,b,vg
! locals
  real, dimension(n) :: b1,flux1,vh
  real w,wp,wm,dw,v
  
  ! unlike the B field, the flux lives on the right cell boundary
  vh=(vg+cshift(vg,1))/2
  where(vh>0)
     flux1=b*vg
  elsewhere
     flux1=cshift(b*vg,1)
  end where
  b1=b-(flux1-cshift(flux1,-1))*dt/2
  do i=1,n
     ip=mod(i,n)+1
     ipp=mod(ip,n)+1
     im=mod(i+n-2,n)+1
     v=vh(i)
     if (v>0) then
        w=vg(i)*b1(i)
        wp=(vg(ip)*b1(ip)-w)/2
        wm=(w-vg(im)*b1(im))/2
     else
        w=vg(ip)*b1(ip)
        wp=(w-vg(ipp)*b1(ipp))/2
        wm=(vg(i)*b1(i)-w)/2
     end if
     dw=0
     if(wm*wp>0) dw=2*wm*wp/(wm+wp)
     flux(i)=(w+dw)*dt
  end do
  b=b-(flux-cshift(flux,-1))
end subroutine tvdb

recursive subroutine tvd1(u,b,n,dt)
  implicit none
  integer n
  real dt
  real u(5,n),b(3,n)
!locals    
  real, dimension(5,n) :: v, u1,wr,wl,fr,fl,flux,dfrp,dfrm&
       ,dflm,dflp,dfl,dfr
  real c

  call mhdflux(v,c,u,b,n)
  wr=u+v
  wl=u-v
  fr=c*wr
  fl=cshift(c*wl,1,2)
  flux=(fr-fl)/2
  u1=u-(flux-cshift(flux,-1,2))*dt/2
  call mhdflux(v,c,u1,b,n)
  wr=u1+v
  wl=u1-v
  fr=c*wr
  dfrp=(cshift(fr,1,2)-fr)/2
  dfrm=(fr-cshift(fr,-1,2))/2
  dfr=0
  where(dfrp*dfrm>0)
     dfr=2*dfrp*dfrm/(dfrp+dfrm)
  end where
  fl=cshift(c*wl,1,2)
  dflp=(fl-cshift(fl,1,2))/2
  dflm=(cshift(fl,-1,2)-fl)/2
  dfl=0
  where(dflp*dflm>0)
     dfl=2*dflp*dflm/(dflp+dflm)
  end where
  flux=(fr-fl+(dfr-dfl))/2
  u=u-(flux-cshift(flux,-1,2))*dt
  return
end subroutine tvd1



recursive subroutine mhdflux(v,c,u,b,n)
  implicit none
  integer n
  real, dimension(5,n)::v,u
  real b(3,n)
  real, intent(out) :: c
  ! locals
  real, dimension(n) :: vx,ps,p
  real gamma

  gamma=5./3.

  vx=u(2,:)/u(1,:)
  ps=(u(5,:)-sum(u(2:4,:)**2,1)/u(1,:)/2)*(gamma-1)+(2-gamma)*sum(b**2,1)/2
  v(1,:)=u(2,:)
  v(2,:)=u(2,:)*vx+ps-b(1,:)**2
  v(3,:)=u(3,:)*vx-b(2,:)*b(1,:)
  v(4,:)=u(4,:)*vx-b(3,:)*b(1,:)
  v(5,:)=(u(5,:)+ps)*vx-b(1,:)*sum(b*u(2:4,:),1)/u(1,:)
  p=ps-sum(b**2,1)/2
  c=maxval(abs(vx)+sqrt(abs( (sum(b**2,1)+gamma*p)/u(1,:))))
  if (c>0) v=v/c
end subroutine mhdflux



subroutine transpose12(ut,bt,u,b,nx,ny,nz)
  implicit none
  integer nx,ny,nz
  real u(5,nx,ny,nz),b(3,nx,ny,nz)
  real ut(5,ny,nx,nz),bt(3,ny,nx,nz)

  integer i,j,k
  real u2(5,ny,nx),b2(3,ny,nx)

!$omp parallel do default(none) shared(u,b,ut,bt,ny,nx,nz) private(i,j,u2,b2)
  do k=1,nz
     do j=1,ny
        do i=1,nx
           u2(:,j,i)=u((/1,3,2,4,5/),i,j,k)
           b2(:,j,i)=b((/2,1,3/),i,j,k)
        end do
     end do
     ut(:,:,:,k)=u2
     bt(:,:,:,k)=b2
  end do
end subroutine transpose12


subroutine transpose13(ut,bt,u,b,nx,ny,nz)
  implicit none
  integer nx,ny,nz
  real u(5,nx,ny,nz),b(3,nx,ny,nz)
  real ut(5,nz,ny,nx),bt(3,nz,ny,nx)

  integer i,j,k
  real u2(5,nz,nx),b2(3,nz,nx)

  if (nx .eq. nz) then
!$omp parallel do default(none) shared(u,b,ut,bt,ny,nx,nz) private(i,k,u2,b2)
  do j=1,ny
     do k=1,nz
        do i=1,nx
           u2(:,k,i)=u((/1,4,3,2,5/),i,j,k)
           b2(:,k,i)=b((/3,2,1/),i,j,k)
        end do
     end do
     ut(:,:,j,:)=u2
     bt(:,:,j,:)=b2
  end do
  else if (nz .eq. 1) then
     call transpose12(ut,bt,u,b,nx,ny,nz)
     do i=1,nx
        do j=1,ny
           ut(:,1,j,i)=ut((/1,4,2,3,5/),1,j,i)
           bt(:,1,j,i)=bt((/3,1,2/),1,j,i)
        end do
     end do
  else if (nx .eq. 1) then
     call transpose12(ut,bt,u,b,ny,nz,nx)
     do i=1,ny
        do j=1,nz
           ut(:,j,i,1)=ut((/1,4,2,3,5/),j,i,1)
           bt(:,j,i,1)=bt((/3,1,2/),j,i,1)
        end do
     end do
  else
     write(*,*) 'nz<>nx not supported'
     stop
  endif
end subroutine transpose13

