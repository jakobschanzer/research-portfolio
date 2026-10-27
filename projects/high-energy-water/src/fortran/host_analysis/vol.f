c ===== vol.f : per-frame plane-based volume; reads "vol" ================
c Fixed-form F77; each unit declares its own IMPLICIT/PARAMETER/COMMON

c ---------------------- helper: minimal image ---------------------------
      subroutine minimg(dx,dy,dz,bx,by,bz)
      implicit real*8 (a-h,o-z)
      real*8 dx,dy,dz,bx,by,bz
      dx = dx - bx*dnint(dx/bx)
      dy = dy - by*dnint(dy/by)
      dz = dz - bz*dnint(dz/bz)
      return
      end

c ---------------------- init: read "vol" --------------------------------
      subroutine vol_init(dr)
      implicit real*8 (a-h,o-z)
      integer NPLANES
      parameter (NPLANES=6)
c types for COMMON members
      real*8 nrm, poff, scheck
      integer ichk, ipl
      logical have_idx
c COMMON declarations (must appear in each subprogram that uses them)
      common /vol_planes/ nrm(3,NPLANES), poff(NPLANES), scheck(NPLANES)
      common /vol_idx/    ichk(NPLANES), ipl(3,NPLANES), have_idx
c locals
      integer i,j,u
      integer tmp(3*NPLANES)
      character*256 line

c init
      have_idx = .false.
      do i=1,NPLANES
         ichk(i) = 0
         do j=1,3
            ipl(j,i) = 0
         enddo
      enddo

      u = 23
      open(unit=u,file='vol',status='old',err=99)

c find "[ checks ]"
 10   continue
      read(u,'(A)',end=98,err=98) line
      if (index(line,'checks').eq.0) goto 10
      read(u,*,end=98,err=98) (ichk(i), i=1,NPLANES)

c find "[ planes ]"
 20   continue
      read(u,'(A)',end=98,err=98) line
      if (index(line,'planes').eq.0) goto 20

c read 18 integers: 3 per plane
      read(u,*,end=98,err=98) (tmp(i), i=1,3*NPLANES)
      do i=1,NPLANES
         ipl(1,i) = tmp(3*(i-1)+1)
         ipl(2,i) = tmp(3*(i-1)+2)
         ipl(3,i) = tmp(3*(i-1)+3)
      enddo

      close(u)
      have_idx = .true.
      return

 98   continue
      close(u)
 99   continue
c leave have_idx = .false. to disable volume if file missing
      return
      end

c --------------- per-frame update: build planes from coords ------------
      subroutine vol_update_frame(ntot,id,x,y,z,bx,by,bz)
      implicit real*8 (a-h,o-z)
      integer NPLANES
      parameter (NPLANES=6)
c types for COMMON members
      real*8 nrm, poff, scheck
      integer ichk, ipl
      logical have_idx
c COMMON declarations
      common /vol_planes/ nrm(3,NPLANES), poff(NPLANES), scheck(NPLANES)
      common /vol_idx/    ichk(NPLANES), ipl(3,NPLANES), have_idx
c arguments & locals
      integer ntot,id(*)
      real*8 x(*),y(*),z(*), bx,by,bz
      integer i,a1,a2,a3,ac
      real*8 x1(3),x2(3),x3(3), xc(3)
      real*8 v1(3),v2(3), nn(3), mag, p

      if (.not.have_idx) return

      do i=1,NPLANES
         a1 = ipl(1,i)
         a2 = ipl(2,i)
         a3 = ipl(3,i)
         ac = ichk(i)

c fetch positions (1-based IDs)
         x1(1)=x(a1)
         x1(2)=y(a1)
         x1(3)=z(a1)
         x2(1)=x(a2)
         x2(2)=y(a2)
         x2(3)=z(a2)
         x3(1)=x(a3)
         x3(2)=y(a3)
         x3(3)=z(a3)
         xc(1)=x(ac)
         xc(2)=y(ac)
         xc(3)=z(ac)

c minimal-image differences for edges
         v1(1)=x3(1)-x1(1)
         v1(2)=x3(2)-x1(2)
         v1(3)=x3(3)-x1(3)
         v2(1)=x2(1)-x1(1)
         v2(2)=x2(2)-x1(2)
         v2(3)=x2(3)-x1(3)
         call minimg(v1(1),v1(2),v1(3),bx,by,bz)
         call minimg(v2(1),v2(2),v2(3),bx,by,bz)

c normal = cross(v2, v1)  (matches your C++ order)
         nn(1)= v2(2)*v1(3) - v2(3)*v1(2)
         nn(2)= v2(3)*v1(1) - v2(1)*v1(3)
         nn(3)= v2(1)*v1(2) - v2(2)*v1(1)

c normalize
         mag = dsqrt(nn(1)*nn(1)+nn(2)*nn(2)+nn(3)*nn(3))
         if (mag.gt.0.d0) then
            nn(1)=nn(1)/mag
            nn(2)=nn(2)/mag
            nn(3)=nn(3)/mag
         endif

c plane:  n·x + p = 0  through x1
         p = -(nn(1)*x1(1)+nn(2)*x1(2)+nn(3)*x1(3))

         nrm(1,i)=nn(1)
         nrm(2,i)=nn(2)
         nrm(3,i)=nn(3)
         poff(i) = p

c orient with check atom: inside if same signed side as check
         scheck(i) = sign(1.d0,
     &        nn(1)*xc(1)+nn(2)*xc(2)+nn(3)*xc(3)+p)
      enddo

      return
      end

c ---------------------- point membership test --------------------------
      logical function vol_inside(xp,yp,zp, bx,by,bz, xc,yc,zc,
     &                            hvx,hvy,hvz)
      implicit real*8 (a-h,o-z)
      integer NPLANES
      parameter (NPLANES=6)
c types for COMMON members
      real*8 nrm, poff, scheck
      integer ichk, ipl
      logical have_idx
c COMMON declarations
      common /vol_planes/ nrm(3,NPLANES), poff(NPLANES), scheck(NPLANES)
      common /vol_idx/    ichk(NPLANES), ipl(3,NPLANES), have_idx
c locals
      integer i
      real*8 s

      vol_inside = .false.
      if (.not.have_idx) return

c inside if point is on the same side as the check for all planes
      do i=1,NPLANES
         s = nrm(1,i)*xp + nrm(2,i)*yp + nrm(3,i)*zp + poff(i)
         if ( s * scheck(i) .lt. 0.d0 ) return
      enddo
      vol_inside = .true.
      return
      end

c -------------------------- reported volume ----------------------------
      double precision function vol_value()
      implicit real*8 (a-h,o-z)
      integer NPLANES
      parameter (NPLANES=6)
c types for COMMON members (not required here, but consistent)
      real*8 nrm, poff, scheck
      integer ichk, ipl
      logical have_idx
      common /vol_planes/ nrm(3,NPLANES), poff(NPLANES), scheck(NPLANES)
      common /vol_idx/    ichk(NPLANES), ipl(3,NPLANES), have_idx

c Keep -1 to leave janal’s printed vbind unchanged for now.
      vol_value = -1.d0
      return
      end
c =======================================================================
