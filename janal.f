      program HostAnalysis
      implicit real*8 (a-h,o-z)
      character*1 line
      logical rev

c ==== VOL-API (only addition) ==========================================
      external vol_init, vol_inside, vol_value, vol_zone
      logical vol_inside
      integer vol_zone
c --- added so vol.f can update per-frame geometry ----------------------
      external vol_update_frame
c =======================================================================

      dimension x(100000),y(100000),z(100000)
      dimension id(100000),ib(0:40),jb(2,0:20)
      dimension sii(200),eii(200),sjj(20),ejj(20)
      dimension aij(20,200),bij(20,200),rex(200)
      dimension rww(20,200),ess(20,200)
      dimension ced(40,108,54),cd(40,108,54),cdc(108,54)
      dimension vn1(108,54),vn2(108,54)
      dimension ein(40),rtg(0:54),rav(108,54)
      dimension ci(2,40),cii(2,40)
      dimension cei(2,40),ceii(2,40)
      dimension gg(108,54),gh(108,54),gh2(108,54)
      dimension gn(108,54),gn2(108,54)
      dimension mba(0:100),mbd(0:100,2)
      dimension nhb(100000),nhb2(100000)
      dimension non(100000),non2(100000)
      dimension po(0:200),ph(0:200)
      dimension hd(108,54)

      common /ran0/ irn

      open(unit=10,file='drive',status='old')
      open(unit=11,file='temperature',status='old')
      open(unit=20,file='top',status='old')
      open(unit=21,file='top-bind',status='old')
      open(unit=22,file='top-hb',status='old')
      open(unit=30,file='pot',status='old')
      open(unit=31,file='pot-sol',status='old')
      open(unit=40,file='misc.out',status='unknown')
      open(unit=41,file='muzr.out',status='unknown')
      open(unit=42,file='enzr.out',status='unknown')
      open(unit=43,file='dist.out',status='unknown')
      open(unit=44,file='fluc.out',status='unknown')
      open(unit=50,file='data.out',status='unknown')

c - begin: set variables
      pi=2.d0*dasin(1.d0)
      vin=2.d0*(13.5d0**3.d0)*pi
      dr=0.25d0
c ==== VOL-API: initialize volume module with dr =========================
      call vol_init(dr)
c =======================================================================
      vav=0.d0
      anout=0.d0
      avin=0.d0
      a2vin=0.d0
      ahbin=0.d0
      vn1b=0.d0
      vn2b=0.d0
      do i=1,200
        rex(i)=0.d0
      enddo
      do j=1,108
        do k=1,54
          hd(j,k)=0.d0
          vn1(j,k)=0.d0
          vn2(j,k)=0.d0
          cdc(j,k)=0.d0
          gg(j,k)=0.d0
          gn(j,k)=0.d0
          gn2(j,k)=0.d0
          gh(j,k)=0.d0
          gh2(j,k)=0.d0
          do i=1,40
            cd(i,j,k)=0.d0
            ced(i,j,k)=0.d0
          enddo
        enddo
      enddo
      do i=1,40
        cei(1,i)=0.d0
        cei(2,i)=0.d0
      enddo
      rtg(0)=-dr/2.d0
      do i=1,54
        rtg(i)=rtg(i-1)+dr
      enddo
      rtg(0)=0.d0
      do i=0,100
        mba(i)=0
        mbd(i,1)=0
        mbd(i,2)=0
      enddo
      do i=0,200
        po(i)=0.d0
        ph(i)=0.d0
      enddo
c - end: set variables

      read(10,*) nskip
      read(10,*) nconf
      read(10,*) nin1
      read(10,*) nin2
      read(10,*) nin3
      read(10,*) rc
      read(10,*) irn
      read(10,*) kji
      read(10,*) i
      read(10,*) rhb
      read(10,*) hang
      read(10,*) rop
      read(10,*) rf
      hang=dcos(hang*pi/180.d0)
      nins=nin1+nin2+nin3
      nn1=nin1
      nn2=nin1+nin2
      rc2=rc*rc
c     rf2=rf*rf
      rev=i.eq.1

      read(11,*) temp
      temp=temp+273.15d0
      write(*,900) temp

      read(20,*) ntot
      do i=1,ntot
        read(20,*) j,id(i)
      enddo

      read(30,*) ntyp
      do i=1,ntyp
        read(30,*) j,sii(i),eii(i),rex(i)
        rex(i)=(rex(i)+rf)**2.d0
      enddo

      read(31,*) nsol
      if (nsol.gt.10) nsol=10
      do i=1,nsol
        read(31,*) j,sjj(i),ejj(i)
      enddo
      nsol2=2*nsol

      do i=1,nsol
        do j=1,ntyp
          if (kji.eq.1) then
            ss=(sii(j)+sjj(i))/2.d0
          else
            ss=dsqrt(sii(j)*sjj(i))
          endif
          rww(i,j)=ss*ss*(2.d0**(1.d0/3.d0))
          ess(i,j)=dsqrt(eii(j)*ejj(i))/temp
          bij(i,j)=4.d0*ess(i,j)*(ss**6.d0)
          aij(i,j)=bij(i,j)*(ss**6.d0)
        enddo
      enddo

      read(21,*) nbt
      nbt=1
      if (nbt.eq.1) then
        k=0
        do i=1,2
          read(21,*) jb(i,0)
          if (jb(i,0).gt.20) then
            write(*,*) 'reconsider host orientation shit'
            stop
          endif
          do j=1,jb(i,0)
            read(21,*) jb(i,j)
            k=k+1
            ib(k)=jb(i,j)
          enddo
        enddo
        ib(0)=k
        read(21,*) icf
        read(21,*) izl,izh
        read(21,*) irm
c ----- NOTE: The cylinder-derived values are still computed here
c ----- but will be overridden by vol_value()/vol_inside for geometry.
        zz=dr*dble(1+izh-izl)
        rr=dr*dble(irm)
        vbind=pi*rr*rr*zz
        rirm=dr*dble(irm)
        zizl=dr*dble(izl-1)-13.5d0
        zizh=dr*dble(izh)-13.5d0
c     else
      endif

      read(22,*) k
      if (k.gt.0) then
        mba(0)=k
        do i=1,k
          read(22,*) mba(i)
        enddo
      endif
      read(22,*) k
      if (k.gt.0) then
        mbd(0,1)=k
        do i=1,k
          read(22,*) mbd(i,1),mbd(i,2)
        enddo
      endif

      read(*,500) line
      if (nskip.gt.0) then
        do ns=1,nskip
          if (mod(ns,1000).eq.0) write(*,*) ns
          read(*,*)(x(i),y(i),z(i),i=1,ntot),bx,by,bz
        enddo
      endif

c- begin: analyze configurations
      do nc=1,nconf
        write(*,*) nc
        read(*,*)(x(i),y(i),z(i),i=1,ntot),bx,by,bz
        do i=1,ntot
          x(i)=x(i)-bx*dnint(x(i)/bx)
          y(i)=y(i)-by*dnint(y(i)/by)
          z(i)=z(i)-bz*dnint(z(i)/bz)
        enddo
        vol=bx*by*bz
        vav=vav+vol

c- begin: center the host
        rr=1.d20
        do i=1,icf
          xx=bx*(ranu(0)-0.5d0)
          yy=by*(ranu(0)-0.5d0)
          zz=bz*(ranu(0)-0.5d0)
          r2=0.d0
          do j=1,ib(0)
            dx=dabs(xx-x(ib(j)))
            dy=dabs(yy-y(ib(j)))
            dz=dabs(zz-z(ib(j)))
            dx=dmin1(dx,bx-dx)
            dy=dmin1(dy,by-dy)
            dz=dmin1(dz,bz-dz)
            r2=r2+dx*dx+dy*dy+dz*dz
          enddo
          if (r2.lt.rr) then
            rr=r2
            xc=xx
            yc=yy
            zc=zz
          endif
        enddo

        do i=1,ntot
          x(i)=x(i)-xc
          y(i)=y(i)-yc
          z(i)=z(i)-zc
          x(i)=x(i)-bx*dnint(x(i)/bx)
          y(i)=y(i)-by*dnint(y(i)/by)
          z(i)=z(i)-bz*dnint(z(i)/bz)
        enddo
        xc=0.d0
        yc=0.d0
        zc=0.d0
        do i=1,ib(0)
          xc=xc+x(ib(i))
          yc=yc+y(ib(i))
          zc=zc+z(ib(i))
        enddo
        xc=xc/dble(ib(0))
        yc=yc/dble(ib(0))
        zc=zc/dble(ib(0))
        do i=1,ntot
          x(i)=x(i)-xc
          y(i)=y(i)-yc
          z(i)=z(i)-zc
          x(i)=x(i)-bx*dnint(x(i)/bx)
          y(i)=y(i)-by*dnint(y(i)/by)
          z(i)=z(i)-bz*dnint(z(i)/bz)
        enddo
c- end: center the host

c- begin: define host vector (unchanged)
        x1=0.d0
        y1=0.d0
        z1=0.d0
        x2=0.d0
        y2=0.d0
        z2=0.d0
        do i=1,jb(1,0)
          x1=x1+x(jb(1,i))
          y1=y1+y(jb(1,i))
          z1=z1+z(jb(1,i))
        enddo
        do i=1,jb(2,0)
          x2=x2+x(jb(2,i))
          y2=y2+y(jb(2,i))
          z2=z2+z(jb(2,i))
        enddo
        x1=x1/dble(jb(1,0))
        y1=y1/dble(jb(1,0))
        z1=z1/dble(jb(1,0))
        x2=x2/dble(jb(2,0))
        y2=y2/dble(jb(2,0))
        z2=z2/dble(jb(2,0))
        dx=x2-x1
        dy=y2-y1
        dz=z2-z1
        rr=dsqrt(dx*dx+dy*dy+dz*dz)
        hvx=dx/rr
        hvy=dy/rr
        hvz=dz/rr
c- end: define host vector

c ---- added: update per-frame volume geometry from coordinates ----------
        call vol_update_frame(ntot,id,x,y,z,bx,by,bz)
c -----------------------------------------------------------------------

c- begin: water distribution analysis
        do i=1,ntot
          nhb(i)=0
          nhb2(i)=0
          non(i)=0
          non2(i)=0
        enddo
        do i=1,ntot
          if (id(i).ne.1) goto 233
          dx1=x(i+1)-x(i)
          dy1=y(i+1)-y(i)
          dz1=z(i+1)-z(i)
          dx2=x(i+2)-x(i)
          dy2=y(i+2)-y(i)
          dz2=z(i+2)-z(i)
          dx1=dx1-bx*dnint(dx1/bx)
          dy1=dy1-by*dnint(dy1/by)
          dz1=dz1-bz*dnint(dz1/bz)
          dx2=dx2-bx*dnint(dx2/bx)
          dy2=dy2-by*dnint(dy2/by)
          dz2=dz2-bz*dnint(dz2/bz)
          rr=dsqrt(dx1*dx1+dy1*dy1+dz1*dz1)
          dx1=dx1/rr
          dy1=dy1/rr
          dz1=dz1/rr
          rr=dsqrt(dx2*dx2+dy2*dy2+dz2*dz2)
          dx2=dx2/rr
          dy2=dy2/rr
          dz2=dz2/rr
          do j=1,ntot
            if (i.eq.j) goto 235
            if (id(j).ne.1) goto 235
            dx=x(j)-x(i)
            dy=y(j)-y(i)
            dz=z(j)-z(i)
            dx=dx-bx*dnint(dx/bx)
            dy=dy-by*dnint(dy/by)
            dz=dz-bz*dnint(dz/bz)
            rr=dsqrt(dx*dx+dy*dy+dz*dz)
            if (rr.lt.rop) non(i)=non(i)+1
            if (rr.gt.rhb) goto 235
            dx=dx/rr
            dy=dy/rr
            dz=dz/rr
            dot=dx1*dx+dy1*dy+dz1*dz
            if (dot.gt.hang) then
              nhb(i)=nhb(i)+1
              nhb(j)=nhb(j)+1
            endif
            dot=dx2*dx+dy2*dy+dz2*dz
            if (dot.gt.hang) then
              nhb(i)=nhb(i)+1
              nhb(j)=nhb(j)+1
            endif
235         continue
          enddo
          kkk=mba(0)
          if (kkk.ne.0) then
            do jj=1,kkk
              j=mba(jj)
              dx=x(j)-x(i)
              dy=y(j)-y(i)
              dz=z(j)-z(i)
              dx=dx-bx*dnint(dx/bx)
              dy=dy-by*dnint(dy/by)
              dz=dz-bz*dnint(dz/bz)
              rr=dsqrt(dx*dx+dy*dy+dz*dz)
              if (rr.lt.rop) non2(i)=non2(i)+1
              if (rr.gt.rhb) goto 245
              dx=dx/rr
              dy=dy/rr
              dz=dz/rr
              dot=dx*dx1+dy*dy1+dz*dz1
              if (dot.gt.hang) nhb2(i)=nhb2(i)+1
              dot=dx*dx2+dy*dy2+dz*dz2
              if (dot.gt.hang) nhb2(i)=nhb2(i)+1
245           continue
            enddo
          endif
          kkk=mbd(0,1)
          if (kkk.ne.0) then
            do jj=1,kkk
              j1=mbd(jj,1)
              j2=mbd(jj,2)
              dx=x(j2)-x(j1)
              dy=y(j2)-y(j1)
              dz=z(j2)-z(j1)
              dx=dx-bx*dnint(dx/bx)
              dy=dy-by*dnint(dy/by)
              dz=dz-bz*dnint(dz/bz)
              rr=dsqrt(dx*dx+dy*dy+dz*dz)
              dx=dx/rr
              dy=dy/rr
              dz=dz/rr
              dxx=x(i)-x(j1)
              dyy=y(i)-y(j1)
              dzz=z(i)-z(j1)
              dxx=dxx-bx*dnint(dxx/bx)
              dyy=dyy-by*dnint(dyy/by)
              dzz=dzz-bz*dnint(dzz/bz)
              rr=dsqrt(dxx*dxx+dyy*dyy+dzz*dzz)
              if (rr.gt.rhb) goto 255
              dxx=dxx/rr
              dyy=dyy/rr
              dzz=dzz/rr
              dot=dx*dxx+dy*dyy+dz*dzz
              if (dot.gt.hang) nhb2(i)=nhb2(i)+1
255           continue
            enddo
          endif
233       continue
        enddo

        j=0
        k=0
        l=0
        nhh=0
        noo=0
        do i=1,ntot
          if (id(i).ne.1) goto 15
          dd=x(i)*x(i)+y(i)*y(i)+z(i)*z(i)
          zz=hvx*x(i)+hvy*y(i)+hvz*z(i)
          rr=dsqrt(dd-zz*zz)
c ======= VOL-API: use your custom volume membership =====================
c Use zero center since coordinates have been re-centered to COM.
          if (vol_inside(x(i),y(i),z(i), bx,by,bz, 0.d0,0.d0,0.d0,
     &                   hvx,hvy,hvz)) then
            k=k+1
            l=l+nhb(i)
          endif
c =======================================================================
          if ((dabs(zz).gt.13.5d0).or.(rr.gt.13.5d0)) then
            nhh=nhh+dble(nhb(i))
            noo=noo+dble(non(i))
            j=j+1
            goto 15
          endif
          iz=1+idint((zz+13.5d0)/dr)
          ir=1+idint(rr/dr)
          gg(iz,ir)=gg(iz,ir)+1.d0
          gh(iz,ir)=gh(iz,ir)+dble(nhb(i))
          gh2(iz,ir)=gh2(iz,ir)+dble(nhb(i)+nhb2(i))
          gn(iz,ir)=gn(iz,ir)+dble(non(i))
          gn2(iz,ir)=gn2(iz,ir)+dble(non(i)+non2(i))
15        continue
        enddo
        anhh=anhh+dble(nhh)
        anoo=anoo+dble(noo)
        anout=anout+dble(j)
        avin=avin+dble(k)
        a2vin=a2vin+dble(k*k)
        ahbin=ahbin+dble(l)
        po(k)=po(k)+1.d0
        ph(l)=ph(l)+1.d0
c- end: water distribution analysis

c- begin: particle insertions
        do i=1,nsol2
          cii(1,i)=0.d0
          cii(2,i)=0.d0
          ceii(1,i)=0.d0
          ceii(2,i)=0.d0
        enddo
        vn1bt=0.d0
        vn2bt=0.d0
        do ni=1,nins
          izone=1
          if (ni.gt.nn1) izone=2
          if (ni.gt.nn2) izone=3
10        xx=bx*(ranu(0)-0.5d0)
          yy=by*(ranu(0)-0.5d0)
          zz=bz*(ranu(0)-0.5d0)
          dd=xx*xx+yy*yy+zz*zz
          zi=hvx*xx+hvy*yy+hvz*zz
          ri=dsqrt(dd-zi*zi)
c ======= VOL-API: zone selection (binding vs outside) ===================
c Default mid-zone
          j=2
c Outside box -> bulk
          if ((dabs(zi).gt.13.5d0).or.(ri.gt.13.5d0)) j=1
c If inside custom volume, set to binding (3); otherwise leave as 1 or 2.
          if (vol_inside(xx,yy,zz, bx,by,bz, 0.d0,0.d0,0.d0,
     &                   hvx,hvy,hvz)) j=3
c =======================================================================
          if (j.ne.izone) goto 10
          izz=1+idint((zi+13.5d0)/dr)
          irr=1+idint(ri/dr)
          do i=1,nsol2
            ein(i)=0.d0
          enddo
          nwi=0
          iov=1
          do i=1,ntot
            k=id(i)
            if (k.eq.0) goto 20
            dx=dabs(x(i)-xx)
            dy=dabs(y(i)-yy)
            dz=dabs(z(i)-zz)
            dx=dmin1(dx,bx-dx)
            dy=dmin1(dy,by-dy)
            dz=dmin1(dz,bz-dz)
            r2=dx*dx+dy*dy+dz*dz
            if ((r2.lt.rex(1)).and.(k.eq.1)) nwi=nwi+1
            if ((k.gt.1).and.(r2.lt.rex(k))) iov=0
            if (r2.lt.rc2) then
              r6=r2**(-3.d0)
              do j=1,nsol
                eij=r6*(aij(j,k)*r6-bij(j,k))
                ein(j)=ein(j)+eij
                if (r2.lt.rww(j,k)) then
                  eij=eij+ess(j,k)
                  ein(j+nsol)=ein(j+nsol)+eij
                endif
              enddo
            endif
20          continue
          enddo
          k=0
          if (izone.eq.1) k=1
          if (izone.eq.3) k=2
          if (k.eq.0) goto 30
          do i=1,nsol2
            f=vol*dexp(-ein(i))
            cii(k,i)=cii(k,i)+f
            ceii(k,i)=ceii(k,i)+ein(i)*f
          enddo
30        continue
          if (izone.eq.1) then
            vn1bt=vn1bt+dble(nwi)
            vn2bt=vn2bt+dble(nwi*nwi)
            goto 40
          endif
          cdc(izz,irr)=cdc(izz,irr)+1.d0
          hd(izz,irr)=hd(izz,irr)+dble(iov)
          vn1(izz,irr)=vn1(izz,irr)+dble(nwi)
          vn2(izz,irr)=vn2(izz,irr)+dble(nwi*nwi)
          do i=1,nsol2
            f=vol*dexp(-ein(i))
            cd(i,izz,irr)=cd(i,izz,irr)+f
            ced(i,izz,irr)=ced(i,izz,irr)+ein(i)*f
          enddo
40        continue
        enddo

        vn1b=vn1b+vn1bt/dble(nin1)
        vn2b=vn2b+vn2bt/dble(nin1)
        do i=1,nsol2
          cii(1,i)=cii(1,i)/dble(nin1)
          cii(2,i)=cii(2,i)/dble(nin2)
          ceii(1,i)=ceii(1,i)/dble(nin1)
          ceii(2,i)=ceii(2,i)/dble(nin2)
        enddo

        write(50,510) nc,vol
        write(50,520) (cii(1,i),i=1,nsol2)
        write(50,520) (cii(2,i),i=1,nsol2)
        write(50,520) (ceii(1,i),i=1,nsol2)
        write(50,520) (ceii(2,i),i=1,nsol2)

        do i=1,nsol2
          ci(1,i)=ci(1,i)+cii(1,i)
          ci(2,i)=ci(2,i)+cii(2,i)
          cei(1,i)=cei(1,i)+ceii(1,i)
          cei(2,i)=cei(2,i)+ceii(2,i)
        enddo

      enddo
c- end: analyze configurations

c- begin: output results
      do i=1,nsol2
        do j=1,2
          f=ci(j,i)+1.d-20
          cei(j,i)=cei(j,i)/f
          ci(j,i)=ci(j,i)/vav
        enddo
      enddo
      vav=vav/dble(nconf)
      do i=1,nsol2
        do j=1,108
          do k=1,54
            f=cd(i,j,k)+1.d-20
            ced(i,j,k)=ced(i,j,k)/f
            f=cdc(j,k)+1.d-20
            cd(i,j,k)=cd(i,j,k)/(vav*f)
          enddo
        enddo
      enddo
      do i=1,108
        do j=1,54
          f=cdc(i,j)+1.d-20
          hd(i,j)=hd(i,j)/f
          vn1(i,j)=vn1(i,j)/f
          vn2(i,j)=vn2(i,j)/f
          vn2(i,j)=vn2(i,j)-vn1(i,j)*vn1(i,j)
        enddo
      enddo
      vn1b=vn1b/dble(nconf)
      vn2b=vn2b/dble(nconf)
      vn2b=vn2b-vn1b*vn1b

      rr=-dr
      do i=1,54
        rr=rr+dr
        dv=dble(nconf)*pi*((rr+dr)**2.d0-rr*rr)*dr
        do j=1,108
          gh(j,i)=gh(j,i)/(gg(j,i)+1.d-20)
          gh2(j,i)=gh2(j,i)/(gg(j,i)+1.d-20)
          gn(j,i)=gn(j,i)/(gg(j,i)+1.d-20)
          gn2(j,i)=gn2(j,i)/(gg(j,i)+1.d-20)
          gg(j,i)=gg(j,i)/dv
        enddo
      enddo
      s1=0.d0
      s2=0.d0
      do i=0,200
        s1=s1+po(i)
        s2=s2+ph(i)
      enddo
      do i=0,200
        po(i)=po(i)/s1
        ph(i)=ph(i)/s2
      enddo

      avin=avin/dble(nconf)
      a2vin=a2vin/dble(nconf)
      a2vin=dsqrt(a2vin-avin*avin)
      ahbin=ahbin/dble(nconf)
      vout=vav-vin
      anhh=anhh/anout
      anoo=anoo/anout
      anout=anout/dble(nconf)
      rho=anout/vout
      write(40,*) vav,'   average system volume'
c ======= VOL-API: report your volume instead of cylinder ===============
      vbind = vol_value()
c =======================================================================
      write(40,*) vbind,'   volume of binding zone'
      write(40,*) avin,'   average # of waters in binding zone'
      write(40,*) a2vin,'   deviation in # of waters in binding zone'
      write(40,*) anhh,'   bulk number of h-bonds'
      write(40,*) ahbin,'   # of hbonds in binding zone'
      write(40,*) ' '
      do i=1,nsol2
        write(40,510) i,ci(1,i),cei(1,i),ci(2,i),cei(2,i)
      enddo
      write(40,*) ' '
      do i=0,200
        write(40,510) i,po(i),ph(i)
      enddo

      if (rev) then
        do i=1,nsol2
          do j=1,108
            do k=1,54
              rav(j,k)=cd(i,j,k)
            enddo
          enddo
          do j=1,108
            do k=1,54
              cd(i,j,k)=(rav(j,k)+rav((109-j),k))/2.d0
            enddo
          enddo
          do j=1,108
            do k=1,54
              rav(j,k)=ced(i,j,k)
            enddo
          enddo
          do j=1,108
            do k=1,54
              ced(i,j,k)=(rav(j,k)+rav((109-j),k))/2.d0
            enddo
          enddo
        enddo
        do j=1,108
          do k=1,54
            rav(j,k)=gg(j,k)
          enddo
        enddo
        do j=1,108
          do k=1,54
            gg(j,k)=(rav(j,k)+rav((109-j),k))/2.d0
          enddo
        enddo
        do j=1,108
          do k=1,54
            rav(j,k)=gh(j,k)
          enddo
        enddo
        do j=1,108
          do k=1,54
            gh(j,k)=(rav(j,k)+rav((109-j),k))/2.d0
          enddo
        enddo
        do j=1,108
          do k=1,54
            rav(j,k)=gh2(j,k)
          enddo
        enddo
        do j=1,108
          do k=1,54
            gh2(j,k)=(rav(j,k)+rav((109-j),k))/2.d0
          enddo
        enddo
        do j=1,108
          do k=1,54
            rav(j,k)=gn(j,k)
          enddo
        enddo
        do j=1,108
          do k=1,54
            gn(j,k)=(rav(j,k)+rav((109-j),k))/2.d0
          enddo
        enddo
        do j=1,108
          do k=1,54
            rav(j,k)=gn2(j,k)
          enddo
        enddo
        do j=1,108
          do k=1,54
            gn2(j,k)=(rav(j,k)+rav((109-j),k))/2.d0
          enddo
        enddo
        do j=1,108
          do k=1,54
            rav(j,k)=vn1(j,k)
          enddo
        enddo
        do j=1,108
          do k=1,54
            vn1(j,k)=(rav(j,k)+rav((109-j),k))/2.d0
          enddo
        enddo
        do j=1,108
          do k=1,54
            rav(j,k)=vn2(j,k)
          enddo
        enddo
        do j=1,108
          do k=1,54
            vn2(j,k)=(rav(j,k)+rav((109-j),k))/2.d0
          enddo
        enddo
        do j=1,108
          do k=1,54
            rav(j,k)=hd(j,k)
          enddo
        enddo
        do j=1,108
          do k=1,54
            hd(j,k)=(rav(j,k)+rav((109-j),k))/2.d0
          enddo
        enddo
      endif

      do i=1,nsol2
        write(41,510) i
        write(41,*) ' '
        write(42,510) i
        write(42,*) ' '
        write(41,530) (rtg(j),j=0,54)
        write(42,530) (rtg(j),j=0,54)
        zz=13.5d0+dr/2.d0
        do j=1,108
          zz=zz-dr
          write(41,540) zz,(cd(i,j,k),k=1,54)
          write(42,540) zz,(ced(i,j,k),k=1,54)
        enddo
        write(41,*) ' '
        write(41,*) 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
        write(41,*) ' '
        write(42,*) ' '
        write(42,*) 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
        write(42,*) ' '
      enddo

      write(43,*) rho,'   bulk density'
      write(43,*) anhh,'   bulk number of h-bonds'
      write(43,*) anoo,'   bulk coordination number'
      write(43,*) ' '
      write(43,*) 'spatial disribtution function'
      write(43,*) ' '
      write(43,530) (rtg(j),j=0,54)
      zz=13.5d0+dr/2.d0
      do i=1,108
        zz=zz-dr
        write(43,540) zz,((gg(i,j)/rho),j=1,54)
      enddo
      write(43,*) ' '
      write(42,*) 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
      write(43,*) ' '
      write(43,530) (rtg(j),j=0,54)
      zz=13.5d0+dr/2.d0
      do i=1,108
        zz=zz-dr
        write(43,540) zz,(gh(i,j),j=1,54)
      enddo
      write(43,*) ' '
      write(42,*) 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
      write(43,*) ' '
      write(43,530) (rtg(j),j=0,54)
      zz=13.5d0+dr/2.d0
      do i=1,108
        zz=zz-dr
        write(43,540) zz,(gh2(i,j),j=1,54)
      enddo
      write(43,*) ' '
      write(42,*) 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
      write(43,*) ' '
      write(43,530) (rtg(j),j=0,54)
      zz=13.5d0+dr/2.d0
      do i=1,108
        zz=zz-dr
        write(43,540) zz,(gn(i,j),j=1,54)
      enddo
      write(43,*) ' '
      write(42,*) 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
      write(43,*) ' '
      write(43,530) (rtg(j),j=0,54)
      zz=13.5d0+dr/2.d0
      do i=1,108
        zz=zz-dr
        write(43,540) zz,(gn2(i,j),j=1,54)
      enddo
      write(43,*) ' '
      write(42,*) 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
      write(43,*) ' '
      write(43,530) (rtg(j),j=0,54)
      zz=13.5d0+dr/2.d0
      do i=1,108
        zz=zz-dr
        write(43,540) zz,((cdc(i,j)/dble(nconf)),j=1,54)
      enddo

      write(44,*) vn1b,'   bulk # of wat'
      write(44,*) vn2b,'   bulk variance # of wat'
      write(44,*) ' '
      write(44,530) (rtg(j),j=0,54)
      zz=13.5d0+dr/2.d0
      do i=1,108
        zz=zz-dr
        write(44,540) zz,(vn1(i,j),j=1,54)
      enddo
      write(44,*) ' '
      write(44,*) 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
      write(44,*) ' '
      write(44,530) (rtg(j),j=0,54)
      zz=13.5d0+dr/2.d0
      do i=1,108
        zz=zz-dr
        write(44,540) zz,(vn2(i,j),j=1,54)
      enddo
      write(44,*) ' '
      write(44,*) 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
      write(44,*) ' '
      write(44,530) (rtg(j),j=0,54)
      zz=13.5d0+dr/2.d0
      do i=1,108
        zz=zz-dr
        write(44,540) zz,(hd(i,j),j=1,54)
      enddo
c- end: output results

500   format(a1)
510   format(i8,4f25.10)
520   format(40e20.12)
525   format(10f20.10)
530   format(f10.3,54f25.3)
540   format(f10.3,54f25.10)
550   format(f15.10)
900   format('The temperature is ',f7.2,' K')

      stop
      end


      double precision function ranu(iqfd)
      implicit real*8 (a-h,o-z)
      common /ran0/ irn
      data icnt /0/,ichg /1167/,ichg0 /1167/
      data imul /1173/,icon /458753759/
      data imod /1048573/,icn0 /458753759/
      data jmul /1161/, jcon /458716759/
      data jmod /1048573/,jrn /124690/
      icnt = icnt + 1
      if (icnt.eq.ichg) then
        jrn = jrn*jmul+jcon
        jrn = mod(jrn,jmod)
        rnj = float(jrn)/float(jmod)
        if (rnj.gt.0.5) then
          fac = 1.+0.5*rnj
        else
          fac = 1.-0.5*rnj
        endif
        fac = float(icn0)*fac
        icon = int(fac)
        jrn = jrn*jmul+jcon
        jrn = mod(jrn,jmod)
        rnj = float(jrn)/float(jmod)
        if (rnj.gt.0.5) then
          fac = 1.+0.5*rnj
        else
          fac = 1.-0.5*rnj
        end if
        fac = float(ichg0)*fac
        ichg = int(fac)
        icnt = 0
      endif
      irn = irn*imul+icon
      irn = mod(irn,imod)
      ranu = float(irn)/float(imod)
      return
      end
