      program ringstats
c --- Fortran 95; compile with:  f95 ring_angle.f -o ring_angle
c --- Reads:
c ---   drive : ntot, nframes, [nskip], [stride]
c ---   vol   : NMON and RING indices (6 heavy ring atoms per monomer)
c ---   stdin : prd.crdbox (TITLE line, then frames of x y z ... bx by bz)
c --- Writes:
c ---   cos_theta_hist.out         (center, P, count, first_frame)
c ---   cos_theta_summary.out      (mean, std, count)
c ---   d12_hist.out .. d15_hist.out (center, P, count, first_frame)
c ---   d_summary.out              (per-k mean/std/count)

      implicit real*8 (a-h,o-z)

      integer MAXATOM, MAXRING, MAXMON, KMAX, NBCTH, NBD
      parameter (MAXATOM=100000, MAXRING=2048, MAXMON=64)
      parameter (KMAX=5)
      parameter (NBCTH=100)          ! cos(theta) bins on [-1,1]
      parameter (NBD=600)            ! distances 0..30 A in 0.05 A bins

c --- Drive params
      integer ntot, nframes, nskip, stride, ios

c --- Vol (ring) indexing
      integer nmon, chunk, nring
      integer iring(MAXRING)

c --- Coords & box
      real*8 x(MAXATOM), y(MAXATOM), z(MAXATOM)
      real*8 bx, by, bz

c --- Unwrapped ring coords (all monomers)
      real*8 xr(MAXRING), yr(MAXRING), zr(MAXRING)

c --- Per-monomer COM & normal, macro COM
      real*8 comx(MAXMON), comy(MAXMON), comz(MAXMON)
      real*8 nxm(MAXMON), nym(MAXMON), nzm(MAXMON)
      real*8 comx_mac, comy_mac, comz_mac

c --- Histograms & stats
      real*8 hcth(NBCTH), cth_min, cth_max, cth_w
      integer first_cth(NBCTH)
      integer hD(KMAX, NBD)
      integer first_d(KMAX, NBD)
      real*8 d_min, d_max, d_w
      real*8 sum_cth, sum2_cth
      integer ns_cth
      real*8 sumd(KMAX), sumd2(KMAX)
      integer nsd(KMAX)

c --- Scratch
      integer i, j, k, f, m, ii, jj, kk
      character*256 line

c =======================
c  Read drive parameters
c =======================
      nskip  = 0
      stride = 1
      open(99,file='drive',status='old',iostat=ios)
      if (ios.ne.0) then
        write(*,*) 'ERROR: cannot open drive'
        stop
      endif
      read(99,*,iostat=ios) ntot
      if (ios.ne.0) stop 'drive: missing ntot'
      read(99,*,iostat=ios) nframes
      if (ios.ne.0) stop 'drive: missing nframes'
      read(99,*,iostat=ios) nskip
      if (ios.ne.0) nskip = 0
      read(99,*,iostat=ios) stride
      if (ios.ne.0) stride = 1
      close(99)

c =================
c  Read VOL (RING)
c =================
      call vol_read_ring('vol', nmon, chunk, iring, nring)

      if (nring.gt.MAXRING) stop 'RING too large for MAXRING'
      if (nmon.gt.MAXMON)   stop 'NMON too large for MAXMON'
      if (chunk.ne.6) then
        write(*,*) 'Note: CHUNK_RING=',chunk,' (expected 6)'
      endif
      if (nring.ne.nmon*chunk) then
        write(*,*) 'ERROR: RING count=',nring,
     &             ' != NMON*CHUNK=',nmon*chunk
        stop
      endif

c ==========================
c  Init histograms & stats
c ==========================
      do i=1,NBCTH
        hcth(i)=0.d0
        first_cth(i)=0
      enddo
      do k=1,KMAX
        do i=1,NBD
          hD(k,i)=0
          first_d(k,i)=0
        enddo
        sumd(k)=0.d0
        sumd2(k)=0.d0
        nsd(k)=0
      enddo
      sum_cth=0.d0
      sum2_cth=0.d0
      ns_cth=0

      cth_min = -1.d0
      cth_max =  1.d0
      cth_w   = (cth_max-cth_min)/dble(NBCTH)

      d_min = 0.d0
      d_max = 30.d0
      d_w   = (d_max-d_min)/dble(NBD)

c ==================
c  Read TITLE line
c ==================
      read(*,'(A)',iostat=ios) line
      if (ios.ne.0) then
        write(*,*) 'ERROR: no TITLE line / bad input'
        stop
      endif

c ===========================
c  Skip initial frames (nskip)
c ===========================
      if (nskip.gt.0) then
        do f=1,nskip
          read(*,*,iostat=ios) (x(i),y(i),z(i),i=1,ntot), bx,by,bz
          if (ios.ne.0) then
            write(*,*) 'WARNING: EOF while skipping frames at f=',f
            goto 900
          endif
        enddo
      endif

c =====================
c  Main frame loop
c =====================
      kk = 0
      f  = 0
 100  continue
      read(*,*,iostat=ios) (x(i),y(i),z(i),i=1,ntot), bx,by,bz
      if (ios.ne.0) goto 900
      f = f + 1

c ---- stride filter (retain absolute frame number f for first_seen)
      if (mod(f-1, stride).ne.0) goto 100

      kk = kk + 1
      if (kk.gt.nframes) goto 900

c ---- unwrap all RING atoms relative to first ring atom
      call unwrap_ring(ntot, x,y,z, bx,by,bz, iring, nring,
     &                 xr,yr,zr)

c ---- macro COM over all ring atoms
      call com_block(xr,yr,zr, 1, nring, comx_mac, comy_mac, comz_mac)

c ---- per-monomer COM and normals (Newell); SIGNED cosθ (no reorient)
      do m=1,nmon
        ii = (m-1)*chunk + 1
        jj = ii + chunk - 1
        call com_block(xr,yr,zr, ii, jj, comx(m), comy(m), comz(m))
        call newell_normal(xr(ii),yr(ii),zr(ii), chunk,
     &                     nxm(m),nym(m),nzm(m))

        call add_costheta_with_first(nxm(m),nym(m),nzm(m),
     &                    comx(m)-comx_mac,
     &                    comy(m)-comy_mac,
     &                    comz(m)-comz_mac,
     &                    cth_min, cth_w, NBCTH,
     &                    hcth, first_cth, f,
     &                    sum_cth, sum2_cth, ns_cth)
      enddo

c ---- distances d12..d15 (k=1..5) between monomer COMs
      do k=1,KMAX
        if (k.ge.nmon) cycle
        do m=1,nmon
          j = m + k
          if (j.gt.nmon) j = j - nmon
          call add_distance_with_first(comx(m),comy(m),comz(m),
     &                      comx(j),comy(j),comz(j),
     &                      d_min, d_w, NBD, KMAX, k,
     &                      hD, first_d, f,
     &                      sumd, sumd2, nsd)
        enddo
      enddo

      goto 100

c =====================
c  Write output files
c =====================
 900  continue
      call write_cos_hist_ff('cos_theta_hist.out',
     &                    cth_min, cth_w, NBCTH, hcth, first_cth)
      call write_cos_summary('cos_theta_summary.out',
     &                       sum_cth, sum2_cth, ns_cth)

      call write_d_hist_ff('d12_hist.out', 1, d_min, d_w,
     &                     NBD, KMAX, hD, first_d)
      call write_d_hist_ff('d13_hist.out', 2, d_min, d_w,
     &                     NBD, KMAX, hD, first_d)
      call write_d_hist_ff('d14_hist.out', 3, d_min, d_w,
     &                     NBD, KMAX, hD, first_d)
      call write_d_hist_ff('d15_hist.out', 4, d_min, d_w,
     &                     NBD, KMAX, hD, first_d)

      call write_d_summary('d_summary.out', KMAX, sumd, sumd2, nsd)

      stop
      end

c ============================
c  ---- subroutines below ----
c ============================

      subroutine vol_read_ring(fname, nmon, chunk, iring, nring)
      implicit real*8 (a-h,o-z)
      character*(*) fname
      integer nmon, chunk, iring(*), nring
      character*256 ln, up
      integer iu, ios, pos, state, val
      integer i, nlen

      nmon = -1
      chunk = 6
      nring = 0
      state = 0

      open(newunit=iu, file=fname, status='old', iostat=ios)
      if (ios.ne.0) then
        write(*,*) 'ERROR: cannot open ',fname
        stop
      endif

 10   continue
      read(iu,'(A)',iostat=ios) ln
      if (ios.ne.0) goto 99

c --- strip comments after '#'
      nlen = len_trim(ln)
      do i=1,nlen
        if (ln(i:i).eq.'#') then
          nlen = i-1
          exit
        endif
      enddo
      if (nlen.le.0) goto 10
      up = ln
      call toupper(up)

c --- NMON
      pos = index(up(1:nlen),'NMON')
      if (pos.gt.0) then
        read(up(pos+4:nlen),*,iostat=ios) nmon
        goto 10
      endif

c --- CHUNK_RING (optional)
      pos = index(up(1:nlen),'CHUNK_RING')
      if (pos.gt.0) then
        read(up(pos+10:nlen),*,iostat=ios) chunk
        goto 10
      endif

c --- RING block
      if (state.eq.0) then
        pos = index(up(1:nlen),'RING')
        if (pos.gt.0) then
          state = 1
        endif
        goto 10
      else
c ----- reading ring indices; expect one integer per line
        read(up(1:nlen),*,iostat=ios) val
        if (ios.eq.0) then
          nring = nring + 1
          iring(nring) = val
        else
c ------- non-integer => end of RING block
          state = 0
        endif
        goto 10
      endif

 99   continue
      close(iu)

      if (nmon.le.0) then
        write(*,*) 'ERROR: NMON not found in vol'
        stop
      endif

      return
      end

      subroutine toupper(s)
      implicit real*8 (a-h,o-z)
      character*(*) s
      integer i, c
      do i=1,len_trim(s)
        c = ichar(s(i:i))
        if (c.ge.97 .and. c.le.122) s(i:i)=achar(c-32)
      enddo
      end

      subroutine unwrap_ring(ntot,x,y,z,bx,by,bz,iring,nring,
     &                       xr,yr,zr)
      implicit real*8 (a-h,o-z)
      integer ntot, nring, iring(*)
      real*8 x(*),y(*),z(*), bx,by,bz, xr(*),yr(*),zr(*)
      real*8 xref,yref,zref, dx,dy,dz
      integer t, ia

c --- anchor = first ring atom
      ia = iring(1)
      xref = x(ia)
      yref = y(ia)
      zref = z(ia)

      do t=1,nring
        ia = iring(t)
        dx = x(ia)-xref
        dy = y(ia)-yref
        dz = z(ia)-zref
        dx = dx - bx*dnint(dx/bx)
        dy = dy - by*dnint(dy/by)
        dz = dz - bz*dnint(dz/bz)
        xr(t) = xref + dx
        yr(t) = yref + dy
        zr(t) = zref + dz
      enddo
      return
      end

      subroutine com_block(x,y,z, i1,i2, cx,cy,cz)
      implicit real*8 (a-h,o-z)
      real*8 x(*),y(*),z(*), cx,cy,cz
      integer i1,i2, i
      cx=0.d0
      cy=0.d0
      cz=0.d0
      do i=i1,i2
        cx=cx+x(i)
        cy=cy+y(i)
        cz=cz+z(i)
      enddo
      cx=cx/dble(i2-i1+1)
      cy=cy/dble(i2-i1+1)
      cz=cz/dble(i2-i1+1)
      return
      end

      subroutine newell_normal(x,y,z, n, nx,ny,nz)
      implicit real*8 (a-h,o-z)
      real*8 x(*),y(*),z(*), nx,ny,nz
      integer n, i, j
      real*8 ax,ay,az, norm
      ax=0.d0
      ay=0.d0
      az=0.d0
      do i=1,n
        j = i+1
        if (j.gt.n) j=1
        ax = ax + (y(i)-y(j))*(z(i)+z(j))
        ay = ay + (z(i)-z(j))*(x(i)+x(j))
        az = az + (x(i)-x(j))*(y(i)+y(j))
      enddo
      norm = dsqrt(ax*ax+ay*ay+az*az)
      if (norm.gt.1.d-12) then
        nx = ax/norm
        ny = ay/norm
        nz = az/norm
      else
        nx = 0.d0
        ny = 0.d0
        nz = 1.d0
      endif
      return
      end

      subroutine add_costheta_with_first(nx,ny,nz, rx,ry,rz,
     &                        cmin, cbinw, NB, hist, first_seen, fr,
     &                        s1, s2, ns)
      implicit real*8 (a-h,o-z)
      integer NB, ib, fr
      real*8 nx,ny,nz, rx,ry,rz, r2, r, cth
      real*8 cmin, cbinw, hist(*), s1, s2
      integer ns, first_seen(*)
      r2 = rx*rx + ry*ry + rz*rz
      if (r2.le.1.d-12) return
      r = dsqrt(r2)
      cth = (nx*rx + ny*ry + nz*rz)/r
      if (cth.lt.-1.d0) cth=-1.d0
      if (cth.gt. 1.d0) cth= 1.d0
      ib = 1 + int( (cth - cmin)/cbinw )
      if (ib.lt.1) ib=1
      if (ib.gt.NB) ib=NB
      hist(ib) = hist(ib) + 1.d0
      if (first_seen(ib).eq.0) first_seen(ib)=fr
      s1 = s1 + cth
      s2 = s2 + cth*cth
      ns = ns + 1
      return
      end

      subroutine add_distance_with_first(x1,y1,z1, x2,y2,z2,
     &                        dmin, dbinw, NB, KX, kidx,
     &                        HD, firstD, fr, s1, s2, ns)
      implicit real*8 (a-h,o-z)
c     NB  = number of distance bins
c     KX  = max k-separations tracked
      integer NB, KX, kidx, ib, fr
      real*8 x1,y1,z1, x2,y2,z2, dx,dy,dz, d, dmin, dbinw
      integer HD(KX, NB), firstD(KX, NB)
      real*8  s1(KX), s2(KX)
      integer ns(KX)

      dx = x2-x1
      dy = y2-y1
      dz = z2-z1
      d  = dsqrt(dx*dx + dy*dy + dz*dz)

      ib = 1 + int( (d - dmin)/dbinw )
      if (ib.lt.1) ib=1
      if (ib.gt.NB) ib=NB

      HD(kidx, ib) = HD(kidx, ib) + 1
      if (firstD(kidx,ib).eq.0) firstD(kidx,ib)=fr
      s1(kidx)     = s1(kidx) + d
      s2(kidx)     = s2(kidx) + d*d
      ns(kidx)     = ns(kidx) + 1
      return
      end

      subroutine write_cos_hist_ff(fname, cmin, cbinw, NB, hist, first)
      implicit real*8 (a-h,o-z)
      character*(*) fname
      integer NB, iu, i
      real*8 cmin, cbinw, hist(*), tot, center
      integer first(*)
      open(newunit=iu, file=fname, status='unknown')
      tot=0.d0
      do i=1,NB
        tot = tot + hist(i)
      enddo
      if (tot.le.0.d0) tot=1.d0
      do i=1,NB
        center = cmin + (dble(i)-0.5d0)*cbinw
        write(iu,'(3f16.6,1x,i10)')
     &    center, hist(i)/tot, hist(i), first(i)
      enddo
      close(iu)
      return
      end

      subroutine write_cos_summary(fname, s1,s2,ns)
      implicit real*8 (a-h,o-z)
      character*(*) fname
      real*8 s1,s2, mu, sig2
      integer ns, iu
      open(newunit=iu, file=fname, status='unknown')
      if (ns.gt.0) then
        mu = s1/dble(ns)
        sig2 = s2/dble(ns) - mu*mu
        if (sig2.lt.0.d0) sig2=0.d0
        write(iu,'(A,1x,f12.6)') 'mean', mu
        write(iu,'(A,1x,f12.6)') 'std', dsqrt(sig2)
        write(iu,'(A,1x,i12)')   'count', ns
      else
        write(iu,*) 'no samples'
      endif
      close(iu)
      return
      end

      subroutine write_d_hist_ff(fname, kidx, dmin, dbinw,
     &                           NB, KX, HD, firstD)
      implicit real*8 (a-h,o-z)
      character*(*) fname
      integer kidx, NB, KX, iu, i
      integer HD(KX, NB), firstD(KX, NB)
      real*8 dmin, dbinw, tot, center
      open(newunit=iu, file=fname, status='unknown')
      tot=0.d0
      do i=1,NB
        tot = tot + dble(HD(kidx,i))
      enddo
      if (tot.le.0.d0) tot=1.d0
      do i=1,NB
        center = dmin + (dble(i)-0.5d0)*dbinw
        write(iu,'(f12.6,1x,f12.6,1x,i12,1x,i10)')
     &    center, dble(HD(kidx,i))/tot, HD(kidx,i), firstD(kidx,i)
      enddo
      close(iu)
      return
      end

      subroutine write_d_summary(fname, KX, s1, s2, ns)
      implicit real*8 (a-h,o-z)
      character*(*) fname
      real*8 s1(*), s2(*), mu, sig2
      integer ns(*), iu, k, KX
      open(newunit=iu, file=fname, status='unknown')
      write(iu,'(A)') '# k  mean(A)   std(A)   count'
      do k=1,KX
        if (ns(k).gt.0) then
          mu = s1(k)/dble(ns(k))
          sig2 = s2(k)/dble(ns(k)) - mu*mu
          if (sig2.lt.0.d0) sig2=0.d0
          write(iu,'(i3,1x,f10.4,1x,f10.4,1x,i12)')
     &      k, mu, dsqrt(sig2), ns(k)
        else
          write(iu,'(i3,1x,A)') k, 'no samples'
        endif
      enddo
      close(iu)
      return
      end
