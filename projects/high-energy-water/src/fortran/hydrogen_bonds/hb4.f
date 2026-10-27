      program sineater_oa
      implicit real*8 (a-h,o-z)
      character*1 line

      dimension x(10000), y(10000), z(10000)
      dimension mactop(100), macbot(100)

      open(unit=1, file='drive', status='old')
      open(unit=2, file='top', status='old')
      open(unit=10, file='dist.out', status='unknown')

      read(*,500) line

      pi = 2.d0 * dasin(1.d0)

      read(1,*) nconf
      read(1,*) nout
      read(1,*) roo
      read(1,*) ang

      read(2,*) nwat
      read(2,*) nhst
      read(2,*) nactop
      do i = 1, nactop
         read(2,*) mactop(i)
      enddo
      read(2,*) nacbot
      do i = 1, nacbot
         read(2,*) macbot(i)
      enddo

      nac = nactop + nacbot

      ntot = 4 * nwat + nhst + 8  ! +8 for sodium ions

      write(10,*) nwat, nhst, ntot

      roo = roo * roo
      cang = dcos(pi * ang / 180.d0)

      hbw = 0.d0
      hbt = 0.d0
      hbb = 0.d0
      hbw2 = 0.d0
      hbt2 = 0.d0
      hbb2 = 0.d0

      do nc = 1, nconf
         read(*,*) (x(i), y(i), z(i), i=1, ntot), bx, by, bz

         do n = 1, nhst
            x(n) = x(n) - bx * dnint(x(n)/bx)
            y(n) = y(n) - by * dnint(y(n)/by)
            z(n) = z(n) - bz * dnint(z(n)/bz)
         enddo

         do n = 1, nwat
            i = 4 * n - 3 + nhst
            if (i + 3 > ntot - 8) cycle
            xx = x(i)
            yy = y(i)
            zz = z(i)
            x(i) = x(i) - bx * dnint(x(i)/bx)
            y(i) = y(i) - by * dnint(y(i)/by)
            z(i) = z(i) - bz * dnint(z(i)/bz)
            do j = 1, 3
               dx = x(i + j) - xx
               dy = y(i + j) - yy
               dz = z(i + j) - zz
               dx = dx - bx * dnint(dx/bx)
               dy = dy - by * dnint(dy/by)
               dz = dz - bz * dnint(dz/bz)
               x(i + j) = x(i) + dx
               y(i + j) = y(i) + dy
               z(i + j) = z(i) + dz
            enddo
         enddo

         nwd = 0
         ntopd = 0
         nbotd = 0

         ! --- Top acceptors loop
         do ii = 1, nwat
            i = 4 * ii - 3 + nhst
            if (i + 3 > ntot - 8) cycle
            do jj = 1, nactop
               j = mactop(jj)
               dx = x(j) - x(i)
               dy = y(j) - y(i)
               dz = z(j) - z(i)
               dx = dx - bx * dnint(dx/bx)
               dy = dy - by * dnint(dy/by)
               dz = dz - bz * dnint(dz/bz)
               r2 = dx*dx + dy*dy + dz*dz
               if (r2 .gt. roo) goto 111
               r2 = dsqrt(r2)
               dx = dx / r2
               dy = dy / r2
               dz = dz / r2
               do k = 1, 2
                  dxh = x(i + k) - x(i)
                  dyh = y(i + k) - y(i)
                  dzh = z(i + k) - z(i)
                  r3 = dsqrt(dxh*dxh + dyh*dyh + dzh*dzh)
                  dxh = dxh / r3
                  dyh = dyh / r3
                  dzh = dzh / r3
                  dot = dx * dxh + dy * dyh + dz * dzh
                  if (dot .gt. cang) ntopd = ntopd + 1
               enddo
111            continue
            enddo
         enddo

         ! --- Bottom acceptors loop
         do ii = 1, nwat
            i = 4 * ii - 3 + nhst
            if (i + 3 > ntot - 8) cycle
            do jj = 1, nacbot
               j = macbot(jj)
               dx = x(j) - x(i)
               dy = y(j) - y(i)
               dz = z(j) - z(i)
               dx = dx - bx * dnint(dx/bx)
               dy = dy - by * dnint(dy/by)
               dz = dz - bz * dnint(dz/bz)
               r2 = dx*dx + dy*dy + dz*dz
               if (r2 .gt. roo) goto 112
               r2 = dsqrt(r2)
               dx = dx / r2
               dy = dy / r2
               dz = dz / r2
               do k = 1, 2
                  dxh = x(i + k) - x(i)
                  dyh = y(i + k) - y(i)
                  dzh = z(i + k) - z(i)
                  r3 = dsqrt(dxh*dxh + dyh*dyh + dzh*dzh)
                  dxh = dxh / r3
                  dyh = dyh / r3
                  dzh = dzh / r3
                  dot = dx * dxh + dy * dyh + dz * dzh
                  if (dot .gt. cang) nbotd = nbotd + 1
               enddo
112            continue
            enddo
         enddo

         ! --- Water-water loop
         do ii = 1, nwat
            i = 4 * ii - 3 + nhst
            if (i + 3 > ntot - 8) cycle
            do jj = 1, nwat
               j = 4 * jj - 3 + nhst
               if (i == j .or. j + 3 > ntot - 8) goto 10
               dx = x(j) - x(i)
               dy = y(j) - y(i)
               dz = z(j) - z(i)
               dx = dx - bx * dnint(dx/bx)
               dy = dy - by * dnint(dy/by)
               dz = dz - bz * dnint(dz/bz)
               r2 = dx*dx + dy*dy + dz*dz
               if (r2 .gt. roo) goto 10
               r2 = dsqrt(r2)
               dx = dx / r2
               dy = dy / r2
               dz = dz / r2
               do k = 1, 2
                  dxh = x(i + k) - x(i)
                  dyh = y(i + k) - y(i)
                  dzh = z(i + k) - z(i)
                  r3 = dsqrt(dxh*dxh + dyh*dyh + dzh*dzh)
                  dxh = dxh / r3
                  dyh = dyh / r3
                  dzh = dzh / r3
                  dot = dx * dxh + dy * dyh + dz * dzh
                  if (dot .gt. cang) nwd = nwd + 1
               enddo
10             continue
            enddo
         enddo

         ! --- Accumulate sums
         hbw = hbw + dble(nwd)
         hbt = hbt + dble(ntopd)
         hbb = hbb + dble(nbotd)
         hbw2 = hbw2 + dble(nwd)
         hbt2 = hbt2 + dble(ntopd)
         hbb2 = hbb2 + dble(nbotd)

         write(*,*) nc, nwd, ntopd, nbotd

         if (mod(nc, nout) .eq. 0) then
            hbw2 = hbw2 / dble(nout)
            hbt2 = hbt2 / dble(nout)
            hbb2 = hbb2 / dble(nout)
            write(10,*) nc, hbw2, hbt2, hbb2
            write(10,*) '------------------'
            hbw2 = 0.d0
            hbt2 = 0.d0
            hbb2 = 0.d0
         endif

      enddo

      hbw = hbw / dble(nconf)
      hbt = hbt / dble(nconf)
      hbb = hbb / dble(nconf)

      write(10,*) nconf, hbw, hbt, hbb

500   format(a1)

      stop
      end
