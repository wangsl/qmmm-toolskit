 
C     $Id: etortor2.F 2 2008-02-28 22:52:52Z wangsl $
c
c
c     #############################################################
c     ##  COPYRIGHT (C) 2003 by Pengyu Ren & Jay William Ponder  ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine etortor2  --  atomwise torsion-torsion Hessian  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "etortor2" calculates the torsion-torsion potential energy
c     second derivatives with respect to Cartesian coordinates
c
c
      subroutine etortor2 (i)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bitor.i'
      include 'bound.i'
      include 'group.i'
      include 'hessn.i'
      include 'ktrtor.i'
      include 'math.i'
      include 'torpot.i'
      include 'tortor.i'
      include 'units.i'
      integer i,j,k,itortor
      integer pos1,pos2
      integer ia,ib,ic,id,ie
      integer nlo,nhi,nt,xlo,ylo
      real*8 e,fgrp,sign
      real*8 angle1,angle2
      real*8 value1,value2
      real*8 cosine1,cosine2
      real*8 dedphi,dedpsi,d12
      real*8 d2edphi2,d2edpsi2
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xid,yid,zid
      real*8 xie,yie,zie
      real*8 xba,yba,zba
      real*8 xcb,ycb,zcb
      real*8 xdc,ydc,zdc
      real*8 xca,yca,zca
      real*8 xdb,ydb,zdb
      real*8 xed,yed,zed
      real*8 xec,yec,zec
      real*8 xt,yt,zt
      real*8 xu,yu,zu
      real*8 xv,yv,zv
      real*8 xtu,ytu,ztu
      real*8 xuv,yuv,zuv
      real*8 xh,yh,x1l,x1u,y1l,y1u
      real*8 rt2,ru2,rv2
      real*8 rtru,rcb,rurv,rdc
      real*8 dphidxt,dphidyt,dphidzt
      real*8 dphidxu,dphidyu,dphidzu
      real*8 dphidxia,dphidyia,dphidzia
      real*8 dphidxib,dphidyib,dphidzib
      real*8 dphidxic,dphidyic,dphidzic
      real*8 dphidxid,dphidyid,dphidzid
      real*8 dpsidxv,dpsidyv,dpsidzv
      real*8 dpsidxu,dpsidyu,dpsidzu
      real*8 dpsidxie,dpsidyie,dpsidzie
      real*8 dpsidxib,dpsidyib,dpsidzib
      real*8 dpsidxic,dpsidyic,dpsidzic
      real*8 dpsidxid,dpsidyid,dpsidzid
      real*8 xycb2,xzcb2,yzcb2
      real*8 xydc2,xzdc2,yzdc2
      real*8 rcbxt,rcbyt,rcbzt,rcbt2
      real*8 rcbxu,rcbyu,rcbzu,rcbu2
      real*8 rdcxv,rdcyv,rdczv,rdcv2
      real*8 rdcxu,rdcyu,rdczu,rdcu2
      real*8 dphidxibt,dphidyibt,dphidzibt
      real*8 dphidxibu,dphidyibu,dphidzibu
      real*8 dphidxict,dphidyict,dphidzict
      real*8 dphidxicu,dphidyicu,dphidzicu
      real*8 dpsidxidu,dpsidyidu,dpsidzidu
      real*8 dpsidxidv,dpsidyidv,dpsidzidv
      real*8 dpsidxicu,dpsidyicu,dpsidzicu
      real*8 dpsidxicv,dpsidyicv,dpsidzicv
      real*8 dxiaxia,dyiayia,dziazia
      real*8 dxibxib,dyibyib,dzibzib
      real*8 dxicxic,dyicyic,dziczic
      real*8 dxidxid,dyidyid,dzidzid
      real*8 dxiayia,dxiazia,dyiazia
      real*8 dxibyib,dxibzib,dyibzib
      real*8 dxicyic,dxiczic,dyiczic
      real*8 dxidyid,dxidzid,dyidzid
      real*8 dxiaxib,dxiayib,dxiazib
      real*8 dyiaxib,dyiayib,dyiazib
      real*8 dziaxib,dziayib,dziazib
      real*8 dxiaxic,dxiayic,dxiazic
      real*8 dyiaxic,dyiayic,dyiazic
      real*8 dziaxic,dziayic,dziazic
      real*8 dxiaxid,dxiayid,dxiazid
      real*8 dyiaxid,dyiayid,dyiazid
      real*8 dziaxid,dziayid,dziazid
      real*8 dxibxic,dxibyic,dxibzic
      real*8 dyibxic,dyibyic,dyibzic
      real*8 dzibxic,dzibyic,dzibzic
      real*8 dxibxid,dxibyid,dxibzid
      real*8 dyibxid,dyibyid,dyibzid
      real*8 dzibxid,dzibyid,dzibzid
      real*8 dxicxid,dxicyid,dxiczid
      real*8 dyicxid,dyicyid,dyiczid
      real*8 dzicxid,dzicyid,dziczid
      real*8 dxibxib2,dyibyib2,dzibzib2
      real*8 dxicxic2,dyicyic2,dziczic2
      real*8 dxidxid2,dyidyid2,dzidzid2
      real*8 dxiexie2,dyieyie2,dziezie2
      real*8 dxibyib2,dxibzib2,dyibzib2
      real*8 dxicyic2,dxiczic2,dyiczic2
      real*8 dxidyid2,dxidzid2,dyidzid2
      real*8 dxieyie2,dxiezie2,dyiezie2
      real*8 dxibxic2,dxibyic2,dxibzic2
      real*8 dyibxic2,dyibyic2,dyibzic2
      real*8 dzibxic2,dzibyic2,dzibzic2
      real*8 dxibxid2,dxibyid2,dxibzid2
      real*8 dyibxid2,dyibyid2,dyibzid2
      real*8 dzibxid2,dzibyid2,dzibzid2
      real*8 dxibxie2,dxibyie2,dxibzie2
      real*8 dyibxie2,dyibyie2,dyibzie2
      real*8 dzibxie2,dzibyie2,dzibzie2
      real*8 dxicxid2,dxicyid2,dxiczid2
      real*8 dyicxid2,dyicyid2,dyiczid2
      real*8 dzicxid2,dzicyid2,dziczid2
      real*8 dxicxie2,dxicyie2,dxiczie2
      real*8 dyicxie2,dyicyie2,dyiczie2
      real*8 dzicxie2,dzicyie2,dziczie2
      real*8 dxidxie2,dxidyie2,dxidzie2
      real*8 dyidxie2,dyidyie2,dyidzie2
      real*8 dzidxie2,dzidyie2,dzidzie2
      real*8 ftt(4),ft12(4)
      real*8 ft1(4),ft2(4)
      logical proceed
c
c
c     compute the Hessian elements of the torsion-torsions
c
      do itortor = 1, ntortor
         j = itt(1,itortor)
         k = itt(2,itortor)
         if (itt(3,itortor) .eq. 1) then
            ia = ibitor(1,j)
            ib = ibitor(2,j)
            ic = ibitor(3,j)
            id = ibitor(4,j)
            ie = ibitor(5,j)
         else
            ia = ibitor(5,j)
            ib = ibitor(4,j)
            ic = ibitor(3,j)
            id = ibitor(2,j)
            ie = ibitor(1,j)
         end if
c
c     decide whether to compute the current interaction
c
         proceed = (i.eq.ia .or. i.eq.ib .or. i.eq.ic
     &                 .or. i.eq.id .or. i.eq.ie)
         if (proceed .and. use_group)
     &      call groups (proceed,fgrp,ia,ib,ic,id,ie,0)
c
c     compute the values of the torsional angles
c
         if (proceed) then
            xia = x(ia)
            yia = y(ia)
            zia = z(ia)
            xib = x(ib)
            yib = y(ib)
            zib = z(ib)
            xic = x(ic)
            yic = y(ic)
            zic = z(ic)
            xid = x(id)
            yid = y(id)
            zid = z(id)
            xie = x(ie)
            yie = y(ie)
            zie = z(ie)
            xba = xib - xia
            yba = yib - yia
            zba = zib - zia
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            xdc = xid - xic
            ydc = yid - yic
            zdc = zid - zic
            xed = xie - xid
            yed = yie - yid
            zed = zie - zid
            if (use_polymer) then
               call image (xba,yba,zba,0)
               call image (xcb,ycb,zcb,0)
               call image (xdc,ydc,zdc,0)
               call image (xed,yed,zed,0)
            end if
            xt = yba*zcb - ycb*zba
            yt = zba*xcb - zcb*xba
            zt = xba*ycb - xcb*yba
            xu = ycb*zdc - ydc*zcb
            yu = zcb*xdc - zdc*xcb
            zu = xcb*ydc - xdc*ycb
            xtu = yt*zu - yu*zt
            ytu = zt*xu - zu*xt
            ztu = xt*yu - xu*yt
            rt2 = xt*xt + yt*yt + zt*zt
            ru2 = xu*xu + yu*yu + zu*zu          
            rtru = sqrt(rt2 * ru2)
            xv = ydc*zed - yed*zdc
            yv = zdc*xed - zed*xdc
            zv = xdc*yed - xed*ydc
            xuv = yu*zv - yv*zu
            yuv = zu*xv - zv*xu
            zuv = xu*yv - xv*yu
            rv2 = xv*xv + yv*yv + zv*zv
            rurv = sqrt(ru2 * rv2)
            if (rtru .ne. 0.0d0 .and. rurv.ne.0.0d0) then
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
               cosine1 = (xt*xu + yt*yu + zt*zu) / rtru
               cosine1 = min(1.0d0,max(-1.0d0,cosine1))
               angle1 = radian * acos(cosine1)
               sign = xba*xu + yba*yu + zba*zu
               if (sign .lt. 0.0d0)  angle1 = -angle1
               value1 = angle1
               rdc = sqrt(xdc*xdc + ydc*ydc + zdc*zdc)
               cosine2 = (xu*xv + yu*yv + zu*zv) / rurv
               cosine2 = min(1.0d0,max(-1.0d0,cosine2))
               angle2 = radian * acos(cosine2)
               sign = xcb*xv + ycb*yv + zcb*zv
               if (sign .lt. 0.0d0)  angle2 = -angle2
               value2 = angle2
c
c     check for inverted chirality at the central atom
c
               call chkttor (ib,ic,id,sign,value1,value2)
c
c     use bicubic interpolation to compute spline values
c
               nlo = 1
               nhi = tnx(k) 
               dowhile (nhi-nlo .gt. 1) 
                  nt = (nhi+nlo) / 2
                  if (ttx(nt,k) .gt. value1) then
                     nhi = nt
                  else  
                     nlo = nt
                  end if
               end do
               xlo = nlo
               nlo = 1
               nhi = tny(k) 
               dowhile (nhi-nlo .gt. 1) 
                  nt = (nhi + nlo)/2
                  if (tty(nt,k) .gt. value2) then
                     nhi = nt
                  else  
                     nlo = nt
                  end if
               end do
               ylo = nlo
               xh = ttx(xlo+1,k) - ttx(xlo,k)
               yh = tty(ylo+1,k) - tty(ylo,k)
               x1l = ttx(xlo,k)
               x1u = ttx(xlo+1,k)
               y1l = tty(ylo,k)
               y1u = tty(ylo+1,k)
               pos2 = ylo*tnx(k) + xlo
               pos1 = pos2 - tnx(k)
               ftt(1) = tbf(pos1,k)
               ftt(2) = tbf(pos1+1,k)
               ftt(3) = tbf(pos2+1,k)
               ftt(4) = tbf(pos2,k)
               ft1(1) = tbx(pos1,k)
               ft1(2) = tbx(pos1+1,k)
               ft1(3) = tbx(pos2+1,k)
               ft1(4) = tbx(pos2,k)
               ft2(1) = tby(pos1,k)
               ft2(2) = tby(pos1+1,k)
               ft2(3) = tby(pos2+1,k)
               ft2(4) = tby(pos2,k)
               ft12(1) = tbxy(pos1,k)
               ft12(2) = tbxy(pos1+1,k)
               ft12(3) = tbxy(pos2+1,k)
               ft12(4) = tbxy(pos2,k)
               call bcuint2 (ftt,ft1,ft2,ft12,x1l,x1u,y1l,y1u,
     &                       value1,value2,e,dedphi,dedpsi,
     &                       d12,d2edphi2,d2edpsi2)
               dedphi = sign * ttorunit * radian * dedphi
               dedpsi = sign * ttorunit * radian * dedpsi
               d2edphi2 = ttorunit * radian**2 * d2edphi2
               d2edpsi2 = ttorunit * radian**2 * d2edpsi2
               d12 = ttorunit * radian**2 * d12
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  dedphi = dedphi * fgrp
                  dedpsi = dedpsi * fgrp
                  d2edphi2 = d2edphi2 * fgrp
                  d2edpsi2 = d2edpsi2 * fgrp
                  d12 = d12 * fgrp
               end if
c
c     abbreviations for first derivative chain rule terms
c
               xca = xic - xia
               yca = yic - yia
               zca = zic - zia
               xdb = xid - xib
               ydb = yid - yib
               zdb = zid - zib
               if (use_polymer) then
                  call image (xca,yca,zca,0)
                  call image (xdb,ydb,zdb,0)
               end if
               dphidxt = (yt*zcb - ycb*zt) / (rt2*rcb)
               dphidyt = (zt*xcb - zcb*xt) / (rt2*rcb)
               dphidzt = (xt*ycb - xcb*yt) / (rt2*rcb)
               dphidxu = -(yu*zcb - ycb*zu) / (ru2*rcb)
               dphidyu = -(zu*xcb - zcb*xu) / (ru2*rcb)
               dphidzu = -(xu*ycb - xcb*yu) / (ru2*rcb)
c
c     abbreviations for second derivative chain rule terms
c
               xycb2 = xcb*xcb + ycb*ycb
               xzcb2 = xcb*xcb + zcb*zcb
               yzcb2 = ycb*ycb + zcb*zcb
               rcbxt = -2.0d0 * rcb * dphidxt
               rcbyt = -2.0d0 * rcb * dphidyt
               rcbzt = -2.0d0 * rcb * dphidzt
               rcbt2 = rcb * rt2
               rcbxu = 2.0d0 * rcb * dphidxu
               rcbyu = 2.0d0 * rcb * dphidyu
               rcbzu = 2.0d0 * rcb * dphidzu
               rcbu2 = rcb * ru2
               dphidxibt = yca*dphidzt - zca*dphidyt
               dphidxibu = zdc*dphidyu - ydc*dphidzu
               dphidyibt = zca*dphidxt - xca*dphidzt
               dphidyibu = xdc*dphidzu - zdc*dphidxu
               dphidzibt = xca*dphidyt - yca*dphidxt
               dphidzibu = ydc*dphidxu - xdc*dphidyu
               dphidxict = zba*dphidyt - yba*dphidzt
               dphidxicu = ydb*dphidzu - zdb*dphidyu
               dphidyict = xba*dphidzt - zba*dphidxt
               dphidyicu = zdb*dphidxu - xdb*dphidzu
               dphidzict = yba*dphidxt - xba*dphidyt
               dphidzicu = xdb*dphidyu - ydb*dphidxu
c
c     chain rule terms for first derivative components
c
               dphidxia = zcb*dphidyt - ycb*dphidzt
               dphidyia = xcb*dphidzt - zcb*dphidxt
               dphidzia = ycb*dphidxt - xcb*dphidyt
               dphidxib = dphidxibt + dphidxibu
               dphidyib = dphidyibt + dphidyibu
               dphidzib = dphidzibt + dphidzibu
               dphidxic = dphidxict + dphidxicu
               dphidyic = dphidyict + dphidyicu
               dphidzic = dphidzict + dphidzicu
               dphidxid = zcb*dphidyu - ycb*dphidzu
               dphidyid = xcb*dphidzu - zcb*dphidxu
               dphidzid = ycb*dphidxu - xcb*dphidyu
c
c     chain rule terms for second derivative components
c
               dxiaxia = rcbxt*dphidxia
               dxiayia = rcbxt*dphidyia - zcb*rcb/rt2
               dxiazia = rcbxt*dphidzia + ycb*rcb/rt2
               dxiaxic = rcbxt*dphidxict + xcb*xt/rcbt2
               dxiayic = rcbxt*dphidyict - dphidzt
     &                      - (xba*zcb*xcb+zba*yzcb2)/rcbt2
               dxiazic = rcbxt*dphidzict + dphidyt
     &                      + (xba*ycb*xcb+yba*yzcb2)/rcbt2
               dxiaxid = 0.0d0
               dxiayid = 0.0d0
               dxiazid = 0.0d0
               dyiayia = rcbyt*dphidyia
               dyiazia = rcbyt*dphidzia - xcb*rcb/rt2
               dyiaxib = rcbyt*dphidxibt - dphidzt
     &                      - (yca*zcb*ycb+zca*xzcb2)/rcbt2
               dyiaxic = rcbyt*dphidxict + dphidzt
     &                      + (yba*zcb*ycb+zba*xzcb2)/rcbt2
               dyiayic = rcbyt*dphidyict + ycb*yt/rcbt2
               dyiazic = rcbyt*dphidzict - dphidxt
     &                      - (yba*xcb*ycb+xba*xzcb2)/rcbt2
               dyiaxid = 0.0d0
               dyiayid = 0.0d0
               dyiazid = 0.0d0
               dziazia = rcbzt*dphidzia
               dziaxib = rcbzt*dphidxibt + dphidyt
     &                      + (zca*ycb*zcb+yca*xycb2)/rcbt2
               dziayib = rcbzt*dphidyibt - dphidxt
     &                      - (zca*xcb*zcb+xca*xycb2)/rcbt2
               dziaxic = rcbzt*dphidxict - dphidyt
     &                      - (zba*ycb*zcb+yba*xycb2)/rcbt2
               dziayic = rcbzt*dphidyict + dphidxt
     &                      + (zba*xcb*zcb+xba*xycb2)/rcbt2
               dziazic = rcbzt*dphidzict + zcb*zt/rcbt2
               dziaxid = 0.0d0
               dziayid = 0.0d0
               dziazid = 0.0d0
               dxibxic = -xcb*dphidxib/(rcb*rcb)
     &             - (yca*(zba*xcb+yt)-zca*(yba*xcb-zt))/rcbt2
     &             - 2.0d0*(yt*zba-yba*zt)*dphidxibt/rt2
     &             - (zdc*(ydb*xcb+zu)-ydc*(zdb*xcb-yu))/rcbu2
     &             + 2.0d0*(yu*zdb-ydb*zu)*dphidxibu/ru2
               dxibyic = -ycb*dphidxib/(rcb*rcb) + dphidzt + dphidzu
     &             - (yca*(zba*ycb-xt)+zca*(xba*xcb+zcb*zba))/rcbt2
     &             - 2.0d0*(zt*xba-zba*xt)*dphidxibt/rt2
     &             + (zdc*(xdb*xcb+zcb*zdb)+ydc*(zdb*ycb+xu))/rcbu2
     &             + 2.0d0*(zu*xdb-zdb*xu)*dphidxibu/ru2
               dxibxid = rcbxu*dphidxibu + xcb*xu/rcbu2
               dxibyid = rcbyu*dphidxibu - dphidzu
     &                      - (ydc*zcb*ycb+zdc*xzcb2)/rcbu2
               dxibzid = rcbzu*dphidxibu + dphidyu
     &                      + (zdc*ycb*zcb+ydc*xycb2)/rcbu2
               dyibzib = ycb*dphidzib/(rcb*rcb)
     &             - (xca*(xca*xcb+zcb*zca)+yca*(ycb*xca+zt))/rcbt2
     &             - 2.0d0*(xt*zca-xca*zt)*dphidzibt/rt2
     &             + (ydc*(xdc*ycb-zu)+xdc*(xdc*xcb+zcb*zdc))/rcbu2
     &             + 2.0d0*(xu*zdc-xdc*zu)*dphidzibu/ru2
               dyibxic = -xcb*dphidyib/(rcb*rcb) - dphidzt - dphidzu
     &             + (xca*(zba*xcb+yt)+zca*(zba*zcb+ycb*yba))/rcbt2
     &             - 2.0d0*(yt*zba-yba*zt)*dphidyibt/rt2
     &             - (zdc*(zdb*zcb+ycb*ydb)+xdc*(zdb*xcb-yu))/rcbu2
     &             + 2.0d0*(yu*zdb-ydb*zu)*dphidyibu/ru2
               dyibyic = -ycb*dphidyib/(rcb*rcb)
     &             - (zca*(xba*ycb+zt)-xca*(zba*ycb-xt))/rcbt2
     &             - 2.0d0*(zt*xba-zba*xt)*dphidyibt/rt2
     &             - (xdc*(zdb*ycb+xu)-zdc*(xdb*ycb-zu))/rcbu2
     &             + 2.0d0*(zu*xdb-zdb*xu)*dphidyibu/ru2
               dyibxid = rcbxu*dphidyibu + dphidzu
     &                      + (xdc*zcb*xcb+zdc*yzcb2)/rcbu2
               dyibyid = rcbyu*dphidyibu + ycb*yu/rcbu2
               dyibzid = rcbzu*dphidyibu - dphidxu
     &                      - (zdc*xcb*zcb+xdc*xycb2)/rcbu2
               dzibxic = -xcb*dphidzib/(rcb*rcb) + dphidyt + dphidyu
     &             - (xca*(yba*xcb-zt)+yca*(zba*zcb+ycb*yba))/rcbt2
     &             - 2.0d0*(yt*zba-yba*zt)*dphidzibt/rt2
     &             + (ydc*(zdb*zcb+ycb*ydb)+xdc*(ydb*xcb+zu))/rcbu2
     &             + 2.0d0*(yu*zdb-ydb*zu)*dphidzibu/ru2
               dzibzic = -zcb*dphidzib/(rcb*rcb)
     &             - (xca*(yba*zcb+xt)-yca*(xba*zcb-yt))/rcbt2
     &             - 2.0d0*(xt*yba-xba*yt)*dphidzibt/rt2
     &             - (ydc*(xdb*zcb+yu)-xdc*(ydb*zcb-xu))/rcbu2
     &             + 2.0d0*(xu*ydb-xdb*yu)*dphidzibu/ru2
               dzibxid = rcbxu*dphidzibu - dphidyu
     &                      - (xdc*ycb*xcb+ydc*yzcb2)/rcbu2
               dzibyid = rcbyu*dphidzibu + dphidxu
     &                      + (ydc*xcb*ycb+xdc*xzcb2)/rcbu2
               dzibzid = rcbzu*dphidzibu + zcb*zu/rcbu2
               dxicxid = rcbxu*dphidxicu - xcb*(zdb*ycb-ydb*zcb)/rcbu2
               dxicyid = rcbyu*dphidxicu + dphidzu
     &                      + (ydb*zcb*ycb+zdb*xzcb2)/rcbu2
               dxiczid = rcbzu*dphidxicu - dphidyu
     &                      - (zdb*ycb*zcb+ydb*xycb2)/rcbu2
               dyicxid = rcbxu*dphidyicu - dphidzu
     &                      - (xdb*zcb*xcb+zdb*yzcb2)/rcbu2
               dyicyid = rcbyu*dphidyicu - ycb*(xdb*zcb-zdb*xcb)/rcbu2
               dyiczid = rcbzu*dphidyicu + dphidxu
     &                      + (zdb*xcb*zcb+xdb*xycb2)/rcbu2
               dzicxid = rcbxu*dphidzicu + dphidyu
     &                      + (xdb*ycb*xcb+ydb*yzcb2)/rcbu2
               dzicyid = rcbyu*dphidzicu - dphidxu
     &                      - (ydb*xcb*ycb+xdb*xzcb2)/rcbu2
               dziczid = rcbzu*dphidzicu - zcb*(ydb*xcb-xdb*ycb)/rcbu2
               dxidxid = rcbxu*dphidxid
               dxidyid = rcbxu*dphidyid + zcb*rcb/ru2
               dxidzid = rcbxu*dphidzid - ycb*rcb/ru2
               dyidyid = rcbyu*dphidyid
               dyidzid = rcbyu*dphidzid + xcb*rcb/ru2
               dzidzid = rcbzu*dphidzid
c
c     get some second derivative chain rule terms by difference
c
               dxiaxib = -dxiaxia - dxiaxic - dxiaxid
               dxiayib = -dxiayia - dxiayic - dxiayid
               dxiazib = -dxiazia - dxiazic - dxiazid
               dyiayib = -dyiayia - dyiayic - dyiayid
               dyiazib = -dyiazia - dyiazic - dyiazid
               dziazib = -dziazia - dziazic - dziazid
               dxibxib = -dxiaxib - dxibxic - dxibxid
               dxibyib = -dyiaxib - dxibyic - dxibyid
               dxibzib = -dxiazib - dzibxic - dzibxid
               dxibzic = -dziaxib - dxibzib - dxibzid
               dyibyib = -dyiayib - dyibyic - dyibyid
               dyibzic = -dziayib - dyibzib - dyibzid
               dzibzib = -dziazib - dzibzic - dzibzid
               dzibyic = -dyiazib - dyibzib - dzibyid
               dxicxic = -dxiaxic - dxibxic - dxicxid
               dxicyic = -dyiaxic - dyibxic - dxicyid
               dxiczic = -dziaxic - dzibxic - dxiczid
               dyicyic = -dyiayic - dyibyic - dyicyid
               dyiczic = -dziayic - dzibyic - dyiczid
               dziczic = -dziazic - dzibzic - dziczid
c
c     increment diagonal and off-diagonal Hessian elements
c
               if (i .eq. ia) then
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxia
     &                             + d2edphi2*dphidxia*dphidxia
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayia
     &                             + d2edphi2*dphidxia*dphidyia
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazia
     &                             + d2edphi2*dphidxia*dphidzia
                  hessx(2,ia) = hessx(2,ia) + dedphi*dxiayia
     &                             + d2edphi2*dphidxia*dphidyia
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayia
     &                             + d2edphi2*dphidyia*dphidyia
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazia
     &                             + d2edphi2*dphidyia*dphidzia
                  hessx(3,ia) = hessx(3,ia) + dedphi*dxiazia
     &                             + d2edphi2*dphidxia*dphidzia
                  hessy(3,ia) = hessy(3,ia) + dedphi*dyiazia
     &                             + d2edphi2*dphidyia*dphidzia
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazia
     &                             + d2edphi2*dphidzia*dphidzia
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxiaxib
     &                             + d2edphi2*dphidxia*dphidxib
                  hessy(1,ib) = hessy(1,ib) + dedphi*dyiaxib
     &                             + d2edphi2*dphidyia*dphidxib
                  hessz(1,ib) = hessz(1,ib) + dedphi*dziaxib
     &                             + d2edphi2*dphidzia*dphidxib
                  hessx(2,ib) = hessx(2,ib) + dedphi*dxiayib
     &                             + d2edphi2*dphidxia*dphidyib
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyiayib
     &                             + d2edphi2*dphidyia*dphidyib
                  hessz(2,ib) = hessz(2,ib) + dedphi*dziayib
     &                             + d2edphi2*dphidzia*dphidyib
                  hessx(3,ib) = hessx(3,ib) + dedphi*dxiazib
     &                             + d2edphi2*dphidxia*dphidzib
                  hessy(3,ib) = hessy(3,ib) + dedphi*dyiazib
     &                             + d2edphi2*dphidyia*dphidzib
                  hessz(3,ib) = hessz(3,ib) + dedphi*dziazib
     &                             + d2edphi2*dphidzia*dphidzib
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxiaxic
     &                             + d2edphi2*dphidxia*dphidxic
                  hessy(1,ic) = hessy(1,ic) + dedphi*dyiaxic
     &                             + d2edphi2*dphidyia*dphidxic
                  hessz(1,ic) = hessz(1,ic) + dedphi*dziaxic
     &                             + d2edphi2*dphidzia*dphidxic
                  hessx(2,ic) = hessx(2,ic) + dedphi*dxiayic
     &                             + d2edphi2*dphidxia*dphidyic
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyiayic
     &                             + d2edphi2*dphidyia*dphidyic
                  hessz(2,ic) = hessz(2,ic) + dedphi*dziayic
     &                             + d2edphi2*dphidzia*dphidyic
                  hessx(3,ic) = hessx(3,ic) + dedphi*dxiazic
     &                             + d2edphi2*dphidxia*dphidzic
                  hessy(3,ic) = hessy(3,ic) + dedphi*dyiazic
     &                             + d2edphi2*dphidyia*dphidzic
                  hessz(3,ic) = hessz(3,ic) + dedphi*dziazic
     &                             + d2edphi2*dphidzia*dphidzic
                  hessx(1,id) = hessx(1,id) + dedphi*dxiaxid
     &                             + d2edphi2*dphidxia*dphidxid
                  hessy(1,id) = hessy(1,id) + dedphi*dyiaxid
     &                             + d2edphi2*dphidyia*dphidxid
                  hessz(1,id) = hessz(1,id) + dedphi*dziaxid
     &                             + d2edphi2*dphidzia*dphidxid
                  hessx(2,id) = hessx(2,id) + dedphi*dxiayid
     &                             + d2edphi2*dphidxia*dphidyid
                  hessy(2,id) = hessy(2,id) + dedphi*dyiayid
     &                             + d2edphi2*dphidyia*dphidyid
                  hessz(2,id) = hessz(2,id) + dedphi*dziayid
     &                             + d2edphi2*dphidzia*dphidyid
                  hessx(3,id) = hessx(3,id) + dedphi*dxiazid
     &                             + d2edphi2*dphidxia*dphidzid
                  hessy(3,id) = hessy(3,id) + dedphi*dyiazid
     &                             + d2edphi2*dphidyia*dphidzid
                  hessz(3,id) = hessz(3,id) + dedphi*dziazid
     &                             + d2edphi2*dphidzia*dphidzid
               else if (i .eq. ib) then
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxibxib
     &                             + d2edphi2*dphidxib*dphidxib
                  hessy(1,ib) = hessy(1,ib) + dedphi*dxibyib
     &                             + d2edphi2*dphidxib*dphidyib
                  hessz(1,ib) = hessz(1,ib) + dedphi*dxibzib
     &                             + d2edphi2*dphidxib*dphidzib
                  hessx(2,ib) = hessx(2,ib) + dedphi*dxibyib
     &                             + d2edphi2*dphidxib*dphidyib
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyibyib
     &                             + d2edphi2*dphidyib*dphidyib
                  hessz(2,ib) = hessz(2,ib) + dedphi*dyibzib
     &                             + d2edphi2*dphidyib*dphidzib
                  hessx(3,ib) = hessx(3,ib) + dedphi*dxibzib
     &                             + d2edphi2*dphidxib*dphidzib
                  hessy(3,ib) = hessy(3,ib) + dedphi*dyibzib
     &                             + d2edphi2*dphidyib*dphidzib
                  hessz(3,ib) = hessz(3,ib) + dedphi*dzibzib
     &                             + d2edphi2*dphidzib*dphidzib
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxib
     &                             + d2edphi2*dphidxib*dphidxia
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayib
     &                             + d2edphi2*dphidyib*dphidxia
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazib
     &                             + d2edphi2*dphidzib*dphidxia
                  hessx(2,ia) = hessx(2,ia) + dedphi*dyiaxib
     &                             + d2edphi2*dphidxib*dphidyia
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayib
     &                             + d2edphi2*dphidyib*dphidyia
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazib
     &                             + d2edphi2*dphidzib*dphidyia
                  hessx(3,ia) = hessx(3,ia) + dedphi*dziaxib
     &                             + d2edphi2*dphidxib*dphidzia
                  hessy(3,ia) = hessy(3,ia) + dedphi*dziayib
     &                             + d2edphi2*dphidyib*dphidzia
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazib
     &                             + d2edphi2*dphidzib*dphidzia
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxibxic
     &                             + d2edphi2*dphidxib*dphidxic
                  hessy(1,ic) = hessy(1,ic) + dedphi*dyibxic
     &                             + d2edphi2*dphidyib*dphidxic
                  hessz(1,ic) = hessz(1,ic) + dedphi*dzibxic
     &                             + d2edphi2*dphidzib*dphidxic
                  hessx(2,ic) = hessx(2,ic) + dedphi*dxibyic
     &                             + d2edphi2*dphidxib*dphidyic
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyibyic
     &                             + d2edphi2*dphidyib*dphidyic
                  hessz(2,ic) = hessz(2,ic) + dedphi*dzibyic
     &                             + d2edphi2*dphidzib*dphidyic
                  hessx(3,ic) = hessx(3,ic) + dedphi*dxibzic
     &                             + d2edphi2*dphidxib*dphidzic
                  hessy(3,ic) = hessy(3,ic) + dedphi*dyibzic
     &                             + d2edphi2*dphidyib*dphidzic
                  hessz(3,ic) = hessz(3,ic) + dedphi*dzibzic
     &                             + d2edphi2*dphidzib*dphidzic
                  hessx(1,id) = hessx(1,id) + dedphi*dxibxid
     &                             + d2edphi2*dphidxib*dphidxid
                  hessy(1,id) = hessy(1,id) + dedphi*dyibxid
     &                             + d2edphi2*dphidyib*dphidxid
                  hessz(1,id) = hessz(1,id) + dedphi*dzibxid
     &                             + d2edphi2*dphidzib*dphidxid
                  hessx(2,id) = hessx(2,id) + dedphi*dxibyid
     &                             + d2edphi2*dphidxib*dphidyid
                  hessy(2,id) = hessy(2,id) + dedphi*dyibyid
     &                             + d2edphi2*dphidyib*dphidyid
                  hessz(2,id) = hessz(2,id) + dedphi*dzibyid
     &                             + d2edphi2*dphidzib*dphidyid
                  hessx(3,id) = hessx(3,id) + dedphi*dxibzid
     &                             + d2edphi2*dphidxib*dphidzid
                  hessy(3,id) = hessy(3,id) + dedphi*dyibzid
     &                             + d2edphi2*dphidyib*dphidzid
                  hessz(3,id) = hessz(3,id) + dedphi*dzibzid
     &                             + d2edphi2*dphidzib*dphidzid
               else if (i .eq. ic) then
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxicxic
     &                             + d2edphi2*dphidxic*dphidxic
                  hessy(1,ic) = hessy(1,ic) + dedphi*dxicyic
     &                             + d2edphi2*dphidxic*dphidyic
                  hessz(1,ic) = hessz(1,ic) + dedphi*dxiczic
     &                             + d2edphi2*dphidxic*dphidzic
                  hessx(2,ic) = hessx(2,ic) + dedphi*dxicyic
     &                             + d2edphi2*dphidxic*dphidyic
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyicyic
     &                             + d2edphi2*dphidyic*dphidyic
                  hessz(2,ic) = hessz(2,ic) + dedphi*dyiczic
     &                             + d2edphi2*dphidyic*dphidzic
                  hessx(3,ic) = hessx(3,ic) + dedphi*dxiczic
     &                             + d2edphi2*dphidxic*dphidzic
                  hessy(3,ic) = hessy(3,ic) + dedphi*dyiczic
     &                             + d2edphi2*dphidyic*dphidzic
                  hessz(3,ic) = hessz(3,ic) + dedphi*dziczic
     &                             + d2edphi2*dphidzic*dphidzic
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxic
     &                             + d2edphi2*dphidxic*dphidxia
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayic
     &                             + d2edphi2*dphidyic*dphidxia
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazic
     &                             + d2edphi2*dphidzic*dphidxia
                  hessx(2,ia) = hessx(2,ia) + dedphi*dyiaxic
     &                             + d2edphi2*dphidxic*dphidyia
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayic
     &                             + d2edphi2*dphidyic*dphidyia
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazic
     &                             + d2edphi2*dphidzic*dphidyia
                  hessx(3,ia) = hessx(3,ia) + dedphi*dziaxic
     &                             + d2edphi2*dphidxic*dphidzia
                  hessy(3,ia) = hessy(3,ia) + dedphi*dziayic
     &                             + d2edphi2*dphidyic*dphidzia
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazic
     &                             + d2edphi2*dphidzic*dphidzia
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxibxic
     &                             + d2edphi2*dphidxic*dphidxib
                  hessy(1,ib) = hessy(1,ib) + dedphi*dxibyic
     &                             + d2edphi2*dphidyic*dphidxib
                  hessz(1,ib) = hessz(1,ib) + dedphi*dxibzic
     &                             + d2edphi2*dphidzic*dphidxib
                  hessx(2,ib) = hessx(2,ib) + dedphi*dyibxic
     &                             + d2edphi2*dphidxic*dphidyib
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyibyic
     &                             + d2edphi2*dphidyic*dphidyib
                  hessz(2,ib) = hessz(2,ib) + dedphi*dyibzic
     &                             + d2edphi2*dphidzic*dphidyib
                  hessx(3,ib) = hessx(3,ib) + dedphi*dzibxic
     &                             + d2edphi2*dphidxic*dphidzib
                  hessy(3,ib) = hessy(3,ib) + dedphi*dzibyic
     &                             + d2edphi2*dphidyic*dphidzib
                  hessz(3,ib) = hessz(3,ib) + dedphi*dzibzic
     &                             + d2edphi2*dphidzic*dphidzib
                  hessx(1,id) = hessx(1,id) + dedphi*dxicxid
     &                             + d2edphi2*dphidxic*dphidxid
                  hessy(1,id) = hessy(1,id) + dedphi*dyicxid
     &                             + d2edphi2*dphidyic*dphidxid
                  hessz(1,id) = hessz(1,id) + dedphi*dzicxid
     &                             + d2edphi2*dphidzic*dphidxid
                  hessx(2,id) = hessx(2,id) + dedphi*dxicyid
     &                             + d2edphi2*dphidxic*dphidyid
                  hessy(2,id) = hessy(2,id) + dedphi*dyicyid
     &                             + d2edphi2*dphidyic*dphidyid
                  hessz(2,id) = hessz(2,id) + dedphi*dzicyid
     &                             + d2edphi2*dphidzic*dphidyid
                  hessx(3,id) = hessx(3,id) + dedphi*dxiczid
     &                             + d2edphi2*dphidxic*dphidzid
                  hessy(3,id) = hessy(3,id) + dedphi*dyiczid
     &                             + d2edphi2*dphidyic*dphidzid
                  hessz(3,id) = hessz(3,id) + dedphi*dziczid
     &                             + d2edphi2*dphidzic*dphidzid
               else if (i .eq. id) then
                  hessx(1,id) = hessx(1,id) + dedphi*dxidxid
     &                             + d2edphi2*dphidxid*dphidxid
                  hessy(1,id) = hessy(1,id) + dedphi*dxidyid
     &                             + d2edphi2*dphidxid*dphidyid
                  hessz(1,id) = hessz(1,id) + dedphi*dxidzid
     &                             + d2edphi2*dphidxid*dphidzid
                  hessx(2,id) = hessx(2,id) + dedphi*dxidyid
     &                             + d2edphi2*dphidxid*dphidyid
                  hessy(2,id) = hessy(2,id) + dedphi*dyidyid
     &                             + d2edphi2*dphidyid*dphidyid
                  hessz(2,id) = hessz(2,id) + dedphi*dyidzid
     &                             + d2edphi2*dphidyid*dphidzid
                  hessx(3,id) = hessx(3,id) + dedphi*dxidzid
     &                             + d2edphi2*dphidxid*dphidzid
                  hessy(3,id) = hessy(3,id) + dedphi*dyidzid
     &                             + d2edphi2*dphidyid*dphidzid
                  hessz(3,id) = hessz(3,id) + dedphi*dzidzid
     &                             + d2edphi2*dphidzid*dphidzid
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxid
     &                             + d2edphi2*dphidxid*dphidxia
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayid
     &                             + d2edphi2*dphidyid*dphidxia
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazid
     &                             + d2edphi2*dphidzid*dphidxia
                  hessx(2,ia) = hessx(2,ia) + dedphi*dyiaxid
     &                             + d2edphi2*dphidxid*dphidyia
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayid
     &                             + d2edphi2*dphidyid*dphidyia
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazid
     &                             + d2edphi2*dphidzid*dphidyia
                  hessx(3,ia) = hessx(3,ia) + dedphi*dziaxid
     &                             + d2edphi2*dphidxid*dphidzia
                  hessy(3,ia) = hessy(3,ia) + dedphi*dziayid
     &                             + d2edphi2*dphidyid*dphidzia
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazid
     &                             + d2edphi2*dphidzid*dphidzia
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxibxid
     &                             + d2edphi2*dphidxid*dphidxib
                  hessy(1,ib) = hessy(1,ib) + dedphi*dxibyid
     &                             + d2edphi2*dphidyid*dphidxib
                  hessz(1,ib) = hessz(1,ib) + dedphi*dxibzid
     &                             + d2edphi2*dphidzid*dphidxib
                  hessx(2,ib) = hessx(2,ib) + dedphi*dyibxid
     &                             + d2edphi2*dphidxid*dphidyib
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyibyid
     &                             + d2edphi2*dphidyid*dphidyib
                  hessz(2,ib) = hessz(2,ib) + dedphi*dyibzid
     &                             + d2edphi2*dphidzid*dphidyib
                  hessx(3,ib) = hessx(3,ib) + dedphi*dzibxid
     &                             + d2edphi2*dphidxid*dphidzib
                  hessy(3,ib) = hessy(3,ib) + dedphi*dzibyid
     &                             + d2edphi2*dphidyid*dphidzib
                  hessz(3,ib) = hessz(3,ib) + dedphi*dzibzid
     &                             + d2edphi2*dphidzid*dphidzib
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxicxid
     &                             + d2edphi2*dphidxid*dphidxic
                  hessy(1,ic) = hessy(1,ic) + dedphi*dxicyid
     &                             + d2edphi2*dphidyid*dphidxic
                  hessz(1,ic) = hessz(1,ic) + dedphi*dxiczid
     &                             + d2edphi2*dphidzid*dphidxic
                  hessx(2,ic) = hessx(2,ic) + dedphi*dyicxid
     &                             + d2edphi2*dphidxid*dphidyic
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyicyid
     &                             + d2edphi2*dphidyid*dphidyic
                  hessz(2,ic) = hessz(2,ic) + dedphi*dyiczid
     &                             + d2edphi2*dphidzid*dphidyic
                  hessx(3,ic) = hessx(3,ic) + dedphi*dzicxid
     &                             + d2edphi2*dphidxid*dphidzic
                  hessy(3,ic) = hessy(3,ic) + dedphi*dzicyid
     &                             + d2edphi2*dphidyid*dphidzic
                  hessz(3,ic) = hessz(3,ic) + dedphi*dziczid
     &                             + d2edphi2*dphidzid*dphidzic
               end if
c
c     abbreviations for first derivative chain rule terms
c
               xec = xie - xic
               yec = yie - yic
               zec = zie - zic
               if (use_polymer) then
                  call image (xec,yec,zec,0)
               end if
               dpsidxu = (yu*zdc - ydc*zu) / (ru2*rdc)
               dpsidyu = (zu*xdc - zdc*xu) / (ru2*rdc)
               dpsidzu = (xu*ydc - xdc*yu) / (ru2*rdc)
               dpsidxv = -(yv*zdc - ydc*zv) / (rv2*rdc)
               dpsidyv = -(zv*xdc - zdc*xv) / (rv2*rdc)
               dpsidzv = -(xv*ydc - xdc*yv) / (rv2*rdc)
c
c     abbreviations for second derivative chain rule terms
c
               xydc2 = xdc*xdc + ydc*ydc
               xzdc2 = xdc*xdc + zdc*zdc
               yzdc2 = ydc*ydc + zdc*zdc
               rdcxu = -2.0d0 * rdc * dpsidxu
               rdcyu = -2.0d0 * rdc * dpsidyu
               rdczu = -2.0d0 * rdc * dpsidzu
               rdcu2 = rdc * ru2
               rdcxv = 2.0d0 * rdc * dpsidxv
               rdcyv = 2.0d0 * rdc * dpsidyv
               rdczv = 2.0d0 * rdc * dpsidzv
               rdcv2 = rdc * rv2
               dpsidxicu = ydb*dpsidzu - zdb*dpsidyu
               dpsidxicv = zed*dpsidyv - yed*dpsidzv
               dpsidyicu = zdb*dpsidxu - xdb*dpsidzu
               dpsidyicv = xed*dpsidzv - zed*dpsidxv
               dpsidzicu = xdb*dpsidyu - ydb*dpsidxu
               dpsidzicv = yed*dpsidxv - xed*dpsidyv
               dpsidxidu = zcb*dpsidyu - ycb*dpsidzu
               dpsidxidv = yec*dpsidzv - zec*dpsidyv
               dpsidyidu = xcb*dpsidzu - zcb*dpsidxu
               dpsidyidv = zec*dpsidxv - xec*dpsidzv
               dpsidzidu = ycb*dpsidxu - xcb*dpsidyu
               dpsidzidv = xec*dpsidyv - yec*dpsidxv
c
c     chain rule terms for first derivative components
c
               dpsidxib = zdc*dpsidyu - ydc*dpsidzu
               dpsidyib = xdc*dpsidzu - zdc*dpsidxu
               dpsidzib = ydc*dpsidxu - xdc*dpsidyu
               dpsidxic = dpsidxicu + dpsidxicv
               dpsidyic = dpsidyicu + dpsidyicv
               dpsidzic = dpsidzicu + dpsidzicv
               dpsidxid = dpsidxidu + dpsidxidv
               dpsidyid = dpsidyidu + dpsidyidv
               dpsidzid = dpsidzidu + dpsidzidv
               dpsidxie = zdc*dpsidyv - ydc*dpsidzv
               dpsidyie = xdc*dpsidzv - zdc*dpsidxv
               dpsidzie = ydc*dpsidxv - xdc*dpsidyv
c
c     chain rule terms for second derivative components
c
               dxibxib2 = rdcxu*dpsidxib
               dxibyib2 = rdcxu*dpsidyib - zdc*rdc/ru2
               dxibzib2 = rdcxu*dpsidzib + ydc*rdc/ru2
               dxibxid2 = rdcxu*dpsidxidu + xdc*xu/rdcu2
               dxibyid2 = rdcxu*dpsidyidu - dpsidzu
     &                       - (xcb*zdc*xdc+zcb*yzdc2)/rdcu2
               dxibzid2 = rdcxu*dpsidzidu + dpsidyu
     &                       + (xcb*ydc*xdc+ycb*yzdc2)/rdcu2
               dxibxie2 = 0.0d0
               dxibyie2 = 0.0d0
               dxibzie2 = 0.0d0
               dyibyib2 = rdcyu*dpsidyib
               dyibzib2 = rdcyu*dpsidzib - xdc*rdc/ru2
               dyibxic2 = rdcyu*dpsidxicu - dpsidzu
     &                       - (ydb*zdc*ydc+zdb*xzdc2)/rdcu2
               dyibxid2 = rdcyu*dpsidxidu + dpsidzu
     &                       + (ycb*zdc*ydc+zcb*xzdc2)/rdcu2
               dyibyid2 = rdcyu*dpsidyidu + ydc*yu/rdcu2
               dyibzid2 = rdcyu*dpsidzidu - dpsidxu
     &                       - (ycb*xdc*ydc+xcb*xzdc2)/rdcu2
               dyibxie2 = 0.0d0
               dyibyie2 = 0.0d0
               dyibzie2 = 0.0d0
               dzibzib2 = rdczu*dpsidzib
               dzibxic2 = rdczu*dpsidxicu + dpsidyu
     &                       + (zdb*ydc*zdc+ydb*xydc2)/rdcu2
               dzibyic2 = rdczu*dpsidyicu - dpsidxu
     &                       - (zdb*xdc*zdc+xdb*xydc2)/rdcu2
               dzibxid2 = rdczu*dpsidxidu - dpsidyu
     &                       - (zcb*ydc*zdc+ycb*xydc2)/rdcu2
               dzibyid2 = rdczu*dpsidyidu + dpsidxu
     &                       + (zcb*xdc*zdc+xcb*xydc2)/rdcu2
               dzibzid2 = rdczu*dpsidzidu + zdc*zu/rdcu2
               dzibxie2 = 0.0d0
               dzibyie2 = 0.0d0
               dzibzie2 = 0.0d0
               dxicxid2 = -xdc*dpsidxic/(rdc*rdc)
     &             - (ydb*(zcb*xdc+yu)-zdb*(ycb*xdc-zu))/rdcu2
     &             - 2.0d0*(yu*zcb-ycb*zu)*dpsidxicu/ru2
     &             - (zed*(yec*xdc+zv)-yed*(zec*xdc-yv))/rdcv2
     &             + 2.0d0*(yv*zec-yec*zv)*dpsidxicv/rv2
               dxicyid2 = -ydc*dpsidxic/(rdc*rdc) + dpsidzu + dpsidzv
     &             - (ydb*(zcb*ydc-xu)+zdb*(xcb*xdc+zdc*zcb))/rdcu2
     &             - 2.0d0*(zu*xcb-zcb*xu)*dpsidxicu/ru2
     &             + (zed*(xec*xdc+zdc*zec)+yed*(zec*ydc+xv))/rdcv2
     &             + 2.0d0*(zv*xec-zec*xv)*dpsidxicv/rv2
               dxicxie2 = rdcxv*dpsidxicv + xdc*xv/rdcv2
               dxicyie2 = rdcyv*dpsidxicv - dpsidzv
     &                       - (yed*zdc*ydc+zed*xzdc2)/rdcv2
               dxiczie2 = rdczv*dpsidxicv + dpsidyv
     &                       + (zed*ydc*zdc+yed*xydc2)/rdcv2
               dyiczic2 = ydc*dpsidzic/(rdc*rdc)
     &             - (xdb*(xdb*xdc+zdc*zdb)+ydb*(ydc*xdb+zu))/rdcu2
     &             - 2.0d0*(xu*zdb-xdb*zu)*dpsidzicu/ru2
     &             + (yed*(xed*ydc-zv)+xed*(xed*xdc+zdc*zed))/rdcv2
     &             + 2.0d0*(xv*zed-xed*zv)*dpsidzicv/rv2
               dyicxid2 = -xdc*dpsidyic/(rdc*rdc) - dpsidzu - dpsidzv
     &             + (xdb*(zcb*xdc+yu)+zdb*(zcb*zdc+ydc*ycb))/rdcu2
     &             - 2.0d0*(yu*zcb-ycb*zu)*dpsidyicu/ru2
     &             - (zed*(zec*zdc+ydc*yec)+xed*(zec*xdc-yv))/rdcv2
     &             + 2.0d0*(yv*zec-yec*zv)*dpsidyicv/rv2
               dyicyid2 = -ydc*dpsidyic/(rdc*rdc)
     &             - (zdb*(xcb*ydc+zu)-xdb*(zcb*ydc-xu))/rdcu2
     &             - 2.0d0*(zu*xcb-zcb*xu)*dpsidyicu/ru2
     &             - (xed*(zec*ydc+xv)-zed*(xec*ydc-zv))/rdcv2
     &             + 2.0d0*(zv*xec-zec*xv)*dpsidyicv/rv2
               dyicxie2 = rdcxv*dpsidyicv + dpsidzv
     &                       + (xed*zdc*xdc+zed*yzdc2)/rdcv2
               dyicyie2 = rdcyv*dpsidyicv + ydc*yv/rdcv2
               dyiczie2 = rdczv*dpsidyicv - dpsidxv
     &                       - (zed*xdc*zdc+xed*xydc2)/rdcv2
               dzicxid2 = -xdc*dpsidzic/(rdc*rdc) + dpsidyu + dpsidyv
     &             - (xdb*(ycb*xdc-zu)+ydb*(zcb*zdc+ydc*ycb))/rdcu2
     &             - 2.0d0*(yu*zcb-ycb*zu)*dpsidzicu/ru2
     &             + (yed*(zec*zdc+ydc*yec)+xed*(yec*xdc+zv))/rdcv2
     &             + 2.0d0*(yv*zec-yec*zv)*dpsidzicv/rv2
               dziczid2 = -zdc*dpsidzic/(rdc*rdc)
     &             - (xdb*(ycb*zdc+xu)-ydb*(xcb*zdc-yu))/rdcu2
     &             - 2.0d0*(xu*ycb-xcb*yu)*dpsidzicu/ru2
     &             - (yed*(xec*zdc+yv)-xed*(yec*zdc-xv))/rdcv2
     &             + 2.0d0*(xv*yec-xec*yv)*dpsidzicv/rv2
               dzicxie2 = rdcxv*dpsidzicv - dpsidyv
     &                       - (xed*ydc*xdc+yed*yzdc2)/rdcv2
               dzicyie2 = rdcyv*dpsidzicv + dpsidxv
     &                       + (yed*xdc*ydc+xed*xzdc2)/rdcv2
               dziczie2 = rdczv*dpsidzicv + zdc*zv/rdcv2
               dxidxie2 = rdcxv*dpsidxidv - xdc*(zec*ydc-yec*zdc)/rdcv2
               dxidyie2 = rdcyv*dpsidxidv + dpsidzv
     &                       + (yec*zdc*ydc+zec*xzdc2)/rdcv2
               dxidzie2 = rdczv*dpsidxidv - dpsidyv
     &                       - (zec*ydc*zdc+yec*xydc2)/rdcv2
               dyidxie2 = rdcxv*dpsidyidv - dpsidzv
     &                       - (xec*zdc*xdc+zec*yzdc2)/rdcv2
               dyidyie2 = rdcyv*dpsidyidv - ydc*(xec*zdc-zec*xdc)/rdcv2
               dyidzie2 = rdczv*dpsidyidv + dpsidxv
     &                       + (zec*xdc*zdc+xec*xydc2)/rdcv2
               dzidxie2 = rdcxv*dpsidzidv + dpsidyv
     &                       + (xec*ydc*xdc+yec*yzdc2)/rdcv2
               dzidyie2 = rdcyv*dpsidzidv - dpsidxv
     &                       - (yec*xdc*ydc+xec*xzdc2)/rdcv2
               dzidzie2 = rdczv*dpsidzidv - zdc*(yec*xdc-xec*ydc)/rdcv2
               dxiexie2 = rdcxv*dpsidxie
               dxieyie2 = rdcxv*dpsidyie + zdc*rdc/rv2
               dxiezie2 = rdcxv*dpsidzie - ydc*rdc/rv2
               dyieyie2 = rdcyv*dpsidyie
               dyiezie2 = rdcyv*dpsidzie + xdc*rdc/rv2
               dziezie2 = rdczv*dpsidzie
c
c     get some second derivative chain rule terms by difference
c
               dxibxic2 = -dxibxib2 - dxibxid2 - dxibxie2
               dxibyic2 = -dxibyib2 - dxibyid2 - dxibyie2
               dxibzic2 = -dxibzib2 - dxibzid2 - dxibzie2
               dyibyic2 = -dyibyib2 - dyibyid2 - dyibyie2
               dyibzic2 = -dyibzib2 - dyibzid2 - dyibzie2
               dzibzic2 = -dzibzib2 - dzibzid2 - dzibzie2
               dxicxic2 = -dxibxic2 - dxicxid2 - dxicxie2
               dxicyic2 = -dyibxic2 - dxicyid2 - dxicyie2
               dxiczic2 = -dxibzic2 - dzicxid2 - dzicxie2
               dxiczid2 = -dzibxic2 - dxiczic2 - dxiczie2
               dyicyic2 = -dyibyic2 - dyicyid2 - dyicyie2
               dyiczid2 = -dzibyic2 - dyiczic2 - dyiczie2
               dziczic2 = -dzibzic2 - dziczid2 - dziczie2
               dzicyid2 = -dyibzic2 - dyiczic2 - dzicyie2
               dxidxid2 = -dxibxid2 - dxicxid2 - dxidxie2
               dxidyid2 = -dyibxid2 - dyicxid2 - dxidyie2
               dxidzid2 = -dzibxid2 - dzicxid2 - dxidzie2
               dyidyid2 = -dyibyid2 - dyicyid2 - dyidyie2
               dyidzid2 = -dzibyid2 - dzicyid2 - dyidzie2
               dzidzid2 = -dzibzid2 - dziczid2 - dzidzie2
c
c     increment diagonal and off-diagonal Hessian elements
c
               if (i .eq. ia) then
                  hessx(1,ib) = hessx(1,ib) + d12*dphidxia*dpsidxib
                  hessy(1,ib) = hessy(1,ib) + d12*dphidyia*dpsidxib
                  hessz(1,ib) = hessz(1,ib) + d12*dphidzia*dpsidxib
                  hessx(2,ib) = hessx(2,ib) + d12*dphidxia*dpsidyib
                  hessy(2,ib) = hessy(2,ib) + d12*dphidyia*dpsidyib
                  hessz(2,ib) = hessz(2,ib) + d12*dphidzia*dpsidyib
                  hessx(3,ib) = hessx(3,ib) + d12*dphidxia*dpsidzib
                  hessy(3,ib) = hessy(3,ib) + d12*dphidyia*dpsidzib
                  hessz(3,ib) = hessz(3,ib) + d12*dphidzia*dpsidzib
                  hessx(1,ic) = hessx(1,ic) + d12*dphidxia*dpsidxic
                  hessy(1,ic) = hessy(1,ic) + d12*dphidyia*dpsidxic
                  hessz(1,ic) = hessz(1,ic) + d12*dphidzia*dpsidxic
                  hessx(2,ic) = hessx(2,ic) + d12*dphidxia*dpsidyic
                  hessy(2,ic) = hessy(2,ic) + d12*dphidyia*dpsidyic
                  hessz(2,ic) = hessz(2,ic) + d12*dphidzia*dpsidyic
                  hessx(3,ic) = hessx(3,ic) + d12*dphidxia*dpsidzic
                  hessy(3,ic) = hessy(3,ic) + d12*dphidyia*dpsidzic
                  hessz(3,ic) = hessz(3,ic) + d12*dphidzia*dpsidzic
                  hessx(1,id) = hessx(1,id) + d12*dphidxia*dpsidxid
                  hessy(1,id) = hessy(1,id) + d12*dphidyia*dpsidxid
                  hessz(1,id) = hessz(1,id) + d12*dphidzia*dpsidxid
                  hessx(2,id) = hessx(2,id) + d12*dphidxia*dpsidyid
                  hessy(2,id) = hessy(2,id) + d12*dphidyia*dpsidyid
                  hessz(2,id) = hessz(2,id) + d12*dphidzia*dpsidyid
                  hessx(3,id) = hessx(3,id) + d12*dphidxia*dpsidzid
                  hessy(3,id) = hessy(3,id) + d12*dphidyia*dpsidzid
                  hessz(3,id) = hessz(3,id) + d12*dphidzia*dpsidzid
                  hessx(1,ie) = hessx(1,ie) + d12*dphidxia*dpsidxie
                  hessy(1,ie) = hessy(1,ie) + d12*dphidyia*dpsidxie
                  hessz(1,ie) = hessz(1,ie) + d12*dphidzia*dpsidxie
                  hessx(2,ie) = hessx(2,ie) + d12*dphidxia*dpsidyie
                  hessy(2,ie) = hessy(2,ie) + d12*dphidyia*dpsidyie
                  hessz(2,ie) = hessz(2,ie) + d12*dphidzia*dpsidyie
                  hessx(3,ie) = hessx(3,ie) + d12*dphidxia*dpsidzie
                  hessy(3,ie) = hessy(3,ie) + d12*dphidyia*dpsidzie
                  hessz(3,ie) = hessz(3,ie) + d12*dphidzia*dpsidzie
               else if (i .eq. ib) then
                  hessx(1,ib) = hessx(1,ib) + dedpsi*dxibxib2
     &                             + d2edpsi2*dpsidxib*dpsidxib
     &                             + 2.0d0*d12*dphidxib*dpsidxib
                  hessy(1,ib) = hessy(1,ib) + dedpsi*dxibyib2
     &                             + d2edpsi2*dpsidxib*dpsidyib
     &                             + d12*dphidxib*dpsidyib
     &                             + d12*dpsidxib*dphidyib
                  hessz(1,ib) = hessz(1,ib) + dedpsi*dxibzib2
     &                             + d2edpsi2*dpsidxib*dpsidzib
     &                             + d12*dphidxib*dpsidzib
     &                             + d12*dpsidxib*dphidzib
                  hessx(2,ib) = hessx(2,ib) + dedpsi*dxibyib2
     &                             + d2edpsi2*dpsidxib*dpsidyib
     &                             + d12*dphidxib*dpsidyib
     &                             + d12*dpsidxib*dphidyib
                  hessy(2,ib) = hessy(2,ib) + dedpsi*dyibyib2
     &                             + d2edpsi2*dpsidyib*dpsidyib
     &                             + d12*dphidyib*dpsidyib
     &                             + d12*dpsidyib*dphidyib
                  hessz(2,ib) = hessz(2,ib) + dedpsi*dyibzib2
     &                             + d2edpsi2*dpsidyib*dpsidzib
     &                             + d12*dphidyib*dpsidzib
     &                             + d12*dpsidyib*dphidzib
                  hessx(3,ib) = hessx(3,ib) + dedpsi*dxibzib2
     &                             + d2edpsi2*dpsidxib*dpsidzib
     &                             + d12*dphidxib*dpsidzib
     &                             + d12*dpsidxib*dphidzib
                  hessy(3,ib) = hessy(3,ib) + dedpsi*dyibzib2
     &                             + d2edpsi2*dpsidyib*dpsidzib
     &                             + d12*dphidyib*dpsidzib
     &                             + d12*dpsidyib*dphidzib
                  hessz(3,ib) = hessz(3,ib) + dedpsi*dzibzib2
     &                             + d2edpsi2*dpsidzib*dpsidzib
     &                             + d12*dphidzib*dpsidzib
     &                             + d12*dpsidzib*dphidzib
                  hessx(1,ia) = hessx(1,ia) + d12*dpsidxib*dphidxia
                  hessy(1,ia) = hessy(1,ia) + d12*dpsidyib*dphidxia
                  hessz(1,ia) = hessz(1,ia) + d12*dpsidzib*dphidxia
                  hessx(2,ia) = hessx(2,ia) + d12*dpsidxib*dphidyia
                  hessy(2,ia) = hessy(2,ia) + d12*dpsidyib*dphidyia
                  hessz(2,ia) = hessz(2,ia) + d12*dpsidzib*dphidyia
                  hessx(3,ia) = hessx(3,ia) + d12*dpsidxib*dphidzia
                  hessy(3,ia) = hessy(3,ia) + d12*dpsidyib*dphidzia
                  hessz(3,ia) = hessz(3,ia) + d12*dpsidzib*dphidzia
                  hessx(1,ic) = hessx(1,ic) + dedpsi*dxibxic2
     &                             + d2edpsi2*dpsidxib*dpsidxic
     &                             + d12*dphidxib*dpsidxic
     &                             + d12*dpsidxib*dphidxic 
                  hessy(1,ic) = hessy(1,ic) + dedpsi*dyibxic2
     &                             + d2edpsi2*dpsidyib*dpsidxic
     &                             + d12*dphidyib*dpsidxic
     &                             + d12*dpsidyib*dphidxic 
                  hessz(1,ic) = hessz(1,ic) + dedpsi*dzibxic2
     &                             + d2edpsi2*dpsidzib*dpsidxic
     &                             + d12*dphidzib*dpsidxic
     &                             + d12*dpsidzib*dphidxic 
                  hessx(2,ic) = hessx(2,ic) + dedpsi*dxibyic2
     &                             + d2edpsi2*dpsidxib*dpsidyic
     &                             + d12*dphidxib*dpsidyic
     &                             + d12*dpsidxib*dphidyic 
                  hessy(2,ic) = hessy(2,ic) + dedpsi*dyibyic2
     &                             + d2edpsi2*dpsidyib*dpsidyic
     &                             + d12*dphidyib*dpsidyic
     &                             + d12*dpsidyib*dphidyic 
                  hessz(2,ic) = hessz(2,ic) + dedpsi*dzibyic2
     &                             + d2edpsi2*dpsidzib*dpsidyic
     &                             + d12*dphidzib*dpsidyic
     &                             + d12*dpsidzib*dphidyic 
                  hessx(3,ic) = hessx(3,ic) + dedpsi*dxibzic2
     &                             + d2edpsi2*dpsidxib*dpsidzic
     &                             + d12*dphidxib*dpsidzic
     &                             + d12*dpsidxib*dphidzic 
                  hessy(3,ic) = hessy(3,ic) + dedpsi*dyibzic2
     &                             + d2edpsi2*dpsidyib*dpsidzic
     &                             + d12*dphidyib*dpsidzic
     &                             + d12*dpsidyib*dphidzic 
                  hessz(3,ic) = hessz(3,ic) + dedpsi*dzibzic2
     &                             + d2edpsi2*dpsidzib*dpsidzic
     &                             + d12*dphidzib*dpsidzic
     &                             + d12*dpsidzib*dphidzic 
                  hessx(1,id) = hessx(1,id) + dedpsi*dxibxid2
     &                             + d2edpsi2*dpsidxib*dpsidxid
     &                             + d12*dphidxib*dpsidxid
     &                             + d12*dpsidxib*dphidxid 
                  hessy(1,id) = hessy(1,id) + dedpsi*dyibxid2
     &                             + d2edpsi2*dpsidyib*dpsidxid
     &                             + d12*dphidyib*dpsidxid
     &                             + d12*dpsidyib*dphidxid 
                  hessz(1,id) = hessz(1,id) + dedpsi*dzibxid2
     &                             + d2edpsi2*dpsidzib*dpsidxid
     &                             + d12*dphidzib*dpsidxid
     &                             + d12*dpsidzib*dphidxid 
                  hessx(2,id) = hessx(2,id) + dedpsi*dxibyid2
     &                             + d2edpsi2*dpsidxib*dpsidyid
     &                             + d12*dphidxib*dpsidyid
     &                             + d12*dpsidxib*dphidyid 
                  hessy(2,id) = hessy(2,id) + dedpsi*dyibyid2
     &                             + d2edpsi2*dpsidyib*dpsidyid
     &                             + d12*dphidyib*dpsidyid
     &                             + d12*dpsidyib*dphidyid 
                  hessz(2,id) = hessz(2,id) + dedpsi*dzibyid2
     &                             + d2edpsi2*dpsidzib*dpsidyid
     &                             + d12*dphidzib*dpsidyid
     &                             + d12*dpsidzib*dphidyid 
                  hessx(3,id) = hessx(3,id) + dedpsi*dxibzid2
     &                             + d2edpsi2*dpsidxib*dpsidzid
     &                             + d12*dphidxib*dpsidzid
     &                             + d12*dpsidxib*dphidzid 
                  hessy(3,id) = hessy(3,id) + dedpsi*dyibzid2
     &                             + d2edpsi2*dpsidyib*dpsidzid
     &                             + d12*dphidyib*dpsidzid
     &                             + d12*dpsidyib*dphidzid 
                  hessz(3,id) = hessz(3,id) + dedpsi*dzibzid2
     &                             + d2edpsi2*dpsidzib*dpsidzid
     &                             + d12*dphidzib*dpsidzid
     &                             + d12*dpsidzib*dphidzid 
                  hessx(1,ie) = hessx(1,ie) + dedpsi*dxibxie2
     &                             + d2edpsi2*dpsidxib*dpsidxie
     &                             + d12*dphidxib*dpsidxie 
                  hessy(1,ie) = hessy(1,ie) + dedpsi*dyibxie2
     &                             + d2edpsi2*dpsidyib*dpsidxie
     &                             + d12*dphidyib*dpsidxie 
                  hessz(1,ie) = hessz(1,ie) + dedpsi*dzibxie2
     &                             + d2edpsi2*dpsidzib*dpsidxie
     &                             + d12*dphidzib*dpsidxie 
                  hessx(2,ie) = hessx(2,ie) + dedpsi*dxibyie2
     &                             + d2edpsi2*dpsidxib*dpsidyie
     &                             + d12*dphidxib*dpsidyie 
                  hessy(2,ie) = hessy(2,ie) + dedpsi*dyibyie2
     &                             + d2edpsi2*dpsidyib*dpsidyie
     &                             + d12*dphidyib*dpsidyie 
                  hessz(2,ie) = hessz(2,ie) + dedpsi*dzibyie2
     &                             + d2edpsi2*dpsidzib*dpsidyie
     &                             + d12*dphidzib*dpsidyie 
                  hessx(3,ie) = hessx(3,ie) + dedpsi*dxibzie2
     &                             + d2edpsi2*dpsidxib*dpsidzie
     &                             + d12*dphidxib*dpsidzie 
                  hessy(3,ie) = hessy(3,ie) + dedpsi*dyibzie2
     &                             + d2edpsi2*dpsidyib*dpsidzie
     &                             + d12*dphidyib*dpsidzie 
                  hessz(3,ie) = hessz(3,ie) + dedpsi*dzibzie2
     &                             + d2edpsi2*dpsidzib*dpsidzie
     &                             + d12*dphidzib*dpsidzie 
               else if (i .eq. ic) then
                  hessx(1,ic) = hessx(1,ic) + dedpsi*dxicxic2
     &                             + d2edpsi2*dpsidxic*dpsidxic
     &                             + d12*dphidxic*dpsidxic
     &                             + d12*dpsidxic*dphidxic 
                  hessy(1,ic) = hessy(1,ic) + dedpsi*dxicyic2
     &                             + d2edpsi2*dpsidxic*dpsidyic
     &                             + d12*dphidxic*dpsidyic
     &                             + d12*dpsidxic*dphidyic 
                  hessz(1,ic) = hessz(1,ic) + dedpsi*dxiczic2
     &                             + d2edpsi2*dpsidxic*dpsidzic
     &                             + d12*dphidxic*dpsidzic
     &                             + d12*dpsidxic*dphidzic 
                  hessx(2,ic) = hessx(2,ic) + dedpsi*dxicyic2
     &                             + d2edpsi2*dpsidxic*dpsidyic
     &                             + d12*dphidxic*dpsidyic
     &                             + d12*dpsidxic*dphidyic 
                  hessy(2,ic) = hessy(2,ic) + dedpsi*dyicyic2
     &                             + d2edpsi2*dpsidyic*dpsidyic
     &                             + d12*dphidyic*dpsidyic
     &                             + d12*dpsidyic*dphidyic 
                  hessz(2,ic) = hessz(2,ic) + dedpsi*dyiczic2
     &                             + d2edpsi2*dpsidyic*dpsidzic
     &                             + d12*dphidyic*dpsidzic
     &                             + d12*dpsidyic*dphidzic 
                  hessx(3,ic) = hessx(3,ic) + dedpsi*dxiczic2
     &                             + d2edpsi2*dpsidxic*dpsidzic
     &                             + d12*dphidxic*dpsidzic
     &                             + d12*dpsidxic*dphidzic 
                  hessy(3,ic) = hessy(3,ic) + dedpsi*dyiczic2
     &                             + d2edpsi2*dpsidyic*dpsidzic
     &                             + d12*dphidyic*dpsidzic
     &                             + d12*dpsidyic*dphidzic 
                  hessz(3,ic) = hessz(3,ic) + dedpsi*dziczic2
     &                             + d2edpsi2*dpsidzic*dpsidzic
     &                             + d12*dphidzic*dpsidzic
     &                             + d12*dpsidzic*dphidzic 
                  hessx(1,ia) = hessx(1,ia) + d12*dpsidxic*dphidxia 
                  hessy(1,ia) = hessy(1,ia) + d12*dpsidyic*dphidxia 
                  hessz(1,ia) = hessz(1,ia) + d12*dpsidzic*dphidxia 
                  hessx(2,ia) = hessx(2,ia) + d12*dpsidxic*dphidyia 
                  hessy(2,ia) = hessy(2,ia) + d12*dpsidyic*dphidyia 
                  hessz(2,ia) = hessz(2,ia) + d12*dpsidzic*dphidyia 
                  hessx(3,ia) = hessx(3,ia) + d12*dpsidxic*dphidzia 
                  hessy(3,ia) = hessy(3,ia) + d12*dpsidyic*dphidzia 
                  hessz(3,ia) = hessz(3,ia) + d12*dpsidzic*dphidzia 
                  hessx(1,ib) = hessx(1,ib) + dedpsi*dxibxic2
     &                             + d2edpsi2*dpsidxic*dpsidxib
     &                             + d12*dphidxic*dpsidxib
     &                             + d12*dpsidxic*dphidxib 
                  hessy(1,ib) = hessy(1,ib) + dedpsi*dxibyic2
     &                             + d2edpsi2*dpsidyic*dpsidxib
     &                             + d12*dphidyic*dpsidxib
     &                             + d12*dpsidyic*dphidxib 
                  hessz(1,ib) = hessz(1,ib) + dedpsi*dxibzic2
     &                             + d2edpsi2*dpsidzic*dpsidxib
     &                             + d12*dphidzic*dpsidxib
     &                             + d12*dpsidzic*dphidxib 
                  hessx(2,ib) = hessx(2,ib) + dedpsi*dyibxic2
     &                             + d2edpsi2*dpsidxic*dpsidyib
     &                             + d12*dphidxic*dpsidyib
     &                             + d12*dpsidxic*dphidyib 
                  hessy(2,ib) = hessy(2,ib) + dedpsi*dyibyic2
     &                             + d2edpsi2*dpsidyic*dpsidyib
     &                             + d12*dphidyic*dpsidyib
     &                             + d12*dpsidyic*dphidyib 
                  hessz(2,ib) = hessz(2,ib) + dedpsi*dyibzic2
     &                             + d2edpsi2*dpsidzic*dpsidyib
     &                             + d12*dphidzic*dpsidyib
     &                             + d12*dpsidzic*dphidyib 
                  hessx(3,ib) = hessx(3,ib) + dedpsi*dzibxic2
     &                             + d2edpsi2*dpsidxic*dpsidzib
     &                             + d12*dphidxic*dpsidzib
     &                             + d12*dpsidxic*dphidzib 
                  hessy(3,ib) = hessy(3,ib) + dedpsi*dzibyic2
     &                             + d2edpsi2*dpsidyic*dpsidzib
     &                             + d12*dphidyic*dpsidzib
     &                             + d12*dpsidyic*dphidzib 
                  hessz(3,ib) = hessz(3,ib) + dedpsi*dzibzic2
     &                             + d2edpsi2*dpsidzic*dpsidzib
     &                             + d12*dphidzic*dpsidzib
     &                             + d12*dpsidzic*dphidzib 
                  hessx(1,id) = hessx(1,id) + dedpsi*dxicxid2
     &                             + d2edpsi2*dpsidxic*dpsidxid
     &                             + d12*dphidxic*dpsidxid
     &                             + d12*dpsidxic*dphidxid 
                  hessy(1,id) = hessy(1,id) + dedpsi*dyicxid2
     &                             + d2edpsi2*dpsidyic*dpsidxid
     &                             + d12*dphidyic*dpsidxid
     &                             + d12*dpsidyic*dphidxid 
                  hessz(1,id) = hessz(1,id) + dedpsi*dzicxid2
     &                             + d2edpsi2*dpsidzic*dpsidxid
     &                             + d12*dphidzic*dpsidxid
     &                             + d12*dpsidzic*dphidxid 
                  hessx(2,id) = hessx(2,id) + dedpsi*dxicyid2
     &                             + d2edpsi2*dpsidxic*dpsidyid
     &                             + d12*dphidxic*dpsidyid
     &                             + d12*dpsidxic*dphidyid 
                  hessy(2,id) = hessy(2,id) + dedpsi*dyicyid2
     &                             + d2edpsi2*dpsidyic*dpsidyid
     &                             + d12*dphidyic*dpsidyid
     &                             + d12*dpsidyic*dphidyid 
                  hessz(2,id) = hessz(2,id) + dedpsi*dzicyid2
     &                             + d2edpsi2*dpsidzic*dpsidyid
     &                             + d12*dphidzic*dpsidyid
     &                             + d12*dpsidzic*dphidyid 
                  hessx(3,id) = hessx(3,id) + dedpsi*dxiczid2
     &                             + d2edpsi2*dpsidxic*dpsidzid
     &                             + d12*dphidxic*dpsidzid
     &                             + d12*dpsidxic*dphidzid 
                  hessy(3,id) = hessy(3,id) + dedpsi*dyiczid2
     &                             + d2edpsi2*dpsidyic*dpsidzid
     &                             + d12*dphidyic*dpsidzid
     &                             + d12*dpsidyic*dphidzid 
                  hessz(3,id) = hessz(3,id) + dedpsi*dziczid2
     &                             + d2edpsi2*dpsidzic*dpsidzid
     &                             + d12*dphidzic*dpsidzid
     &                             + d12*dpsidzic*dphidzid 
                  hessx(1,ie) = hessx(1,ie) + dedpsi*dxicxie2
     &                             + d2edpsi2*dpsidxic*dpsidxie
     &                             + d12*dphidxic*dpsidxie 
                  hessy(1,ie) = hessy(1,ie) + dedpsi*dyicxie2
     &                             + d2edpsi2*dpsidyic*dpsidxie
     &                             + d12*dphidyic*dpsidxie 
                  hessz(1,ie) = hessz(1,ie) + dedpsi*dzicxie2
     &                             + d2edpsi2*dpsidzic*dpsidxie
     &                             + d12*dphidzic*dpsidxie 
                  hessx(2,ie) = hessx(2,ie) + dedpsi*dxicyie2
     &                             + d2edpsi2*dpsidxic*dpsidyie
     &                             + d12*dphidxic*dpsidyie 
                  hessy(2,ie) = hessy(2,ie) + dedpsi*dyicyie2
     &                             + d2edpsi2*dpsidyic*dpsidyie
     &                             + d12*dphidyic*dpsidyie 
                  hessz(2,ie) = hessz(2,ie) + dedpsi*dzicyie2
     &                             + d2edpsi2*dpsidzic*dpsidyie
     &                             + d12*dphidzic*dpsidyie 
                  hessx(3,ie) = hessx(3,ie) + dedpsi*dxiczie2
     &                             + d2edpsi2*dpsidxic*dpsidzie
     &                             + d12*dphidxic*dpsidzie 
                  hessy(3,ie) = hessy(3,ie) + dedpsi*dyiczie2
     &                             + d2edpsi2*dpsidyic*dpsidzie
     &                             + d12*dphidyic*dpsidzie 
                  hessz(3,ie) = hessz(3,ie) + dedpsi*dziczie2
     &                             + d2edpsi2*dpsidzic*dpsidzie
     &                             + d12*dphidzic*dpsidzie 
               else if (i .eq. id) then
                  hessx(1,id) = hessx(1,id) + dedpsi*dxidxid2
     &                             + d2edpsi2*dpsidxid*dpsidxid
     &                             + d12*dphidxid*dpsidxid
     &                             + d12*dpsidxid*dphidxid 
                  hessy(1,id) = hessy(1,id) + dedpsi*dxidyid2
     &                             + d2edpsi2*dpsidxid*dpsidyid
     &                             + d12*dphidxid*dpsidyid
     &                             + d12*dpsidxid*dphidyid 
                  hessz(1,id) = hessz(1,id) + dedpsi*dxidzid2
     &                             + d2edpsi2*dpsidxid*dpsidzid
     &                             + d12*dphidxid*dpsidzid
     &                             + d12*dpsidxid*dphidzid 
                  hessx(2,id) = hessx(2,id) + dedpsi*dxidyid2
     &                             + d2edpsi2*dpsidxid*dpsidyid
     &                             + d12*dphidxid*dpsidyid
     &                             + d12*dpsidxid*dphidyid 
                  hessy(2,id) = hessy(2,id) + dedpsi*dyidyid2
     &                             + d2edpsi2*dpsidyid*dpsidyid
     &                             + d12*dphidyid*dpsidyid
     &                             + d12*dpsidyid*dphidyid 
                  hessz(2,id) = hessz(2,id) + dedpsi*dyidzid2
     &                             + d2edpsi2*dpsidyid*dpsidzid
     &                             + d12*dphidyid*dpsidzid
     &                             + d12*dpsidyid*dphidzid 
                  hessx(3,id) = hessx(3,id) + dedpsi*dxidzid2
     &                             + d2edpsi2*dpsidxid*dpsidzid
     &                             + d12*dphidxid*dpsidzid
     &                             + d12*dpsidxid*dphidzid 
                  hessy(3,id) = hessy(3,id) + dedpsi*dyidzid2
     &                             + d2edpsi2*dpsidyid*dpsidzid
     &                             + d12*dphidyid*dpsidzid
     &                             + d12*dpsidyid*dphidzid 
                  hessz(3,id) = hessz(3,id) + dedpsi*dzidzid2
     &                             + d2edpsi2*dpsidzid*dpsidzid
     &                             + d12*dphidzid*dpsidzid
     &                             + d12*dpsidzid*dphidzid 
                  hessx(1,ia) = hessx(1,ia) +  d12*dpsidxid*dphidxia
                  hessy(1,ia) = hessy(1,ia) +  d12*dpsidyid*dphidxia  
                  hessz(1,ia) = hessz(1,ia) +  d12*dpsidzid*dphidxia
                  hessx(2,ia) = hessx(2,ia) +  d12*dpsidxid*dphidyia
                  hessy(2,ia) = hessy(2,ia) +  d12*dpsidyid*dphidyia
                  hessz(2,ia) = hessz(2,ia) +  d12*dpsidzid*dphidyia
                  hessx(3,ia) = hessx(3,ia) +  d12*dpsidxid*dphidzia
                  hessy(3,ia) = hessy(3,ia) +  d12*dpsidyid*dphidzia
                  hessz(3,ia) = hessz(3,ia) +  d12*dpsidzid*dphidzia
                  hessx(1,ib) = hessx(1,ib) + dedpsi*dxibxid2
     &                             + d2edpsi2*dpsidxid*dpsidxib
     &                             + d12*dphidxid*dpsidxib
     &                             + d12*dpsidxid*dphidxib 
                  hessy(1,ib) = hessy(1,ib) + dedpsi*dxibyid2
     &                             + d2edpsi2*dpsidyid*dpsidxib
     &                             + d12*dphidyid*dpsidxib
     &                             + d12*dpsidyid*dphidxib 
                  hessz(1,ib) = hessz(1,ib) + dedpsi*dxibzid2
     &                             + d2edpsi2*dpsidzid*dpsidxib
     &                             + d12*dphidzid*dpsidxib
     &                             + d12*dpsidzid*dphidxib 
                  hessx(2,ib) = hessx(2,ib) + dedpsi*dyibxid2
     &                             + d2edpsi2*dpsidxid*dpsidyib
     &                             + d12*dphidxid*dpsidyib
     &                             + d12*dpsidxid*dphidyib 
                  hessy(2,ib) = hessy(2,ib) + dedpsi*dyibyid2
     &                             + d2edpsi2*dpsidyid*dpsidyib
     &                             + d12*dphidyid*dpsidyib
     &                             + d12*dpsidyid*dphidyib 
                  hessz(2,ib) = hessz(2,ib) + dedpsi*dyibzid2
     &                             + d2edpsi2*dpsidzid*dpsidyib
     &                             + d12*dphidzid*dpsidyib
     &                             + d12*dpsidzid*dphidyib 
                  hessx(3,ib) = hessx(3,ib) + dedpsi*dzibxid2
     &                             + d2edpsi2*dpsidxid*dpsidzib
     &                             + d12*dphidxid*dpsidzib
     &                             + d12*dpsidxid*dphidzib 
                  hessy(3,ib) = hessy(3,ib) + dedpsi*dzibyid2
     &                             + d2edpsi2*dpsidyid*dpsidzib
     &                             + d12*dphidyid*dpsidzib
     &                             + d12*dpsidyid*dphidzib 
                  hessz(3,ib) = hessz(3,ib) + dedpsi*dzibzid2
     &                             + d2edpsi2*dpsidzid*dpsidzib
     &                             + d12*dphidzid*dpsidzib
     &                             + d12*dpsidzid*dphidzib 
                  hessx(1,ic) = hessx(1,ic) + dedpsi*dxicxid2
     &                             + d2edpsi2*dpsidxid*dpsidxic
     &                             + d12*dphidxid*dpsidxic
     &                             + d12*dpsidxid*dphidxic 
                  hessy(1,ic) = hessy(1,ic) + dedpsi*dxicyid2
     &                             + d2edpsi2*dpsidyid*dpsidxic
     &                             + d12*dphidyid*dpsidxic
     &                             + d12*dpsidyid*dphidxic 
                  hessz(1,ic) = hessz(1,ic) + dedpsi*dxiczid2
     &                             + d2edpsi2*dpsidzid*dpsidxic
     &                             + d12*dphidzid*dpsidxic
     &                             + d12*dpsidzid*dphidxic 
                  hessx(2,ic) = hessx(2,ic) + dedpsi*dyicxid2
     &                             + d2edpsi2*dpsidxid*dpsidyic
     &                             + d12*dphidxid*dpsidyic
     &                             + d12*dpsidxid*dphidyic 
                  hessy(2,ic) = hessy(2,ic) + dedpsi*dyicyid2
     &                             + d2edpsi2*dpsidyid*dpsidyic
     &                             + d12*dphidyid*dpsidyic
     &                             + d12*dpsidyid*dphidyic 
                  hessz(2,ic) = hessz(2,ic) + dedpsi*dyiczid2
     &                             + d2edpsi2*dpsidzid*dpsidyic
     &                             + d12*dphidzid*dpsidyic
     &                             + d12*dpsidzid*dphidyic 
                  hessx(3,ic) = hessx(3,ic) + dedpsi*dzicxid2
     &                             + d2edpsi2*dpsidxid*dpsidzic
     &                             + d12*dphidxid*dpsidzic
     &                             + d12*dpsidxid*dphidzic 
                  hessy(3,ic) = hessy(3,ic) + dedpsi*dzicyid2
     &                             + d2edpsi2*dpsidyid*dpsidzic
     &                             + d12*dphidyid*dpsidzic
     &                             + d12*dpsidyid*dphidzic 
                  hessz(3,ic) = hessz(3,ic) + dedpsi*dziczid2
     &                             + d2edpsi2*dpsidzid*dpsidzic
     &                             + d12*dphidzid*dpsidzic
     &                             + d12*dpsidzid*dphidzic 
                  hessx(1,ie) = hessx(1,ie) + dedpsi*dxidxie2
     &                             + d2edpsi2*dpsidxid*dpsidxie
     &                             + d12*dphidxid*dpsidxie 
                  hessy(1,ie) = hessy(1,ie) + dedpsi*dyidxie2
     &                             + d2edpsi2*dpsidyid*dpsidxie
     &                             + d12*dphidyid*dpsidxie 
                  hessz(1,ie) = hessz(1,ie) + dedpsi*dzidxie2
     &                             + d2edpsi2*dpsidzid*dpsidxie
     &                             + d12*dphidzid*dpsidxie 
                  hessx(2,ie) = hessx(2,ie) + dedpsi*dxidyie2
     &                             + d2edpsi2*dpsidxid*dpsidyie
     &                             + d12*dphidxid*dpsidyie 
                  hessy(2,ie) = hessy(2,ie) + dedpsi*dyidyie2
     &                             + d2edpsi2*dpsidyid*dpsidyie
     &                             + d12*dphidyid*dpsidyie 
                  hessz(2,ie) = hessz(2,ie) + dedpsi*dzidyie2
     &                             + d2edpsi2*dpsidzid*dpsidyie
     &                             + d12*dphidzid*dpsidyie 
                  hessx(3,ie) = hessx(3,ie) + dedpsi*dxidzie2
     &                             + d2edpsi2*dpsidxid*dpsidzie
     &                             + d12*dphidxid*dpsidzie 
                  hessy(3,ie) = hessy(3,ie) + dedpsi*dyidzie2
     &                             + d2edpsi2*dpsidyid*dpsidzie
     &                             + d12*dphidyid*dpsidzie 
                  hessz(3,ie) = hessz(3,ie) + dedpsi*dzidzie2
     &                             + d2edpsi2*dpsidzid*dpsidzie
     &                             + d12*dphidzid*dpsidzie
               else if (i .eq. ie) then
                  hessx(1,ie) = hessx(1,ie) + dedpsi*dxiexie2
     &                             + d2edpsi2*dpsidxie*dpsidxie
                  hessy(1,ie) = hessy(1,ie) + dedpsi*dxieyie2
     &                             + d2edpsi2*dpsidxie*dpsidyie
                  hessz(1,ie) = hessz(1,ie) + dedpsi*dxiezie2
     &                             + d2edpsi2*dpsidxie*dpsidzie
                  hessx(2,ie) = hessx(2,ie) + dedpsi*dxieyie2
     &                             + d2edpsi2*dpsidxie*dpsidyie
                  hessy(2,ie) = hessy(2,ie) + dedpsi*dyieyie2
     &                             + d2edpsi2*dpsidyie*dpsidyie
                  hessz(2,ie) = hessz(2,ie) + dedpsi*dyiezie2
     &                             + d2edpsi2*dpsidyie*dpsidzie
                  hessx(3,ie) = hessx(3,ie) + dedpsi*dxiezie2
     &                             + d2edpsi2*dpsidxie*dpsidzie
                  hessy(3,ie) = hessy(3,ie) + dedpsi*dyiezie2
     &                             + d2edpsi2*dpsidyie*dpsidzie
                  hessz(3,ie) = hessz(3,ie) + dedpsi*dziezie2
     &                             + d2edpsi2*dpsidzie*dpsidzie
                  hessx(1,ia) = hessx(1,ia) + d12*dpsidxie*dphidxia
                  hessy(1,ia) = hessy(1,ia) + d12*dpsidyie*dphidxia  
                  hessz(1,ia) = hessz(1,ia) + d12*dpsidzie*dphidxia
                  hessx(2,ia) = hessx(2,ia) + d12*dpsidxie*dphidyia
                  hessy(2,ia) = hessy(2,ia) + d12*dpsidyie*dphidyia
                  hessz(2,ia) = hessz(2,ia) + d12*dpsidzie*dphidyia
                  hessx(3,ia) = hessx(3,ia) + d12*dpsidxie*dphidzia
                  hessy(3,ia) = hessy(3,ia) + d12*dpsidyie*dphidzia
                  hessz(3,ia) = hessz(3,ia) + d12*dpsidzie*dphidzia
                  hessx(1,ib) = hessx(1,ib) + dedpsi*dxibxie2
     &                             + d2edpsi2*dpsidxie*dpsidxib
     &                             + d12*dpsidxie*dphidxib 
                  hessy(1,ib) = hessy(1,ib) + dedpsi*dxibyie2
     &                             + d2edpsi2*dpsidyie*dpsidxib
     &                             + d12*dpsidyie*dphidxib 
                  hessz(1,ib) = hessz(1,ib) + dedpsi*dxibzie2
     &                             + d2edpsi2*dpsidzie*dpsidxib
     &                             + d12*dpsidzie*dphidxib 
                  hessx(2,ib) = hessx(2,ib) + dedpsi*dyibxie2
     &                             + d2edpsi2*dpsidxie*dpsidyib
     &                             + d12*dpsidxie*dphidyib 
                  hessy(2,ib) = hessy(2,ib) + dedpsi*dyibyie2
     &                             + d2edpsi2*dpsidyie*dpsidyib
     &                             + d12*dpsidyie*dphidyib 
                  hessz(2,ib) = hessz(2,ib) + dedpsi*dyibzie2
     &                             + d2edpsi2*dpsidzie*dpsidyib
     &                             + d12*dpsidzie*dphidyib 
                  hessx(3,ib) = hessx(3,ib) + dedpsi*dzibxie2
     &                             + d2edpsi2*dpsidxie*dpsidzib
     &                             + d12*dpsidxie*dphidzib 
                  hessy(3,ib) = hessy(3,ib) + dedpsi*dzibyie2
     &                             + d2edpsi2*dpsidyie*dpsidzib
     &                             + d12*dpsidyie*dphidzib 
                  hessz(3,ib) = hessz(3,ib) + dedpsi*dzibzie2
     &                             + d2edpsi2*dpsidzie*dpsidzib
     &                             + d12*dpsidzie*dphidzib 
                  hessx(1,ic) = hessx(1,ic) + dedpsi*dxicxie2
     &                             + d2edpsi2*dpsidxie*dpsidxic
     &                             + d12*dpsidxie*dphidxic 
                  hessy(1,ic) = hessy(1,ic) + dedpsi*dxicyie2
     &                             + d2edpsi2*dpsidyie*dpsidxic
     &                             + d12*dpsidyie*dphidxic 
                  hessz(1,ic) = hessz(1,ic) + dedpsi*dxiczie2
     &                             + d2edpsi2*dpsidzie*dpsidxic
     &                             + d12*dpsidzie*dphidxic 
                  hessx(2,ic) = hessx(2,ic) + dedpsi*dyicxie2
     &                             + d2edpsi2*dpsidxie*dpsidyic
     &                             + d12*dpsidxie*dphidyic 
                  hessy(2,ic) = hessy(2,ic) + dedpsi*dyicyie2
     &                             + d2edpsi2*dpsidyie*dpsidyic
     &                             + d12*dpsidyie*dphidyic 
                  hessz(2,ic) = hessz(2,ic) + dedpsi*dyiczie2
     &                             + d2edpsi2*dpsidzie*dpsidyic
     &                             + d12*dpsidzie*dphidyic 
                  hessx(3,ic) = hessx(3,ic) + dedpsi*dzicxie2
     &                             + d2edpsi2*dpsidxie*dpsidzic
     &                             + d12*dpsidxie*dphidzic 
                  hessy(3,ic) = hessy(3,ic) + dedpsi*dzicyie2
     &                             + d2edpsi2*dpsidyie*dpsidzic
     &                             + d12*dpsidyie*dphidzic 
                  hessz(3,ic) = hessz(3,ic) + dedpsi*dziczie2
     &                             + d2edpsi2*dpsidzie*dpsidzic
     &                             + d12*dpsidzie*dphidzic 
                  hessx(1,id) = hessx(1,id) + dedpsi*dxidxie2
     &                             + d2edpsi2*dpsidxid*dpsidxie
     &                             + d12*dphidxid*dpsidxie 
                  hessy(1,id) = hessy(1,id) + dedpsi*dxidyie2
     &                             + d2edpsi2*dpsidxid*dpsidyie
     &                             + d12*dphidxid*dpsidyie 
                  hessz(1,id) = hessz(1,id) + dedpsi*dxidzie2
     &                             + d2edpsi2*dpsidxid*dpsidzie
     &                             + d12*dphidxid*dpsidzie 
                  hessx(2,id) = hessx(2,id) + dedpsi*dyidxie2
     &                             + d2edpsi2*dpsidyid*dpsidxie
     &                             + d12*dphidyid*dpsidxie 
                  hessy(2,id) = hessy(2,id) + dedpsi*dyidyie2
     &                             + d2edpsi2*dpsidyid*dpsidyie
     &                             + d12*dphidyid*dpsidyie 
                  hessz(2,id) = hessz(2,id) + dedpsi*dyidzie2
     &                             + d2edpsi2*dpsidyid*dpsidzie
     &                             + d12*dphidyid*dpsidzie 
                  hessx(3,id) = hessx(3,id) + dedpsi*dzidxie2
     &                             + d2edpsi2*dpsidzid*dpsidxie
     &                             + d12*dphidzid*dpsidxie 
                  hessy(3,id) = hessy(3,id) + dedpsi*dzidyie2
     &                             + d2edpsi2*dpsidzid*dpsidyie
     &                             + d12*dphidzid*dpsidyie 
                  hessz(3,id) = hessz(3,id) + dedpsi*dzidzie2
     &                             + d2edpsi2*dpsidzid*dpsidzie
     &                             + d12*dphidzid*dpsidzie 
               end if
            end if
         end if
      end do
      return
      end
