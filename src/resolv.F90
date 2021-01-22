      subroutine resolv(c, cold, mdim, eig, nocc) 
      USE vast_kind_param, ONLY:  double 
      use common_arrays_C, only : nfirst, nlast
      use molkst_C, only : norbs, numat
      implicit none
      integer , intent(in) :: mdim 
      integer , intent(in) :: nocc 
      real(double) , intent(inout) :: c(mdim,mdim) 
      real(double) , intent(in) :: cold(mdim,mdim) 
      real(double) , intent(in) :: eig(nocc) 
!
      integer , dimension(4) :: idegen 
      integer :: loop, j, k, nsec, i, ij, ii, jj, l, mo1, mo2, mo3, mo4 
      real(double), dimension(10) :: sec 
      real(double), dimension(16) :: vec 
      real(double), dimension(4) :: eigs 
      real(double) :: thresh, half, suml, sumi, sum, coi, coj, sum1, sum2, sum3, sum4 
!***********************************************************************
!
!   RESOLVE removes any ill-definition in the LMOs.
!
!   If two or more LMOs have the same atomic contributions, then
!   any linear combination of these LMOs is an acceptable solution.
!   This is undesirable, in that the energies of these LMOs is therefore
!   not defined.  RESOLVE will identify such sets of LMOs, and resolve
!   them so that they have a zero energy interaction, that is, the
!   integral <psi(1)|F|psi(2)> is zero.
!
!***********************************************************************
      thresh = 1.D-3 
      half = 0.5D0 - thresh 
      l160: do loop = 1, nocc 
!
!  Is the LMO a potential candidate for degeneracy?
!
        do j = 1, numat 
          if (nlast(j) == nfirst(j)) cycle  
          suml = 0.D0 
          do k = nfirst(j), nlast(j) 
            suml = suml + c(k,loop)**2 
          end do 
!
!  Only LMOs that are 50% or 100% on an atom are potential candidates.
!
          if (suml<=thresh .or. suml>=half) cycle  
          cycle  l160 
        end do 
!
!   LMO 'LOOP' is a candidate.  Now identify any related LMOs
!
        nsec = 1 
        idegen(1) = loop 
       loop1: do i = loop+1,nocc
         !C
         !C   First do a quick check.
         !C
         !            SUMI=0.D0
         !            DO 30 K=NFIRST(L),NLAST(L)
         !   30       SUMI=SUMI+C(K,I)**2
         !            IF(ABS(SUML-SUMI) .GT. THRESH) GOTO 60
         !C
         !C   LMO 'I' is not disqualified.  Now do a complete check.
         !C
         do j = 1,numat
            suml = 0.d0
            sumi = 0.d0
            do k = nfirst(j),nlast(j)
               suml = suml + c(k,loop)**2
               sumi = sumi + c(k,i)**2
            end do
            if (abs(suml-sumi) > thresh) cycle loop1
         end do
         !
         !   LMO 'I' is similar to LMO LOOP
         !
         nsec = nsec + 1
         idegen(nsec) = i
      end do loop1
        if (nsec /= 1) then 
!
!   Build small secular determinant.
!
          ij = 0 
          do ii = 1, nsec 
            i = idegen(ii) 
            do jj = 1, ii 
              j = idegen(jj) 
              sum = 0.D0 
              do l = 1, nocc 
                coi = 0.D0 
                coj = 0.D0 
                do k = 1, norbs 
                  coi = coi + cold(k,l)*c(k,i) 
                  coj = coj + cold(k,l)*c(k,j) 
                end do 
                sum = sum + coi*eig(l)*coj 
              end do 
              ij = ij + 1 
              sec(ij) = sum 
            end do 
          end do 
!
!   Diagonalize, to identify LCMO
!
          call rsp (sec, nsec, eigs, vec) 
          sum = eigs(1) ! dummy use of eigs
!
!    Crude, but fast, way of rotating LMOs
!
          select case (nsec)  
          case default 
            mo1 = idegen(1) 
            mo2 = idegen(2) 
            do i = 1, norbs 
              sum1 = vec(1)*c(i,mo1) + vec(2)*c(i,mo2) 
              sum2 = vec(3)*c(i,mo1) + vec(4)*c(i,mo2) 
              c(i,mo1) = sum1 
              c(i,mo2) = sum2 
            end do 
          case (3)  
            mo1 = idegen(1) 
            mo2 = idegen(2) 
            mo3 = idegen(3) 
            do i = 1, norbs 
              sum1 = vec(1)*c(i,mo1) + vec(2)*c(i,mo2) + vec(3)*c(i,mo3) 
              sum2 = vec(4)*c(i,mo1) + vec(5)*c(i,mo2) + vec(6)*c(i,mo3) 
              sum3 = vec(7)*c(i,mo1) + vec(8)*c(i,mo2) + vec(9)*c(i,mo3) 
              c(i,mo1) = sum1 
              c(i,mo2) = sum2 
              c(i,mo3) = sum3 
            end do 
          case (4)  
            mo1 = idegen(1) 
            mo2 = idegen(2) 
            mo3 = idegen(3) 
            mo4 = idegen(4) 
            do i = 1, norbs 
              sum1 = vec(1)*c(i,mo1) + vec(2)*c(i,mo2) + vec(3)*c(i,mo3) + vec(&
                4)*c(i,mo4) 
              sum2 = vec(5)*c(i,mo1) + vec(6)*c(i,mo2) + vec(7)*c(i,mo3) + vec(&
                8)*c(i,mo4) 
              sum3 = vec(9)*c(i,mo1) + vec(10)*c(i,mo2) + vec(11)*c(i,mo3) + &
                vec(12)*c(i,mo4) 
              sum4 = vec(13)*c(i,mo1) + vec(14)*c(i,mo2) + vec(15)*c(i,mo3) + &
                vec(16)*c(i,mo4) 
              c(i,mo1) = sum1 
              c(i,mo2) = sum2 
              c(i,mo3) = sum3 
              c(i,mo4) = sum4 
            end do 
          end select 
        endif 
      end do l160 
      return  
      end subroutine resolv 
