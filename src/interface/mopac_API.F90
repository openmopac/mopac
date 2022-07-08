module mopac_api
  implicit none

  ! unit number for dummy mopac input file; used only to initialize mopac
  integer :: MOPAC_IF = 1001

  ! keywords specified by API
  character :: keywords*3000

  contains
  
    subroutine mopac_initialize(Natom, XYZ, labels)
      use mopac_interface_flags, only : use_api, reset_mopac_initialize
      implicit none
      ! variables passed to API
      integer :: Natom
      double precision :: XYZ(3, Natom)
      character(3) :: labels(Natom)
      character :: tmpkeywords*3000

      ! local variables
      integer :: i,j

      ! perform checks on the supplied keywords to ensure calculation is compatible with API
      tmpkeywords=trim(keywords)
      call upcase(tmpkeywords, len(trim(tmpkeywords)))
      if (index(tmpkeywords,'1SCF') .eq. 0) then
         write(*,"('API CURRENTLY ONLY SUPPORTS 1SCF CALCULATIONS')")
         stop
      end if

      if (index(tmpkeywords,'1SCF') .eq. 0) then
         write(*,"('API CURRENTLY REQUIRES THE XYZ KEYWORD; IT HAS BEEN ADDED')")
         keywords=trim(keywords)//" XYZ"
      end if
      
      if (index(tmpkeywords,'GEO-OK') .eq. 0) then
         write(*,"('WARNING:')")
         write(*,"('IT IS RECOMMENDED TO RUN THE API WITH THE GEO-OK KEYWORD')")
      end if

      if (index(tmpkeywords,'NOSYM') .eq. 0) then
         write(*,"('WARNING:')")
         write(*,"('IF THE SYMMETRY WILL CHANGE DURING YOUR CALCULATIONS, RUN WITH')")
         write(*,"('KEYWORD NOSYM')")
      end if
      
      ! make sure that mopac knows that the API is being used.
      use_api = .true.

      ! check if a previous system have been initialized through the API
      if (reset_mopac_initialize) then
         ! mopac has been called already, clean up arrays
         call setup_mopac_arrays(0,0)
      else
         ! make sure the API knows that it has previously been initialized
         reset_mopac_initialize = .true.
      end if

      ! create mopac input file for initialization
      open(MOPAC_IF, file='MOPAC input', status='unknown')
      write(MOPAC_IF, "(A)") trim(keywords)
      write(MOPAC_IF, "(A)") "Dummy MOPAC input file for API"
      write(MOPAC_IF, *)
      do i = 1, Natom
         write(MOPAC_IF, "(XA10,f18.8, ' 1', f18.8, ' 1', f18.8, ' 1')") labels(i), (XYZ(j,i),j=1,3)
      end do
      write(MOPAC_IF, *)
      close(MOPAC_IF)

      ! let default mopac routine setup storage, etc for calculations
      call run_mopac

      return

    end subroutine mopac_initialize

    ! driver routine to obtain the SCF energy and gradient
    subroutine mopac_return_energy_and_grad(Natom, XYZ, energy, dXYZ, reset_density, get_grad)
      use common_arrays_C, only : grad, xparam
      use molkst_C, only : escf, moperr, nvar
      use mopac_interface_flags, only : reset_iter_L, SCF_avail
      implicit none
      ! variables passed to API
      integer, intent(in) :: Natom
      double precision, intent(in) :: XYZ(3, Natom)
      double precision, intent(out) :: energy, dXYZ(3, Natom)
      logical, intent(in), optional :: reset_density, get_grad

      ! local variables
      integer :: i, j, k
      logical :: calc_grad
      
      ! map API passed geometry into xparam
      k = 0
      do i = 1,Natom
         do j = 1, 3
            k = k + 1
            xparam(k) = XYZ(j,i)
         end do
      end do

      ! zero gradient
      grad(:nvar) = 0.d0

      ! determine if gradient will be calculated
      if (present(get_grad)) then
         if (get_grad) then
            calc_grad = .true.
         else
            calc_grad = .false.
         end if
      else
         calc_grad = .true.
      end if
      
      ! optionally reset density
      if (present(reset_density)) then
         if (reset_density) reset_iter_L = .true.
      endif
      
      ! obtain new SCF
      call compfg(xparam, .true., escf, .true., grad, calc_grad)
      if (moperr) then
         write(*,"('MOPAC FAILED. STOPPING')")
         SCF_avail = .false.
         stop
      else
         SCF_avail = .true.
      end if

      ! map mopac variables into API variables
      energy = escf
      if (calc_grad) then
         do i = 1, Natom
            k = 3*(i-1)
            do j = 1, 3
               dXYZ(j,i)=grad(k+j)
            end do
         end do
      end if

      return
    end subroutine mopac_return_energy_and_grad

    ! API routine to obtain just the SCF energy
    subroutine mopac_return_energy(Natom, XYZ, energy, reset_density)
      use common_arrays_C, only : grad, xparam
      use molkst_C, only : escf, moperr, nvar
      use mopac_interface_flags, only : reset_iter_L, SCF_avail
      implicit none
      ! variables passed to API
      integer, intent(in) :: Natom
      double precision, intent(in) :: XYZ(3, Natom)
      logical, intent(in), optional :: reset_density
      double precision, intent(out) :: energy

      ! local variables
      integer :: i, j, k
      double precision :: dXYZ(3, Natom)

      if (present(reset_density)) then
         call mopac_return_energy_and_grad(Natom, XYZ, energy, dXYZ, &
              reset_density = reset_density, get_grad = .false.)
      else
         call mopac_return_energy_and_grad(Natom, XYZ, energy, dXYZ, &
              get_grad = .false.)
      end if

      return
    end subroutine mopac_return_energy

    ! driver routine to obtain atomic partial charges
    subroutine mopac_return_charges(Natom, charges)
      use common_arrays_C, only: p, labels
      use parameters_C, only: tore
      use mopac_interface_flags, only : SCF_avail
      implicit none
      integer, intent(in) :: Natom
      double precision, intent(out) :: charges(Natom)
      integer :: i, l

      ! Make sure that an SCF is available
      if (.not. SCF_avail) then
         write(*,*) 'Obtain an SCF first'
         stop
      endif
      ! Obtain the charges
      call CHRGE(P, Charges)
      ! subtract off the number of protons
      do i=1,Natom
         l = labels(i)
         Charges(i) = tore(l) - Charges(i)
      end do

    end subroutine mopac_return_charges

    ! driver routine to obtain bond order matrix
    subroutine mopac_return_bonds(Natom, BO)
      use common_arrays_C, only : p, labels, bondab
      use mopac_interface_flags, only : SCF_avail
      implicit none
      integer, intent(in) :: Natom
      double precision, intent(out) :: BO(Natom, Natom)
      integer i, j, k

      if (.not. SCF_avail) then
         write(*,*) 'Obtain an SCF first'
         stop
      endif
      ! get bond order matrix
      call BONDS(P)
      k = 0
      do i= 1,Natom
         do j=1, i
            k = k + 1
            BO(i,j) = bondab(k)
            BO(j,i) = BO(i,j)
         end do
      end do
      return
    end subroutine mopac_return_bonds
  end module mopac_api
