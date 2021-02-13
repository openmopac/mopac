      MODULE dmecip_I   
      INTERFACE
      subroutine dmecip (coeffs, deltap, delta, eig, vectci, conf)
      use molkst_C, only: norbs
      use meci_C, only : lab, nmos
      implicit none
      double precision, dimension (lab), intent (inout) :: eig, vectci
      double precision, dimension (lab**2), intent (in) :: conf
      double precision, dimension (norbs, norbs), intent (in) :: coeffs
      double precision, dimension (norbs, nmos), intent (inout) :: delta
      double precision, dimension (nmos, nmos), intent (inout) :: deltap
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
