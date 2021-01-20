double precision function parips (eigs, loop)
    use molkst_C, only : norbs, nopen, nclose, nalpha
    use param_global_C, only : lions
    implicit none
    integer, intent(in) :: loop
    double precision, dimension (norbs), intent (in) :: eigs
    if(nclose /= 0) then
!
!  Only allow for higher IPs for closed-shell systems
!
      parips = -eigs(lions(loop))
      if (nopen /= nclose) then
           parips = Min (parips,-eigs(nopen))
      end if
    elseif (nalpha > 0) then
      parips = -eigs(nalpha)
    else
      parips = 0.d0
    end if
end function parips
