      MODULE write_trajectory_I   
        INTERFACE
          subroutine write_trajectory(xyz, mode, charge, escf, ekin, time, xtot)
          use molkst_C, only : numat
          double precision, optional :: escf, ekin, time, xtot
          double precision, intent (in) :: xyz(3,numat)
          double precision, optional :: charge(numat)
          integer, intent (in) :: mode
          end subroutine
        END INTERFACE 
      END MODULE 
