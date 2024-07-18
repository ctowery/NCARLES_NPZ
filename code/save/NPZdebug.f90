SUBROUTINE NPZdebug(Phyto, Zoo, Nutrients, iz, code_location)
    
    USE pars
    USE fields
    USE con_data
    USE con_stats
    
    INCLUDE 'mpif.h'
    INTEGER :: iz
    CHARACTER(len=*) :: code_location

    REAL :: Phyto, Zoo, Nutrients
    CHARACTER(len=100) :: output_filename

   

    ! Write averages to file
    output_filename = 'debugNPZ.txt'
    open(10, file=output_filename, access='append')
    write(10, '(A, /)') 'Code: ' // trim(code_location)
    write(10, '(A, I10, /)') 'iz: ', iz
    write(10, '(A, F10.6, /)') 'P: ', Phyto
    write(10, '(A, F10.6, /)') 'Z: ', Zoo
    write(10, '(A, F10.6, /)') 'N: ', Nutrients
    write(10, '(A, /)') '------'
    close(10)

    RETURN
END SUBROUTINE NPZdebug
