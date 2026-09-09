! -*-f90-*-
module inputparser
  implicit none

contains

  subroutine get_string_parameter_base(fn,parname,par,found)

    implicit none
    character(*),intent(in) :: fn
    character(*),intent(in) ::  parname
    character(*),intent(out) ::  par
    logical,intent(out) :: found

    character*(200) line_string
    integer i,j,l,ll
    integer isokay
    character*(200) temp_string

    open(unit=27,file=fn,status="old",action="read")

10  continue
    read(27,'(a)',end=19) line_string
    ! separator is an equal sign '=', # is comment
    i = index(line_string,'=')
    j = index(line_string,'#')

    if (i.eq.0.or.j.eq.1) goto 10
    !   if the whole line is a comment or there is no
    !   equal sign, then go on to the next line    

    if(j.gt.0.and.j.lt.i) goto 10
    !   if there is an equal sign, but it is in a comment
    !   then go on to the next line


    ! is this the right parameter? If not, cycle
    temp_string=trim(adjustl(line_string(1:i-1)))
    l=len(parname)
    if(parname.ne.temp_string(1:l)) goto 10

    !  If there is a comment in the line, exclude it!
    l = len(line_string)
    if (j.gt.0) l = j - 1

    par = line_string(i+1:l)
    ! now remove potential crap!
    do ll=1,len(par)
       if(par(ll:ll).eq.'\t') par(ll:ll) = ' '
       if(par(ll:ll).eq.'"') par(ll:ll) = ' '
       if(par(ll:ll).eq."'") par(ll:ll) = ' '
    enddo
    do ll=1,len(par) 
       if(isokay(par(ll:ll)).ne.1) then
          par(ll:ll) = ' '
       endif
    enddo

    ! adjust left...
    ! ll = len(par)
    ! do while (par(1:1).eq.' '.or.par(1:1).eq.'\t') 
    !    par(1:ll-1) = par(2:ll)
    ! end do
    ! get rid of trailing blanks
    ! j = index(par," ")
    ! par = par(1:j-1)

    
    par = adjustl(par)
    par = trim(par)

    ! now look for " or ' and remove them
    j=index(par,'"')
    if(j.ne.0) stop "No quotes in my strings, please!"

    j=index(par,"'")
    if(j.ne.0) stop "No quotes in my strings, please!"

    found = .true.

    close(27)
    return


19  continue

    found = .false.
    par = ""
    close(27)
    return

  end subroutine get_string_parameter_base

  
  subroutine get_string_parameter(fn,parname,par,par_default)

    implicit none
    character(*),intent(in) :: fn
    character(*),intent(in) ::  parname
    character(*),intent(out) ::  par
    character(*),intent(in),optional :: par_default

    character*256 line_string
    logical :: found
    integer :: ios

    call get_string_parameter_base(fn,parname,line_string,found)

    if(found)then
       !read(line_string,*,iostat=ios) par
       par = trim(line_string)
       ! if(ios /= 0) then
       !    write(6,*) "Uh. Bad string parameter ",trim(parname)
       !    write(6,*) "Please check input file!"
       !    call flush(6)
       !    stop
       ! endif
    else
       if(present(par_default))then
          par = par_default
       else
          write(6,*) "could not be read, and no default given!"
          call flush(6)
          stop
       endif
       
    endif

  end subroutine get_string_parameter

  subroutine get_double_parameter(fn,parname,par,par_default)

    implicit none
    character(*),intent(in) :: fn

    character(*),intent(in) :: parname
    character*256 :: line_string
    real*8,intent(out) :: par
    real*8,intent(in),optional :: par_default
    logical :: found
    integer :: ios

    call get_string_parameter_base(fn,parname,line_string,found)

    if(found)then
       read(line_string,*,iostat=ios) par

       if(ios /= 0) then
          write(6,*) "Uh. Bad double parameter ",trim(parname)
          write(6,*) "Please check input file!"
          call flush(6)
          stop
       endif
    else
       if(present(par_default))then
          par = par_default
       else
          write(6,*) "could not be read, and no default given!"
          call flush(6)
          stop
       endif
    endif

  end subroutine get_double_parameter

  subroutine get_integer_parameter(fn,parname,par,par_default)

    implicit none
    character(*),intent(in) :: fn

    character(*),intent(in) :: parname
    
    character*256 :: line_string
    integer,intent(out) :: par
    integer,intent(in),optional :: par_default
    logical :: found
    integer :: ios

    call get_string_parameter_base(fn,parname,line_string,found)

    if(found)then
       read(line_string,*,iostat=ios) par

       if(ios /= 0) then
          write(6,*) "Uh. Bad integer parameter ",trim(parname)
          write(6,*) "Please check input file!"
          call flush(6)
          stop
       endif
    else
       if(present(par_default))then
          par = par_default
       else
          write(6,*) "could not be read, and no default given!"
          call flush(6)
          stop
       endif
    endif

  end subroutine get_integer_parameter

  subroutine get_logical_parameter(fn,parname,par,par_default)

    implicit none
    character(*),intent(in) :: fn

    character(*),intent(in) :: parname
    
    character*256 :: line_string
    logical,intent(out) :: par
    logical,intent(in),optional :: par_default
    logical :: found
    integer :: ios

    call get_string_parameter_base(fn,parname,line_string,found)
    
    if(found)then
       read(line_string,*,iostat=ios) par

       if(ios /= 0) then
          write(6,*) "Uh. Bad logical parameter ",trim(parname)
          write(6,*) "Please check input file!"
          call flush(6)
          stop
       endif
    else
       if(present(par_default))then
          par = par_default
       else
          write(6,*) "could not be read, and no default given!"
          call flush(6)
          stop
       endif
    endif
    
  end subroutine get_logical_parameter

end module inputparser

function isokay(checkme)
  implicit none
  integer isokay
  character(1) checkme
  character(len=67) okaychars
  integer i
  isokay = 0
  okaychars = "1234567890.-+/_abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
  do i=1,len(okaychars)
     if(checkme(1:1) .eq. okaychars(i:i)) then
        isokay = 1
        exit
     endif
  enddo

end function isokay
