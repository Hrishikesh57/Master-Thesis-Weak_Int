module weak_angular_labels
    implicit none
    private
    public :: p_two_j, is_p12_state, j_state_factor
contains
    integer function p_two_j(state_name) result(two_j)
        character(len=*), intent(in) :: state_name
        character(len=len(state_name)) :: label
        integer :: i, j, code, ppos
        label = ''; j = 0; two_j = 0
        do i=1,len_trim(state_name)
            code = iachar(state_name(i:i))
            if (code>=iachar('A') .and. code<=iachar('Z')) code=code+32
            if ((code>=iachar('a') .and. code<=iachar('z')) .or. &
                (code>=iachar('0') .and. code<=iachar('9'))) then
                j=j+1; label(j:j)=achar(code)
            end if
        end do
        ppos=index(label,'p')
        if (ppos==0) return
        do i=1,ppos-1
            if (label(i:i)<'0' .or. label(i:i)>'9') return
        end do
        ! Preserve the historical bare np shorthand for np1/2. Never inspect n
        ! when choosing j. Braces, underscores and / in formatted labels are removed.
        select case(trim(label(ppos+1:)))
        case ('', '1', '12')
            two_j=1
        case ('3', '32')
            two_j=3
        end select
    end function p_two_j

    logical function is_p12_state(state_name)
        character(len=*), intent(in) :: state_name
        is_p12_state=p_two_j(state_name)==1
    end function is_p12_state

    double precision function j_state_factor(state_name)
        character(len=*), intent(in) :: state_name
        select case(p_two_j(state_name))
        case(1)
            j_state_factor=1d0
        case(3)
            j_state_factor=sqrt(2d0)
        case default
            error stop 'Unsupported p-state angular label'
        end select
    end function j_state_factor
end module weak_angular_labels
