module bspqed_function
    implicit none
    private

    public :: function_type

    type, abstract :: function_type
    contains
        procedure(function_evaluation), deferred, pass(self) :: eval
    end type function_type

    abstract interface
        double precision function function_evaluation(x, self)
            import :: function_type
            double precision, intent(in) :: x
            class(function_type), intent(in) :: self
        end function function_evaluation
    end interface

end module bspqed_function
