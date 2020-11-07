! Fortan only uses exclamation points for comments. There are no block comments.
! I am less familiar with Fortran, so I'll just go with the basics

! Let's create a module that contains all info for our integral methods
module Integration_Class
    ! We must include 'implicit none' so that the Fortran compiler doesn't use implicit types
    implicit none
    ! Public means that all methods are accessible
    public

    ! let's create a class (by creating a type in fortran)
    type, public :: Integration
        ! Real is float, kind=8 specifies that it is 8 bytes
        real(kind=8)    ::  low
        real(kind=8)    ::  high
        integer(kind=8) ::  steps

        contains
        ! Our class consists of one method: integrate()
        ! This is linked to a function defined in the module called simps_integrate
            procedure   ::  integrate   =>  simps_integrate
        ! and a function to check if steps is even
            procedure   ::  even        =>  check_if_even
    end type Integration

contains
    ! The function takes argument this (a reference to class Integration) and returns integral_result
    function simps_integrate(this, f) result(int_result)
        class(Integration), intent(in)  ::  this

        ! Create abstract interface so we can point to function
        abstract interface
            function func(z)
                real(kind=8)                :: func
                real(kind=8), intent(in)    :: z
            end function func
        end interface

        ! func is pointer to a function
        procedure(func), pointer, intent(in)    ::  f
        real(kind=8)                            ::  int_result

        real(kind=8)                        ::  dx
        ! Here we create two arrays, type double, length = # of steps
        real(kind=8), dimension(this%steps) ::  x, fs
        integer(kind=8)                     ::  i

        dx = (this%high - this%low) / this%steps
        do i=1, this%steps+1
            x(i) = this%low + ((i-1) * dx)
            fs(i) = f(x(i))
        end do

        do i=1, this%steps+1
            if (i == 1 .or. i == this%steps + 1) then
                int_result = int_result + fs(i)
            else if (modulo(i, 2) == 0) then
                int_result = int_result + (4 * fs(i))
            else
                int_result = int_result + (2 * fs(i))
            end if
        end do

        int_result = int_result * (dx / 3)
    end function

    ! We use subroutine instead of function here because we aren't returning anything
    subroutine check_if_even(this)
        class(Integration), intent(inout)  ::  this

        if (modulo(this%steps, 2) /= 0) then
            this%steps = this%steps + 1
        end if
    end subroutine
end module

module test_functions
    implicit none
contains
    ! Create test functions
    function func1(x)
        real(kind=8)                :: func1
        real(kind=8), intent(in)    :: x
        func1 = x**2
        return
    end function

    function func2(x)
        real(kind=8)                :: func2
        real(kind=8), intent(in)    :: x
        func2 = cosh(x)
        return
    end function
end module


program test_integrations
    use Integration_Class
    use test_functions
    implicit none

    real(kind=8)        ::  result1, result2
    type(Integration)   ::  integral1, integral2

    abstract interface
            function func(x)
                real(kind=8)                :: func
                real(kind=8), intent(in)    :: x
            end function func
    end interface

    procedure(func), pointer :: f_ptr

    ! instantiate and check if steps is even (we must call subroutines)
    integral1 = Integration(0, 5, 100)
    call integral1%even

    integral2 = Integration(0, 10, 500)
    call integral2%even

    ! Point our pointer toward the first function
    f_ptr => func1
    result1 = integral1%integrate(f_ptr)

    ! Now points toward second function
    f_ptr => func2
    result2 = integral2%integrate(f_ptr)

    ! Print results to screen
    write (*, *) "Integral 1: x**2 from 0 to 5: ", result1
    write (*, *) "Integral 2: cosh(x) from 0 to 10: ", result2
end program test_integrations
