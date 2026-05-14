! **************************************************** !
! ******************* helpers.f90 ******************** !
! **************************************************** !
!
!    Leonardo Dagdug, Ivan Pompa-García, Jason Peña
!
! From the book:
! **************************************************** !
!      Diffusion Under Confinement:
!               A Journey Through Counterintuition    
! **************************************************** !

! **************************************************** !
! This file is used to define a general module
! that stores definitions, constants, functions,
! and procedures that may be helpful in the study
! of diffusion.
! **************************************************** !

module helpers

! From the intrinsics of Fortran, we import certain
! constants to make our data types portable between 
! systems.
! The «iso_fortran_env» intrinsic is available for 
! Fortran >= 2003, but the needed data type constants
! were introduced in Fortran >= 2008
!
!   sp -> Single Precision
!   dp -> Double Precision
!   qp -> Quadruple Precision
!
! There are many more functions, procedures, and
! definitions inside «iso_fortran_env». Using the
! «only» keyword, we only import the ones we need.
use, intrinsic :: iso_fortran_env, &
      only: sp => real32, dp => real64, qp => real128, &
          i8 => int8, i16 => int16, i32 => int32, &
          i64 => int64

! Mandatory declaration of data type variables
! and constants.
implicit none

! Define a generic name for the number to string
! converting functions. This is how we avoid having to
! select the specific function; the compiler will do it
! for us.
interface nstr
  module procedure i8str, i16str, i32str, i64str, &
                  spstr, dpstr, qpstr
end interface

! **************************************************** !
!               Declaration of constants
! **************************************************** !

! pi/4
real(kind=dp), parameter :: PI4 = atan(1.0_dp)
! pi/2
real(kind=dp), parameter :: PI2 = 2.0_dp * PI4
! pi
real(kind=dp), parameter :: PI = 4.0_dp * PI4

! Tabulator
character, parameter :: TAB = char(9)
! New line
character, parameter :: TNL = new_line('l')

! Default decimal digits when converting to char
integer, parameter :: DDECN = 6


! **************************************************** !
!               Functions and Subroutines
! **************************************************** !
contains

! To make «repeatable» simulations, we need to set a
! seed before generating our PRNs. Fortran requires
! an array of length «n», which length depends
! directly on the system and the compiler. But setting
! the seed helps us to make sure that we can obtain
! repeatable results, at least in a machine with the
! same version of the compiler.
! This function sets a seed in a very basic manner
! by fixing all the numbers in the seed array to
! the same integer.
subroutine setseed(nseed)
  implicit none

  ! Input parameter
  integer, intent(in) :: nseed
  ! Variable to hold the size of the needed array
  integer :: n
  ! Array to set the seed
  integer, allocatable :: seed(:)

  ! Read the size of the needed array
  call random_seed(size=n)

  ! Allocate the array
  allocate( seed(n) )
  ! Set all the elements in the array to «nseed»
  seed = nseed

  ! Set the seed
  call random_seed(put=seed)
end subroutine setseed

! Fortran contains random_number procedure to
! generate uniformly distributed pseudorandom numbers
! within the interval 0 <= x < 1. Occasionally it is
! also advisable to use it as a function instead.
! Also, returning 1.0_dp - x gives us the possibility
! to obtain 1.0 as a value, and the uniform distribution
! will hold. This is the interval 0 < x <= 1.
function usrand() result(fn_res)
  implicit none

  ! Output value
  real(kind=dp) :: fn_res

  ! This gives us: 0 <= x < 1
  call random_number(fn_res)

  ! Now, the interval is: 0 < x <= 1
  fn_res = 1.0_dp - fn_res
end function usrand

! Sometimes we will need to get a pseudorandom number
! within an interval (a,b] other than (0,1]. So, we made
! a simple transformation.
! We are asumming that a < b.
function urand(a, b) result(fn_res)
  implicit none

  ! Input parameters
  real(kind=dp), intent(in) :: a, b
  ! Output value
  real(kind=dp) :: fn_res

  fn_res = a + usrand() * (b - a)
end function urand

! For the specific case of Brownian motion, the steps 
! are normally distributed with parameters
!   mu = 0 (mean)
!   sigma = \sqrt{2 D dt} (standard deviation)
! By default we want to generate normally distributed
! pseudorandom numbers with
!   mu = 0 (mean)
!   sigma = 1 (standard deviation)
!
! The transformation from U~(0,1) to N(0,1) is made via
! the Box-Müller transformation.
!
! Original Box & Müller paper:
!       https://doi.org/10.1214%2Faoms%2F1177706645
! 
function nrand(imu, isigma) result (fn_res)
  implicit none

  ! Input parameters
  real(kind=dp), intent(in), optional :: imu, isigma
  ! Output value
  real(kind=dp) :: fn_res
  ! Parameters for calculation
  real(kind=dp) :: mu, sigma
  ! Variables to hold normally distributed PRNs
  real(kind=dp) :: u1, u2

  ! Variable to save the non-used PRN
  real(kind=dp), save :: saved_val
  ! Flag to mark if we have a previously
  ! calculated value
  logical, save :: saved = .false.

  ! If there is a previously generated value ready
  ! to be used...
  if(saved) then
    ! Set the return variable to the saved value
    fn_res = saved_val
    ! We must rememeber that the value was used. Then
    ! in the next function call, all calculations
    ! must be performed again.
    saved  = .false.

    ! Exit from this function, no calculation needed.
    return
  end if

  ! Check for parameter values. If some of them 
  ! are missing, then the default values will be used.
  if( present(imu) ) then
    mu = imu
  else
    mu = 0.0_dp
  end if

  if( present(isigma) ) then
    sigma = isigma
  else
    sigma = 1.0_dp
  end if

  ! Generates two numbers uniformly distributed within
  ! the interval (0,1), required by the
  ! Box-Müller transformation.
  u1 = usrand()
  u2 = usrand()

  ! The Box-Müller transformation.
  ! Both numbers, u1 and u2 are used. If we change the
  ! cosine <-> sine, we can generate another normally 
  ! distributed number independent from the first one.
  ! The election of one of the trigonometric functions
  ! to return one PRN is completely arbitrary.
  ! The second generated value is held and used
  ! in the next function call to avoid a new
  ! calculation.

  ! One of the generated numbers is saved to be used
  ! in the next function call without calculations.
  saved_val = sigma * sqrt( -2.0_dp * log(u1) ) &
                          * sin( 2.0_dp * PI * u2 ) + mu
  ! Flag to inform that one value is saved.  
  saved = .true.

  ! The number to be returned in this function call.
  fn_res = sigma * sqrt( -2.0_dp * log(u1) ) &
                          * cos( 2.0_dp * PI * u2 ) + mu
end function nrand

! This family of functions i*str(num) makes it possible
! to print out an integer concatenated with a string,
! avoiding the odd space characters introduced by
! Fortran.
function i8str(num) result(fn_res)
    implicit none

    ! Input integer value
    integer(kind=i8), intent(in) :: num
    ! Temporary string of 'sufficient' length
    character(len=256) :: stmp
    ! Dynamic length variable to hold the final string
    character(len=:), allocatable :: fn_res

    ! Copy the integer into the temporary character 
    ! variable.
    ! This copies the «num» integer value into the
    ! «stmp» character variable using the I0 format,
    ! which means that there will be no leading or
    ! padding zeros. This is their most compact form.
    write(stmp, '(I0)') num
    ! Move the string to the left
    stmp = adjustl(stmp)

    ! Cut the moved string from the start to the 
    ! spaces at the right of it.
    fn_res = stmp( 1:len_trim(stmp) )
end function i8str

function i16str(num) result(fn_res)
  implicit none

  ! Input integer value
  integer(kind=i16), intent(in) :: num
  ! Temporary string of 'sufficient' length
  character(len=256) :: stmp
  ! Dynamic length variable to hold the final string
  character(len=:), allocatable :: fn_res

  ! Copy the integer into the temporary character 
  ! variable
  write(stmp, '(I0)') num
  ! Move the string to the left
  stmp = adjustl(stmp)

  ! Cut the moved string from the start to the 
  ! spaces at the right of it.
  fn_res = stmp( 1:len_trim(stmp) )
end function i16str

function i32str(num) result(fn_res)
  implicit none

  ! Input integer value
  integer(kind=i32), intent(in) :: num
  ! Temporary string of 'sufficient' length
  character(len=256) :: stmp
  ! Dynamic length variable to hold the final string
  character(len=:), allocatable :: fn_res

  ! Copy the integer into the temporary character 
  ! variable
  write(stmp, '(I0)') num
  ! Move the string to the left
  stmp = adjustl(stmp)

  ! Cut the moved string from the start to the 
  ! spaces at the right of it.
  fn_res = stmp( 1:len_trim(stmp) )
end function i32str

function i64str(num) result(fn_res)
  implicit none

  ! Input integer value
  integer(kind=i64), intent(in) :: num
  ! Temporary string of 'sufficient' length
  character(len=256) :: stmp
  ! Dynamic length variable to hold the final string
  character(len=:), allocatable :: fn_res

  ! Copy the integer into the temporary character 
  ! variable
  write(stmp, '(I0)') num
  ! Move the string to the left
  stmp = adjustl(stmp)

  ! Cut the moved string from the start to the 
  ! spaces at the right of it.
  fn_res = stmp( 1:len_trim(stmp) )
end function i64str

! This family of functions *pstr(num) makes it possible
! to print out a real number concatenated with a string,
! avoiding the odd space characters introduced by
! Fortran.
function spstr(num, idigits) result(fn_res)
  implicit none

  ! Input real single-precision value
  real(kind=sp), intent(in) :: num
  ! Digits after decimal point to be printed
  integer, intent(in), optional :: idigits
  ! The actual value used
  integer :: digits
  ! Temporary string of 'sufficient' length
  character(len=256) :: stmp
  ! Dynamic length variable to hold the final string
  character(len=:), allocatable :: fn_res

  ! Check for parameter values. If some of them 
  ! are missing, then the default values will be used.
  if( present(idigits) ) then
    digits = idigits
  else
    digits = DDECN
  end if

  ! Copy the real_sp into the temporary character 
  ! variable
  write(stmp, '(F100.' // nstr(digits) // ')') num
  ! Move the string to the left
  stmp = adjustl(stmp)

  ! Cut the moved string from the start to the 
  ! spaces to the right of it.
  fn_res = stmp( 1:len_trim(stmp) )
end function spstr

function dpstr(num, idigits) result(fn_res)
  implicit none

  ! Input real double-precision value
  real(kind=dp), intent(in) :: num
  ! Digits after decimal point to be printed
  integer, intent(in), optional :: idigits
  ! The actual used value
  integer :: digits
  ! Temporary string of 'sufficient' length
  character(len=256) :: stmp
  ! Dynamic length variable to hold the final string
  character(len=:), allocatable :: fn_res

  ! Check for parameter values. If some of them 
  ! are missing, then the default values will be used.
  if( present(idigits) ) then
    digits = idigits
  else
    digits = DDECN
  end if

  ! Copy the real_sp into the temporary character 
  ! variable
  write(stmp, '(F100.' // nstr(digits) // ')') num
  ! Move the string to the left
  stmp = adjustl(stmp)

  ! Cut the moved string from the start to the 
  ! spaces at the right of it.
  fn_res = stmp( 1:len_trim(stmp) )
end function dpstr

function qpstr(num, idigits) result(fn_res)
  implicit none

  ! Input real quadruple-precision value
  real(kind=qp), intent(in) :: num
  ! Digits after decimal point to be printed
  integer, intent(in), optional :: idigits
  ! The actual value used
  integer :: digits
  ! Temporary string of 'sufficient' length
  character(len=256) :: stmp
  ! Dynamic length variable to hold the final string
  character(len=:), allocatable :: fn_res

  ! Check for parameter values. If some of them 
  ! are missing, then the default values will be used.
  if( present(idigits) ) then
    digits = idigits
  else
    digits = DDECN
  end if

  ! Copy the real_sp into the temporary character 
  ! variable
  write(stmp, '(F100.' // nstr(digits) // ')') num
  ! Move the string to the left
  stmp = adjustl(stmp)

  ! Cut the moved string from the start to the 
  ! spaces at the right of it.
  fn_res = stmp( 1:len_trim(stmp) )
end function qpstr

end module helpers