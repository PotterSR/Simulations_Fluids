module parameters
    use helpers
    implicit none 

    !
    ! Definimos el tamaño del sistema
    !
    integer(kind=i64), parameter :: N = 600
    real(kind=dp),     parameter :: sigma = 1.0_dp, L = 10.0_dp
    real(kind=dp),     parameter :: skin = 0.5_dp * sigma    ! = 0.5 para el valor de los vecinos
    integer(kind=i64), parameter :: n_cycle = 100000_i64
    integer(kind=i64), parameter :: mxnb = 50
    real(kind=dp)                :: delta = 0.2_dp
    !
    !
    ! Calculamos el packing fraction del sistema
    ! en 3D
    real(kind=dp),     parameter :: eta = real(N, kind=dp) * (4.0_dp/3.0_dp) * PI *(sigma/2.0_dp)**3 / L**3

end module parameters