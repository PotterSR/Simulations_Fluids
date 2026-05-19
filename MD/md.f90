Program MolecularDynamics
    !
    use omp_lib
    use helpers
    use functions
    !
    use parameters, only  : N
    use parameters, only  : eta, L
    use parameters, only  : sigma
    use parameters, only  : n_steps
    use parameters, only  : mxnb

    Implicit none
    !
    !
    !Declaramos la variables a usar
    !
    !Iteradores y tamaño del sistema
    !
    !
    ! integer(kind=i64) :: ii, jj
    real(kind=dp)     :: xx(N), yy(N), zz(N)
    real(kind=dp)     :: vx(N), vy(N), vz(N)
    integer(kind=i64) :: n_neigh(N) 
    integer(kind=i64) :: neigh_list(mxnb, N)

    ! Generamos las posiciones de las partículas
    ! en forma de grid

    write(*, *) "Packing Fraction eta = ", eta

    call init_condition_fcc(xx, yy, zz)

    call md(xx, yy, zz, vx, vz, vy, n_neigh, neigh_list)

    call save_positions(xx, yy, zz)




end program MolecularDynamics