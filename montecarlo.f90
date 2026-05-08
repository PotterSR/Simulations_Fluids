Program Montecarlo
    !
    use omp_lib
    use helpers
    use functions
    !
    use parameters, only  : N
    use parameters, only  : eta, L
    use parameters, only  : sigma
    use parameters, only  : n_cycle
    use parameters, only  : delta
    use parameters, only  : mxnb

    Implicit none
    !
    !Declaramos la variables a usar
    !
    !Iteradores y tamaño del sistema
    !
    !
    ! integer(kind=i64) :: ii, jj
    real(kind=dp)     :: xx(N), yy(N), zz(N)
    integer(kind=i64) :: n_neigh(N) 
    integer(kind=i64) :: neigh_list(N, mxnb)
    !
    !
    write(*, *) "Packing Fraction eta = ", eta
    !
    !
    !
    ! Generamos las posiciones de las partículas
    ! en forma de grid
    call init_codition(xx, yy, zz)
    ! Corremos el código de Montecarlo para esferas duras
    call montecarlo_hs(xx, yy, zz, n_neigh, neigh_list)
    ! Salvamos las posiciones un .xyz
    call save_positions(xx, yy, zz)


end Program Montecarlo