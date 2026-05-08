module functions

    use helpers
    use parameters

contains 

    subroutine init_codition(xx, yy, zz)
        !
        ! Nos genera una condición inicial en forma de grid
        !
        real(kind=dp), intent(inout) :: xx(N), yy(N), zz(N)
        real(kind=dp)              :: a
        integer(kind=i64)          :: ii, jj, ll, idx, n_side

        n_side = ceiling( real(N, kind=dp)**(1.0_dp/3.0_dp) )
        a      = L / real(n_side, kind=dp)
        idx = 1
        do jj = 0, n_side - 1
            do ii = 0, n_side - 1
                do ll = 0, n_side -1
                    if (idx > N) exit
                    xx(idx) = (real(ii, kind=dp) + 0.5_dp) * a
                    yy(idx) = (real(jj, kind=dp) + 0.5_dp) * a
                    zz(idx) = (real(ll, kind=dp) + 0.5_dp) * a
                    idx = idx + 1
                end do
                if (idx > N) exit
            end do
            if (idx > N) exit
        end do

    end subroutine init_codition

    subroutine save_positions(xx, yy, zz)
        real(kind=dp), intent(in) :: xx(N), yy(N), zz(N)
        integer(kind=i64)          :: ii

        open(unit=10, file="positions.xyz", status="replace", action="write")
        write(10, '(I0)') N
        write(10, '(A, F0.6, A, F0.6, A, F0.6, A)') &
        'Lattice="', L, ' 0.0 0.0 0.0 ', L, ' 0.0 0.0 0.0 ', L, '"'

        do ii = 1, N
            write(10, *) "A", xx(ii), yy(ii), zz(ii)
        end do
        close(10)

    end subroutine save_positions

    subroutine montecarlo_hs(xx, yy, zz)
    !***********************************************
    ! Implementación de Montecarlo para esferas duras
    !   
    ! Se rechaza el movimiento si algún traslape con
    ! alguna otra partícula
    !
    !***********************************************
        real(kind=dp), intent(inout) :: xx(N), yy(N), zz(N)
        integer(kind=i64) :: o, n_accept, n_total, ii, jj
        logical           :: overlap
        real(kind=dp)     :: x_trial, y_trial, z_trial
        real(kind=dp)     :: r2, dx, dy, dz

        overlap = .false.

        n_accept = 0
        n_total  = 0
        write(*, *) delta
        do jj = 1, n_cycle
            !
            !
            overlap = .false.
            ! Tomamos un un valor aletorio para las particulas
            ! 
            o = int( urand(0.0_dp, real(N, kind=dp)) , kind=i64) + 1
            !
            !
            ! Generamos el desplazamiento para ver si lo aceptamos 
            x_trial = xx(o) + delta*(urand(-0.5_dp, 0.5_dp))
            y_trial = yy(o) + delta*(urand(-0.5_dp, 0.5_dp))
            z_trial = zz(o) + delta*(urand(-0.5_dp, 0.5_dp))
            !
            ! Periodic Boundary Conditions
            x_trial = x_trial - L*floor(x_trial/L)
            y_trial = y_trial - L*floor(y_trial/L)
            z_trial = z_trial - L*floor(z_trial/L)
            ! Calculamos la energía de la partícula

            do ii = 1, N 
                if(ii == o) cycle 

                dx = x_trial - xx(ii)
                dy = y_trial - yy(ii)
                dz = z_trial - zz(ii)
                !
                ! Condiciones de imagen mínima
                !
                dx = dx - L*nint(dx/L)
                dy = dy - L*nint(dy/L)
                dz = dz - L*nint(dz/L)

                r2 = dx*dx + dy*dy  + dz*dz           
                if(r2 < (sigma)**2) then
                    overlap = .true.
                    exit
                end if
            end do 

            if(.not. overlap) then
                xx(o) = x_trial
                yy(o) = y_trial
                zz(o) = z_trial
                n_accept = n_accept + 1
            end if

            !if (mod(jj, 10000) == 0) then
            !    if (real(n_accept, kind=dp) / 1000.0_dp > 0.2_dp) then
            !        delta = delta * 1.05_dp  
            !    else
            !        delta = delta * 0.95_dp    
            !    end if
            !    n_total  = n_total + n_accept
            !    n_accept = 0
            !end if
        end do

        write(*, *) "Montecarlo Done"
        write(*, *) "*******************************"
        write(*, *) "Ratio de Aceptados: "
        write(*, *) real(n_accept, kind=dp) / real(n_cycle, kind=dp)
        write(*, *) delta
        write(*, *) "*******************************"
    end subroutine montecarlo_hs

    subroutine build_neighbor_list(xx, yy, zz, n_neigh, neigh_list)
        real(kind=dp), intent(in) :: xx(N), yy(N), zz(N)
        integer(kind=i64),intent(out) :: n_neigh(N)
        integer(kind=i64),intent(out) :: neigh_list(N, mxnb)

        integer(kind=i64) :: ii, jj
        real(kind=dp)     :: rc, rc2
        real(kind=dp)     :: dx, dy, dz, r2
        !
        ! Damos los valores para el rc de la lista
        rc  = sigma + skin
        rc2 = rc**2
        !
        ! Inicializamos los vectores 
        neigh_list(:,:) = 0_i64
        n_neigh(:)      = 0_i64
        !
        do ii = 1, N 
            do jj = ii+1, N 
                dx = xx(ii) - xx(jj)
                dy = yy(ii) - yy(jj)
                dz = zz(ii) - zz(jj)
                dx = dx - L*nint(dx/L)
                dy = dy - L*nint(dy/L)
                dz = dz - L*nint(dz/L)
                r2 = dx*dx + dy*dy + dz*dz

                if(r2 < rc2) then
                    ! Cuantos vecinos tiene
                    n_neigh(ii) =  n_neigh(ii) + 1
                    n_neigh(jj) =  n_neigh(jj) + 1
                    ! Estamos checando si la lista no excede el límite que le dimos
                    if (n_neigh(ii) > MXNB .or. n_neigh(jj) > MXNB) then
                        write(0,*) "ERROR: MXNB too small, increase it in params.f90"
                        stop
                    end if
                    ! Añadimos que partícula es el vecino
                    neigh_list(ii, n_neigh(ii)) = jj
                    neigh_list(jj, n_neigh(jj)) = ii
                end if 
            end do
        end do
    end subroutine build_neighbor_list

end module functions