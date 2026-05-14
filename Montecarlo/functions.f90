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

    subroutine build_neighbor_list(xx, yy, zz, n_neigh, neigh_list)
        real(kind=dp), intent(in) :: xx(N), yy(N), zz(N)
        integer(kind=i64),intent(out) :: n_neigh(N)
        integer(kind=i64),intent(out) :: neigh_list(mxnb, N)

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
        do ii = 1, N-1 ! ii = 1 -> jj =2 ... ii = N-1, jj = N
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
                    neigh_list(n_neigh(ii), ii) = jj
                    neigh_list(n_neigh(jj), jj) = ii
                end if 
            end do
        end do
    end subroutine build_neighbor_list

    subroutine radial_distribution(xx, yy, zz, rdf)
        real(kind=dp), intent(in)    :: xx(N), yy(N), zz(N)
        real(kind=dp), intent(inout)   :: rdf(nbins)
        integer(kind=i64) :: ii, jj, bin
        real(kind=dp)     :: dbin, r2, dx, dy, dz, L2, r
        real(kind=dp)     :: volr, nid, r_in, r_out

        L2 = (L/2.0_dp)**2

        dbin = L/(2.0_dp * real(nbins, kind=dp)) ! Tamaño de los bins
        rdf(:) = 0.0_dp ! Inicializamos el rdf

        do ii = 1, N-1
            do jj = ii+1, N 
                dx = xx(ii) - xx(jj)
                dy = yy(ii) - yy(jj)
                dz = zz(ii) - zz(jj)

                dx = dx - L*nint(dx/L)
                dy = dy - L*nint(dy/L)
                dz = dz - L*nint(dz/L)
                
                r2 = dx**2 + dy**2 + dz**2
                ! Solamente aceptamos si esta a menos de L/2
                if(r2 < L2) then
                    r = sqrt(r2)
                    bin = int(r/dbin) + 1_i64
                    rdf(bin) = rdf(bin) + 2.0_dp
                end if
            end do 
        end do

        ! Normalización
        do ii = 1, nbins
            r_in  = real(ii-1, kind=dp) * dbin
            r_out = real(ii, kind=dp) * dbin 
            volr  = (4.0_dp/3.0_dp) * PI * (r_out**3 - r_in**3)
            nid   = (real(N, kind=dp) / L**3) * volr
            
            rdf(ii) = rdf(ii) / (nid * real(N, kind=dp))
        end do

        !open(unit=10, file="rdf.dat", status="replace", action="write")
        !
!
        !do ii = 1, nbins
        !    r = dbin*(real(ii, kind=dp) + 0.5_dp)
        !    write(10, *) r, rdf(ii)
        !end do 
!
        !close(10)
        

    end subroutine radial_distribution



    subroutine montecarlo_hs(xx, yy, zz, n_neigh, neigh_list)
    !***********************************************
    ! Implementación de Montecarlo para esferas duras
    !   
    ! Se rechaza el movimiento si hay algún traslape con
    ! alguna otra partícula
    !
    !***********************************************
        real(kind=dp), intent(inout) :: xx(N), yy(N), zz(N)
        integer(kind=i64), intent(inout)    :: n_neigh(N), neigh_list(mxnb, N)
        integer(kind=i64) :: o, n_accept, n_total, ii, jj, kk, n_frame
        logical           :: overlap
        real(kind=dp)     :: x_trial, y_trial, z_trial
        real(kind=dp)     :: xxo(N), yyo(N), zzo(N)
        real(kind=dp)     :: r2, dx, dy, dz, max_disp, r
        real(kind=dp)     :: rdf_frame(nbins), rdf(nbins)
        real(kind=dp)     :: dbin
        real(kind=dp)     :: max_disp2, d2
        character(len=50) :: filename
        !
        !
        !
        ! Inicializamos las variables
        overlap = .false.
        n_accept = 0
        n_total  = 0
        max_disp = 0.0_dp
        !
        !
        ! Inicializmos la variables antes del Montecarlo
        !
        ! Inicializamos las la lista de vecinos
        call build_neighbor_list(xx, yy, zz, n_neigh, neigh_list)
        xxo = xx
        yyo = yy 
        zzo = zz
        max_disp2 = 0.02_dp
        !
        ! Termalización del sistema
        !
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

            do ii = 1, n_neigh(o)
                
                kk = neigh_list(ii, o)

                dx = x_trial - xx(kk)
                dy = y_trial - yy(kk)
                dz = z_trial - zz(kk)
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
                d2 = (xx(o)-xxo(o))**2 + (yy(o)-yyo(o))**2 + (zz(o)-zzo(o))**2
                max_disp2 = max(max_disp2, d2)
            end if

            if (max_disp2 > (sigma+skin)**2) then
                call build_neighbor_list(xx, yy, zz, n_neigh, neigh_list)
                xxo = xx
                yyo = yy
                zzo = zz
                max_disp2 = 0.0_dp
            end if

            if (mod(jj, 1000) == 0) then
                if (real(n_accept, kind=dp) / 1000.0_dp > 0.5_dp) then
                    delta = delta * 1.05_dp  
                else
                    delta = delta * 0.95_dp    
                end if
                n_accept = 0
            end if
        end do

        !
        write(*, *) "*********************************************"
        write(*, *) "Termalización Lista"
        write(*, *) delta
        write(*, *) "*********************************************"

        !Ciclos de Producción
        n_frame  = 0
        n_accept = 0
        max_disp = 0.0_dp
        rdf(:)       = 0.0_dp
        rdf_frame(:) = 0.0_dp
        !
        !
        do jj = 1, 10_i64*n_cycle
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

            do ii = 1, n_neigh(o)
                
                kk = neigh_list(ii, o)

                dx = x_trial - xx(kk)
                dy = y_trial - yy(kk)
                dz = z_trial - zz(kk)
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
                d2 = (xx(o)-xxo(o))**2 + (yy(o)-yyo(o))**2 + (zz(o)-zzo(o))**2
                max_disp2 = max(max_disp2, d2)
            end if

            if (max_disp2 > (0.5_dp*sigma)**2) then
                call build_neighbor_list(xx, yy, zz, n_neigh, neigh_list)
                xxo = xx
                yyo = yy
                zzo = zz
                max_disp2 = 0.0_dp
            end if

            !
            ! Aquí vamos acumulando los observables del sistema
            !
            if(mod(jj, n_sample) == 0) then
                call radial_distribution(xx, yy, zz, rdf_frame)
                
                do ii = 1, nbins 
                    rdf(ii) = rdf(ii) + rdf_frame(ii) 
                end do

                n_frame = n_frame + 1

            end if

        end do

        !Tomamos el valor promedio de la función de distrubución radial
        ! 
        do ii = 1, nbins
            rdf(ii) = rdf(ii)/n_frame
        end do 

        dbin = L/(2.0_dp * real(nbins, kind=dp))

        
        write(filename, '(A, F5.3, A)') 'rdf_eta', eta, '.dat'
        open(unit=10, file=trim(filename), status="replace", action="write")
        do ii = 1, nbins
            r = dbin*(real(ii, kind=dp) + 0.5_dp)
            write(10, *) r, rdf(ii)
        end do 

        close(10)


        !write(*, *) "Montecarlo Done"
        !write(*, *) "*******************************"
        !write(*, *) "Ratio de Aceptados: "
        !write(*, *) real(n_accept, kind=dp) / real(n_cycle, kind=dp)
        !write(*, *) delta
        !write(*, *) "*******************************"
    end subroutine montecarlo_hs





end module functions