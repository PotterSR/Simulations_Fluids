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


    subroutine md(xx, yy, zz, vx, vy, vz)
    !***************************************************
    ! Implementación de Dinámica Molecular usando Lenard-Jones
    ! Usando algoritmo de Verlet para integral posiciones
    
        real(kind=dp), intent(inout)  :: xx(N), yy(N), zz(N), vx(N), vy(N), vz(N)
        real(kind=dp)                 :: xxold(N), yyold(N), zzold(N)
        real(kind=dp)                 :: xnew, ynew, znew
        real(kind=dp)                 :: sumvsq, sumvx, sumvy, sumvz
        real(kind=dp)                 :: ax(N), ay(N), az(N)
        
        integer(kind=i64) :: n_neigh(N)
        integer(kind=i64) :: neigh_list(N, mxnb)
        integer(kind=i64) :: ii, jj, kk

        real(kind=dp)     :: dt2 
        
        
        dt2 = dt**2
        ! Inicializamos las variables
        ax(:)  = 0.0_dp
        ay(:)  = 0.0_dp
        az(:)  = 0.0_dp
        sumvx  = 0.0_dp
        sumvy  = 0.0_dp
        sumvz  = 0.0_dp
        sumvsq = 0.0_dp

        ! Inicializamos la velocidades
        do ii = 1, N 
            vx(ii) = urand(-1.0_dp, 1.0_dp)
            vy(ii) = urand(-1.0_dp, 1.0_dp)
            vz(ii) = urand(-1.0_dp, 1.0_dp)
        end do
        ! Restamos la velocidad del centro de masa para evitar el drift

        do ii = 1, N 
            sumvx = sumvx + vx(ii)
            sumvy = sumvy + vy(ii)
            sumvz = sumvz + vz(ii)
        end do

        ! Normalizamos
        sumvx = sumvx / real(N, kind=dp)
        sumvy = sumvy / real(N, kind=dp)
        sumvz = sumvz / real(N, kind=dp)

        do ii = 1, n
            vx(ii) = vx(ii) - sumvx
            vy(ii) = vy(ii) - sumvy
            vz(ii) = vz(ii) - sumvz
        end do

        call build_neighbor_list(xx, yy, zz, n_neigh, neigh_list)

        call force(xx, yy, zz, ax, ay, az, n_neigh, neigh_list)

        do ii = 1, N 
            xxold(ii) = xx(ii) - vx(ii)*dt + 0.5_dp*ax(ii)*dt2   
            yyold(ii) = yy(ii) - vy(ii)*dt + 0.5_dp*ay(ii)*dt2    
            zzold(ii) = zz(ii) - vz(ii)*dt + 0.5_dp*az(ii)*dt2    
        end do 


        do ii = 1, n_steps

            sumvx = 0.0_dp
            sumvy = 0.0_dp
            sumvz = 0.0_dp

            call  force(xx, yy, zz, ax, ay, az, n_neigh, neigh_list)

            do jj = 1, N 
                xnew = 2.0_dp * xx(jj) - xxold(jj) + dt2 * ax(jj)
                ynew = 2.0_dp * yy(jj) - yyold(jj) + dt2 * ay(jj)
                znew = 2.0_dp * zz(jj) - zzold(jj) + dt2 * az(jj)
                ! Periodic Boundary Conditions
                xnew = xnew - L*floor(xnew/L)
                ynew = ynew - L*floor(ynew/L)
                znew = znew - L*floor(znew/L)

                vx(jj) = (xnew - xxold(jj))/ (2.0_dp * dt)
                vy(jj) = (ynew - yyold(jj))/ (2.0_dp * dt)
                vz(jj) = (znew - zzold(jj))/ (2.0_dp * dt)
                
                xxold(jj) = xx(jj)
                yyold(jj) = yy(jj)
                zzold(jj) = zz(jj)

                xx(jj) = xnew
                yy(jj) = ynew
                zz(jj) = znew


            end do 

        end do


    end subroutine md

    subroutine force(xx, yy, zz, ax, ay, az, n_neigh, neigh_list)
        real(kind=dp), intent(in)  :: xx(N), yy(N), zz(N)
        real(kind=dp), intent(out) :: ax(N), ay(N), az(N)
        integer(kind=i64), intent(in) :: n_neigh(N)
        integer(kind=i64), intent(in) :: neigh_list(N, mxnb)

        real(kind=dp), parameter :: eps     = 1.0_dp      ! LJ energy scale
        real(kind=dp), parameter :: r_cut   = 2.5_dp * sigma  ! standard LJ cutoff
        real(kind=dp), parameter :: r_cut2  = r_cut**2

        real(kind=dp)            :: r2 
        real(kind=dp)            :: sr2, sr6, sr12, fij
        real(kind=dp)            :: dx, dy, dz
        integer(kind=i64)        :: ii, jj, kk


        ! Inicializamos las variables en cero
        ax(:) = 0.0_dp
        ay(:) = 0.0_dp
        az(:) = 0.0_dp

        do ii = 1, N-1

            do jj = 1, n_neigh(ii) 
                kk = neigh_list(ii, jj)
                if (kk <= ii) cycle 

                dx = xx(ii) - xx(kk)
                dy = yy(ii) - yy(kk)
                dz = zz(ii) - zz(kk)

                ! Condiciones de imagen minima
                dx = dx - L*nint(dx/L)
                dy = dy - L*nint(dy/L)
                dz = dz - L*nint(dz/L)

                r2 = dx*dx + dy*dy + dz*dz

                if(r2 < r_cut2) then
                    sr2  = (sigma)**2 /r2
                    sr6  = sr2**3
                    sr12 = sr6**2

                    fij = 4.0_dp * eps * (12.0_dp*sr12 - 6.0_dp*sr6) / r2

                    ax(ii) = ax(ii) + fij * dx
                    ay(ii) = ay(ii) + fij * dy
                    az(ii) = az(ii) + fij * dz
                    ! Tercera Ley de Newton
                    ax(kk) = ax(kk) - fij * dx
                    ay(kk) = ay(kk) - fij * dy
                    az(kk) = az(kk) - fij * dz
                end if
            end do
        end do 
    end subroutine force



end module functions