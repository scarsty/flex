module functions_init
    use functions_base
    use parameters
    use parameters2

contains

    subroutine init()
        implicit none
        integer i

        T_beta = 1d0/kB/T
        T_eV = kB*T

        iter_method = 0

        ! I_chi, (nb*nb) order
        I_chi=complex_0
        do i=1,nb*nb
            I_chi(i,i)=complex_1
        enddo

        ! I_G, nb order
        I_G=complex_0
        do i=1,nb
            I_G(i,i)=complex_1
        enddo

        mu_history=0d0
        mu_error=0d0

        mu_A = 0
        do i=1,mu_num
            mu_A(0,i)=-1d0
            mu_A(i,0)=-1d0
        enddo
        mu_b = 0d0
        mu_b(0) = -1d0

        mixer_beta0=mixer_beta

        select case (mixer_method)
            case (2:3)

            case (4)
                allocate(Jacobian(total_grid,total_grid),stat=alloc_error)
                !write(stdout,*) alloc_error
        end select

    end subroutine

    subroutine destroy()
        implicit none

        select case (mixer_method)
            case (2:3)

            case (4)
                deallocate(Jacobian,stat=alloc_error)
        end select
    end subroutine

    subroutine init_Kpoints()
        implicit none
        integer ikx, iky, irx, iry, info
        real(8) rdotk, temp(2)
        complex(8) fac

        ! 计算k点的坐标
        write(stdout, *) "Building k-points grid..."
        !zero_k = 1    ! k原点
        do ikx = 1, nkx
            do iky = 1, nky
                k(ikx, iky, 1)=1d0/nkx*(ikx-1)
                k(ikx, iky, 2)=1d0/nky*(iky-1)
                if (k(ikx, iky, 1)>0.5) k(ikx, iky, 1)=k(ikx, iky, 1)-1
                if (k(ikx, iky, 2)>0.5) k(ikx, iky, 2)=k(ikx, iky, 2)-1
                write(stdout, '(2I3,2F9.4)') ikx, iky, k(ikx,iky,:)
            enddo
        enddo
        write(stdout, *)

        ! 反傅里叶变换h0到k空间
        h0_k = complex_0
        do ikx=1,nkx; do iky=1,nky
            do irx=-rx,rx; do iry=-ry,ry
                temp=[irx,iry]
                rdotk = two_pi*dot_product(k(ikx,iky,:),temp)
                fac=exp(complex_i*rdotk)
                h0_k(:,:,ikx,iky)=h0_k(:,:,ikx,iky)+fac*h0_r(:,:,irx,iry)
            enddo; enddo

            ! 计算特征值和特征向量
            h0_k_=h0_k(:,:,ikx,iky)
            u_h0_k_=h0_k_
            call zheev('V','L',nb,u_h0_k_,nb,ev_h0_k_,ev_h0_k_lwork,nb*nb,ev_h0_k_rwork,info)
            u_h0_k(:,:,ikx,iky)=u_h0_k_
            ev_h0_k(:,ikx,iky)=ev_h0_k_
            ! write(stdout,*) diag_h0_tilde_k_
        enddo; enddo

    end subroutine

    subroutine init_U()
        implicit none
        integer ix, iy
        ! U
        ! 能带下标ab, cd -> (a+(b-1)*nb, c+(d-1)*nb)
        ! real, dimension (nb*nb, nb*nb):: U_s, U_c, U_ud, U_uu
        U_ud = 0d0
        U_uu = 0d0
        do ix=1,nb
            do iy=1,nb
                if (ix==iy) then
                    U_ud(sub_g2chi(ix,iy), sub_g2chi(ix,iy))=h1_U
                else
                    U_ud(sub_g2chi(ix,ix), sub_g2chi(iy,iy))=h1_Up
                    U_ud(sub_g2chi(ix,iy), sub_g2chi(ix,iy))=h1_J
                    U_ud(sub_g2chi(ix,iy), sub_g2chi(iy,ix))=h1_Jp
                    U_uu(sub_g2chi(ix,ix), sub_g2chi(iy,iy))=h1_Up-h1_J
                    U_uu(sub_g2chi(ix,iy), sub_g2chi(ix,iy))=-h1_Up+h1_J
                endif
            enddo
        enddo
        U_s = U_ud-U_uu
        U_c = U_ud+U_uu
    end subroutine

end module
