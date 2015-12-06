module functions_mixer
    use functions_base
    use parameters
    use parameters2

contains
    ! 初始化混合器
    subroutine mixer_init()
        implicit none
        integer i

        mixer_G = complex_0
        mixer_beta=mixer_beta0
        !G_mixer(:,:,:,:,:,1)=G_out
        !error_mixer(:,:,:,:,:,1) = G_out-G0

        mixer_pointer=1
        mixer_pointer2=1
        mixer_order=0
        mixer_A = 0
        do i=1,mix_num
            mixer_A(0,i)=-1d0
            mixer_A(i,0)=-1d0
        enddo
        mixer_b = 0d0
        mixer_b(0) = -1d0

        if (mixer_method==4) then
            Jacobian=complex_0
            !这里我也不知道为什么要用负定阵, 资料上给的好像是正定
            do i=1,total_grid
                Jacobian(i,i)=-complex_1*mixer_beta
            enddo
        endif

    end subroutine


    subroutine mixer_linear()
        implicit none
        integer i
        real(8) beta, min_error, min_beta, cur_error
        !real(8), external :: dznrm2

        if (mixer_beta0==0) then
            call cal_best_mixer_beta()
        else
            G=mixer_beta*G_out+(1-mixer_beta)*G
        endif

    end subroutine

    ! pulay mixer 相关
    ! 移动指针
    function mixerIncPointer(p, n, max_value)
        implicit none
        integer mixerIncPointer, p, n, max_value

        mixerIncPointer = p+n
        mixerIncPointer = mod(mixerIncPointer, max_value)
        if (mixerIncPointer==0) mixerIncPointer=max_value
    end function

    ! G_out是新的, G是上一步
    ! 混合算法
    ! http://vergil.chemistry.gatech.edu/notes/diis/node2.html
    subroutine mixer_Pulay()
        implicit none
        integer n, i, info, next_pointer, prev_pointer
        integer ipiv(mix_num+1)
        real(8), dimension (mix_num*2) :: lwork
        real(8) e, e0, max_error
        logical find_bigger

        ! method.3 - Refined Pulay方法, 保留残差较小的, 实际上高度非线性时没啥作用
        if (mixer_method==3) then
            if (G_iter<=mix_keep) then
                mixer_order=G_iter
            else
                mixer_error_=G_out-G
                e0=real(GProduct(mixer_error_,mixer_error_))
                ! 在误差列表中找最大的取代
                find_bigger=.false.
                max_error=0
                do i=1,mix_keep
                    if ((mixer_A(i,i))>max_error) then
                        max_error = (mixer_A(i,i))
                        mixer_pointer=i
                    endif
                enddo
                if (max_error>e0) then
                    find_bigger=.true.
                endif
                ! 如果没有找到取代位置, 则在后面找一个空位
                if (.not. find_bigger) then
                    mixer_pointer=mixer_pointer2+mix_keep
                    mixer_pointer2 = mixerIncPointer(mixer_pointer2,1,mix_num-mix_keep)
                    mixer_order=min(mixer_order+1,mix_num)
                endif
            endif
        else
            mixer_order=min(mixer_order+1,mix_num)
        endif

        !prev_pointer=mixerIncPointer(mixer_pointer,-1)

        ! 误差是G_out-G, 这个值越小则说明G是一个接近好的解, 而非G_out
        mixer_error(:,:,:,:,:,mixer_pointer)=G-G_out
        if (mixer_beta==0) then
            call cal_best_mixer_beta()
            mixer_G(:,:,:,:,:,mixer_pointer)=G
        else
            mixer_G(:,:,:,:,:,mixer_pointer)=mixer_beta*G_out+(1-mixer_beta)*G
        endif
        !mixer_G(:,:,:,:,:,mixer_pointer)=G_out

        ! A_ij=e_i**H*e_j
        mixer_error_=mixer_error(:,:,:,:,:,mixer_pointer)
        !omp parallel do private(mixer_error2_,e)
        do i=1,mix_num
            mixer_error2_=mixer_error(:,:,:,:,:,i)
            e=real(GProduct(mixer_error_,mixer_error2_))
            mixer_A(mixer_pointer,i)=e
            mixer_A(i,mixer_pointer)=e
            !write(*,*)e
        enddo
        !omp end parallel do

        next_pointer=mixerIncPointer(mixer_pointer,1,mix_num)
        mixer_pointer=next_pointer

        mixer_A1=mixer_A
        mixer_x=mixer_b

        n=mixer_order
        ! write(stderr,*) n
        ! 系数矩阵实际上多一行
        call dsysv('U', n+1, 1, mixer_A1, mix_num+1, ipiv, mixer_x, mix_num+1, lwork, 2*mix_num, info)
        !call zgesv(n, 1, Pulay_A1, mix_num+1, ipiv, Pulay_x, mix_num+1, info)
        G=complex_0
        do i=1,n
            G=G+mixer_G(:,:,:,:,:,i)*(mixer_x(i))
            !write(stderr,*) n,real(mixer_x(i)), real(mixer_A(i,i))
        enddo
        !G=mixer_beta*G_out+(1-mixer_beta)*G
        !call writematrix(Pulay_A,11)
        !if (n==1) stop
    end subroutine

    !未完成
    subroutine mixer_Broyden()
        implicit none
        complex(8) fac, fac2
        complex(8), external :: zdotc
        integer i

        ! delta R
        !        G_error=G_out-G
        !        G_out=G_error-G_error0
        !        fac=1/dznrm2(total_grid,G_out,1)**2
        !        deltaG=G-G_prev
        !        G_prev=G
        !        if (G_iter>1) then
        !            call zgemv('N',total_grid,total_grid,-complex_1,Jacobian,total_grid, &
            !                G_out,1,complex_1,deltaG,1)
        !            call zgerc(total_grid,total_grid,fac,deltaG,1, &
            !                G_out,1,Jacobian,total_grid)
        !        endif
        !        call zgemv('N',total_grid,total_grid,-complex_1,Jacobian,total_grid, &
            !            G_error,1,complex_1,G,1)
        !        G_error0=G_error
        !do i=1,total_grid
        !write(stderr,*)G
        !enddo
        !stop
    end subroutine


    subroutine cal_best_mixer_beta()
        implicit none
        integer i,f_(1),f
        real(8) beta(3), min_error, min_beta, beta1, beta2, error0(3)

        !        if (G_iter>20)then
        !            mixer_beta=0.001
        !            G=mixer_beta*G_out+(1-mixer_beta)*G
        !            mixer_method=1
        !            mixer_beta=0.01
        !            return
        !        endif


        !        G_in0=G
        !        G_out0=G_out
        !        G_error=G_out-G
        !        min_error=1d100
        !        beta(1)=-1
        !        beta(3)=1
        !        beta(2)=(beta(1)+beta(3))/2
        !
        !        G=beta(1)*G_error+G_in0
        !        call cal_G_out()
        !        G=G_out-G
        !        error0(1)=dznrm2(total_grid,G,1)
        !
        !        G=beta(3)*G_error+G_in0
        !        call cal_G_out()
        !        G=G_out-G
        !        error0(3)=dznrm2(total_grid,G,1)
        !
        !        G=beta(2)*G_error+G_in0
        !        call cal_G_out()
        !        G=G_out-G
        !        error0(2)=dznrm2(total_grid,G,1)
        !        i=0
        !        do i=1,10
        !            f_=minloc(error0)
        !            f=f_(1)
        !            if (f==2)then
        !                !中间的最小，截断大的一边
        !                f_=maxloc(error0)
        !                f=f_(1)
        !                beta(f)=beta(2)
        !                error0(f)=error0(2)
        !
        !                beta(2)=(beta(1)+beta(3))/2
        !                G=beta(2)*G_error+G_in0
        !                call cal_G_out()
        !                G=G_out-G
        !                error0(2)=dznrm2(total_grid,G,1)
        !
        !            else
        !                !中间的不是最小，找小的，外推
        !                f_=minloc(error0)
        !                f=f_(1)
        !
        !                beta(4-f)=beta(2)
        !                error0(4-f)=error0(2)
        !
        !                beta(2)=beta(f)
        !                error0(2)=error0(f)
        !
        !                beta(f)=beta(2)+(beta(2)-beta(4-f))
        !
        !                G=beta(f)*G_error+G_in0
        !                call cal_G_out()
        !                G=G_out-G
        !                error0(f)=dznrm2(total_grid,G,1)
        !
        !
        !            endif
        !
        !            !write(stdout,*) beta
        !            !write(stdout,*) error0
        !        enddo
        !        f_=minloc(error0)
        !        mixer_beta=beta(f_(1))
        !
        !        if (mixer_beta==0) then
        !            f_=maxloc(error0)
        !            mixer_beta=0.5*(sum(beta)-beta(f_(1)))
        !        endif
        !
        !        G=mixer_beta*G_out0+(1-mixer_beta)*G_in0
    end subroutine

    ! 修改的牛顿迭代, 使用已经得到的结果拟合一个多项式, 求其导数代入牛顿迭代
    ! 最多100次, 100次仍不收敛报错不管了
    subroutine modify_mu()
        implicit none
        integer n, mu_pointer, i, info
        real(8) d
        integer ipiv(mu_num+1)
        real(8), dimension (mu_num*2) :: lwork

        mu_pointer=density_iter-1
        n=density_iter
        ! 误差是cur_density-target_density
        mu_history(mu_pointer)=mu
        mu_b(mu_pointer)=cur_density-target_density

        do i=0,mu_num
            !mu_A(i,mu_pointer)=mu_history(i)**mu_pointer
            mu_A(mu_pointer,i)=mu_history(mu_pointer)**i
        enddo

        mu_A1=mu_A
        mu_x=mu_b

        call dgesv(n,1,mu_A1,mu_num+1,ipiv,mu_x,mu_num+1,info)

        d=0d0
        do i=1,mu_pointer
            d=d+mu_x(i)*i*mu**(i-1)
        enddo

        mu=mu-mu_b(mu_pointer)/d

        if (density_iter==1) then
            mu=maxval(eigen_value)
            !mu=mu+10d0
        endif

    end subroutine

    !使用Pulay混合得到一个新的mu
    !如果对单值函数使用Pulay方法, 得到的系数矩阵很大可能是奇异阵, 会导致结果不正确, 废弃
    subroutine modify_mu_pulay()
        implicit none
        integer n, mu_pointer, i, info
        real(8) e
        integer ipiv(mu_num+1)
        real(8), dimension (mu_num*2) :: lwork

        mu_pointer=mod(density_iter,mu_num)
        if (mu_pointer==0) mu_pointer=mu_num
        n=min(density_iter,mu_num)
        ! 误差是cur_density-target_density
        mu_error(mu_pointer)=cur_density-target_density
        mu_history(mu_pointer)=mu

        ! A_ij=e_i**H*e_j
        !!$omp parallel do private(e)
        do i=1,mu_num
            e=mu_error(i)*mu_error(mu_pointer)
            mu_A(mu_pointer,i)=e
            mu_A(i,mu_pointer)=e
        enddo
        !!$omp end parallel do


        mu_A1=mu_A
        mu_x=mu_b

        call dsysv('U', n+1, 1, mu_A1, mu_num+1, ipiv, mu_x, mu_num+1, lwork, 2*mu_num, info)
        !call zgesv(n, 1, Pulay_A1, mix_num+1, ipiv, Pulay_x, mix_num+1, info)
        mu=0d0
        do i=1,n
            mu=mu+mu_history(i)*mu_x(i)
            !write(stdout,*) mu_history(i),mu_x(i),mu_error(i)
        enddo
        do i=1,n
            !mu=mu+mu_history(i)*mu_x(i)
            !write(stdout,*) mu_history(i),mu_x(i),mu_error(i)
        enddo
        !write(stdout,*) mu
        if (density_iter==1) then
            mu=maxval(eigen_value)
            !mu=eigen_value(5)
        endif

    end subroutine

end module
