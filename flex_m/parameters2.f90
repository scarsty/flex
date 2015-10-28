!����ר��
module parameters2

    use constants
    implicit none

    public
    ! ���ֺ���, �������ֺ���
    ! ���ܺ���, �������ܺ���
    complex(8), dimension (nb, nb, nkx, nky, minomegaf:maxomegaf) :: &
        G, conjgG, G0, G1, sigma, sigma0 ! ���ֺ���������

    complex(8), dimension (nb, nb, nkx, nky, dft_grid) :: &
        r_tau1, r_tau2

    ! ������, susceptibilities, effective interactions
    complex(8), dimension (nb*nb, nb*nb, nkx, nky, minomegab:maxomegab) :: &
        chi_0, V !,V_s, chi_s, chi_c,

    complex(8), dimension (nb*nb, nb*nb, nkx, nky, dft_grid) :: &
        r_tau_sqr


    ! ������, ��λ����
    complex(8), dimension (nb*nb, nb*nb) :: U_s, U_c, U_ud, U_uu, I_chi
    ! dyson���̵ĸ����������
    complex(8), dimension (nb, nb) :: I_G, G_, G0_, sigma_

    ! H0
    ! u_h0_k�Ǽ���G0��ʹ�õĸ���������, G0=u_h0_k**H*diag(1/(i*omage-e+mu))*u_h0_k
    complex(8), dimension (nb, nb, -rx:rx, -ry:ry) :: h0_r
    complex(8), dimension (nb, nb, nkx, nky) :: h0_k
    real(8), dimension (nb, nb, nkx, nky) :: u_tilde_k, h0_tilde_k

    ! k�ռ��Ӧ
    real(8), dimension (nkx, nky, 2) :: k
    ! ����k�Ĳ��Ӧ��k
    integer, dimension (nk, nk) :: k_minus, k_plus

    ! ����g0�ĸ���
    complex(8), dimension (nb, nb, nkx, nky) :: u_h0_k
    complex(8), dimension (nb, nb) :: I_g0, i_plus, i_minus, h0_k_, diag_h0_G0_, u_h0_k_, ev_h0_k_lwork
    real(8), dimension (nb, nb) :: u_g0, u_tilde_k_, h0_tilde_k_, diag_h0_tilde_k_lwork, diag_test, ev_h0_k_rwork
    real(8), dimension (nb, nkx, nky) :: ev_h0_k
    real(8), dimension (nb) :: ev_h0_k_

    ! ����chi�ĸ���
    complex(8), dimension (nb*nb, nb*nb) :: chi_0_, chi_c_, chi_s_, Iminuschi_0_

    ! ����Ҷ�任, ������Ҷ�任�ĸ���
    complex(8), dimension (nkx, nky, dft_grid) :: dft_in
    complex(8), dimension (nkx, nky, dft_grid) :: dft_out

    ! Pulay mixer ���
    integer, parameter :: mix_num  = 50, mix_keep = 40
    integer mixer_pointer
    complex(8), dimension (nb, nb, nkx, nky, minomegaf:maxomegaf, mix_num) :: mixer_G, mixer_error
    complex(8), dimension (nb, nb, nkx, nky, minomegaf:maxomegaf) :: mixer_G_, mixer_error_
    complex(8), dimension (0:mix_num, 0:mix_num) :: mixer_A, mixer_A1
    complex(8), dimension(0:mix_num) :: mixer_x, mixer_b

    ! ռ�������
    real(8) mu_less(3), mu_more(3), density_less(3), density_more(3)
    integer mu_less_count, mu_more_count

    ! ����sigma����, density����
    integer sigma_iter, density_iter, total_iter
    real(8) cur_density, density_base
    logical sigma_conv, density_conv

end module parameters2
