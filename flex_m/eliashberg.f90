subroutine eliashberg()
	use Constants
	use myfunctions
	use parameters
	use parameters2
	implicit none

	! 未完成
	! 厄立希伯格方程, 直接应用上面得到的组合

	! 自旋态不是3就是1
	! 含矩阵乘, 需改写

	integer elia1, elia2
	integer ikk,ikq,iomegak,iomegaq, k_kplusq, omega_kplusq,k_kminusq, omega_kminusq, itau, k_0minusq, k_qminusk
	integer ikk1, ikk2, iomegak1, iomegak2, k_kminusk, omega_kminusk
	integer l1,m1,l2,m2,l3,m3
	! complex, dimension (nb*nb*nk*(nomega*2-1), nb*nb*nk*(nomega*2-1)):: Elishaberg

	! ---------------------------------------------------------------------------------------------

	do ikq=1,nk; do iomegaq=-nomega,nomega
		chi_c_ = chi_0(:, :, ikq, iomegaq)
		chi_s_ = chi_0(:, :, ikq, iomegaq)
		if (spin_state==3) then
			V_s(:, :, ikq, iomegaq) = U_ud - 0.5*ABA(U_s,chi_s_) &
				-0.5*ABA(U_c, chi_c_)
		else
			V_s(:, :, ikq, iomegaq) = U_ud + 1.5*ABA(U_s,chi_s_) &
				-0.5*ABA(U_c, chi_c_)
		endif
	enddo; enddo


	! Elishberg方程所需矩阵
	! 原方程包含负号, 使用减法
	! 必要时改稀疏阵
	! sub_g2e转换坐标
	! Elishaberg=complex_0
	! 再说吧
	do l1=1,nb; do m1=1,nb
		do ikk1=1,nk; do iomegak1=-nomega, nomega
			elia1 = sub_g2e(l1,m1,ikk1,iomegak1)
			do ikk2=1,nk; do iomegak2=-nomega, nomega
				omega_kminusk = iomegak1-iomegak2
				if (abs(omega_kminusk)<=nomega) then
					k_kminusk = k_minus(ikk1, ikk2)
					do l2=1,nb; do m2=1,nb
						elia2 = sub_g2e(l2,m2,ikk2,iomegak2)
						do l3=1,nb; do m3=1,nb

							!Elishaberg(elia1, elia2) = Elishaberg(elia1, elia2) &
								!   - &
								!   V_s(sub_g2chi(l1,l3), sub_g2chi(m3,m1), k_kminusk, omega_kminusk) &
								!   *G(l3,l2,ikk2,iomegak2)*conjg(G(m3,m2,ikk2,iomegak2))

						enddo; enddo
					enddo; enddo
				endif
			enddo; enddo
		enddo; enddo
	enddo; enddo

	! 求特征值和特征向量, 调用数学库, 未完成

	!call ()

end subroutine eliashberg
