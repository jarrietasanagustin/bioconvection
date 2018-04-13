subroutine third_stage(k,k_t)

use definiciones
implicit none
integer(4)::i,j,k,k_t


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!Tercera etapa del esquema RK3
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

do i=2,nx_i

do j=2,nz_i

k_omega(1,1)=-u_old(i,j)*(omega_old(i+1,j)-omega_old(i-1,j))/(2.0d0*delta_x)-&
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w_old(i,j)*(omega_old(i,j+1)-omega_old(i,j-1))/(2.0d0*delta_z)+&
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
(omega_old(i+1,j)-2.0d0*omega_old(i,j)+omega_old(i-1,j))/(sqrt(gr)*delta_x**2)+&
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
(omega_old(i,j+1)-2.0d0*omega_old(i,j)+omega_old(i,j-1))/(sqrt(gr)*delta_z**2)-&
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
(n_old(i+1,j)-n_old(i-1,j))/(2.0d0*delta_x)

k_omega(2,1)=-u_aprox(i,j)*(omega_aprox(i+1,j)-omega_aprox(i-1,j))/(2.0d0*delta_x)-&
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w_aprox(i,j)*(omega_aprox(i,j+1)-omega_aprox(i,j-1))/(2.0d0*delta_z)+&
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
(omega_aprox(i+1,j)-2.0d0*omega_aprox(i,j)+omega_aprox(i-1,j))/(sqrt(gr)*delta_x**2)+&
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
(omega_aprox(i,j+1)-2.0d0*omega_aprox(i,j)+omega_aprox(i,j-1))/(sqrt(gr)*delta_z**2)-&
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
(n_aprox(i+1,j)-n_aprox(i-1,j))/(2.0d0*delta_x)

k_omega(3,1)=-u_aprox_1(i,j)*(omega_aprox_1(i+1,j)-omega_aprox_1(i-1,j))/(2.0d0*delta_x)-&
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w_aprox_1(i,j)*(omega_aprox_1(i,j+1)-omega_aprox_1(i,j-1))/(2.0d0*delta_z)+&
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
(omega_aprox_1(i+1,j)-2.0d0*omega_aprox_1(i,j)+omega_aprox_1(i-1,j))/(sqrt(gr)*delta_x**2)+&
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
(omega_aprox_1(i,j+1)-2.0d0*omega_aprox_1(i,j)+omega_aprox_1(i,j-1))/(sqrt(gr)*delta_z**2)-&
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
(n_aprox_1(i+1,j)-n_aprox_1(i-1,j))/(2.0d0*delta_x)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    k_n(1,1)=-((u_old(i+1,j)+u_photo(i+1,j))*n_old(i+1,j)&
        -(u_old(i-1,j)+u_photo(i-1,j))*n_old(i-1,j))/(2.0d0*delta_x)-&
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
((w_old(i,j+1)+w_photo(i,j+1))*n_old(i,j+1)-&
       (w_old(i,j-1)+w_photo(i,j-1))*n_old(i,j-1))/(2.0d0*delta_z)+&
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
(n_old(i+1,j)-2.0d0*n_old(i,j)+n_old(i-1,j))/(sc*sqrt(gr)*(delta_x**2))+&
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        (n_old(i,j+1)-2.0d0*n_old(i,j)+n_old(i,j-1))/(sc*sqrt(gr)*(delta_z**2))
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    k_n(2,1)=-((u_aprox(i+1,j)+u_photo(i+1,j))*n_aprox(i+1,j)&
        -(u_aprox(i-1,j)+u_photo(i-1,j))*n_aprox(i-1,j))/(2.0d0*delta_x)-&
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
((w_aprox(i,j+1)+w_photo(i,j+1))*n_aprox(i,j+1)-&
       (w_aprox(i,j-1)+w_photo(i,j-1))*n_aprox(i,j-1))/(2.0d0*delta_z)+&
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
(n_aprox(i+1,j)-2.0d0*n_aprox(i,j)+n_aprox(i-1,j))/(sc*sqrt(gr)*(delta_x**2))+&
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        (n_aprox(i,j+1)-2.0d0*n_aprox(i,j)+n_aprox(i,j-1))/(sc*sqrt(gr)*(delta_z**2))
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    k_n(3,1)=-((u_aprox_1(i+1,j)+u_photo(i+1,j))*n_aprox_1(i+1,j)&
        -(u_aprox_1(i-1,j)+u_photo(i-1,j))*n_aprox_1(i-1,j))/(2.0d0*delta_x)-&
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
((w_aprox_1(i,j+1)+w_photo(i,j+1))*n_aprox_1(i,j+1)-&
       (w_aprox_1(i,j-1)+w_photo(i,j-1))*n_aprox_1(i,j-1))/(2.0d0*delta_z)+&
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
(n_aprox_1(i+1,j)-2.0d0*n_aprox_1(i,j)+n_aprox_1(i-1,j))/(sc*sqrt(gr)*(delta_x**2))+&
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        (n_aprox_1(i,j+1)-2.0d0*n_aprox_1(i,j)+n_aprox_1(i,j-1))/(sc*sqrt(gr)*(delta_z**2))
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


k_n_p(1,1)=-((u_old(i+1,j))*n_passive(i+1,j)&
        -(u_old(i-1,j))*n_passive(i-1,j))/(2.0d0*delta_x)-&
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
((w_old(i,j+1))*n_passive(i,j+1)-&
       (w_old(i,j-1))*n_passive(i,j-1))/(2.0d0*delta_z)+&
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
(n_passive(i+1,j)-2.0d0*n_passive(i,j)+n_passive(i-1,j))/(pe*(delta_x**2))+&
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        (n_passive(i,j+1)-2.0d0*n_passive(i,j)+n_passive(i,j-1))/(pe*(delta_z**2))
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


k_n_p(2,1)=-((u_aprox(i+1,j))*n_p_aprox(i+1,j)&
        -(u_aprox(i-1,j))*n_p_aprox(i-1,j))/(2.0d0*delta_x)-&
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
((w_aprox(i,j+1))*n_p_aprox(i,j+1)-&
       (w_aprox(i,j-1))*n_p_aprox(i,j-1))/(2.0d0*delta_z)+&
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
(n_p_aprox(i+1,j)-2.0d0*n_p_aprox(i,j)+n_p_aprox(i-1,j))/(pe*(delta_x**2))+&
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        (n_p_aprox(i,j+1)-2.0d0*n_p_aprox(i,j)+n_p_aprox(i,j-1))/(pe*(delta_z**2))
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k_n_p(3,1)=-((u_aprox_1(i+1,j))*n_p_aprox_2(i+1,j)&
        -(u_aprox_1(i-1,j))*n_p_aprox_2(i-1,j))/(2.0d0*delta_x)-&
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
((w_aprox_1(i,j+1))*n_p_aprox_2(i,j+1)-&
       (w_aprox_1(i,j-1))*n_p_aprox_2(i,j-1))/(2.0d0*delta_z)+&
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
(n_p_aprox_2(i+1,j)-2.0d0*n_p_aprox_2(i,j)+n_p_aprox_2(i-1,j))/(pe*(delta_x**2))+&
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        (n_p_aprox_2(i,j+1)-2.0d0*n_p_aprox_2(i,j)+n_p_aprox_2(i,j-1))/(pe*(delta_z**2))
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




omega_new(i,j)=omega_old(i,j)+delta_t/6.0d0*(k_omega(1,1)+4.0d0*k_omega(2,1)+k_omega(3,1))
n_new(i,j)=n_old(i,j)+delta_t/6.0d0*(k_n(1,1)+4.0d0*k_n(2,1)+k_n(3,1))
n_p_new(i,j)=n_passive(i,j)+delta_t/6.0d0*(k_n_p(1,1)+4.0d0*k_n_p(2,1)+k_n_p(3,1))
end do

end do



end subroutine third_stage
