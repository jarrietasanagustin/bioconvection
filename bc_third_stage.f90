subroutine bc_third_stage(k,k_t)

use definiciones
implicit none

integer(4)::iter,i,j,k,k_t
real(8)::a_top,b_top,c_top,a_bot,b_bot,c_bot
real(8)::a_left,b_left,c_left,a_right,b_right,c_right


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!Calculo de velocidades y aplicacion de condiciones de
!contorno
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!----------------------------------------
!Condiciones de contorno para la vorticidad
!---------------------------------------

do j=1,nz_i+1
i=1
!omega_new(i,j)=(psi_new(i+2,j)*(x(i+1)-x(i))**3-psi_new(i+1,j)*(x(i+2)-x(i))**3)/&
!(0.5d0*((x(i+2)-x(i))**3*(x(i+1)-x(i))**2-(x(i+1)-x(i))**3*(x(i+2)-x(i))**2))
!omega_aprox_2(i,j)=-2.0d0*psi_new(i+1,j)*(h_x**2)/delta_x**2
!omega_aprox_2(i,j)=(psi_new(i+2,j)-8.0d0*psi_new(i+1,j))/(2.0d0*delta_x**2)
omega_new(i,j)=(psi_new(i+2,j)-8.0d0*psi_new(i+1,j))/(2.0d0*delta_x**2)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=nx_i+1
!omega_new(i,j)=(psi_new(i-2,j)*(x(i)-x(i-1))**3-psi_new(i-1,j)*(x(i)-x(i-2))**3)/&
!(0.5d0*((x(i)-x(i-2))**3*(x(i)-x(i-1))**2-(x(i)-x(i-1))**3*(x(i)-x(i-2))**2))
!omega_aprox_2(i,j)=-2.0d0*psi_new(i-1,j)*(h_x**2)/delta_x**2
!omega_aprox_2(i,j)=(psi_new(i-2,j)-8.0d0*psi_new(i-1,j))/(2.0d0*delta_x**2)
omega_new(i,j)=(psi_new(i-2,j)-8.0d0*psi_new(i-1,j))/(2.0d0*delta_x**2)
end do

do i=1,nx_i+1
j=1
!omega_new(i,j)=(psi_new(i,j+2)*(z(j+1)-z(j))**3-psi_new(i,j+1)*(z(j+2)-z(j))**3)/&
!(0.5d0*((z(j+2)-z(j))**3*(z(j+1)-z(j))**2-(z(j+1)-z(j))**3*(z(j+2)-z(j))**2))
!omega_aprox_2(i,j)=-2.0d0*psi_new(i,j+1)*(h_z**2)/delta_z**2
!omega_aprox_2(i,j)=(psi_new(i,j+2)-8.0d0*psi_new(i,j+1))/(2.0d0*delta_z**2)
omega_new(i,j)=(psi_new(i,j+2)-8.0d0*psi_new(i,j+1))/(2.0d0*delta_z**2)
j=nz_i+1
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!omega_new(i,j)=(psi_new(i,j-2)*(z(j)-z(j-1))**3-psi_new(i,j-1)*(z(j)-z(j-2))**3)/&
!(0.5d0*((z(j)-z(j-2))**3*(z(j)-z(j-1))**2-(z(j)-z(j-1))**3*(z(j)-z(j-2))**2))
!omega_aprox_2(i,j)=-2.0d0*psi_new(i,j-1)*(h_z**2)/delta_z**2
!omega_aprox_2(i,j)=(psi_new(i,j-2)-8.0d0*psi_new(i,j-1))/(2.0d0*delta_z**2)
omega_new(i,j)=(psi_new(i,j-2)-8.0d0*psi_new(i,j-1))/(2.0d0*delta_z**2)
end do
do j=1,nz_i
i=1
!u_aprox_2(i,j)=0.0d0
!w_aprox_2(i,j)=0.0d0
u_new(i,j)=0.0d0
w_new(i,j)=0.0d0
i=nx_i+1
!u_aprox_2(i,j)=0.0d0
!w_aprox_2(i,j)=0.0d0
u_new(i,j)=0.0d0
w_new(i,j)=0.0d0
end do
do i=1,nx_i+1
j=1
!u_aprox_2(i,j)=0.0d0
!w_aprox_2(i,j)=0.0d0
u_new(i,j)=0.0d0
w_new(i,j)=0.0d0
j=nz_i+1
!u_aprox_2(i,j)=0.0d0
!w_aprox_2(i,j)=0.0d0
u_new(i,j)=0.0d0
w_new(i,j)=0.0d0
end do

!--------------------------------------
!Velocidad en el interior del dominio
!--------------------------------------

do i=2,nx_i
do j=2,nz_i
!u_aprox_2(i,j)=(psi_new(i,j+1)-psi_new(i,j-1))/(2.0d0*delta_z)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!w_aprox_2(i,j)=-(psi_new(i+1,j)-psi_new(i-1,j))/(2.0d0*delta_x)
u_new(i,j)=(psi_new(i,j+1)-psi_new(i,j-1))/(2.0d0*delta_z)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w_new(i,j)=-(psi_new(i+1,j)-psi_new(i-1,j))/(2.0d0*delta_x)
end do
end do

do i=1,nx_i+1
j=1
a_bot=(z(j+1)-z(j))/((z(j+2)-z(j))*(z(j+1)-z(j))-(z(j+2)-z(j))**2)
b_bot=-(z(j+2)-z(j))/((z(j+1)-z(j))**2-(z(j+2)-z(j))*(z(j+1)-z(j)))
c_bot=((z(j+2)-z(j))**2-(z(j+1)-z(j))**2)/((z(j+2)-z(j))*(z(j+1)-z(j))**2-(z(j+2)-z(j))**2*(z(j+1)-z(j)))
!c_bot=-1.0d0/(z(j+1)-z(j))
!b_bot=1.0d0/(z(j+1)-z(j))
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!n_new(i,j)=(b_bot*n_new(i,j+1)+a_bot*n_new(i,j+2))/(sc*sqrt(gr)*w_photo(i,j)-c_bot)
n_new(i,j)=(2.0d0*n_new(i,j+1)-0.5d0*n_new(i,j+2))/(sc*sqrt(gr)*w_photo(i,j)*delta_z+1.5d0)
n_p_new(i,j)=(2.0d0*n_p_new(i,j+1)-0.5d0*n_p_new(i,j+2))/(1.5d0)
!n_aprox_2(i,j)=(b_bot*n_aprox_2(i,j+1)+a_bot*n_aprox_2(i,j+2))/(pe*w_photo(i,j)-c_bot)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

j=nz_i+1


a_top=(-(z(j-2)-z(j))**2+(z(j-1)-z(j))**2)/(-(z(j-2)-z(j))*(z(j-1)-z(j))**2+(z(j-2)-z(j))**2*(z(j-1)-z(j)))
b_top=(z(j-2)-z(j))/(-(z(j-1)-z(j))**2+(z(j-2)-z(j))*(z(j-1)-z(j)))
c_top=-(z(j-1)-z(j))/(-(z(j-2)-z(j))*(z(j-1)-z(j))+(z(j-2)-z(j))**2)
!n_aprox_2(i,j)=(b_top*n_aprox_2(i,j-1)+c_top*n_aprox_2(i,j-2))/(pe*w_photo(i,j)-a_top)
!n_new(i,j)=(b_top*n_new(i,j-1)+c_top*n_new(i,j-2))/(sc*sqrt(gr)*w_photo(i,j)-a_top)
n_new(i,j)=(2.0d0*n_new(i,j-1)-0.5d0*n_new(i,j-2))/(-sc*sqrt(gr)*w_photo(i,j)*delta_z+1.5d0)
n_p_new(i,j)=(2.0d0*n_p_new(i,j-1)-0.5d0*n_p_new(i,j-2))/(1.5d0)
end do

i=1
a_left=(x(i+1)-x(i))/((x(i+2)-x(i))*(x(i+1)-x(i))-(x(i+2)-x(i))**2)
b_left=-(x(i+2)-x(i))/((x(i+1)-x(i))**2-(x(i+2)-x(i))*(x(i+1)-x(i)))
c_left=((x(i+2)-x(i))**2-(x(i+1)-x(i))**2)/((x(i+2)-x(i))*(x(i+1)-x(i))**2-(x(i+2)-x(i))**2*(x(i+1)-x(i)))
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=nx_i+1
a_right=(-(x(i-2)-x(i))**2+(x(i-1)-x(i))**2)/(-(x(i-2)-x(i))*(x(i-1)-x(i))**2+(x(i-2)-x(i))**2*(x(i-1)-x(i)))
b_right=(x(i-2)-x(i))/(-(x(i-1)-x(i))**2+(x(i-2)-x(i))*(x(i-1)-x(i)))
c_right=-(x(i-1)-x(i))/(-(x(i-2)-x(i))*(x(i-1)-x(i))+(x(i-2)-x(i))**2)


do j=2,nz_i
i=1
!n_aprox_2(i,j)=(b_left*n_aprox_2(i+1,j)+a_left*n_aprox_2(i+2,j))/(pe*u_photo(i,j)-c_left)
n_new(i,j)=(b_left*n_new(i+1,j)+a_left*n_new(i+2,j))/(sc*sqrt(gr)*u_photo(i,j)-c_left)
n_p_new(i,j)=(b_left*n_p_new(i+1,j)+a_left*n_p_new(i+2,j))/(-c_left)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=nx_i+1
!n_aprox_2(i,j)=(b_right*n_aprox_2(i-1,j)+c_right*n_aprox_2(i-2,j))/(pe*u_photo(i,j)-a_right)
n_new(i,j)=(b_right*n_new(i-1,j)+c_right*n_new(i-2,j))/(sc*sqrt(gr)*u_photo(i,j)-a_right)
n_p_new(i,j)=(b_right*n_p_new(i-1,j)+c_right*n_p_new(i-2,j))/(-a_right)
end do

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!Fin condiciones de contorno
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end subroutine bc_third_stage
