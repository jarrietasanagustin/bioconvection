module definiciones
implicit none

integer(4)::nz_i,nx_i,nt_i,ne_i
real(8)::nz,nx,nt,delta_t
real(8)::eps
real(8)::t_max,h_x,h_z,x_foco_1,x_foco_2
real(8)::sc,gr,re,ra,pr,alpha_phot,v_fall

real(8)::beta_z,tau,b_mesh,b_mesh_z,delta_z,delta_x,pe

integer(4)::i_write,num_files
real(8),dimension(:),allocatable::time,z,x

real(8),dimension(4,1)::k_omega,k_n,k_n_p




real(8),dimension(:,:),allocatable::psi_new,u_new,w_new,omega_new,n_new,psi_old,omega_old,n_old,u_old,w_old,n_passive,n_p_aprox&
                    ,n_p_aprox_2,n_p_new
real(8),dimension(:,:),allocatable::u_aprox,w_aprox,omega_aprox,n_aprox,u_aprox_1,w_aprox_1,&
                                    n_aprox_1,omega_aprox_1,omega_aprox_2,n_aprox_2,u_aprox_2,w_aprox_2

real(8),dimension(:,:),allocatable::w_photo,u_photo,intensity,u_photo_1,w_photo_1,intensity_1,intensity_2,u_photo_2,w_photo_2

real(8),parameter::relax_f=1.0d0,relax_gauss=0.5d0
real(8),parameter::relax_psi=1.95d0
real(8)::x_foco,z_foco,i_max

end module definiciones

