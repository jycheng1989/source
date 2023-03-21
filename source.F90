!Author Junyi Cheng Oct.6 2022
!2D (radial and energy) source with the OBR5 formula (Physics of Plasmas 15, 052308 (2008))
!Currently focus on the ion. Will update to electron later.
module source_2d
  use initial_perturbation 
  integer :: src_2d_period=1 !source period
  integer :: num_psi_src, num_energy_src ! the grid numer for psi and energy space
  integer :: length
  real (8) :: denergy_src, cut_energy_src ! grid size and cutoff of energy
  real (8) :: gamma_src ! frequency of the source rate 
  real (8), allocatable :: df_src(:,:),data_source(:,:),g_local_rev(:,:),g_local(:,:), vol_local(:,:), vol_local_space(:), f_src(:,:), f0_grid_src(:,:) ! perturbed f in (psi, energy) space
  real (8), allocatable :: dn_src(:),n_profile_list(:),t_profile_list(:) ! velocity-space integral of perturbed f in (psi, energy) space
  real (8), allocatable :: n_src(:) ! velocity-space intergral of f0 in (psi, energy) space

  type source_ptl_type
    real (8), allocatable :: phase(:,:)
    integer (8), allocatable :: gid(:)
  end type source_ptl_type

  type(source_ptl_type), allocatable :: source_ptl  
contains
  
  !Allocate memories and initialize paramers
  subroutine source_2d_init(grid)
    use initial_perturbation
    use grid_class
    use sml_module
    implicit none
    type(grid_type) :: grid

    num_energy_src=esrc_in
    gamma_src=gamma_in
    cut_energy_src=vcut_in 
    num_psi_src=grid%npsi_surf2 
    denergy_src=cut_energy_src/real(num_energy_src-1)
    allocate(df_src(num_psi_src,num_energy_src),f_src(num_psi_src,num_energy_src),f0_grid_src(num_psi_src,num_energy_src),g_local_rev(num_psi_src,num_energy_src),g_local(num_psi_src,num_energy_src),vol_local(num_psi_src,num_energy_src))
    allocate(dn_src(num_psi_src),vol_local_space(num_psi_src))
    allocate(n_src(num_psi_src))
    allocate(n_profile_list(num_psi_src),t_profile_list(num_psi_src))

    allocate(source_ptl)

   !vol_local_space=diag_1d_vol
   ! do m=1,num_psi_src
   !   psi=psi_tmp(m)
   !   ter=eq_ftn(psi,r,z,eq_temp(sp_type))
   !   do j=1,num_energy_src
   !      vfrac1=j*denergy_src*ter
   !      vfrac=(j+0.5)*denergy_src*ter
   !      vfrac2=(j+1)*denergy_src*ter
   !      v=sqrt(2.0*vfrac)
   !      dv=sqrt(2.0*vfrac2)-sqrt(2.0*vfrac1)
   !      vol_local(m,j)=1.0/(vol_local_space(m)*4*sml_pi*dv*v**2)
   !   enddo
   !  enddo

  end subroutine

  !Add 2D source in weight
  subroutine source_2d(grid,sp)
    use grid_class
    use ptl_module
    use sml_module
    implicit none
    integer :: np
    type(grid_type) :: grid
    type(species_type) :: sp
    logical, save :: first=.true.
    integer, parameter :: ict1=ptl_nphase+1
    integer, parameter :: ict2=ptl_nphase+ptl_nconst

    if(first)then
      call source_2d_init(grid)
      first=.false.
    endif

    np=sp%num
    allocate(source_ptl%phase(ict2,np),source_ptl%gid(np)) 
    df_src=0D0
    f_src=0D0
    f0_grid_src=0D0
    dn_src=0D0
    n_src=0D0
    g_local_rev=0D0
    n_profile_list=0D0
    t_profile_list=0D0
    g_local=0D0
    if(sml_mype==0)then
      if (sp%type==1) then ! ions
         print *, 'ions'
      elseif (sp%type==0) then ! electrons
         print *, 'electrons'
         stop
      endif
    endif
    !generate df_src, dn_src, n_src
    call generate_src(grid,sp)
    !apply source to weight
    !call apply_src(grid,sp)

  end subroutine source_2d

  !geneate nesseary source terms: df_src, dn_src, n_src
  subroutine generate_src(grid,sp)
    use grid_class
    use ptl_module
    use eq_module
    use diag_module
    use sml_module
    use mpi
    implicit none
    type(grid_type), intent(in) :: grid
    type(ptl_type_aosoa) :: ptli
    type(species_type) :: sp
    integer :: i,j,m,ip,ie,isize,ierr
    real (8) :: r,z,phi,B,psi
    real (8) :: wpsi,we,pn,deltaf,mu,rho,ter,en_nor,en_ev,w00,w01,w10,w11,g_rev
    real (8) :: n_profile,dn_tmp,n_tmp,f_tmp,g_tmp,df_tmp,vfrac1,vfrac,vfrac2,v,dv,d_tmp
    real (8), dimension(:), allocatable :: psi_tmp
    real (8), dimension(:,:), allocatable :: df_src_tmp,f_src_tmp,g_src_tmp,g_local_rev_tmp,information_save
    integer, dimension(:,:), allocatable :: slice_save
    logical, dimension(:), allocatable :: particle_save
    integer, parameter :: sp_type=1 !ion
    integer, parameter :: ict1=ptl_nphase+1
    integer, parameter :: ict2=ptl_nphase+ptl_nconst
    real (8), external :: b_interpol,psi_interpol

    allocate(df_src_tmp(num_psi_src,num_energy_src),f_src_tmp(num_psi_src,num_energy_src))
    allocate(g_src_tmp(num_psi_src,num_energy_src))
    allocate(g_local_rev_tmp(num_psi_src,num_energy_src))
    length=sp%num
    !allocate(data_source(1:sp%num,26))
    allocate(information_save(1:sp%num,4))
    allocate(slice_save(1:sp%num,2))
    allocate(particle_save(1:sp%num))
    information_save(:,:)=0D0
    slice_save(:,:)=-1
    particle_save(:)=.false.
    !data_source(:,:)=0D0

    allocate(psi_tmp(num_psi_src))
    psi_tmp=grid%psi_surf2

    df_src_tmp(:,:)=0D0
    f_src_tmp(:,:)=0D0
    g_src_tmp(:,:)=0D0
    g_local_rev_tmp(:,:)=0D0
    vol_local(:,:)=0.0

    vol_local_space=diag_1d_vol
    do m=1,num_psi_src
      psi=psi_tmp(m)
      ter=eq_ftn(psi,r,z,eq_temp(sp_type))
      do j=1,num_energy_src
         vfrac1=j*denergy_src*ter
         vfrac=(j+0.5)*denergy_src*ter
         vfrac2=(j+1)*denergy_src*ter
         v=sqrt(2.0*vfrac)
         dv=sqrt(2.0*vfrac2)-sqrt(2.0*vfrac1)
         !v=sqrt(2.0*vfrac*sml_ev2j/ptl_mass(sp_type))
         !dv=sqrt(2.0*vfrac2*sml_ev2j/ptl_mass(sp_type))-sqrt(2.0*vfrac1*sml_ev2j/ptl_mass(sp_type))
         vol_local(m,j)=1.0/(vol_local_space(m)*4*sml_pi*dv*v**2)*(2*sml_pi)**1.5!*(1/(2*sml_pi*vfrac)**1.5)
      enddo
     enddo

    !calculate delta f in psi and energy space
    !!$omp parallel do
    do i=1, sp%num
       ptli=sp%ptl(SIND(i))
       if(ptli%gid(AIND(i))<=0) cycle

       r=ptli%ph(AIND(i),pir) !R
       z=ptli%ph(AIND(i),piz) !Z
       phi=ptli%ph(AIND(i),pip) !phi
       psi=psi_interpol(r,z,0,0)
       if(psi<0)then
         print *,'print psi in source', psi
         cycle
       endif
       B=b_interpol(r,z,phi) !B field

       !psi=psi_interpol(r,z,0,0) !psi
       
       if(.not. is_rgn12(r,z,psi) .or. z < eq_x_z .or. z > eq_x2_z) cycle ! skip for private region of lower divertor
       call search_psi(grid,psi,wpsi,ip)
       if(ip <0 .or. wpsi <0) cycle
       !CONVERT_GRID2 is asumed, but kept the other for future extension.
!#ifdef CONVERT_GRID2
!       call search_ptl_1d(grid,psi,wpsi,ip)
!       if (ip .lt. 1) cycle
!#else
!       pn=(psi-diag_1d_pin)*diag_1d_dp_inv
!       ip=floor(pn)+1
!       if(ip <1 .or. diag_1d_npsi <= ip) cycle    ! EXIT for out of diag_1d_pin/pout range
!       wpsi=1D0 - pn + real(ip-1,8)
!#endif
       
       mu=ptli%ct(AIND(i),pim)
       rho=ptli%ph(AIND(i),pirho)

       en_ev=(abs(mu*B) + 0.5D0*(ptl_charge(sp_type)*rho*B)**2/ptl_mass(sp_type))*sml_j2ev
       ter=eq_ftn(psi,r,z,eq_temp(sp_type))
       en_nor=en_ev/ter

       !data_source(i,1)=r
       !data_source(i,2)=z
       !data_source(i,3)=phi
       !data_source(i,4)=psi
       !data_source(i,5)=B
       !data_source(i,6)=rho
       !data_source(i,7)=mu
       !data_source(i,8)=ter
       !data_source(i,9)=en_nor
       !data_source(i,10)=en_ev
       !data_source(i,11)=ptli%ph(AIND(i),piw1)
       !data_source(i,12)=ptli%ph(AIND(i),pif0)
       !data_source(i,13)=ptli%ph(AIND(i),piw0)
       !data_source(i,19)=en_nor/denergy_src
       !data_source(i,20)=denergy_src
       ie=max(floor(en_nor/denergy_src),0)
       !data_source(i,17)=real(ie)
       !data_source(i,16)=we
       if(ie>(num_energy_src-2) .or. ie <0)cycle
       if(en_nor .gt. cut_energy_src) cycle
       
       we=1D0-(en_nor/denergy_src-real(ie))
       ie=ie+1
       if(we<0 .or. wpsi <0)cycle
       !data_source(i,15)=wpsi
       !data_source(i,18)=real(ip)

       !data_source(i,21)=ptli%ph(AIND(i),piw1)*ptli%ph(AIND(i),pif0)/ptli%ph(AIND(i),piw0)
       if(ie .le. 0)cycle
       if(ip .le. 0)cycle
       if(ip+1 .gt. num_psi_src)cycle       

       particle_save(i)=.true.
       slice_save(i,1)=ip
       slice_save(i,2)=ie
       information_save(i,1)=wpsi
       information_save(i,2)=we
       information_save(i,3)=ter
       information_save(i,4)=en_nor 
       w00=wpsi*we
       w01=wpsi*(1-we)
       w10=(1-wpsi)*we
       w11=(1-wpsi)*(1-we)
      
       !deltaf/g
       deltaf=ptli%ph(AIND(i),piw1)*ptli%ct(AIND(i),piw0)
       !if(abs(deltaf)>1)cycle
       !data_source(i,14)=deltaf
       !!$omp atomic update 
       df_src_tmp(ip,ie)=df_src_tmp(ip,ie)+deltaf*w00
       !!$omp atomic update
       df_src_tmp(ip,ie+1)=df_src_tmp(ip,ie+1)+deltaf*w01
       !!$omp atomic update 
       df_src_tmp(ip+1,ie)=df_src_tmp(ip+1,ie)+deltaf*w10
       !!$omp atomic update
       df_src_tmp(ip+1,ie+1)=df_src_tmp(ip+1,ie+1)+deltaf*w11


       !f_0/g
       deltaf=ptli%ct(AIND(i),piw0)
       f_src_tmp(ip,ie)=f_src_tmp(ip,ie)+deltaf*w00
       !!$omp atomic update
       f_src_tmp(ip,ie+1)=f_src_tmp(ip,ie+1)+deltaf*w01
       !!$omp atomic update 
       f_src_tmp(ip+1,ie)=f_src_tmp(ip+1,ie)+deltaf*w10
       !!$omp atomic update
       f_src_tmp(ip+1,ie+1)=f_src_tmp(ip+1,ie+1)+deltaf*w11

       !!$omp atomic update 
       !g_src_tmp(ip,ie)=g_src_tmp(ip,ie)+w00
       !!$omp atomic update
       !g_src_tmp(ip,ie+1)=g_src_tmp(ip,ie+1)+w01
       !!$omp atomic update
       !g_src_tmp(ip+1,ie)=g_src_tmp(ip+1,ie)+w10
       !!$omp atomic update
       !g_src_tmp(ip+1,ie+1)=g_src_tmp(ip+1,ie+1)+w11

       !g_rev=ptli%ph(AIND(i),pif0)/ptli%ph(AIND(i),piw0)
       !g_local_rev_tmp(ip,ie)=g_local_rev_tmp(ip,ie)+g_rev*w00
       !g_local_rev_tmp(ip,ie+1)=g_local_rev_tmp(ip,ie+1)+g_rev*w01
       !g_local_rev_tmp(ip+1,ie)=g_local_rev_tmp(ip+1,ie)+g_rev*w10
       !g_local_rev_tmp(ip+1,ie+1)=g_local_rev_tmp(ip+1,ie+1)+g_rev*w11
    enddo

    !reduction of df_src_tmp and g_src_tmp
    isize=num_psi_src*num_energy_src
    call mpi_allreduce(MPI_IN_PLACE,df_src_tmp,isize,MPI_REAL8,MPI_SUM,sml_comm,ierr)
    call mpi_allreduce(MPI_IN_PLACE,f_src_tmp,isize,MPI_REAL8,MPI_SUM,sml_comm,ierr)
    !call mpi_allreduce(MPI_IN_PLACE,g_src_tmp,isize,MPI_REAL8,MPI_SUM,sml_comm,ierr)
    !call mpi_allreduce(MPI_IN_PLACE,g_local_rev_tmp,isize,MPI_REAL8,MPI_SUM,sml_comm,ierr)

    !calculate df_src
    !!$omp parallel do collapse(2) 
    do i=1,num_psi_src
       do j=1,num_energy_src
          !if(g_src_tmp(i,j)>0.0)then
            df_src(i,j)=df_src_tmp(i,j)*vol_local(i,j)!/g_src_tmp(i,j)
            f_src(i,j)=f_src_tmp(i,j)*vol_local(i,j)
            !g_local_rev(i,j)=g_local_rev_tmp(i,j)/g_src_tmp(i,j)
          !endif
       enddo
    enddo
    !g_local=g_src_tmp
    !allocate(psi_tmp(num_psi_src))
    !psi_tmp=grid%psi_surf2
    r=0D0
    z=0D0
    !!$omp parallel do
    do m=1,num_psi_src
      psi=psi_tmp(m)
      ter=eq_ftn(psi,r,z,eq_temp(sp_type))
      n_profile=eq_ftn(psi,r,z,eq_den(sp_type))
      n_profile_list(m)=n_profile
      t_profile_list(m)=ter
      dn_tmp=0.0
      n_tmp=0.0
      do j=1,num_energy_src
         vfrac1=j*denergy_src*ter
         vfrac=(j+0.5)*denergy_src*ter
         vfrac2=(j+1)*denergy_src*ter
         !v=sqrt(2.0*vfrac*sml_ev2j/ptl_mass(sp_type))
         !dv=sqrt(2.0*vfrac2*sml_ev2j/ptl_mass(sp_type))-sqrt(2.0*vfrac1*sml_ev2j/ptl_mass(sp_type))
         v=sqrt(2.0*vfrac)
         dv=sqrt(2.0*vfrac2)-sqrt(2.0*vfrac1)
         f_tmp=n_profile*sqrt(1D0/ter**3)*exp(-vfrac/ter)
         f0_grid_src(m,j)=f_tmp
         !if(g_local_rev(m,j)==0)cycle
         dn_tmp=dn_tmp+df_src(m,j)*4*sml_pi*dv*v**2/(2*sml_pi)**1.5
         n_tmp=n_tmp+f_tmp*4*sml_pi*dv*v**2/(2*sml_pi)**1.5
      enddo
      dn_src(m)=dn_tmp
      n_src(m)=n_tmp 
    enddo

    !!$omp parallel do
    do i=1, sp%num
       ptli=sp%ptl(SIND(i))
       if(ptli%gid(AIND(i))<=0 .or. particle_save(i)==.false.) cycle
       
       ip=slice_save(i,1)
       ie=slice_save(i,2)
       
       wpsi=information_save(i,1)
       we=information_save(i,2)
       ter=information_save(i,3)
       en_nor=information_save(i,4)

       w00=wpsi*we
       w01=wpsi*(1-we)
       w10=(1-wpsi)*we
       w11=(1-wpsi)*(1-we)
       n_profile=eq_ftn(psi,r,z,eq_den(sp_type))

       df_tmp=df_src(ip,ie)*w00+df_src(ip+1,ie)*w10+df_src(ip,ie+1)*w01+df_src(ip+1,ie+1)*w11
       g_rev=g_local_rev(ip,ie)*w00+g_local_rev(ip+1,ie)*w10+g_local_rev(ip,ie+1)*w01+g_local_rev(ip+1,ie+1)*w11
       if(g_rev==0)cycle
       dn_tmp=dn_src(ip)*wpsi+dn_src(ip+1)*(1-wpsi)
       n_tmp=n_src(ip)*wpsi+n_src(ip+1)*(1-wpsi)
       f_tmp=n_profile*sqrt(1D0/ter**3)*exp(-en_nor)
       !data_source(i,22)=ptli%ph(AIND(i),piw1)
       !data_source(i,23)=gamma_src*(df_tmp-f_tmp*dn_tmp/n_tmp)/(ptli%ph(AIND(i),pif0)/ptli%ph(AIND(i),piw0))
       sp%ptl(SIND(i))%ph(AIND(i),piw1)=sp%ptl(SIND(i))%ph(AIND(i),piw1)-gamma_src*(df_tmp/(ptli%ct(AIND(i),pif0))-f_tmp/ptli%ct(AIND(i),pif0)*dn_tmp/n_tmp)
       !sp%ptl(SIND(i))%ph(AIND(i),piw2)=0.0
       !data_source(i,24)=ptli%ph(AIND(i),piw1)
       !data_source(i,25)=ptli%ph(AIND(i),piw2)
       !data_source(i,26)=gamma_src*(df_tmp-f_tmp*dn_tmp/n_tmp)
    enddo

    call output_src

    do i=1, sp%num
       source_ptl%phase(1:ptl_nphase,i)=sp%ptl(SIND(i))%ph(AIND(i),:)
       source_ptl%phase(ict1:ict2,i)=sp%ptl(SIND(i))%ct(AIND(i),:)
       source_ptl%gid(i)=sp%ptl(SIND(i))%gid(AIND(i))
    enddo       
    deallocate(psi_tmp)
    deallocate(df_src_tmp)
    deallocate(g_src_tmp)
    deallocate(g_local_rev_tmp)
    !deallocate(data_source)
    deallocate(information_save)
    deallocate(slice_save)
    deallocate(particle_save) 
        
  end subroutine generate_src

  subroutine transfer_source_ptl(isp) bind(C,name="transfer_source_ptl")
    use ptl_module
    implicit none
    integer, intent(in), value :: isp
    integer :: i
    integer, parameter :: ict1=ptl_nphase+1
    integer, parameter :: ict2=ptl_nphase+ptl_nconst
    do i=1, spall_global(isp)%num
      spall_global(isp)%ptl(SIND(i))%ph(AIND(i),:)=source_ptl%phase(1:ptl_nphase,i)
      spall_global(isp)%ptl(SIND(i))%ct(AIND(i),:)=source_ptl%phase(ict1:ict2   ,i)
      spall_global(isp)%ptl(SIND(i))%gid(AIND(i))=source_ptl%gid(i)
    enddo
    deallocate(source_ptl%phase)
    deallocate(source_ptl%gid)
  end subroutine transfer_source_ptl
 
  subroutine output_src
    use sml_module
    use adios2_comm_module
    implicit none

    integer :: ierr

    type(adios2_io), save :: io
    type(adios2_engine), save :: engine
    type(adios2_variable) :: varid
    logical, save :: init=.true.
    integer (kind=8), dimension(2) :: gdims_2d, goffset_2d, ldims_2d
    integer (kind=8), dimension(1) :: gdims_1d, goffset_1d, ldims_1d
    character (len=256) :: filename
    if(sml_mype==0)then

      write(filename,'("source_2d/xgc.source",".",i5.5,".bp")') sml_gstep
      if(init)then
        init=.false.
        gdims_2d(1)=1_8*num_psi_src
        gdims_2d(2)=1_8*num_energy_src
        goffset_2d(1)=0_8
        goffset_2d(2)=0_8
        ldims_2d(1)=1_8*num_psi_src
        ldims_2d(2)=1_8*num_energy_src

        call adios2_declare_io(io,adios2obj,'source_2d',ierr)
        call adios2_define_variable(varid,io,'df_src',adios2_type_dp,2,gdims_2d,goffset_2d,ldims_2d,adios2_constant_dims,ierr)
        call adios2_define_variable(varid,io,'g_local_rev',adios2_type_dp,2,gdims_2d,goffset_2d,ldims_2d,adios2_constant_dims,ierr)
        call adios2_define_variable(varid,io,'f0_grid_src',adios2_type_dp,2,gdims_2d,goffset_2d,ldims_2d,adios2_constant_dims,ierr)
        call adios2_define_variable(varid,io,'f_src',adios2_type_dp,2,gdims_2d,goffset_2d,ldims_2d,adios2_constant_dims,ierr)

        gdims_1d(1)=1_8*num_psi_src
        goffset_1d(1)=0_8
        ldims_1d(1)=1_8*num_psi_src
        
        call adios2_define_variable(varid,io,'dn_src',adios2_type_dp,1,gdims_1d,goffset_1d,ldims_1d,adios2_constant_dims,ierr)
        call adios2_define_variable(varid,io,'n_src',adios2_type_dp,1,gdims_1d,goffset_1d,ldims_1d,adios2_constant_dims,ierr)
        call adios2_define_variable(varid,io,'n_profile_list',adios2_type_dp,1,gdims_1d,goffset_1d,ldims_1d,adios2_constant_dims,ierr)
        call adios2_define_variable(varid,io,'t_profile_list',adios2_type_dp,1,gdims_1d,goffset_1d,ldims_1d,adios2_constant_dims,ierr)
        !write data
      endif
  
        !gdims_2d(1)=1_8*length
        !gdims_2d(2)=26_8
        !goffset_2d(1)=0_8
        !goffset_2d(2)=0_8
        !ldims_2d(1)=1_8*length
        !ldims_2d(2)=26_8
        !call adios2_define_variable(varid,io,'data_source',adios2_type_dp,2,gdims_2d,goffset_2d,ldims_2d,adios2_constant_dims,ierr)
 
      call adios2_open(engine,io,filename,adios2_mode_write,sml_comm_self,ierr)

      call adios2_begin_step(engine, adios2_step_mode_append, ierr)

      call adios2_put(engine, "df_src", df_src, ierr)
      call adios2_put(engine, "g_local_rev", g_local_rev, ierr)
      call adios2_put(engine, "f0_grid_src", f0_grid_src, ierr)
      call adios2_put(engine, "f_src", f_src, ierr)
      !call adios2_put(engine, "data_source", data_source, ierr)
      call adios2_put(engine, "dn_src", dn_src, ierr)
      call adios2_put(engine, "n_src", n_src, ierr)
      call adios2_put(engine, "n_profile_list", n_profile_list, ierr)
      call adios2_put(engine, "t_profile_list", t_profile_list, ierr)

      call adios2_end_step(engine, ierr)
      call adios2_close(engine,ierr)
    endif
  end subroutine output_src 

  subroutine search_psi(grid,psi,wpsi,ip)
    use eq_module
    use grid_class
    use sml_module
    implicit none
    type(grid_type), intent(in) :: grid
    real (kind=8), intent(in) :: psi
    real (kind=8), intent(out) :: wpsi
    integer, intent(out) :: ip
    real (kind=8) :: psi_list(num_psi_src)
    integer :: i
    
    ip=-1
    wpsi=-1D0
    psi_list=grid%psi_surf2
    do i=1,num_psi_src-1
       if(psi>sml_inpsi .and. psi<sml_outpsi)then
         if(psi>=psi_list(i) .and. psi<psi_list(i+1))then
           ip=i
           wpsi=(psi_list(i+1)-psi)/(psi_list(i+1)-psi_list(i)) 
         endif
       endif
    enddo 
    
  end subroutine search_psi
end module source_2d
