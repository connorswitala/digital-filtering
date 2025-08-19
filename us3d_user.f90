
      module userdata

      implicit none
      !save

      real(8), parameter :: pi = 4*atan(1.0d0)
      real(8), parameter :: r0 = 8.8035d-2
      real(8), parameter :: t0 = 216.65d0
      real(8), parameter :: u0 = 2950.7d0
      real(8), parameter :: amp = 1.0d-3
      real(8), parameter :: dx0 = 0.0001d0 !!<<<<<<CHECK
      real(8), parameter :: mindiss   = 0.0d0
      real(8), parameter :: x_sponge  = 0.95d0
      real(8), parameter :: theta = -45.0d0 * pi/180.0d0
      real(8) :: p0, a0

      contains

      !---------------------------------------
      subroutine my_user_init(cmod,iop,ier)
      use sizing, only : net
      use geometry, only : xce
      use simvars
      use restartio
      use us3d_core_module, only: update_dfromp
      character(LEN=*), intent(IN) :: cmod
      integer, intent(IN) :: iop
      integer, intent(OUT) :: ier
      integer :: i
      
      p0 = r0*t0*gsp(1)
      a0 = sqrt(1.4*p0/r0)

      do i=1,net
            u(i) = u0
            v(i) = 0.0d0  
            w(i) = 0.0d0
            t(i) = t0
            r(i) = r0
            cs(1,i) = 1.0d0
            p(i) = p0
      enddo

      call update_dfromp(1,net)
      call dataio("w", ier)
      return
      end subroutine my_user_init

      !---------------------------------------
      subroutine my_user_main_pre(ier)
         use geometry, only : xce, sx
         use simvars
         use sizing
         use models
         use connect
         use us3d_core_module, only: update_dfromp
         implicit none 
         integer, intent(OUT) :: ier
         real*8 :: p_pert,r_pert,u_pert,T_pert
         integer :: i,ii, j,k,ibc,bcn,zn,j1,j2
         real(8) :: omega, lambda, fmax, yc, ramp, width
         real(8) :: ffreq,per,tol,pw         
 
         ier = 0
         fmax = 718.94d3 !frequency of unstable wavelength
         lambda = (u0-sqrt(1.4d0*p0/r0))/fmax
         !print *, "Lambda", lambda
         omega  = 2.0d0*pi * fmax

         yc = 0.0145256 !center of pulse 
         width = 1e-1 !width of 99% drop in magnitude
         ffreq = 1.0d3 !forcing frequency
         per = 1.0d0/ffreq
         tol = 1.0d-8
         pw = 2.0d-7 !pulse width
        
         !print *,"Freq", u0/lambda

         do k=1,nft-1   
            ibc = jbcmap(k)
            bcn = bcs%bcn(ibc)
            zn = bcs%zn(ibc)
            ! Check for Inflow conditions
            if ((bcn.eq.10).or.(bcn.eq.20).or.(bcn.eq.30)) then !inflow
               
               ! provide the current rank               
               j1 = jmaster(k,1)    ! Starting and ending face range
               j2 = jmaster(k,2)
               do j=j1,j2
                  i     = ife(j,1) !inner cell
                  ii    = ife(j,2) !ghost cell
                  if (MODULO(sim_time,per)<(tol+pw)) then
                  ramp = amp * exp(-1.0*((xce(2,ii) - yc)**2)/ &
                                   (2*(width/3.03485)**2))
                  else
                  ramp = 0.0d0
                  endif 
                  p_pert =ramp*sin(2*pi/lambda*(xce(1,ii)*cos(theta) + &
                          xce(2,ii)*sin(theta)) - omega*sim_time)
                  !print *, p_pert, ramp, xce(2,ii)! pi, lambda, omega, theta
                  r_pert = (1/a0)**2 * p_pert
                  u_pert = -1.0d0/r0/a0 * p_pert
                  T_pert = (1.4-1.0)*t0/r0/a0/a0 * p_pert
                 ! if (xce(2,ii).gt.0.003.and.xce(2,ii).lt.0.007) then
                  r(ii) = r0 + r_pert
                  !print *,"r=", r(ii)
                  u(ii) = u0 + u_pert
                  !print *,"u=", u(ii)
                  v(ii) = 0.0d0
                  w(ii) = 0.0d0
                  p(ii) = p0 + p_pert
                  !print *,"p=", p(ii)
                  t(ii) = t0 + T_pert
                  !print *,"t=", t(ii)
                 ! endif

               enddo
            endif
         enddo
 
         ! open(100,file='data.dat',status='unknown')
         ! do i=1,nel
         !    write(100,'(99e20.8)') xce(1,i),r(i),u(i),v(i),w(i),t(i),p(i)
         ! enddo  
         ! close(100)

 
         return
      end subroutine my_user_main_pre

      subroutine my_user_dswitch(idgo,ibc,ibt,j,ml,mr,divl,divr,dfac)
         use geometry, only : xcf,voli
         use models, only: ns
         use simvars
         use mpivars
         use sizing
         use connect
 
         implicit none
         integer, intent(IN)  :: idgo,ibc,ibt,j
         Real(8), intent(IN)  :: ml,mr,divl,divr
         Real(8), dimension(:), intent(OUT) :: dfac
         Real(8) :: div,vort(3),alpha,eps_sw,beta,vmag,div2,vmag2,al,ar
         integer :: i,ii 
           
         eps_sw = 0.1d0 
         i = ife(j,1)
         ii = ife(j,2) 
       
         !left
         div = grad(1,ns+1,i) + grad(2,ns+2,i) + grad(3,ns+3,i)
         div2=div*div
         vort(1) = grad(2,ns+3,i) - grad(3,ns+2,i) ! dw/dy - dv/dz
         vort(2) = grad(3,ns+1,i) - grad(1,ns+3,i) ! du/dz - dw/dx
         vort(3) = grad(1,ns+2,i) - grad(2,ns+1,i) ! dv/dx - du/dy
         vmag2 = dot_product(vort,vort)
            
         vmag = sqrt(u(i)**2 + v(i)**2 + w(i)**2)
         beta = div2/(div2+(0.005*vmag*voli(i))**2 + 0.01)
         al = -div/(sqrt(vmag2) + eps_sw * u0/dx0)

         !right
         div = grad(1,ns+1,ii) + grad(2,ns+2,ii) + grad(3,ns+3,ii)
         div2=div*div
         vort(1) = grad(2,ns+3,ii) - grad(3,ns+2,ii) ! dw/dy - dv/dz
         vort(2) = grad(3,ns+1,ii) - grad(1,ns+3,ii) ! du/dz - dw/dx
         vort(3) = grad(1,ns+2,ii) - grad(2,ns+1,ii) ! dv/dx - du/dy
         vmag2 = dot_product(vort,vort)

         vmag = sqrt(u(ii)**2 + v(ii)**2 + w(ii)**2)
         beta = div2/(div2+(0.005*vmag*voli(ii))**2 + 0.01)
         ar = -div/(sqrt(vmag2) + eps_sw * u0/dx0)

         alpha = max(al,ar,0.0d0)
         dfac(:) = min(alpha*beta, 1.0d0)  + mindiss
         
         return
       end subroutine my_user_dswitch

      end module userdata
 
      !-------------------------
      subroutine user_initialize(component,debug,ier)
      use us3d_user_module, only : user_init, &
                                  user_main_pre, user_dswitch
      use userdata

      implicit none
      character(LEN=*), intent(IN) :: component
      logical, intent(IN) :: debug
      integer, intent(OUT) :: ier

      ier = 0
      !call get_ppw(ppw)
      user_init      => my_user_init
      user_main_pre  => my_user_main_pre
      user_dswitch   => my_user_dswitch

      return
      end subroutine user_initialize
