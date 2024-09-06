program main
  use real_precision
  use utils
  use diaglib
!
! simple minded test program to call the lobpcg routine
!
  implicit none
  integer  :: n, itmax, m_max, n_want, iwhat
  real(dp) :: tol
!
! initialize:
!
  n      = 1000 
  n_want = 20  
  tol    = 1.0e-6_dp
  itmax  = 1000
  m_max  = 20
  nmult  = 0
  tdscf  = .false.
  i_alg  = 0
!
  if ((m_max*n_want) .gt. n) then
       write(6,*)
       write(6,*) 'WARNING'
       write(6,*) 'Rouché-Capelli is encountered, please change your parameters'
       write(6,*)
       stop
  end if
!
  call print_header
! 
  write(6,*)
  1000 format(t3,' simple-minded test driver. select ',/, &
              t3,'   1 for symmetric eigenvalue problems ',/, &
              t3,'   2 for symmetric generalized eigenvalue problems ',/, &
              t3,'   3 for linear-response equations (SCF-like)',/, &
              t3,'   4 for linear-response equations (CASSCF-like).',/, &
              t3,'   5 for standard-response equations (CASSCF-like).')
  write(6,1000)
!  read(5,*) iwhat
  iwhat = 5
  write(6,*)
!
  if (iwhat.eq.1) then 
    call test_symm(.true.,n,n_want,tol,itmax,m_max)
  else if (iwhat.eq.2) then 
    call test_geneig(.true.,n,n_want,tol,itmax,m_max)
  else if (iwhat.eq.3) then 
    call test_scflr(.true.,n,n_want,tol,itmax,m_max)
  else if (iwhat.eq.4) then 
    call test_caslr(.false.,n,n_want,tol,itmax,m_max)
  else if (iwhat.eq.5) then
    call test_caslr_std(.true.,n,n_want,tol,itmax,m_max)  
  else 
    write(6,*) ' invalid selection. aborting ...'
  end if
!
end program main
!
  subroutine print_header
    implicit none
!
    1000 format( &
  ' diaglib - a fortran library of matrix-free iterative algorithms to ',/, &
  ' compute a few eigenvalues and eigenvectors of large matrices.',/, &
  ' ==================================================================',/, &
  /, &
  ' Implementation by',/, &
  /, &
  '   Ivan Giannì, Tommaso Nottoli, Federica Pes, Antoine Levitt ',/, &
  '   and Filippo Lipparini',/, &
  '   MoLECoLab Pisa',/, &
  '   Department of Chemistry and Industrial Chemistry',/, &
  '   University of Pisa',/, &
  '   Via G. Moruzzi 13, I-56124, Pisa, Italy')
    write(6,1000)
    return
  end subroutine print_header
!
  subroutine mmult(n,m,x,ax)
    use utils
    implicit none
!
!   simple-minded matrix-vector multiplication routine, that needs to be passed as 
!   an argument to lobpcg
!
    integer,                   intent(in)    :: n, m
    real(dp),  dimension(n,m), intent(in)    :: x
    real(dp),  dimension(n,m), intent(inout) :: ax
!
    integer :: icol
!
    nmult = nmult + m
    do icol = 1, m
      ax(:,icol) = matmul(a,x(:,icol))
    end do
    return
  end subroutine mmult
!
  subroutine smult(n,m,x,sx)
    use utils
    implicit none
!
!   simple-minded matrix-vector multiplication routine, that needs to be passed as 
!   an argument to lobpcg
!
    integer,                   intent(in)    :: n, m
    real(dp),  dimension(n,m), intent(in)    :: x
    real(dp),  dimension(n,m), intent(inout) :: sx
!
    integer :: icol
!
!
!
    if (tdscf) then
!
!     nothing to do, as s is the identity. not the most elegant way of doing this,
!     but it works.
!
      sx = x
      return
    end if
!
    nmult = nmult + m
    do icol = 1, m
      sx(:,icol) = matmul(s,x(:,icol))
    end do
    return
  end subroutine smult
!
  subroutine mprec(n,m,fac,x,px)
    use utils
    implicit none
!
!   simple-minded shift-and-invert preconditioner routine, that needs to be passed as 
!   an argument to lobpcg
!
    integer,                   intent(in)    :: n, m
    real(dp),                  intent(in)    :: fac
    real(dp),  dimension(n,m), intent(in)    :: x
    real(dp),  dimension(n,m), intent(inout) :: px
!
    integer             :: i, icol
    real(dp), parameter :: tol = 1.0d-5
!
    do icol = 1, m
      do i = 1, n
        if (abs(a(i,i)+fac).gt.tol) px(i,icol) = x(i,icol)/(a(i,i) + fac)
      end do
    end do
    return
  end subroutine mprec
!
  subroutine apbvec(n,m,x,y)
    use utils
    implicit none
    integer,                  intent(in)    :: n, m
    real(dp), dimension(n,m), intent(in)    :: x
    real(dp), dimension(n,m), intent(inout) :: y
!
    integer icol
!
    do icol = 1, m
      y(:,icol) = matmul(apb,x(:,icol))
    end do
!
    return
  end subroutine apbvec
  subroutine ambvec(n,m,x,y)
    use utils
    implicit none
    integer,                  intent(in)    :: n, m
    real(dp), dimension(n,m), intent(in)    :: x
    real(dp), dimension(n,m), intent(inout) :: y
!
    integer icol
!
    do icol = 1, m
      y(:,icol) = matmul(amb,x(:,icol))
    end do
!
    return
  end subroutine ambvec
  subroutine spdvec(n,m,x,y)
    use utils
    implicit none
    integer,                  intent(in)    :: n, m
    real(dp), dimension(n,m), intent(in)    :: x
    real(dp), dimension(n,m), intent(inout) :: y
!
    integer icol
!
    do icol = 1, m
      y(:,icol) = matmul(spd,x(:,icol))
    end do
!
    return
  end subroutine spdvec
  subroutine smdvec(n,m,x,y)
    use utils
    implicit none
    integer,                  intent(in)    :: n, m
    real(dp), dimension(n,m), intent(in)    :: x
    real(dp), dimension(n,m), intent(inout) :: y
!
    integer icol
!
    do icol = 1, m
      y(:,icol) = matmul(smd,x(:,icol))
    end do
!
    return
  end subroutine smdvec
!
  subroutine lrprec_1(n,m,fac,xp,xm,yp,ym)
    use utils
    implicit none
    integer,                  intent(in)    :: n, m
    real(dp),                 intent(in)    :: fac
    real(dp), dimension(n,m), intent(in)    :: xp, xm
    real(dp), dimension(n,m), intent(inout) :: yp, ym
!
    integer :: i, icol
!
    yp = 0.0d0
    ym = 0.0d0
!
    do icol = 1, m
      do i = 1, n
        yp(i,icol) = - (1.0d0/(aa(i,i)**2 - fac**2 * sigma(i,i)**2)) * (aa(i,i) * xp(i,icol) + fac * sigma(i,i) * xm(i,icol))
        ym(i,icol) = - (1.0d0/(aa(i,i)**2 - fac**2 * sigma(i,i)**2)) * (aa(i,i) * xm(i,icol) + fac * sigma(i,i) * xp(i,icol))
      end do
    end do
!
    return
  end subroutine lrprec_1
!
  subroutine lrprec_2(n,m,fac,xp,xm,yp,ym)
    use utils
    implicit none
    integer,                  intent(in)    :: n, m
    real(dp),                 intent(in)    :: fac
    real(dp), dimension(n,m), intent(in)    :: xp, xm
    real(dp), dimension(n,m), intent(inout) :: yp, ym
!
    integer  :: i, icol
    real(dp) :: denom
!
    yp = 0.0d0
    ym = 0.0d0
!
    do icol = 1, m
      do i = 1, n
        denom = fac * fac * aa(i,i) * aa(i,i) - sigma(i,i) * sigma(i,i)
        denom = 1.0d0/denom
        yp(i,icol) = denom * (fac * aa(i,i) * xp(i,icol) + sigma(i,i) * xm(i,icol))
        ym(i,icol) = denom * (fac * aa(i,i) * xm(i,icol) + sigma(i,i) * xp(i,icol))
      end do
    end do
!
    return
  end subroutine lrprec_2
!
  subroutine test_symm(check_lapack,n,n_want,tol,itmax,m_max)
    use real_precision
    use utils
    use diaglib, only : lobpcg_driver, davidson_driver
    implicit none
    logical,  intent(in) :: check_lapack
    integer,  intent(in) :: n, n_want, itmax, m_max
    real(dp), intent(in) :: tol
!
!   test iterative solvers for a standard symmetric eigenvalue problem.
!
    logical  :: ok
    integer  :: i, j, n_eig, info, lwork
    real(dp) :: lw(1)
    real(dp), allocatable :: diagonal(:), eig(:), evec(:,:)
    real(dp), allocatable :: work(:), w(:), a_copy(:,:)
!
    external :: mmult, mprec, smult
!
    1000 format(t3,' eigenvalue # ',i6,': ',f12.6,/,t3,' eigenvector: ')
!
!
!   allocate memory for the a matrix 
!
    allocate (a(n,n))
!
!   build a symmetric matrix:
!
    do i = 1, n
      a(i,i) = real(i,kind=dp) + 1.0_dp
      do j = 1, i - 1
        a(j,i) = 1.0_dp/real(i+j,kind=dp)
        a(i,j) = a(j,i)
      end do
    end do
!
!   if required, solve the problem with a dense lapack routine:
!
    if (check_lapack) then
      allocate (a_copy(n,n), w(n))
      a_copy = a
      call dsyev('v','l',n,a_copy,n,w,lw,-1,info)
      lwork = int(lw(1))
      allocate (work(lwork))
      call dsyev('v','l',n,a_copy,n,w,work,lwork,info)
!
!     write the results on file for comparison:
!
      open (unit = 10, file = 'lapack.txt', form = 'formatted', access = 'sequential')
      do i = 1, n_want
        write(10,1000) i, w(i)
!
!       fix the phase so that the first element of the eigenvector is positive.
!
        if (a_copy(1,i) .lt. 0.0_dp) a_copy(:,i) = - a_copy(:,i)
        write(10,'(10f12.6)') a_copy(:,i)
        write(10,*)
      end do
      close (10)
    end if
!
!   allocate and gather the diagonal:
!
    allocate (diagonal(n))
    do i = 1, n
      diagonal(i) = a(i,i)
    end do
!
!   for better convergence, we seek more eigenpairs and stop the iterations when the
!   required ones are converged.
!
    n_eig = min(2*n_want, n_want + 5)
!
!   allocate memory for the eigenvalues and eigenvectors:
!
    allocate (eig(n), evec(n,n_eig))
!
!   compute a guess for the eigenvector (see guess_evec for more information)
!
    call guess_evec(4,n,n_eig,diagonal,evec)
!
!   call lobpcg:
!
    call lobpcg_driver(.true.,.false.,n,n_want,n_eig,itmax,tol,0.0_dp,mmult, &
                       mprec,smult,eig,evec,ok)
!
!   write the converged results on file for comparison:
!
    open (unit = 10, file = 'lobcpg.txt', form = 'formatted', access = 'sequential')
    do i = 1, n_want
      write(10,1000) i, eig(i)
      if (evec(1,i) .lt. 0.0d0) evec(:,i) = - evec(:,i)
      write(10,'(10f12.6)') evec(:,i)
      write(10,*)
    end do
    close (10)
!
!   make a guess for the eigenvector
!
    call guess_evec(4,n,n_eig,diagonal,evec)
!
!   call the davidson driver:
!
    call davidson_driver(.true.,n,n_want,n_eig,itmax,tol,m_max,0.0_dp,mmult, &
                         mprec,eig,evec,ok)
!
!   write the converged results on file for comparison:
!
    open (unit = 10, file = 'davidson.txt', form = 'formatted', access = 'sequential')
    do i = 1, n_want
      write(10,1000) i, eig(i)
      if (evec(1,i) .lt. 0.0d0) evec(:,i) = - evec(:,i)
      write(10,'(10f12.6)') evec(:,i)
      write(10,*)
    end do
    close (10)
!
    return
  end subroutine test_symm
!
  subroutine test_geneig(check_lapack,n,n_want,tol,itmax,m_max)
    use real_precision
    use utils
    use diaglib, only : lobpcg_driver, gen_david_driver
    implicit none
    logical,  intent(in) :: check_lapack
    integer,  intent(in) :: n, n_want, itmax, m_max
    real(dp), intent(in) :: tol
!
    logical  :: ok
    integer  :: i, j, n_eig, info, lwork
    real(dp) :: lw(1)
    real(dp), allocatable :: diagonal(:), eig(:), evec(:,:)
    real(dp), allocatable :: work(:), w(:), a_copy(:,:), s_copy(:,:)
!
    external :: mmult, mprec, smult
!
    1000 format(t3,' eigenvalue # ',i6,': ',f12.6,/,t3,' eigenvector: ')
!
!
!   allocate memory for the a and s matrices
!
    allocate (a(n,n), s(n,n))
!
!   build a symmetric, positive definite metric matrix 
!
    call random_number(a)
    s = matmul(transpose(a),a)
!
!   build a symmetric matrix:
!
    a = 0.0_dp
    do i = 1, n
      a(i,i) = real(i,kind=dp) + 1.0_dp
      do j = 1, i - 1
        a(j,i) = 1.0_dp/real(i+j,kind=dp)
        a(i,j) = a(j,i)
      end do
    end do
!
!   if required, solve the problem with a dense lapack routine:
!
    if (check_lapack) then
      allocate (a_copy(n,n), s_copy(n,n), w(n))
      a_copy = a
      s_copy = s
      call dsygv(1,'v','l',n,a_copy,n,s_copy,n,w,lw,-1,info)
      lwork = int(lw(1))
      allocate (work(lwork))
      call dsygv(1,'v','l',n,a_copy,n,s_copy,n,w,work,lwork,info)
!
!     write the results on file for comparison:
!
      open (unit = 10, file = 'lapack.txt', form = 'formatted', access = 'sequential')
      do i = 1, n_want
        write(10,1000) i, w(i)
!
!       fix the phase so that the first element of the eigenvector is positive.
!
        if (a_copy(1,i) .lt. 0.0_dp) a_copy(:,i) = - a_copy(:,i)
        write(10,'(10f12.6)') a_copy(:,i)
        write(10,*)
      end do
      close (10)
    end if
!
!   allocate and gather the diagonal:
!
    allocate (diagonal(n))
    do i = 1, n
      diagonal(i) = a(i,i) - s(i,i)
    end do
!
!   for better convergence, we seek more eigenpairs and stop the iterations when the
!   required ones are converged.
!
    n_eig = min(2*n_want, n_want + 5)
!
!   allocate memory for the eigenvalues and eigenvectors:
!
    allocate (eig(n), evec(n,n_eig))
!
!   compute a guess for the eigenvector (see guess_evec for more information)
!
    call guess_evec(4,n,n_eig,diagonal,evec)
!
!   call lobpcg:
!
    call lobpcg_driver(.true.,.true.,n,n_want,n_eig,itmax,tol,0.0_dp,mmult, &
                       mprec,smult,eig,evec,ok)
!
!   write the converged results on file for comparison:
!
    open (unit = 10, file = 'lobcpg.txt', form = 'formatted', access = 'sequential')
    do i = 1, n_want
      write(10,1000) i, eig(i)
      if (evec(1,i) .lt. 0.0d0) evec(:,i) = - evec(:,i)
      write(10,'(10f12.6)') evec(:,i)
      write(10,*)
    end do
    close (10)
!
!   make a guess for the eigenvector
!
    call guess_evec(4,n,n_eig,diagonal,evec)
!
!   call the davidson driver:
!
    call gen_david_driver(.true.,n,n_want,n_eig,itmax,tol,m_max,0.0_dp,mmult, &
                         mprec,smult,eig,evec,ok)
!
!   write the converged results on file for comparison:
!
    open (unit = 10, file = 'davidson.txt', form = 'formatted', access = 'sequential')
    do i = 1, n_want
      write(10,1000) i, eig(i)
      if (evec(1,i) .lt. 0.0d0) evec(:,i) = - evec(:,i)
      write(10,'(10f12.6)') evec(:,i)
      write(10,*)
    end do
    close (10)
!
    return
  end subroutine test_geneig
!
  subroutine test_caslr(check_lapack,n,n_want,tol,itmax,m_max)
    use real_precision
    use utils
    use diaglib, only : caslr_driver, caslr_eff_driver, prt_generic
    implicit none
    logical,  intent(in) :: check_lapack
    integer,  intent(in) :: n, n_want, itmax, m_max
    real(dp), intent(in) :: tol
!
!   this subroutine builds a generalized eigenvalue problem such as the one encountered
!   in casscf linear response theory, and then uses lapack and iterative routines to
!   solve it.
!
    logical  :: ok
    integer  :: i, j, lwork, info
    integer  :: n2, n_eig
    real(dp) :: sqrttwo, lw(1)
    real(dp), allocatable :: evec(:,:), eig(:), diagonal(:)
    real(dp), allocatable :: work(:), w(:)
!
    external apbvec, ambvec, spdvec, smdvec, lrprec_1, lrprec_2
!
    integer :: n_seed, k
    integer, dimension(:), allocatable :: iseed
    integer, dimension(:), allocatable :: seed_array
    integer :: t
    1000 format(t3,' eigenvalue # ',i6,': ',f12.6,/,t3,' eigenvector: ')
!
!   the actual size of the problem is 2*n:
!
!   Initialize seeds with current time
    call random_seed(size = n_seed)
    allocate(iseed(n_seed))
!!    call system_clock(t)
!!    iseed = t + (/(k, k = 0, n_seed-1)/)
    iseed = 10   
    call random_seed(put = iseed)
    deallocate(iseed)
!
    n2 = 2 * n
!
!   allocate memory for the a, b, apb, amb, sigma, delta, spd, smd matrices:
!
    allocate (aa(n,n), bb(n,n), apb(n,n), amb(n,n), sigma(n,n), delta(n,n), spd(n,n), smd(n,n))
!
!   build a positive definite, symmetric apb, amb, sigma matrices:
!  
    do i = 1, n
      apb(i,i) = 5.0_dp + real(i,kind=dp)
      do j = 1, i - 1
        apb(j,i) = 1.0_dp/real(i+j,kind=dp)
        apb(i,j) = apb(j,i)
      end do 
    end do
    do i = 1, n
      amb(i,i) =2.0_dp + dble(i)
      do j = 1, i - 1
        apb(j,i) = 0.2_dp/dble(i+j)
        apb(i,j) = apb(j,i)
      end do 
    end do
!  
!    call random_number(sigma)
!    print *, "SIGMA"
!    call prt_generic(n,n,n,n,sigma,2)
!!    
    delta = matmul(transpose(sigma),sigma)
    sigma = delta
    do i = 1, n
      sigma(i,i) = sigma(i,i) + 1.0_dp
    end do
!  
!   build antisymmetric delta:
!  
    call random_number(delta)
    delta = delta - transpose(delta)
!
!    print *, "DELTA"
!    call prt_generic(n,n,n,n,delta,2)
!   build a and b:
!  
    aa = 0.5_dp * (apb + amb)
    bb = 0.5_dp * (apb - amb)
!  
!   build sigma + delta and sigma - delta
!  
    spd = sigma + delta
    smd = sigma - delta 
!
    if (check_lapack) then
!
!     build the complete matrices:
!    
      allocate (a(n2,n2), s(n2,n2), w(n2))
      a(1:n,   1:n)    = aa
      a(n+1:n2,n+1:n2) = aa
      a(1:n,   n+1:n2) = bb
      a(n+1:n2,1:n)    = bb
      s(1:n,   1:n)    = sigma
      s(n+1:n2,n+1:n2) = - sigma
      s(1:n,   n+1:n2) = delta
      s(n+1:n2,1:n)    = - delta
!
!   if required, solve the generalized eigenvalue problem with a dense
!   lapack routine:
!
      call dsygv(1,'V','L',n2,s,n2,a,n2,w,lw,-1,info)
      lwork = int(lw(1))
      allocate (work(lwork))
      call dsygv(1,'V','L',n2,s,n2,a,n2,w,work,lwork,info)
!
!     write the results on a file for comparison.
!
      open (unit = 10, file = 'lapack.txt', form = 'formatted', access = 'sequential')
      do i = 1, n_want
        write(10,1000) i, 1.0_dp/w(n2-i+1)
!
!       fix the phase so that the first element of the eigenvector is positive.
!
        if (s(1,n2-i+1) .lt. 0.0_dp) s(:,n2-i+1) = - s(:,n2-i+1)
        write(10,'(10f12.6)') s(:,n2-i+1)
        write(10,*)
      end do
      close (10)
!
!     free the memory:
!
      deallocate (a, s, w, work)
    end if
!
!   allocate space for the diagonal:
!
    allocate (diagonal(n))
!
!   gather diag(a) - diag(sigma):
!
    do i = 1, n
      diagonal(i) = aa(i,i) - sigma(i,i)
    enddo
!
!   for better convergence, we seek more eigenpairs and stop the iterations when the
!   required ones are converged.
!
    n_eig = min(2*n_want, n_want + 5)
!
!   allocate memory for the eigenvectors and eigenvalues:
!
    allocate (evec(n2,n_eig), eig(n_eig))
!
!   make a guess for the eigenvector (see guess_evec for more information)
!
    call guess_evec(4,n2,n_eig,diagonal,evec)
!
!   call the traditional solver:
!
    write(6,*) ' traditional implementation'
    write(6,*)
    call caslr_driver(.true.,n,n2,n_want,n_eig,itmax,tol,m_max,apbvec,ambvec, &
                      spdvec,smdvec,lrprec_1,eig,evec,ok)
!
!   write the converged results on file for comparison:
!
    sqrttwo = sqrt(2.0_dp)
    open (unit = 10, file = 'caslr.txt', form = 'formatted', access = 'sequential')
    do i = 1, n_want
      write(10,1000) i, eig(i)
      if (evec(1,i) .lt. 0.0d0) evec(:,i) = - evec(:,i)
      write(10,'(10f12.6)') evec(:,i)/sqrttwo
      write(10,*)
    end do
    close (10)
!
!   make a guess for the eigenvector (see guess_evec for more information)
!
    call guess_evec(4,n2,n_eig,diagonal,evec)
!
!   call the traditional solver using the helmich-paris algorithm:
!
    i_alg = 1
    write(6,*) ' helmich-paris implementation'
    write(6,*)
    call caslr_driver(.true.,n,n2,n_want,n_eig,itmax,tol,m_max,apbvec,ambvec, &
                      spdvec,smdvec,lrprec_1,eig,evec,ok)
    i_alg = 0
!
!   write the converged results on file for comparison:
!
    sqrttwo = sqrt(2.0_dp)
    open (unit = 10, file = 'cashp.txt', form = 'formatted', access = 'sequential')
    do i = 1, n_want
      write(10,1000) i, eig(i)
      if (evec(1,i) .lt. 0.0d0) evec(:,i) = - evec(:,i)
      write(10,'(10f12.6)') evec(:,i)/sqrttwo
      write(10,*)
    end do
    close (10)
!
!   make a guess for the eigenvector (see guess_evec for more information)
!
    call guess_evec(4,n2,n_eig,diagonal,evec)
!
!   call the modified solver:
!
    write(6,*) ' new implementation'
    write(6,*)
    call caslr_eff_driver(.true.,n,n2,n_want,n_eig,itmax,tol,m_max,apbvec,ambvec, &
                          spdvec,smdvec,lrprec_2,eig,evec,ok)
!
!   write the converged results on file for comparison:
!
    open (unit = 10, file = 'caslr_eff.txt', form = 'formatted', access = 'sequential')
    do i = 1, n_want
      write(10,1000) i, eig(i)
      if (evec(1,i) .lt. 0.0d0) evec(:,i) = - evec(:,i)
      write(10,'(10f12.6)') evec(:,i)/2.0_dp
      write(10,*)
    end do
    close (10)
!
    return
  end subroutine test_caslr
!
  subroutine test_caslr_std(check_lapack,n,n_want,tol,itmax,m_max)
    use real_precision
    use utils
    use diaglib, only : caslr_std, caslr_eff_std, prt_generic
    implicit none
    logical,  intent(in) :: check_lapack
    integer,  intent(in) :: n, n_want, itmax, m_max
    real(dp), intent(in) :: tol
!
!   this subroutine builds a set of linear equations for a real symmetric matrix,  such as the one encountered
!   in casscf linear response theory, and then uses lapack and iterative routines to
!   solve it.
!
    logical               :: ok, imag, verbose
    integer               :: i, j, info
    integer               :: n2, n_eig
    real(dp), allocatable :: ipiv(:) 
    real(dp)              :: omega
    real(dp)              :: sqrttwo, lw(1)
    real(dp)              :: zero = 0.0d0
    real(dp), allocatable :: vec(:,:), diagonal(:), eig(:), g(:,:), g_half(:,:)
!
    external apbvec, ambvec, spdvec, smdvec, lrprec_1, lrprec_2
!
    integer :: n_seed, k
    integer, dimension(:), allocatable :: iseed
    integer, dimension(:), allocatable :: seed_array
    integer :: t
!
!   set a couple of logicals
!
    verbose = .true.
    imag    = .false.
!
!   Initialize seeds with current time
    call random_seed(size = n_seed)
    allocate(iseed(n_seed))
    call system_clock(t)
    iseed = t + (/(k, k = 0, n_seed-1)/)
!    iseed = 10   
    call random_seed(put = iseed)
    deallocate(iseed)
!
!   the actual size of the problem is 2*n:
!
    n2 = 2 * n
!
!   enter the value of omega
!
    omega = 1.0_dp
!
!   allocate memory for the a, b, apb, amb, sigma, delta, spd, smd and  g matrices:
!
    allocate (aa(n,n), bb(n,n), apb(n,n), amb(n,n), sigma(n,n), delta(n,n), spd(n,n), smd(n,n), &
             g(n2,n_want), g_half(n,n_want))
!    
!   build a positive definite, symmetric apb, amb, sigma matrices:
!  
    do i = 1, n
      apb(i,i) = 5.0_dp + real(i,kind=dp)
      do j = 1, i - 1
        apb(j,i) = 1.0_dp/real(i+j,kind=dp)
        apb(i,j) = apb(j,i)
      end do 
    end do
    do i = 1, n
      amb(i,i) = 2.0_dp + dble(i)
      do j = 1, i - 1
        apb(j,i) = 0.2_dp/dble(i+j)
        apb(i,j) = apb(j,i)
      end do 
    end do
    delta = matmul(transpose(sigma),sigma)
    sigma = delta
    do i = 1, n
      sigma(i,i) = sigma(i,i) + 1.0_dp
    end do
!  
!   build antisymmetric delta:
!  
    call random_number(delta)
    delta = delta - transpose(delta)
!
    aa = 0.5_dp * (apb + amb)
    bb = 0.5_dp * (apb - amb)
!  
!   build sigma + delta and sigma - delta
!  
    spd = sigma + delta
    smd = sigma - delta
! 
!   make a distinction between a real and an imaginary perturbation for the costrunction of the gradient vector
!
    if (imag) then
      call random_number(g(1:n,:))
      g(n+1:n2,:) = -g(1:n,:) 
    else
      call random_number(g(1:n,:))
      g(n+1:n2,:) = g(1:n,:)
    end if
    g_half(1:n,:) = g(1:n,:)
!
!   allocate space for dgesv
    if (check_lapack) then
!
!     build the complete matrices:
!    
      allocate (a(n2,n2), s(n2,n2))
!
      a(1:n,   1:n)    =    aa
      a(n+1:n2,n+1:n2) =    aa
      a(1:n,   n+1:n2) =    bb
      a(n+1:n2,1:n)    =    bb
      s(1:n,   1:n)    =   sigma
      s(n+1:n2,n+1:n2) = - sigma
      s(1:n,   n+1:n2) =   delta
      s(n+1:n2,1:n)    = - delta
!
!   multiply directly omega to s
!
    s = omega * s
!
!   compute a-s
!
    a = a - s
!
    allocate(ipiv(n2))
!    
!   solve the linear response equation system using a dense linear algebra routine 
!
    call dgesv(n2,n_want,a,n2,ipiv,g,n2,info)
    deallocate(ipiv)
!
!   verify the result
! 
    open (unit = 10, file = 'test_caslr_std.txt', form = 'formatted', access = 'sequential')
    do i = 1, n2
      write(10,'(10f12.6)') g(i,:)
      write(10,*)
    end do
    close (10)
!
!     free the memory:
      deallocate(a,s)
    end if
      deallocate (g)
!
!   allocate space for the diagonal:
!
    allocate (diagonal(n))
!
!   gather diag(a) - diag(sigma):
!
    do i = 1, n
      diagonal(i) = aa(i,i) - sigma(i,i)
    end do
!
!   for better convergence, we seek more eigenpairs and stop the iterations when the
!   required ones are converged.
!
    n_eig = n_want
!
!   allocate memory for the solution vectors:
!
    allocate (vec(n2,n_eig))
!
!   make a guess for the solution vector (see guess_evec for more information)
!
    vec = 0.0d0
!    call guess_evec(3,n2,n_eig,diagonal,vec)
!
!   call the traditional solver:
!
!    write(6,*) ' traditional implementation'
!    write(6,*)
!    call caslr_std(verbose,n,n2,n_want,n_eig,itmax,tol,m_max,apbvec,ambvec, &
!                      spdvec,smdvec,lrprec_1,vec,ok,omega,g_half,imag)

!   write the converged results on file for comparison:
!
!    open (unit = 20, file = 'caslr_std.txt', status = 'replace', form = 'formatted', access = 'sequential')
!    do i = 1, n2
!      write(20,'(10f12.6)') vec(i,:)
!      write(20,*)
!    end do
!    close (20)
!
!   make a guess for the solution vector (see guess_evec for more information)
!
    call guess_evec(1,n2,n_eig,diagonal,vec)
!
!   call the modified solver:
!
    write(6,*) ' new implementation'
    write(6,*)
    call caslr_eff_std(verbose,n,n2,n_want,n_eig,itmax,tol,m_max,apbvec,ambvec, &
                          spdvec,smdvec,lrprec_1,vec,ok,omega,g_half,imag)
!
!   write the converged results on file for comparison:
!
    open (unit = 70, file = 'caslr_eff_std.txt', form = 'formatted', access = 'sequential')
    do i = 1, n2
      write(70,'(10f12.6)')  vec(i,:)
      write(70,*)
    end do
    close (70)
    return
  end subroutine test_caslr_std
!  
  subroutine test_scflr(check_lapack,n,n_want,tol,itmax,m_max)
    use real_precision
    use utils
    use diaglib, only : caslr_driver, caslr_eff_driver
    implicit none
    logical,  intent(in) :: check_lapack
    integer,  intent(in) :: n, n_want, itmax, m_max
    real(dp), intent(in) :: tol
!
!   this subroutine builds a generalized eigenvalue problem such as the one encountered
!   in hf/dft linear response theory, and then uses lapack and iterative routines to
!   solve it.
!
    logical  :: ok
    integer  :: i, j, lwork, info
    integer  :: n2, n_eig
    real(dp) :: sqrttwo, lw(1)
    real(dp), allocatable :: evec(:,:), eig(:), diagonal(:)
    real(dp), allocatable :: work(:), w(:)
!
    external apbvec, ambvec, spdvec, smdvec, lrprec_1, lrprec_2
!
    1000 format(t3,' eigenvalue # ',i6,': ',f12.6,/,t3,' eigenvector: ')
!
!   the actual size of the problem is 2*n:
!
    n2 = 2 * n
!
!   set the flag that avoids multiplications with the metric
!
    tdscf = .true.
!
!   allocate memory for the a, b, apb, amb, sigma, delta, spd, smd matrices:
!
    allocate (aa(n,n), bb(n,n), apb(n,n), amb(n,n), sigma(n,n), delta(n,n), spd(n,n), smd(n,n))
!
!   build a positive definite, symmetric apb, amb matrices:
!  
    do i = 1, n
      apb(i,i) = 5.0_dp + real(i,kind=dp)
      do j = 1, i - 1
        apb(j,i) = 1.0_dp/real(i+j,kind=dp)
        apb(i,j) = apb(j,i)
      end do 
    end do
    do i = 1, n
      amb(i,i) =2.0_dp + dble(i)
      do j = 1, i - 1
        apb(j,i) = 0.2_dp/dble(i+j)
        apb(i,j) = apb(j,i)
      end do 
    end do
!  
!   set sigma to the identity:
!
    sigma = 0.0_dp
    do i = 1, n
      sigma(i,i) = 1.0_dp
    end do
!  
!   set delta to zero:
!  
    delta = 0.0_dp
!  
!   build a and b:
!  
    aa = 0.5_dp * (apb + amb)
    bb = 0.5_dp * (apb - amb)
!  
!   build sigma + delta and sigma - delta
!  
    spd = sigma
    smd = sigma
!
    if (check_lapack) then
!
!     build the complete matrices:
!    
      allocate (a(n2,n2), s(n2,n2), w(n2))
      a(1:n,   1:n)    = aa
      a(n+1:n2,n+1:n2) = aa
      a(1:n,   n+1:n2) = bb
      a(n+1:n2,1:n)    = bb
      s(1:n,   1:n)    = sigma
      s(n+1:n2,n+1:n2) = - sigma
      s(1:n,   n+1:n2) = delta
      s(n+1:n2,1:n)    = - delta
!
!   if required, solve the generalized eigenvalue problem with a dense
!   lapack routine:
!
      call dsygv(1,'V','L',n2,s,n2,a,n2,w,lw,-1,info)
      lwork = int(lw(1))
      allocate (work(lwork))
      call dsygv(1,'V','L',n2,s,n2,a,n2,w,work,lwork,info)
!
!     write the results on a file for comparison.
!
      open (unit = 10, file = 'lapack.txt', form = 'formatted', access = 'sequential')
      do i = 1, n_want
        write(10,1000) i, 1.0_dp/w(n2-i+1)
!
!       fix the phase so that the first element of the eigenvector is positive.
!
        if (s(1,n2-i+1) .lt. 0.0_dp) s(:,n2-i+1) = - s(:,n2-i+1)
        write(10,'(10f12.6)') s(:,n2-i+1)
        write(10,*)
      end do
      close (10)
!
!     free the memory:
!
      deallocate (a, s, w, work)
    end if
!
!   allocate space for the diagonal:
!
    allocate (diagonal(n))
!
!   gather diag(a) - diag(sigma):
!
    do i = 1, n
      diagonal(i) = aa(i,i) - sigma(i,i)
    enddo
!
!   for better convergence, we seek more eigenpairs and stop the iterations when the
!   required ones are converged.
!
    n_eig = min(2*n_want, n_want + 5)
!
!   allocate memory for the eigenvectors and eigenvalues:
!
    allocate (evec(n2,n_eig), eig(n_eig))
!
!   make a guess for the eigenvector (see guess_evec for more information)
!
    call guess_evec(4,n2,n_eig,diagonal,evec)
!
!   call the traditional solver:
!
    call caslr_driver(.true.,n,n2,n_want,n_eig,itmax,tol,m_max,apbvec,ambvec, &
                      spdvec,smdvec,lrprec_1,eig,evec,ok)
!
!   write the converged results on file for comparison:
!
    sqrttwo = sqrt(2.0_dp)
    open (unit = 10, file = 'caslr.txt', form = 'formatted', access = 'sequential')
    do i = 1, n_want
      write(10,1000) i, eig(i)
      if (evec(1,i) .lt. 0.0d0) evec(:,i) = - evec(:,i)
      write(10,'(10f12.6)') evec(:,i)/sqrttwo
      write(10,*)
    end do
    close (10)
!
!   make a guess for the eigenvector (see guess_evec for more information)
!
    call guess_evec(4,n2,n_eig,diagonal,evec)
!
!   call the modified solver:
!
    call caslr_eff_driver(.true.,n,n2,n_want,n_eig,itmax,tol,m_max,apbvec,ambvec, &
                          spdvec,smdvec,lrprec_2,eig,evec,ok)
!
!   write the converged results on file for comparison:
!
    open (unit = 10, file = 'caslr_eff.txt', form = 'formatted', access = 'sequential')
    do i = 1, n_want
      write(10,1000) i, eig(i)
      if (evec(1,i) .lt. 0.0d0) evec(:,i) = - evec(:,i)
      write(10,'(10f12.6)') evec(:,i)/2.0_dp
      write(10,*)
    end do
    close (10)
!
    return
  end subroutine test_scflr
!
  subroutine guess_evec(iwhat,n,m,diagonal,evec)
    use real_precision
    implicit none
    integer,                  intent(in)    :: iwhat, n, m
    real(dp), dimension(n),   intent(in)    :: diagonal
    real(dp), dimension(n,m), intent(inout) :: evec
!
!   guess the eigenvector.
!
    integer              :: i, ipos, n_seed, k, t
    integer, allocatable :: iseed(:)
    logical, allocatable :: mask(:)
!
!   initialize a random number generator in a predictible way.
!
    call random_seed(size=n_seed)
    allocate (iseed(n_seed))
    iseed = 1
!    call system_clock(t)
!    iseed = t + (/(k, k = 0, n_seed-1)/)
    call random_seed(put=iseed)
    deallocate (iseed)
!
    evec = 0.0_dp
!
    allocate (mask(n))
!
    if (iwhat.eq.1) then
!
!     get the minimum element of the diagonal
!
      mask = .true.
      do i = 1, m
        ipos = minloc(diagonal,dim=1,mask=mask)
        mask(ipos) = .false.   
        evec(ipos,i) = 1.0d0
      enddo
    else if (iwhat.eq.2) then
!
!     get the maximum element of the diagonal
!
      mask = .true.
      do i = 1, m
        ipos = maxloc(diagonal,dim=1,mask=mask)
        mask(ipos) = .false.   
        evec(ipos,i) = 1.0d0
      enddo
    else if (iwhat.eq.3) then
!
!     random vector between 0 and 1
!
      call random_number(evec)
    else if (iwhat.eq.4) then
!
!     random vector between -0.5 and 0.5
!
      call random_number(evec)
      evec = evec - 0.50_dp
    else if (iwhat.eq.5) then
!
      call random_number(evec)
      evec = evec * 0.01_dp
!
!     get the maximum element of the diagonal
!
      mask = .true.
      do i = 1, m
        ipos = maxloc(diagonal,dim=1,mask=mask)
        mask(ipos) = .false.   
        evec(ipos,i) = evec(ipos,i) + 1.0d0
      enddo

    else if (iwhat.eq.6) then
!
      call random_number(evec)
      evec = evec * 0.01_dp
!
!     get the minimum element of the diagonal
!
      mask = .true.
      do i = 1, m
        ipos = minloc(diagonal,dim=1,mask=mask)
        mask(ipos) = .false.   
        evec(ipos,i) = evec(ipos,i) + 1.0d0
      enddo
    end if
    return
  end subroutine guess_evec
