module diaglib
  use real_precision
  implicit none
!
! diaglib - a fortran library of matrix-free iterative algorithms to
! compute a few eigenvalues and eigenvectors of large matrices.
! ==================================================================
!
! Implementation by
!
!   Ivan Gianni', Tommaso Nottoli, Federica Pes, Antoine Levitt, and Filippo Lipparini
!   MoLECoLab Pisa
!   Department of Chemistry and Industrial Chemistry
!   University of Pisa
!   Via G. Moruzzi 13, I-56124, Pisa, Italy
!
! Pisa, november 2022
!
!                                                                     
!                                       mm                            
!                                    mMMm                             
!                                  mMMMMm         m                   
!                                 mMMMMm          mMm                 
!                                 mMMMMm          mMm                 
!                                 mMMMMMm        mMMm                 
!                                 MMMMMMMMMMMMMMMMMMm                 
!                                mMMMMMMMMMMMMMMMMMm                  
!       __  ___      __    ____________      __MMMm     __            
!      /  |/  /___  / /   / ____/ ____/___  / /  ____ _/ /_           
!     / /|_/ / __ \/ /   / __/ / /   / __ \/ /  / __ `/ __ \          
!    / /  / / /_/ / /___/ /___/ /___/ /_/ / /__/ /_/ / /_/ /          
!   /_/  /_/\__________/_____/\____/_____/_____|__,_/_.___/           
!           /_  __/ __ \/ __ \/ /  / ___/                             
!            / / / / / / / / / /   \__ \                              
!           / / / /_/ / /_/ / /___ __/ /                              
!          /_/  \____/\____/_____/____/                               
!            mMMMMMMMMMMMMMMMm                                        
!          mMMMMMMMMMMMMMMMm                                          
!        mMMMMMMMMMMMMMMMMM   + ------------------------------------ +
!       mMMMMMMMMMMMMMMMMm    |            D I A G L I B             |
!      mMMMMMMMMMMMMMMMMMm    + ------------------------------------ +
!      mMMMMm       mMMMMMm   | I. Gianni', T. Nottoli, F. Lipparini |
!      mMMMm       mMMMMMMm   |                                      |
!       mMm       mMMMMMMm    |                            ver 1.0   |
!        m       mMMMMMMm     |              molecolab.dcci.unipi.it |
!               mMMMMMm       + ------------------------------------ +
!                                                                     
! description of the library:
! ===========================
!
! diaglib provides an implementation of two matrix-free algorithms to
! compute a few eigenvalues and eigenvectors of a large, possibly sparse
! matrix.
!
! the available algorithms are
!
! 1) locally optimal block preconditioned conjugate gradient 
!
! 2) davidson-liu
!
! both algorithms require two user-provided routines to apply the matrx
! and a suitable preconditioner to a set of vectors.
! such routines have the following interface:
!
!   subroutine matvec(n,m,x,ax)
!   subroutine precnd(n,m,shift,x,ax)
!
! where n,m are integers and x(n,m) and ax(n,m) are double precision
! arrays.
! as using the first eigenvalue in a shift-and-invert spirit is very 
! common, a double precision scalar shift is also passed to precnd.
!
! both implementations favor numerical stability over efficiency and are
! targeted at applications in molecular quantum chemistry, such as in
! (full) ci or augmented hessian calculations, where typically m << n.
!
! list of provided routines:
! ==========================
!
! lobpcg_driver:   main lobpcg driver.
!
! davidson_driver: main davidson-liu driver
!
! ortho_vs_x:      subroutine to orthogonalize a set of vectors w against
!                  another set of orthonormal vectors v.
!                  the resulting w vectors are then also orthonormalized.
!                  the linear transformation applied to w can also be 
!                  applied to a second set of vectors, tipically, aw.
!
! ortho:           subroutine to perferm the orthonormalization of a set
!                  of vectors v using QR decomposition.
!                  the linear transformation applied to v can also be 
!                  applied to a second set of vectors, tipically, av.
!
! ortho_cd:        subroutine to perferm the orthonormalization of a set
!                  of vectors v using cholesky decomposition with level
!                  shifting and iterative refinement.
!                  the linear transformation applied to v can also be 
!                  applied to a second set of vectors, tipically, av.
!
! b_ortho_vs_x:    same as ortho_vs_x, but vectors are b-orthogonalized
!                  against the existing ones, where b is the metric in the
!                  generalized eigenvalues problem.
!
! b_ortho:         same as ortho_cd, but the orthogonalization is performed
!                  with respect to the metric b.
!
! check_guess:     routine to check whether a guess is provided by the user
!                  and whether the provided vectors are already orthonormal.
!                  if no guess is provided, a random guess is generated.
!
! get_coeffs:      utility to compute the expasion coefficient for the p
!                  vectors. used in lobpcg.
!
! diag_shift:      utility to apply a diagonal shift to the diagonal of a
!                  given matrix. 
!                  used in ortho_cd.
!
! get_mem_lapack:  utility to compute how much scratch memory is required by
!                  the various lapack routines called in the code.
!
! check_mem:       utility to check whether memory was allocated or deallocated
!                  successfully.
!
! prtmat:          utility for debug printing.
!
! get_time:        subroutine to get cpu and wall time.
!
! please, see the specific routines for a more detailed description.
!
! dependencies:
! =============
!
! the following blas/lapack routines are used:
!
!   blas:          daxpy, dcopy, dgemm, dnrm2
!   lapack:        dgeqrf, dpotrf, dsyev, dtrsm 
!
! shared variables:
! =================
!
  private
!
! useful constants
!
  real(dp), parameter    :: zero = 0.0_dp, one = 1.0_dp, two = 2.0_dp, ten = 10.0_dp
!
! convergence thresholds for orthogonalization
!
!  real(dp), parameter    :: tol_ortho = 1.0e-13_dp
  real(dp), parameter    :: tol_ortho = two * epsilon(one)
! 
! memory and info for lapack routines
!
  integer                :: lwork, info
  real(dp), allocatable  :: work(:), tau(:)
!
! timings:
!
  real(dp)               :: t1(2), t2(2), t_diag(2), t_ortho(2), &
                            t_mv(2), t_tot(2), t_ls(2)
!
! subroutines:
! ============
!
  public :: lobpcg_driver, davidson_driver, gen_david_driver, caslr_driver, caslr_eff_driver, &
            ortho, b_ortho, ortho_cd, ortho_vs_x, b_ortho_vs_x, caslr_std, prt_generic, &
            caslr_eff_std   
!
  contains
!
  subroutine lobpcg_driver(verbose,gen_eig,n,n_targ,n_max,max_iter,tol, &
                           shift,matvec,precnd,bvec,eig,evec,ok)
!
!   main driver for lobpcg.
!
!   input variables:
!   ================
!
!   verbose:  logical, whether to print various information at each 
!             iteration (eigenvalues, residuals...).
!
!   gen_eig:  logical, whether a generalized eigenvalue problem has
!             to be solved. 
!
!   n:        integer, size of the matrix to be diagonalized.
!
!   n_targ:   integer, number of required eigenpairs.
!
!   n_max:    integer, maximum size of the search space. should be 
!             >= n_targ. note that eig and evec should be allocated 
!             n_max and (n,n_max) rather than n_targ and (n,n_targ). 
!             for better convergence, a value larger than n_targ (eg.,
!             n_targ + 10) is recommended.
!   
!   max_iter: integer, maximum allowed number of iterations.
!
!   tol:      double precision real, the convergence threshold.
!
!   shift:    double precision real, a diagonal level shifting parameter
!
!   matvec:   external subroutine that performs the matrix-vector
!             multiplication
!
!   precnd:   external subroutine that applies a preconditioner.
!
!   bvec:     external subroutine that applies the metric to a vector.
!             only referenced is gen is true.
!
!   output variables:
!   =================
!
!   eig:      double precision array of size n_max. if ok is true, 
!             the computed eigenvalues in asceding order.
!
!   evec:     double precision array of size (n,n_max). 
!             in input, a guess for the eigenvectors.
!             if ok is true, in output the computed eigenvectors.
!
!   ok:       logical, true if lobpcg converged.
!
    logical,                      intent(in)    :: verbose, gen_eig
    integer,                      intent(in)    :: n, n_targ, n_max
    integer,                      intent(in)    :: max_iter
    real(dp),                     intent(in)    :: tol, shift
    real(dp), dimension(n_max),   intent(inout) :: eig
    real(dp), dimension(n,n_max), intent(inout) :: evec
    logical,                      intent(inout) :: ok
    external                                    :: matvec, precnd, bvec
!
!   local variables:
!   ================
!
    integer               :: it, i_eig, n_act, ind_x, ind_w, ind_p, &
                             len_a, len_u, n_max_l
    integer               :: istat
    real(dp)              :: sqrtn, tol_rms, tol_max, xx(1)
    real(dp), allocatable :: space(:,:), aspace(:,:), bspace(:,:), &
                             a_red(:,:), e_red(:), r(:,:), r_norm(:,:)
    real(dp), allocatable :: u_x(:,:), u_p(:,:), x_new(:,:), ax_new(:,:), &
                             bx_new(:,:)
    logical,  allocatable :: done(:)
!
!   external functions:
!   ===================
!
    real(dp)              :: dnrm2
    external              :: dnrm2, daxpy, dsyev, dcopy
!
!   start by allocating memory for the various lapack routines
!
    lwork = get_mem_lapack(n,n_max)
    allocate (work(lwork), tau(2*n_max), stat=istat)
    call check_mem(istat)
!
!   allocate memory for the expansion space, the corresponding 
!   matrix-multiplied vectors and the residual:
!
    len_a = 3*n_max
    allocate (space(n,len_a), aspace(n,len_a), bspace(n,len_a), &
              r(n,n_max), stat=istat)
    call check_mem(istat)
!
!   allocate memory for the reduced matrix and its eigenvalues:
!
    allocate (a_red(len_a,len_a), e_red(len_a), stat=istat)
    call check_mem(istat)
!
!   allocate memory for temporary copies of x, ax, and bx:
!
    allocate (x_new(n,n_max), ax_new(n,n_max), bx_new(n,n_max), stat=istat)
    call check_mem(istat)
!
!   allocate memory for convergence check
!
    allocate (done(n_max), r_norm(2,n_max), stat=istat)
    call check_mem(istat)
!
!   clean out:
!
    t_diag   = zero
    t_ortho  = zero
    t_mv     = zero
    t_tot    = zero
    space    = zero
    bspace   = zero
    aspace   = zero
    a_red    = zero
    n_max_l  = n_max
!
    call get_time(t_tot)
!
!   check whether we have a guess for the eigenvectors in evec, and
!   whether it is orthonormal.
!   if evec is zero, create a random guess.
!
    call check_guess(n,n_max,evec)
!
!   if required, compute b*evec and b-orthogonalize the guess
!
    if (gen_eig) then 
      call bvec(n,n_max,evec,bx_new)
      call b_ortho(n,n_max,evec,bx_new)
    end if
!
!   compute the first eigenpairs by diagonalizing the reduced matrix:
!
    call dcopy(n*n_max,evec,1,space,1)
    if (gen_eig) call dcopy(n*n_max,bx_new,1,bspace,1)
    call get_time(t1)
    call matvec(n,n_max,space,aspace)
    call get_time(t2)
    t_mv = t_mv + t2 - t1
    if (shift.ne.zero) call daxpy(n*n_max,shift,space,1,aspace,1)
    call dgemm('t','n',n_max,n_max,n,one,space,n,aspace,n,zero,a_red,len_a)
    call get_time(t1)
    call dsyev('v','l',n_max,a_red,len_a,e_red,work,lwork,info)
    call get_time(t2)
    t_diag = t_diag + t2 - t1
    eig = e_red(1:n_max)
!
!   get the ritz vectors:
!
    call dgemm('n','n',n,n_max,n_max,one,space,n,a_red,len_a,zero,evec,n)
    call dcopy(n*n_max,evec,1,space,1)
    call dgemm('n','n',n,n_max,n_max,one,aspace,n,a_red,len_a,zero,evec,n)
    call dcopy(n*n_max,evec,1,aspace,1)
!
!   if required, also get b times the ritz vector:
!
    if (gen_eig) then 
      call dgemm('n','n',n,n_max,n_max,one,bspace,n,a_red,len_a,zero,evec,n)
      call dcopy(n*n_max,evec,1,bspace,1)
    end if
!
!   do the first iteration explicitly. 
!   build the residuals:
!
    call dcopy(n*n_max,aspace,1,r,1)
    if (gen_eig) then 
      do i_eig = 1, n_max
        call daxpy(n,-eig(i_eig),bspace(:,i_eig),1,r(:,i_eig),1)
      end do
    else
      do i_eig = 1, n_max
        call daxpy(n,-eig(i_eig),space(:,i_eig),1,r(:,i_eig),1)
      end do
    end if
!
!   compute the preconditioned residuals:
!
    ind_x = 1
    ind_w = ind_x + n_max 
    call precnd(n,n_max,shift-eig(ind_x),r(1,ind_x),space(1,ind_w))
!
!   orthogonalize:
!
    call get_time(t1)
    if (gen_eig) then
      call b_ortho_vs_x(n,n_max_l,n_max_l,space,bspace,space(1,ind_w))
!
!     after b_ortho, w is b-orthogonal to x, and orthonormal. 
!     compute the application of b to w, and b-orthonormalize it.
!
      call bvec(n,n_max,space(1,ind_w),bspace(1,ind_w))
      call b_ortho(n,n_max,space(1,ind_w),bspace(1,ind_w))
    else
      call ortho_vs_x(n,n_max_l,n_max_l,space,space(1,ind_w),xx,xx)
    end if
    call get_time(t2)
    t_ortho = t_ortho + t2 - t1
!
!   we are now ready to start the main loop.
!   initialize a few parameters
!
    tol_rms = tol
    tol_max = ten*tol
    sqrtn   = sqrt(real(n,dp))
    ok      = .false.
    done    = .false.
    n_act   = n_max
!
    1030 format(t5,'LOBPCG iterations (tol=',d10.2,'):',/, &
                t5,'------------------------------------------------------------------',/, &
                t7,'  iter  root              eigenvalue','         rms         max ok',/, &
                t5,'------------------------------------------------------------------')
    1040 format(t9,i4,2x,i4,f24.12,2d12.4,l3)
!
    if (verbose) write(6,1030) tol
!
    do it = 1, max_iter
!
!     perform the matrix-vector multiplication for this iteration:
!
      call get_time(t1)
      call matvec(n,n_act,space(1,ind_w),aspace(1,ind_w))
      call get_time(t2)
      t_mv = t_mv + t2 - t1
      if (shift.ne.zero) call daxpy(n*n_act,shift,space(1,ind_w),1,aspace(1,ind_w),1)
!
!     build the reduced matrix and diagonalize it:
!
      len_u = n_max + 2*n_act
      if (it.eq.1) len_u = 2*n_max
      call dgemm('t','n',len_u,len_u,n,one,space,n,aspace,n,zero,a_red,len_a)
!
      call get_time(t1)
      call dsyev('v','l',len_u,a_red,len_a,e_red,work,lwork,info)
      call get_time(t2)
      t_diag = t_diag + t2 - t1
!
!     if dsyev failed, print an error message and abort (this should not happen)
!
      if (info.ne.0) then
        write(6,'(t3,a,i6)') 'dsyev failed. info = ',info
      end if
      eig = e_red(1:n_max)
!
!     update x and ax, and, if required, bx:
!
      call dgemm('n','n',n,n_max,len_u,one,space,n,a_red,len_a,zero,x_new,n)
      call dgemm('n','n',n,n_max,len_u,one,aspace,n,a_red,len_a,zero,ax_new,n)
      if (gen_eig) then
        call dgemm('n','n',n,n_max,len_u,one,bspace,n,a_red,len_a,zero,bx_new,n)
      end if
!
!     compute the residuals and their rms and sup norms:
!
      call dcopy(n*n_max,ax_new,1,r,1)
      do i_eig = 1, n_max
!
!       if the eigenvalue is already converged, skip it.
!
        if (done(i_eig)) cycle
!
        if (gen_eig) then
          call daxpy(n,-eig(i_eig),bx_new(:,i_eig),1,r(:,i_eig),1)
        else
          call daxpy(n,-eig(i_eig),x_new(:,i_eig),1,r(:,i_eig),1)
        end if
        r_norm(1,i_eig) = dnrm2(n,r(:,i_eig),1)/sqrtn
        r_norm(2,i_eig) = maxval(abs(r(:,i_eig)))
      end do
!
!     only lock the first converged eigenvalues/vectors.
!
      do i_eig = 1, n_max
        if (done(i_eig)) cycle
        done(i_eig)     = r_norm(1,i_eig).lt.tol_rms .and. &
                          r_norm(2,i_eig).lt.tol_max .and. &
                          it.gt.1
        if (.not. done(i_eig)) then
          done(i_eig+1:n_max) = .false.
          exit
        end if
      end do
!
!     print some information and check for convergence:
!
      if (verbose) then
        do i_eig = 1, n_targ
          write(6,1040) it, i_eig, eig(i_eig) - shift, r_norm(:,i_eig), done(i_eig)
        end do
        write(6,*)
      end if
      if (all(done(1:n_targ))) then
        call dcopy(n*n_max,x_new,1,evec,1)
        ok = .true.
        exit
      end if
!
!     compute the number of active eigenvalues. 
!     converged eigenvalues and eigenvectors will be locked and kept 
!     for orthogonalization purposes.
!
      n_act = n_max - count(done)
      ind_x = n_max - n_act + 1
      ind_p = ind_x + n_act
      ind_w = ind_p + n_act
!
!     compute the new p and ap vectors only for the active eigenvectors.
!     this is done by computing the expansion coefficients u_p of x_new
!     -x in the basis of (x,p,w), and then by orthogonalizing then to
!     the coefficients u_x of x_new. 
!
      allocate (u_x(len_u,n_max), u_p(len_u,n_act), stat = istat)
      call check_mem(istat)
!
      call get_coeffs(len_a,len_u,n_max,n_act,a_red,u_x,u_p)
!
!     p  = space  * u_p
!     ap = aspace * u_p
!     bp = bspace * u_p
!     note that this is numerically safe, as u_p is orthogonal.
!
      call dgemm('n','n',n,n_act,len_u,one,space,n,u_p,len_u,zero,evec,n)
      call dcopy(n_act*n,evec,1,space(1,ind_p),1)
      call dgemm('n','n',n,n_act,len_u,one,aspace,n,u_p,len_u,zero,evec,n)
      call dcopy(n_act*n,evec,1,aspace(1,ind_p),1)
!
      if (gen_eig) then
        call dgemm('n','n',n,n_act,len_u,one,bspace,n,u_p,len_u,zero,evec,n)
        call dcopy(n_act*n,evec,1,bspace(1,ind_p),1)
      end if
!
      deallocate(u_x, u_p, stat = istat)
      call check_mem(istat)
!
!     now, move x_new and ax_new into space and aspace.
!
      call dcopy(n*n_max,x_new,1,space,1)
      call dcopy(n*n_max,ax_new,1,aspace,1)
      if (gen_eig) then 
        call dcopy(n*n_max,bx_new,1,bspace,1)
      end if
!
!     compute the preconditioned residuals w:
!
      call precnd(n,n_act,shift-eig(1),r(1,ind_x),space(1,ind_w))
!
!     orthogonalize w against x and p, and then orthonormalize it:
!
      call get_time(t1)
      if (gen_eig) then 
        call b_ortho_vs_x(n,n_max_l+n_act,n_act,space,bspace,space(1,ind_w))
        call bvec(n,n_act,space(1,ind_w),bspace(1,ind_w))
        call b_ortho(n,n_act,space(1,ind_w),bspace(1,ind_w))
      else
        call ortho_vs_x(n,n_max_l+n_act,n_act,space,space(1,ind_w),xx,xx)
      end if
      call get_time(t2)
      t_ortho = t_ortho + t2 - t1
!
    end do
!
    call get_time(t2)
    t_tot = t2 - t_tot
!
!   if required, print timings
!
    1000 format(t3,'timings for lobpcg (cpu/wall):   ',/, &
                t3,'  matrix-vector multiplications: ',2f12.4,/, &
                t3,'  diagonalization:               ',2f12.4,/, &
                t3,'  orthogonalization:             ',2f12.4,/, &
                t3,'                                 ',24('='),/,  &
                t3,'  total:                         ',2f12.4)
    if (verbose) write(6,1000) t_mv, t_diag, t_ortho, t_tot
!
!   deallocate memory and return.
!
    deallocate (work, tau, space, aspace, bspace, r, a_red, e_red, & 
                x_new, ax_new, bx_new, done, r_norm, stat = istat)
    call check_mem(istat)
!
    return

  end subroutine lobpcg_driver
!
  subroutine caslr_driver(verbose,n,n2,n_targ,n_max,max_iter,tol,max_dav, &
                          apbmul,ambmul,spdmul,smdmul,lrprec,eig,evec,ok)
    use utils
    implicit none
    logical, intent(in)                          :: verbose
    integer,                       intent(in)    :: n, n2, n_targ, n_max
    integer,                       intent(in)    :: max_iter, max_dav
    real(dp),                      intent(in)    :: tol
    real(dp), dimension(n_max),    intent(inout) :: eig
    real(dp), dimension(n2,n_max), intent(inout) :: evec
    logical,                       intent(inout) :: ok
    external                                     :: apbmul, ambmul, spdmul, smdmul, &
                                                    lrprec
!
!   local variables:
!   ================
!
    integer, parameter    :: min_dav = 10
    integer               :: istat
!
!   actual expansion space size and total dimension
!
    integer               :: dim_dav, lda, lda2
!
!   number of active vectors at a given iteration, and indices to access them
!
    integer               :: n_act, ind, i_beg
!
!   current size and total dimension of the expansion space
!
    integer               :: m_dim, ldu
!
!   number of frozen (i.e. converged) vectors
!
    integer               :: n_frozen
!
    integer               :: i, j, k, it, i_eig
!
    integer               :: lwork_svd
!
    real(dp)              :: sqrtn, tol_rms, tol_max, growth
    real(dp)              :: xx(1), lw_svd(1)
!
!   arrays to control convergence and orthogonalization
!
    logical,  allocatable :: done(:)
!
!   expansion spaces, residuals and their norms.
!
    real(dp), allocatable :: vp(:,:), vm(:,:), lvp(:,:), lvm(:,:), bvp(:,:), bvm(:,:)
    real(dp), allocatable :: rp(:,:), rm(:,:), rr(:,:), r_norm(:,:)
!
!   eigenvectors of the reduced problem and components of the ritz vectors:
!
    real(dp), allocatable :: up(:,:), um(:,:), eigp(:,:), eigm(:,:), bp(:,:), bm(:,:)
!
!   subspace matrix and eigenvalues.
!
    real(dp), allocatable :: a_red(:,:), a_copy(:,:), s_red(:,:), s_copy(:,:), e_red(:)
    real(dp), allocatable :: epmat(:,:), emmat(:,:), smat(:,:)
!
!   auxiliary quantities for the helmich-paris algorithm
!
    real(dp), allocatable :: u_1(:,:), sv_1(:), vt_1(:,:), u_2(:,:), sv_2(:), vt_2(:,:), &
                             ept(:,:), emt(:,:), cmat(:,:), xpt(:,:), xmt(:,:), scr(:,:)
    real(dp), allocatable :: work_svd(:)
!
!   restarting variables
!
    integer               :: n_rst
    logical               :: restart
!
!   external functions:
!   ===================
!
    real(dp)              :: dnrm2
    external              :: dcopy, dnrm2, dgemm, dsyev
!
!   compute the actual size of the expansion space, checking that
!   the input makes sense.
!   no expansion space smaller than max_dav = 10 is deemed acceptable.
!
    dim_dav = max(min_dav,max_dav)
    lda     = dim_dav*n_max
    lda2    = 2 * lda
!
!   start by allocating memory for the various lapack routines
!
    lwork = get_mem_lapack(n,n_max)
    allocate (work(lwork), tau(n_max), stat=istat)
    call check_mem(istat)
!
!   allocate memory for the expansion space, the corresponding
!   matrix-multiplied vectors and the residual:
!
    allocate (vp(n,lda), vm(n,lda), lvp(n,lda), lvm(n,lda), bvp(n,lda), bvm(n,lda), &
              rp(n,n_max), rm(n,n_max), rr(n,n_max), stat = istat)
    call check_mem(istat)
!
!   allocate memory for convergence check
!
    allocate (done(n_max), r_norm(2,n_max), stat=istat)
    call check_mem(istat)
!
!   allocate memory for the reduced matrix and its eigenvalues:
!
    allocate (a_red(lda2,lda2), a_copy(lda2,lda2), s_red(lda2,lda2), s_copy(lda2,lda2), &
      e_red(lda2), epmat(lda,lda), emmat(lda,lda), smat(lda,lda), stat=istat)
    call check_mem(istat)
!
!   allocate memory for the plus and minus eigenvector components:
!
    allocate (up(lda,n_max), um(lda,n_max), eigp(n,n_max), eigm(n,n_max), bp(n,n_max), bm(n,n_max), stat = istat)
    call check_mem(istat)
!
!   allocate additional memory for the helmich-paris algorithm:
!
    if (i_alg.eq.1) then
      allocate (u_1(lda,lda), sv_1(lda), vt_1(lda,lda), u_2(lda,lda), sv_2(lda), vt_2(lda,lda), &
                ept(lda,lda), emt(lda,lda), cmat(lda,lda), xpt(lda,lda), xmt(lda,lda), scr(lda,lda), stat = istat)
      call check_mem(istat)
!
      call dgesvd('a','a',lda,lda,u_1,lda,sv_1,u_1,lda,vt_1,lda,lw_svd,-1,info)
      lwork_svd = int(lw_svd(1))
      allocate (work_svd(lwork_svd))
    end if
!
!   set the tolerances and compute a useful constant to compute rms norms:
!
    sqrtn   = sqrt(real(n,dp))
    tol_rms = tol
    tol_max = 10.0_dp * tol
!
!   clean out various quantities
!
    t_diag  = zero
    t_ortho = zero
    t_mv    = zero
    t_tot   = zero
    vp      = zero
    vm      = zero
    bvp     = zero
    bvm     = zero
    lvp     = zero
    lvm     = zero
    a_red   = zero
    ok      = .false.
    done    = .false.
!
    call get_time(t_tot)
!
!   move the guess into the expansion space.
!
    do i_eig = 1, n_max
      vp(:,i_eig) = evec(1:n,i_eig) + evec(n+1:n2,i_eig)
      vm(:,i_eig) = evec(1:n,i_eig) - evec(n+1:n2,i_eig)
    end do
!
    call ortho_cd(n,n_max,vp,growth,ok)
    call ortho_cd(n,n_max,vm,growth,ok)
!
    n_act = n_max
    ind   = 1
    i_beg = 1
!
!   initialize the counter for the expansion of the subspace
!
    m_dim = 1
    ldu   = 0
!
!   initialize to false the restart
!
    restart = .false.
!
!   main loop:
!
    1030 format(t5,'Davidson-Liu iterations (tol=',d10.2,'):',/, &
                t5,'------------------------------------------------------------------',/, &
                t7,'  iter  root              eigenvalue','         rms         max ok',/, &
                t5,'------------------------------------------------------------------')
    1040 format(t9,i4,2x,i4,f24.12,2d12.4,l3)
!
    if (verbose) write(6,1030) tol
!
    n_rst   = 0
    do it = 1, max_iter
!
!     update the size of the expansion space.
!
      ldu = ldu + n_act
!
!     perform this iteration's matrix-vector multiplication:
!
      call get_time(t1)
      call apbmul(n,n_act,vp(1,i_beg),lvp(1,i_beg))
      call ambmul(n,n_act,vm(1,i_beg),lvm(1,i_beg))
      call spdmul(n,n_act,vp(1,i_beg),bvm(1,i_beg))
      call smdmul(n,n_act,vm(1,i_beg),bvp(1,i_beg))
      call get_time(t2)
      t_mv = t_mv + t2 - t1
!
!     update the reduced matrix
!
      call dgemm('t','n',ldu,ldu,n,one,vp,n,lvp,n,zero,epmat,lda)
      call dgemm('t','n',ldu,ldu,n,one,vm,n,lvm,n,zero,emmat,lda)
      call dgemm('t','n',ldu,ldu,n,one,vm,n,bvm,n,zero,smat,lda)
!
      a_red = zero
      s_red = zero
      a_red(1:ldu,1:ldu)             = epmat(1:ldu,1:ldu)
      a_red(ldu+1:2*ldu,ldu+1:2*ldu) = emmat(1:ldu,1:ldu)
      s_red(1:ldu,ldu+1:2*ldu)       = transpose(smat(1:ldu,1:ldu))
      s_red(ldu+1:2*ldu,1:ldu)       = smat(1:ldu,1:ldu)
!
      call get_time(t1)
      if (i_alg.eq.0) then
!
!       default algorithm: solve the 2n-dimensional inverse problem.
!
        call dsygv(1,'v','l',2*ldu,s_red,2*lda,a_red,2*lda,e_red,work,lwork,info)
!
!       extract the eigenvalues and compute the ritz approximation to the
!       eigenvectors
!
        do i_eig = 1, n_max
          eig(i_eig)      = one/e_red(2*ldu - i_eig + 1)
          up(1:ldu,i_eig) = s_red(1:ldu,2*ldu - i_eig + 1)
          um(1:ldu,i_eig) = s_red(ldu+1:2*ldu,2*ldu - i_eig + 1)
        end do
!
      else if (i_alg.eq.1) then
!
!       helmich-paris algorithm
!
        scr = smat
!
!       compute the svd of smat
!
        call dgesvd('a','a',ldu,ldu,scr,lda,sv_1,u_1,lda,vt_1,lda,work_svd,lwork_svd,info)
!
!       divide the columns of u_1 and the rows of v_1 by the square root of the singular values
!
        do i = 1, ldu
          u_1(:,i) = u_1(:,i) / sqrt(sv_1(i))
          vt_1(i,:) = vt_1(i,:) / sqrt(sv_1(i))
        end do
!
!       for the projected ep and em matrices:
!
        call dgemm('n','t',ldu,ldu,ldu,one,epmat,lda,vt_1,lda,zero,scr,lda)
        call dgemm('n','n',ldu,ldu,ldu,one,vt_1,lda,scr,lda,zero,ept,lda)
        call dgemm('n','n',ldu,ldu,ldu,one,emmat,lda,u_1,lda,zero,scr,lda)
        call dgemm('t','n',ldu,ldu,ldu,one,u_1,lda,scr,lda,zero,emt,lda)
!
!       compute their cholesky decomposition:
!
        call dpotrf('l',ldu,ept,lda,info)
        call dpotrf('l',ldu,emt,lda,info)
!
!       assemble c = (l^-)^t l^+
!
        do j = 1, ldu
          do i = 1, ldu
            cmat(i,j) = zero
            do k = max(i,j), ldu
              cmat(i,j) = cmat(i,j) + emt(k,i) * ept(k,j)
            end do
          end do
        end do
!
!       compute its svd:
!
        call dgesvd('a','a',ldu,ldu,cmat,lda,sv_2,u_2,lda,vt_2,lda,work_svd,lwork_svd,info)
!
!       assemble xpt and xmt
!
        do i = 2, ldu
          do j = 1, i-1
            ept(j,i) = zero
            emt(j,i) = zero
          end do
        end do
        call dgemm('n','n',ldu,ldu,ldu,one,emt,lda,u_2,lda,zero,scr,lda)
        call dgemm('t','n',ldu,ldu,ldu,one,vt_1,lda,scr,lda,zero,xpt,lda)
        call dgemm('n','t',ldu,ldu,ldu,one,ept,lda,vt_2,lda,zero,scr,lda)
        call dgemm('n','n',ldu,ldu,ldu,one,u_1,lda,scr,lda,zero,xmt,lda)
!
!       normalize and gather the eigenvalues and eigenvectors:
!
        do i_eig = 1, n_max
          eig(i_eig)  = sv_2(ldu - i_eig + 1)
          up(1:ldu,i_eig) = xpt(1:ldu,ldu - i_eig + 1)/(sqrt(two) * sv_2(ldu - i_eig + 1))
          um(1:ldu,i_eig) = xmt(1:ldu,ldu - i_eig + 1)/(sqrt(two) * sv_2(ldu - i_eig + 1))
        end do
!
!       gather the eigenvalues
!
      end if
      call get_time(t2)
      t_diag = t_diag + t2 - t1
      call dgemm('n','n',n,n_max,ldu,one,vp,n,up,lda,zero,eigp,n)
      call dgemm('n','n',n,n_max,ldu,one,vm,n,um,lda,zero,eigm,n)
!
      do i_eig = 1, n_max
        evec(1:n,i_eig)    = eigp(:,i_eig) + eigm(:,i_eig)
        evec(n+1:n2,i_eig) = eigp(:,i_eig) - eigm(:,i_eig)
      end do
!
!     compute the residuals, and their rms and sup norms:
!
      call dgemm('n','n',n,n_max,ldu,one,lvp,n,up,lda,zero,rp,n)
      call dgemm('n','n',n,n_max,ldu,one,lvm,n,um,lda,zero,rm,n)
      call dgemm('n','n',n,n_max,ldu,one,bvp,n,um,lda,zero,bp,n)
      call dgemm('n','n',n,n_max,ldu,one,bvm,n,up,lda,zero,bm,n)

      do i_eig = 1, n_targ
!
!       if the eigenvalue is already converged, skip it.
!
        if (done(i_eig)) cycle
!
        call daxpy(n,-eig(i_eig),bp(:,i_eig),1,rp(:,i_eig),1)
        call daxpy(n,-eig(i_eig),bm(:,i_eig),1,rm(:,i_eig),1)
        r_norm(1,i_eig) = dnrm2(n,rp(:,i_eig),1)/sqrtn + dnrm2(n,rm(:,i_eig),1)/sqrtn
        r_norm(2,i_eig) = maxval(abs(rp(:,i_eig))) + maxval(abs(rm(:,i_eig)))
      end do
!
!     check convergence. lock the first contiguous converged eigenvalues
!     by setting the logical array "done" to true.
!
      do i_eig = 1, n_targ
        if (done(i_eig)) cycle
        done(i_eig)     = r_norm(1,i_eig).lt.tol_rms .and. &
                          r_norm(2,i_eig).lt.tol_max .and. &
                          it.gt.1
        if (.not.done(i_eig)) then
          done(i_eig+1:n_max) = .false.
          exit
        end if
      end do
!
!     print some information:
!
      if (verbose) then
        do i_eig = 1, n_targ
          write(6,1040) it, i_eig, eig(i_eig), r_norm(:,i_eig), done(i_eig)
        end do
        write(6,*)
      end if
!
      if (all(done(1:n_targ))) then
        ok = .true.
        exit
      end if
!
!     check whether an update is required.
!     if not, perform a davidson restart.
!
      if (m_dim .lt. dim_dav) then
!
!       compute the preconditioned residuals using davidson's procedure
!       note that this is done with a user-supplied subroutine, that can
!       be generalized to experiment with fancy preconditioners that may
!       be more effective than the diagonal one, as in the original
!       algorithm.
!
        m_dim = m_dim + 1
        i_beg = i_beg + n_act
        n_act = n_max
        n_frozen = 0
        do i_eig = 1, n_targ
          if (done(i_eig)) then
            n_act = n_act - 1
            n_frozen = n_frozen + 1
          else
            exit
          end if
        end do
        ind   = n_max - n_act + 1
        call lrprec(n,n_act,eig(ind),rp(1,ind),rm(1,ind),vp(1,i_beg),vm(1,i_beg))
!
!       orthogonalize the new vectors to the existing ones and then
!       orthonormalize them.
!
        call get_time(t1)
!        call ortho_vs_x(.false.,n,ldu,n_act,vp,vp(1,i_beg),xx,xx)
        call ortho_vs_x(n,ldu,n_act,vp,vp(1,i_beg),xx,xx)
!        call ortho_vs_x(.false.,n,ldu,n_act,vm,vm(1,i_beg),xx,xx)
        call ortho_vs_x(n,ldu,n_act,vm,vm(1,i_beg),xx,xx)
        call get_time(t2)
        t_ortho = t_ortho + t2 - t1
      else
        if (verbose) write(6,'(t7,a)') 'Restarting davidson.'
        restart = .true.
!
!       initialize indexes back to their starting values
!
        ldu   = 0
        i_beg = 1
        m_dim = 1
        n_rst = 0
!
        n_act = n_max
        vp = zero
        vm = zero
!
!       put current eigenvectors into the first position of the
!       expansion space
!
        do i_eig = 1, n_max
          vp(:,i_eig) = evec(1:n,i_eig) + evec(n+1:n2,i_eig)
          vm(:,i_eig) = evec(1:n,i_eig) - evec(n+1:n2,i_eig)
        end do
!
!        call ortho_cd(.false.,n,n_max,vp,xx,ok)
        call ortho_cd(n,n_max,vp,growth,ok)
!        call ortho_cd(.false.,n,n_max,vm,xx,ok)
        call ortho_cd(n,n_max,vm,growth,ok)
!
        lvp   = zero
        lvm   = zero
        bvp   = zero
        bvm   = zero
        a_red = zero
        s_red = zero
        epmat = zero
        emmat = zero
        smat  = zero
!
      end if
      if (verbose) write(6,1050) n_targ, n_act, n_frozen
    end do
!
    call get_time(t1)
    t_tot = t1 - t_tot
    1000 format(t3,'timings for caslr (cpu/wall):   ',/, &
                t3,'  matrix-vector multiplications: ',2f12.4,/, &
                t3,'  diagonalization:               ',2f12.4,/, &
                t3,'  orthogonalization:             ',2f12.4,/, &
                t3,'                                 ',24('='),/,  &
                t3,'  total:                         ',2f12.4)
    if (verbose) write(6,1000) t_mv, t_diag, t_ortho, t_tot
    deallocate(work,tau,vp,vm,lvp,lvm,bvp,bvm,rp,rm,rr,done,r_norm,a_red,a_copy,s_red,s_copy,e_red, &
               epmat,emmat,smat,up,um,eigp,eigm,bp,bm)
!
1050 format(t5,'----------------------------------------',/,&
            t7,'# target vectors:    ',i4,/,&
            t7,'# new vectors added: ',i4,/,&
            t7,'# converged vectors: ',i4,/,&
            t5,'----------------------------------------')
    return
  end subroutine caslr_driver
!
  subroutine caslr_std(verbose,n,n2,n_targ,n_max,max_iter,tol,max_dav, &
                          apbmul,ambmul,spdmul,smdmul,lrprec,vec,ok,omega,g_half,imag)
    use utils
    implicit none
    logical,                       intent(in)    :: verbose, imag
    integer,                       intent(in)    :: n, n2, n_targ, n_max
    integer,                       intent(in)    :: max_iter, max_dav
    real(dp),                      intent(in)    :: omega 
    real(dp),                      intent(in)    :: tol
    real(dp), dimension(n2,n_max), intent(inout) :: vec
    real(dp), dimension(n,n_max),  intent(in)    :: g_half
    logical,                       intent(inout) :: ok
    external                                     :: apbmul, ambmul, spdmul, smdmul, &
                                                    lrprec
!
!   local variables:
!   ================
!
    integer, parameter    :: min_dav = 10
    integer               :: istat
!
!   actual expansion space size and total dimension
!
    integer               :: dim_dav, lda, lda2
!
!   number of active vectors at a given iteration, and indices to access them
!
    integer               :: n_act, ind, i_beg
!
!   current size and total dimension of the expansion space
!
    integer               :: m_dim, ldu
!
!   number of frozen (i.e. converged) vectors
!
    integer               :: n_frozen
!
    integer               :: i, j, k, it, i_eig
!
!
    real(dp)              :: sqrtn, tol_rms, tol_max, growth
    real(dp)              :: xx(1), lw_svd(1)
!
!   arrays to control convergence and orthogonalization:
!
    logical,  allocatable :: done(:)
!
!   expansion spaces, residuals and their norms:
!
    real(dp), allocatable :: vp(:,:), vm(:,:), lvp(:,:), lvm(:,:), bvp(:,:), bvm(:,:)
    real(dp), allocatable :: rp(:,:), rm(:,:), rr(:,:), r_norm(:,:)
!
!   solution vectors of the reduced problem and components of the ritz vectors:
!
    real(dp), allocatable :: up(:,:), um(:,:), bp(:,:), bm(:,:), vecp(:,:), vecm(:,:)
!
!   subspace matrices:
!
    real(dp), allocatable :: epmat(:,:), emmat(:,:), smat(:,:), m_red(:,:)
!
!   subspace gradient vectors:
    real(dp), allocatable :: gp(:,:), gm(:,:), gpm(:,:)
!
!   restarting variables
!
    integer               :: n_rst
    logical               :: restart
!
!   external functions:
!   ===================
!
    real(dp)              :: dnrm2
    external              :: dcopy, dnrm2, dgemm, dgesv
!
!   external arrays and integers for Lapacks
!
    integer               :: info
    real(dp), allocatable :: ipiv(:)
    real(dp)              :: orto
!   compute the actual size of the expansion space, checking that
!   the input makes sense.
!   no expansion space smaller than max_dav = 10 is deemed acceptable.
!
    dim_dav = max(min_dav,max_dav)
    lda     = dim_dav*n_max
    lda2    = 2 * lda
!
!   start by allocating memory for the various lapack routines
!
    allocate (ipiv(n2), tau(n_max), stat=istat)
    call check_mem(istat)
!
!   allocate memory for the expansion space, the corresponding 
!   matrix-multiplied vectors and the residual:
!
    allocate (vp(n,lda), vm(n,lda), lvp(n,lda), lvm(n,lda), bvp(n,lda), bvm(n,lda), &
              rp(n,n_max), rm(n,n_max), rr(n,n_max), stat = istat)
!
!   allocate memory for the reduced arrays:
    allocate (gp(lda,n_max), gm(lda,n_max), gpm(lda2,n_max), m_red(lda2,lda2), &
              epmat(lda,lda), emmat(lda,lda), smat(lda,lda), stat=istat)
!
    call check_mem(istat)
!
!   allocate memory for convergence check
!
    allocate (done(n_max), r_norm(2,n_max), stat=istat)
    call check_mem(istat)
!
!   allocate memory for the plus and minus vector components:
!
    allocate (up(lda,n_max), um(lda,n_max), vecp(n,n_max), vecm(n,n_max), bp(n,n_max), bm(n,n_max), stat = istat)
    call check_mem(istat)
!
!   set the tolerances and compute a useful constant to compute rms norms:
!
    sqrtn   = sqrt(real(n,dp))
    tol_rms = tol
    tol_max = 10.0_dp * tol
!
!   clean out various quantities
!
    t_ortho = zero
    t_mv    = zero
    t_ls    = zero
    t_tot   = zero
    vp      = zero
    vm      = zero
    bvp     = zero
    bvm     = zero
    lvp     = zero
    lvm     = zero
    epmat   = zero
    smat    = zero
    emmat   = zero
    m_red   = zero
    gpm     = zero
    gp      = zero
    gm      = zero
    up      = zero
    um      = zero
    ok      = .false.
    done    = .false.
!
    call get_time(t_tot)
!
!   move the guess into the expansion space.
!
    do i_eig = 1, n_max
      vp(:,i_eig) = vec(1:n,i_eig) + vec(n+1:n2,i_eig)
      vm(:,i_eig) = vec(1:n,i_eig) - vec(n+1:n2,i_eig)
    end do
!
    call get_time(t1)
    call ortho_cd(n,n_max,vp,growth,ok)
    call ortho_cd(n,n_max,vm,growth,ok)
    call get_time(t2)
    t_ortho = t_ortho + t2 - t1
!
    n_act = n_max
    ind   = 1
    i_beg = 1
!
!   initialize the counter for the expansion of the subspace
!
    m_dim = 1
    ldu   = 0
!
!   initialize to false the restart
!
    restart = .false.
!
!   main loop:
!
    1030 format(t5,'Linear system iterations (tol=',d10.2,'):',/, &
                t5,'------------------------------------------------------------------',/, &
                t7,'  iter  vec      rms          max     ok'   ,/, &
                t5,'------------------------------------------------------------------')
    1040 format(t9,i4,2x,i4,2d12.4,l3)
!
!
    if (verbose) write(6,1030) tol
!
    n_rst   = 0
    do it = 1, max_iter
!
!     update the size of the expansion space.
!
      ldu = ldu + n_act
!
!     perform this iteration's matrix-vector multiplication:
!
      call get_time(t1)
      call apbmul(n,n_act,vp(1,i_beg),lvp(1,i_beg))
      call ambmul(n,n_act,vm(1,i_beg),lvm(1,i_beg))
      call spdmul(n,n_act,vp(1,i_beg),bvm(1,i_beg))
      call smdmul(n,n_act,vm(1,i_beg),bvp(1,i_beg))
      call get_time(t2)
      t_mv = t_mv + t2 - t1
!
!     update the reduced matrix 
!
      epmat = zero
      emmat = zero
      smat  = zero
      call dgemm('t','n',ldu,ldu,n,one,vp,n,lvp,n,zero,epmat,lda)
      call dgemm('t','n',ldu,ldu,n,one,vm,n,lvm,n,zero,emmat,lda)
      call dgemm('t','n',ldu,ldu,n,one,vm,n,bvm,n,zero,smat,lda)
!
!     because of the symmetry of g and depending on the nature of the perturbation
!     it is either gp = 0 or gm = 0
!
      gp = zero
      gm = zero
      if (imag) then
        call dgemm('t','n',ldu,n_max,n,one,vm,n,g_half,n,zero,gm,lda)
      else
        call dgemm('t','n',ldu,n_max,n,one,vp,n,g_half,n,zero,gp,lda)
      end if
!
!     build the 2n-dimensional matrix
!     
      m_red = zero
      m_red(1:ldu,1:ldu)             = epmat(1:ldu,1:ldu)
      m_red(ldu+1:2*ldu,ldu+1:2*ldu) = emmat(1:ldu,1:ldu)
      m_red(1:ldu,ldu+1:2*ldu)       =-omega*transpose(smat(1:ldu,1:ldu))
      m_red(ldu+1:2*ldu,1:ldu)       =-omega*smat(1:ldu,1:ldu)
!
!     build the 2n-dimensional gradient vector
!
      gpm                = zero
      gpm(1:ldu,:)       = gp(1:ldu,:)
      gpm(ldu+1:2*ldu,:) = gm(1:ldu,:)
!
!
!     default algorithm: solve the 2n-dimensional linear system.
!
      call get_time(t1)
        call dgesv(2*ldu,n_max,m_red,lda2,ipiv,gpm,lda2,info)
      call get_time(t2)
!
!     extract the coefficients by column
!
      up          = zero
      um          = zero
      up(1:ldu,:) = gpm(1:ldu,:)
      um(1:ldu,:) = gpm(ldu+1:2*ldu,:)
!
!
      t_ls = t_ls + t2 - t1
!
!     compute the ritz approximation to the solution vectors
!         
      vecp = zero
      vecm = zero
      call dgemm('n','n',n,n_max,ldu,one,vp,n,up,lda,zero,vecp,n)
      call dgemm('n','n',n,n_max,ldu,one,vm,n,um,lda,zero,vecm,n)
!
!     put the expansion vectors in the expansion subspace
!
      vec = zero
      do i = 1, n_max
        vec(1:n,i)    = vecp(:,i) + vecm(:,i)
        vec(n+1:n2,i) = vecp(:,i) - vecm(:,i)
      end do
!
!     compute the residuals, and their rms and sup norms:
!
      rp = zero
      rm = zero
      bp = zero
      bm = zero
      call dgemm('n','n',n,n_max,ldu,one,lvp,n,up,lda,zero,rp,n)
      call dgemm('n','n',n,n_max,ldu,one,lvm,n,um,lda,zero,rm,n)
      call dgemm('n','n',n,n_max,ldu,one,bvp,n,um,lda,zero,bp,n)
      call dgemm('n','n',n,n_max,ldu,one,bvm,n,up,lda,zero,bm,n)
!     
!     the construction of the residuals is influenced by the type of perturbation
!
      if (imag) then
        rm = rm - g_half
      else
        rp = rp - g_half
      end if
!
      do i_eig = 1, n_targ
        if (done(i_eig)) cycle
!
        call daxpy(n,-omega,bp(:,i_eig),1,rp(:,i_eig),1)
        call daxpy(n,-omega,bm(:,i_eig),1,rm(:,i_eig),1)
        r_norm(1,i_eig) = dnrm2(n,rp(:,i_eig),1)/sqrtn + dnrm2(n,rm(:,i_eig),1)/sqrtn
        r_norm(2,i_eig) = maxval(abs(rp(:,i_eig))) + maxval(abs(rm(:,i_eig)))
      end do
!
! check convergence. lock the first contiguous converged eigenvalues
!     by setting the logical array "done" to true.
!
      do i_eig = 1, n_targ
        if (done(i_eig)) cycle
        done(i_eig)     = r_norm(1,i_eig).lt.tol_rms .and. &
                          r_norm(2,i_eig).lt.tol_max .and. &
                          it.gt.1
        if (.not.done(i_eig)) then
          done(i_eig+1:n_max) = .false.
          exit
        end if
      end do
!
!     print some information:
!
      if (verbose) then
        do i_eig = 1, n_targ
          write(6,1040) it, i_eig, r_norm(:,i_eig), done(i_eig)
        end do
        write(6,*) 
      end if
!
      if (all(done(1:n_targ))) then
        ok = .true.
        exit
      end if
!
!     check whether an update is required. 
!     if not, perform a davidson restart.
!
      if (m_dim .lt. dim_dav) then
!
!       compute the preconditioned residuals using davidson's procedure
!       note that this is done with a user-supplied subroutine, that can
!       be generalized to experiment with fancy preconditioners that may
!       be more effective than the diagonal one, as in the original 
!       algorithm.
!
        m_dim = m_dim + 1
        i_beg = i_beg + n_act
        n_act = n_max
        n_frozen = 0
        do i_eig = 1, n_targ
          if (done(i_eig)) then
            n_act = n_act - 1
            n_frozen = n_frozen + 1
          else
            exit
          end if
        end do
        ind   = n_max - n_act + 1
        call lrprec(n,n_act,omega,rp(1,ind),rm(1,ind),vp(1,i_beg),vm(1,i_beg))
!
!       orthogonalize the new vectors to the existing ones and then
!       orthonormalize them.
!
        call get_time(t1)
        call ortho_vs_x(n,ldu,n_act,vp,vp(1,i_beg),xx,xx,dropping=.true.)
        call ortho_vs_x(n,ldu,n_act,vm,vm(1,i_beg),xx,xx,dropping=.true.)
        call get_time(t2)
        t_ortho = t_ortho + t2 - t1
      else
        if (verbose) write(6,'(t7,a)') 'Restarting davidson.'
        restart = .true.
!
!       initialize indexes back to their starting values 
!
        ldu   = 0
        i_beg = 1
        m_dim = 1
        n_rst = 0
!
        n_act = n_max 
        vp = zero
        vm = zero
!
!      put current vectors into the first position of the 
!       expansion space
!
        do i_eig = 1, n_max
          vp(:,i_eig) = vec(1:n,i_eig) + vec(n+1:n2,i_eig)
          vm(:,i_eig) = vec(1:n,i_eig) - vec(n+1:n2,i_eig)
        end do
!
        call get_time(t1)
        call ortho_cd(n,n_max,vp,growth,ok)
        call ortho_cd(n,n_max,vm,growth,ok)
        call get_time(t2)
        t_ortho = t_ortho + t2 - t1
!
        lvp    = zero
        lvm    = zero
        bvp    = zero
        bvm    = zero
        epmat  = zero
        emmat  = zero
        smat   = zero
        m_red  = zero
        gp     = zero
        gm     = zero
        gpm    = zero
!
      end if
      if (verbose) write(6,1050) n_targ, n_act, n_frozen
    end do
!
    call get_time(t1)
    t_tot = t1 - t_tot
   1000 format(t3,'timings for caslr (cpu/wall):   ',/, &
                t3,'  matrix-vector multiplications: ',2f12.4,/, &
                t3,'  orthogonalization:             ',2f12.4,/, &
                t3,'  reduced system solving:        ',2f12.4,/, &
                t3,'                                 ',24('='),/,  &
                t3,'  total:                         ',2f12.4)
    if (verbose) write(6,1000) t_mv, t_ortho, t_ls, t_tot
    deallocate(ipiv,tau,vp,vm,lvp,lvm,bvp,bvm,rp,rm,rr,done,r_norm, &
               epmat,emmat,smat,up,um,vecp,vecm,bp,bm,gp,gm,gpm,m_red)
!
1050 format(t5,'----------------------------------------',/,&
            t7,'# target vectors:    ',i4,/,&
            t7,'# new vectors added: ',i4,/,&
            t7,'# converged vectors: ',i4,/,&
            t5,'----------------------------------------')
    return
  end subroutine caslr_std
!
  subroutine caslr_eff_std(verbose,n,n2,n_targ,n_max,max_iter,tol,max_dav, &
                           apbmul,ambmul,spdmul,smdmul,lrprec,vec,ok,omega,g_half,imag)
!
!   Main driver for the efficient solution to the following standard response equation:
!
!   / A  B \ / Y \     /  S  D \ / Y \   /  Q  \
!   |      | |   | - w |       | |   | = |     |
!   \ B  A / \ Z /     \ -D -S / \ Z /   \  Q* /
!
!   where A, B, S are symmetric matrices and D is antysimmetric, while Q* is 
!   either Q (real perturbation) or -Q (imaginary perturbation).
!
!   If (w, Y, Z) are a solution, then (-w, Z, Y) is also a solution. 
!   Following J. Chem. Phys., 118, 522 (2003), we enforce this property in 
!   the iterative procedure by expanding the eigenvector as
!
!   (Y, Z) = (b+, b+) + (b-, - b-)
!
!   Since the metric is symmetric and positive definite, it's possible
!   to use expansion vectors that are orthogonal with respect to the
!   dot product defined by the metric, which in turn results in a
!   Rayleigh-Ritz procedure that requires the solution of the reduced problem:
!
!   / 1  -wS^t \ / u+ \    / g+ \
!   |          | |    | =  |    |
!   \ -wS   1  / \ u- /    \ g- /
!
!   which can be reduced to a half-sized eigenvalue problem
!
!   u+ = (1 - w^2 S^t S)^(-1) * (g+ + w S^t g-)
!   u- = g- + w S u+
!
!   input variables
!   ===============
!
!   verbose:  logical, whether to print various information at each 
!             iteration (eigenvalues, residuals...).
!
!   imag:     logical, true if the perturbation is imaginary.
!
!   n:        integer, size of the A, B, S, D matrices
!
!   n2:       integer, n2 = 2*n, size of the generalized eigenvalue
!             problem.
!
!   n_targ:   integer, number of required solutions.
!
!   n_max:    integer, maximum size of the search space. It should be 
!             >= n_targ. 
!   
!   max_iter: integer, maximum allowed number of iterations.
!
!   max_dav:  integer, maximum allowed number of iterations before a
!             Davidson restart is forced. When n_max eigenvalues are 
!             searched, this implies a maximum dimension of the 
!             expansion subspacem of n_max*max_dav.
!
!   omega:    double precision real, frequency in a.u.
!
!   tol:      double precision real, the convergence threshold.
!
!   g_half:   double precision real array of size (n,n_max),
!             represents (half) the gradient.
!
!   apbmul:   external subroutine to apply (A+B) to a vector.
!
!   ambmul:   external subroutine to apply (A-B) to a vector.
!
!   spdmul:   external subroutine to apply (S+D) to a vector.
!
!   smdmul:   external subroutine to apply (S-D) to a vector.
!
!   lrprec:   external subroutine that applies a preconditioner.
!
!   output variables:
!   =================
!
!   vec:     double precision real array of size (n2,n_max).
!            In input, a guess for the solution vectors.
!            If ok is true, in output the computed vectors.
!
!   ok:      logical, true if caslr_eff_std converged.
!
    use utils
    implicit none
    logical, intent(in)                          :: verbose, imag
    integer,                       intent(in)    :: n, n2, n_targ, n_max
    integer,                       intent(in)    :: max_iter, max_dav
    real(dp),                      intent(in)    :: omega
    real(dp),                      intent(in)    :: tol
    real(dp), dimension(n,n_max),  intent(in)    :: g_half
    external                                     :: apbmul, ambmul, spdmul, smdmul, &
                                                    lrprec
    real(dp), dimension(n2,n_max), intent(inout) :: vec
    logical,                       intent(inout) :: ok
!
!   local variables:
!   ================
!
!   min_dav: integer, minimum acceptable dimension of the expansion subspace.
!
!   istat: integer, used as the value of the parameter "stat" to manage
!          allocation errors.
!    
!   dim_dav: integer, actual expansion space size.
!
!   lda,lda2: integers, total dimension of the expansion arrays.
! 
!   n_act, ind, i_beg: integers, number of active vectors at a given iteration, 
!                      and indices to access them.
!
!   n_p, n_m: integers, number of non-dropped vectors in the "symmetric"
!             and "antisymmetric" subspaces.
!
!   m_dim, ldu: integers, respectively the current Davidson iteration and
!               the current dimension.
!
!   n_frozen: integer, number of frozen (i.e. converged) vectors.
!
!   n_drop: integer, number of dropped expansion vectors.
!    
!   it, i_std: integers, used for loops.
!
!   sqrtn: double precision real, used to compute rms norms, coupled 
!          with the intrinsic function "real".
!
!   tol_rms, tol_max: double precision reals, set to tol and 10*tol, respectively;
!                     used to check convergence.
!
!   done: array of logicals used to control convergence 
!         and orthogonalization.
!
!   vp, vm: real double precision arrays, the symmetric (v+) and 
!           antisymmetric (v-) expansion spaces, respectively.
!
!   lvp, lvm: real double precision arrays, obtained by applying (A+B) to vp
!             and (A-B) to vm, respectively.
!
!   bvp, bvm: real double precision arrays, obtained by applying (S+D)
!             to vp and (S-D) to vm. 
! 
!   rp, rm, r_norm: real double precision arrays: the symmetric (rp), antisymmetric
!                   (rm) residuals and their overall norm.
!
!   um: real double precision array, containing the antisymmetric 
!       coefficients of the Ritz vector, obtained by solving 
!       the reduced problem.
!
!   bp, bm: real double precision arrays, obtained by matrix-matrix 
!           multiplication of bvp and vm in the first case, bvm 
!           and vp in the second.
!
!   vecp, vecm: real double precision arrays, the symmetric and antisymmetric
!               components of the Ritz vector. They are the solutions of the 
!               matrix-matrix products vp and up, vm and um, respectively.
!
!   smat: real double precision array, S. Obtained by applying 
!         (vm^T) to bvm. It occupies the top right and the bottom left
!         positions in the reduced problem's matrix.
!
!   ss_mat: real double precision array, the identity matrix minus 
!           (w^2)*(s^t)s, where w is omega.
!
!   gp,gm: real double precision arrays, the subspace gradient vectors; 
!          gp and gm are the projections of the gradient in the symmetric
!          and antisymmetric subspaces, respectively.
!
!   s_up,s_gm: real double precision arrays, used to construct the 
!              coefficients up and um; s_up is obtained by applying s to s_gm 
!              which is initially w*s^t gm + gp, then takes the role of the 
!              symmetric up coefficients, once the reduced linear system
!              has been solved. 
!
!   restart: logical, gives Davidson' routine a restart order. It's set to
!            true when m_dim is equal to m_max.
!
!   cnt: integer, used to count every time the program enters a 
!        particular loop.
!
    integer, parameter    :: min_dav = 10
    integer               :: istat
    integer               :: dim_dav, lda, lda2
    integer               :: n_act, ind, i_beg, n_p, n_m
    integer               :: m_dim, ldu
    integer               :: n_frozen, n_drop
    integer               :: it, i_std
    real(dp)              :: sqrtn, tol_rms, tol_max
    logical,  allocatable :: done(:)
    real(dp), allocatable :: vp(:,:), vm(:,:), lvp(:,:), lvm(:,:), bvp(:,:), bvm(:,:)
    real(dp), allocatable :: rp(:,:), rm(:,:), r_norm(:,:)
    real(dp), allocatable :: um(:,:), bp(:,:), bm(:,:), vecp(:,:), vecm(:,:)
    real(dp), allocatable :: smat(:,:), ss_mat(:,:)
    real(dp), allocatable :: gp(:,:), gm(:,:)
    real(dp), allocatable :: s_up(:,:), s_gm(:,:)
    logical               :: restart
    integer               :: cnt
!
!   External functions:
!   ===================
!
    real(dp)              :: dnrm2
    external              :: dcopy, dnrm2, dgesv, dgetri, dgetrfi, dgemv
!
!   External arrays and integers used for Lapacks.
!   ==============================================
!
    integer               :: info
    real(dp), allocatable :: ipiv(:)
!
!   Compute the actual size of the expansion space, checking that
!   the input makes sense.
!   No expansion space smaller than min_dav = 10 is deemed acceptable.
!
    dim_dav = max(min_dav,max_dav)
    lda     = dim_dav*n_max
    lda2    = 2 * lda
!
!   Start by allocating memory for the various lapack routines:
!
    allocate (tau(n_max), ipiv(n), stat=istat)
    call check_mem(istat)
!
!   Allocate memory for the expansion space, the corresponding 
!   matrix-multiplied vectors and residuals:
!
    allocate (vp(n,lda), vm(n,lda), lvp(n,lda), lvm(n,lda), bvp(n,lda), bvm(n,lda), &
              rp(n,n_max), rm(n,n_max), stat = istat)
    call check_mem(istat)
!
!   Allocate memory for convergence check:
!
    allocate (done(n_max), r_norm(2,n_max), stat=istat)
    call check_mem(istat)
!
!   Allocate memory for the reduced arrays:
!
    allocate (gp(lda,n_max), gm(lda,n_max), smat(lda,lda), stat=istat)
!
!   Allocate memory for the matrix-matrix products and copies of the reduced arrays:
!
    allocate (ss_mat(lda,lda), s_up(lda,n_max), s_gm(lda,n_max), stat=istat)
    call check_mem(istat)
!
    call check_mem(istat)
!
!   Allocate memory for the plus and minus vector components:
!
    allocate (um(lda,n_max), vecp(n,n_max), vecm(n,n_max), bp(n,n_max), bm(n,n_max), stat = istat)
    call check_mem(istat)
!
!   Set the tolerances and compute a useful constant to compute rms norms:
!
    sqrtn   = sqrt(real(n,dp))
    tol_rms = tol
    tol_max = 10.0_dp * tol
!
!   Clean out various quantities. The quantities t_* are used to indicate 
!   orthogonalization (t_ortho), matrix-vector multiplication (t_mv) and
!   linear system solving (t_ls) times.
!
    t_ortho      = zero
    t_mv         = zero
    t_ls         = zero
    t_tot        = zero
    vp           = zero
    vm           = zero
    bvp          = zero
    bvm          = zero
    lvp          = zero
    lvm          = zero
    ok           = .false.
    done         = .false.
!
!   Begin total run time. 
!
    call get_time(t_tot)
!
!   Use the guess to construct the symmetric (Y+Z) and 
!   antisymmetric (Y-Z) expansion spaces for the first
!   iteration.
!
    do i_std = 1, n_max
      vp(:,i_std) = vec(1:n,i_std) + vec(n+1:n2,i_std)
      vm(:,i_std) = vec(1:n,i_std) - vec(n+1:n2,i_std)
    end do
!
!   Lambda-orthogonalize the expansion space to the metric.
!   This is done firstly by constructing lvp and then by
!   orhogonalizing it to vp. The same is done for lvm, vm.
!
    call get_time(t1)
    call apbmul(n,n_max,vp,lvp)
    call get_time(t2)
    t_mv = t_mv + t2 - t1
    call get_time(t1)
    call b_ortho(n,n_max,vp,lvp)
    call get_time(t2)
    t_ortho = t_ortho + t2 - t1
    call get_time(t1)
    call ambmul(n,n_max,vm,lvm)
    call get_time(t2)
    t_mv = t_mv + t2 - t1
    call get_time(t1)
    call b_ortho(n,n_max,vm,lvm)
    call get_time(t2)
    t_ortho = t_ortho + t2 - t1
!   
!   Set the number of active vectors to the maximum value.  
!
    n_act = n_max
    ind   = 1
    i_beg = 1
!
!   Initialize the counter for the expansion of the subspace
!
    m_dim  = 1
    ldu    = 0
    n_p    = 0
    n_m    = 0
    n_drop = 0
!
!   Initialize to false the restart.
!
    restart = .false.
!
!   Main loop:
!
    1030 format(t5,'Linear system iterations (tol=',d10.2,'):',/, &
                t5,'------------------------------------------------------------------',/, &
                t7,'  iter  vec      rms          max     ok'   ,/, &
                t5,'------------------------------------------------------------------')
    1040 format(t9,i4,2x,i4,2d12.4,l3)
!
    if (verbose) write(6,1030) tol
!
    do it = 1, max_iter
!
!     Update the size of the expansion space.
!
      ldu = ldu + n_act
!
!     Perform this iteration's matrix-vector multiplications.
!
      call get_time(t1)
      call spdmul(n,n_act,vp(1,i_beg),bvm(1,i_beg))
      call smdmul(n,n_act,vm(1,i_beg),bvp(1,i_beg))
      call get_time(t2)
      t_mv = t_mv + t2 - t1
!
!     Update the reduced matrix smat (s). 
!
      smat = zero
      call dgemm('t','n',ldu,ldu,n,one,vm,n,bvm,n,zero,smat,lda)
!
!     Because of the symmetry of g and whether the perturbation is complex 
!     or real, discriminate between gp = 0 or gm = 0, respectively.
!
      gp = zero
      gm = zero
      cnt = 0
      if (imag) then
        do i_std = 1, n_max
          if (done(i_std)) cycle
          cnt = cnt + 1
          call dgemv('t',n,ldu,one,vm,n,g_half(:,i_std),1,zero,gm(:,cnt),1)
        end do
      else
        do i_std = 1, n_max
          if (done(i_std)) cycle
          cnt = cnt + 1
          call dgemv('t',n,ldu,one,vp,n,g_half(:,i_std),1,zero,gp(:,cnt),1)
        end do
      end if
!      
!     Assemble Id - (w^2)s^t s = ss_mat.
!
      ss_mat = zero
      call dgemm('t','n',ldu,ldu,ldu,-omega**2,smat,lda,smat,lda,zero,ss_mat,lda)
      do i_std = 1, ldu
        ss_mat(i_std,i_std) = ss_mat(i_std,i_std) + 1.d0
      end do
!
!     Assemble w*s^t gm + gp = s_gm.
!
      s_gm = zero
      if (imag) then
      call dgemm('t','n',ldu,n_act,ldu,omega,smat,lda,gm,lda,zero,s_gm,lda)
      end if
      s_gm(:,1:n_act) = s_gm(:,1:n_act) + gp(:,1:n_act)
!
!     Solve the reduced linear system: up = (1 - w^2 s^t s)^(-1) * (gp + w s^t gm).
!     The routine dgesv collects the solutions, the up coefficients, in s_gm.
!
      call get_time(t1)
      call dgesv(ldu,n_act,ss_mat,lda,ipiv,s_gm,lda,info)
      call get_time(t2)
      t_ls = t_ls + t2 - t1
!
!     Assemble (s)s_gm = s_up.
!
      s_up = zero
      call dgemm('n','n',ldu,n_act,ldu,one,smat,lda,s_gm,lda,zero,s_up,lda)
!
!     Build um = w*s_up + gm.
!
      um = zero
      do i_std = 1, n_act
        um(:,i_std) = gm(:,i_std) + (omega * s_up(:,i_std))
      end do
!
!     Assemble the symmetric (vecp) and antysimmetric (vecm) vectors 
!     as (vp)s_gm and (vm)um, respectively: the solution vectors are
!     the sum of these quantities.
!
      call dgemm('n','n',n,n_act,ldu,one,vp,n,s_gm,lda,zero,vecp,n)
      call dgemm('n','n',n,n_act,ldu,one,vm,n,um,lda,zero,vecm,n)
!
!     Assemble this iteration's Ritz vectors.
!
      cnt = 0
      do i_std = 1, n_max
        if (done(i_std)) cycle
        cnt = cnt + 1
        vec(1:n,i_std)    = vecp(:,cnt) + vecm(:,cnt)
        vec(n+1:n2,i_std) = vecp(:,cnt) - vecm(:,cnt)
      end do
!
!     Begin the costruction of the symmetric and antisymmetric residuals,
!     defined respectively as:
!     rp = (lvp)up - w*(bvp)um; rm = (lvm)um - w*(bvm)up.
!
      cnt = 0
      do i_std = 1, n_max
        if (done(i_std)) cycle
        cnt = cnt + 1
        call dgemv('n',n,ldu,one,lvp,n,s_gm(:,cnt),1,zero,rp(:,i_std),1)
        call dgemv('n',n,ldu,one,lvm,n,um(:,cnt),1,zero,rm(:,i_std),1)
        call dgemv('n',n,ldu,one,bvp,n,um(:,cnt),1,zero,bp(:,i_std),1)
        call dgemv('n',n,ldu,one,bvm,n,s_gm(:,cnt),1,zero,bm(:,i_std),1)
      end do
!
!     The construction of the residuals is again influenced by the type of perturbation.
!     If complex, rm contains the contribution of the gradient. Else, it is rp that 
!     depends on the gradient.
! 
      if (imag) then
        do i_std = 1, n_max
          if (done(i_std)) cycle
          rm(:,i_std) = rm(:,i_std) - g_half(:,i_std)
        end do
      else
        do i_std = 1, n_max
          if (done(i_std)) cycle
          rp(:,i_std) = rp(:,i_std) - g_half(:,i_std)
          end do
      end if
!     
!     Complete the construction of the residuals. If some vectors have already converged,
!     the corresponding residual won't be built. The corresponding dnrm norm and max norm
!     are also computed.
!
      do i_std = 1, n_max
         if (done(i_std)) cycle
!
         call daxpy(n,-omega,bp(:,i_std),1,rp(:,i_std),1)
         call daxpy(n,-omega,bm(:,i_std),1,rm(:,i_std),1)
         r_norm(1,i_std) = (dnrm2(n,rp(:,i_std),1) + dnrm2(n,rm(:,i_std),1))/(sqrt(two)*sqrtn)
         r_norm(2,i_std) = (maxval(abs(rp(:,i_std))) + maxval(abs(rm(:,i_std))))/(sqrt(two))
       end do
!
!     Check convergence by checking the norms. Sets to true the logical done in the position 
!     which corresponds to the converged vectors, in order to skip this check on the next
!     iteration.
!
      do i_std = 1, n_max
        if (done(i_std)) cycle
        done(i_std)     = r_norm(1,i_std).lt.tol_rms .and. &
                          r_norm(2,i_std).lt.tol_max .and. &
                          it.gt.1
      end do
!
!     Print some information.
!
      if (verbose) then
        do i_std = 1, n_max
          write(6,1040) it, i_std, r_norm(:,i_std), done(i_std)
        end do
        write(6,*) 
      end if
!
!     If all vectors have converged, set ok to true and exit.
!
      if (all(done(1:n_max))) then
        ok = .true.
        exit
      end if
!
!     Check whether an update is required. If not, perform a Davidson restart.
!
      if (m_dim .lt. dim_dav) then
!
        m_dim = m_dim + 1
        i_beg = i_beg + n_act
        n_act = n_max 
!
!       Update the number of frozen vectors. If there are, n_act is modified, so
!       that no more operations will include them.
!
        n_frozen = 0
        do i_std = 1, n_max
          if (done(i_std)) then
            n_act = n_act - 1
            n_frozen = n_frozen + 1
          end if
        end do
        ind = n_max - n_act + 1
!        
!       Compute the preconditioned residuals using Davidson's procedure.
!       Note that this is done with a user-supplied subroutine, that can
!       be generalized to experiment with fancy preconditioners that may
!       be more effective than the diagonal one, as in the original 
!       algorithm.
!
        call lrprec(n,n_act,omega,rp(1,ind),rm(1,ind),vp(1,i_beg),vm(1,i_beg))
!
!       B-orthogonalize the new vectors to the existing ones, with respect to
!       the metric, then orthonormalize them.
!
        n_drop = 0
        n_p    = n_act
        n_m    = n_act
        call get_time(t1)
        call b_ortho_vs_x(n,ldu,n_p,vp,lvp,vp(1,i_beg),dropping=.true.,tol_o=tol)
        call get_time(t2)
        t_ortho = t_ortho + t2 - t1
!
        call get_time(t1)
        call apbmul(n,n_p,vp(1,i_beg),lvp(1,i_beg))
        call get_time(t2)
        t_mv = t_mv + t2 - t1
        call get_time(t1)
        call b_ortho(n,n_p,vp(1,i_beg),lvp(1,i_beg))
        call b_ortho_vs_x(n,ldu,n_m,vm,lvm,vm(1,i_beg),dropping=.true.,tol_o=tol)
        call get_time(t2)
        t_ortho = t_ortho + t2 - t1
!
!       Set n_act equal to the smaller of the numbers n_p and n_m, in order
!       to eventually keep the dimension of the two spaces equal.
!
        n_drop = n_act - min(n_p,n_m)
        n_act  = min(n_p,n_m)
!
        call get_time(t1)
        call ambmul(n,n_act,vm(1,i_beg),lvm(1,i_beg))
        call get_time(t2)
        t_mv = t_mv + t2 - t1
        call get_time(t1)
        call b_ortho(n,n_act,vm(1,i_beg),lvm(1,i_beg))
        call get_time(t2)
        t_ortho = t_ortho + t2 - t1
      else
!
!       Once dim_dav has been reached, restart Davidson.
!
        if (verbose) write(6,'(t7,a)') 'Restarting davidson.'
        restart = .true.
!
!       Initialize indexes back to their starting values. 
!
        ldu   = 0
        i_beg = 1
        m_dim = 1
!
        n_frozen = 0
        do i_std = 1, n_max
          if (done(i_std)) then
            n_act = n_act - 1
            n_frozen = n_frozen + 1
          end if
        end do
        n_act = n_max - n_frozen
        print *, "n_act=", n_act
!
!       Put the current solution vectors into the first position of the 
!       expansion spaces.
!
        vp = zero
        vm = zero
        cnt = 0
        do i_std = 1, n_max
          if(done(i_std)) cycle
          cnt = cnt +1 
          vp(:,cnt) = vec(1:n,i_std) + vec(n+1:n2,i_std)
          vm(:,cnt) = vec(1:n,i_std) - vec(n+1:n2,i_std)
        end do
!
        lvp   = zero
        lvm   = zero
!
        call get_time(t1)
        call apbmul(n,n_act,vp,lvp)
        call get_time(t2)
        t_mv = t_mv + t2 - t1
        call get_time(t1)
        call b_ortho(n,n_act,vp,lvp)
        call get_time(t2)
        t_ortho = t_ortho + t2 - t1
        call get_time(t1)
        call ambmul(n,n_act,vm,lvm)
        call get_time(t2)
        t_mv = t_mv + t2 - t1
        call get_time(t1)
        call b_ortho(n,n_act,vm,lvm)
        call get_time(t2)
        t_ortho = t_ortho + t2 - t1
!        
        bvp      = zero
        bvm      = zero
        smat     = zero
        ss_mat   = zero
        gp       = zero
        gm       = zero
        s_up     = zero
        s_gm     = zero
        n_drop   = 0
!
      end if
!      
      if (verbose) write(6,1050) n_max, n_act, n_frozen
!     
    end do
!
    call get_time(t1)
    t_tot = t1 - t_tot
    1000 format(t3,'timings for caslr_eff (cpu/wall):   ',/, &
                t3,'  matrix-vector multiplications: ',2f12.4,/, &
                t3,'  orthogonalization:             ',2f12.4,/, &
                t3,'  reduced system solving:        ',2f12.4,/, &
                t3,'                                 ',24('='),/,  &
                t3,'  total:                         ',2f12.4)
    if (verbose) write(6,1000) t_mv, t_ortho, t_ls, t_tot
    deallocate(ipiv,tau,vp,vm,lvp,lvm,bvp,bvm,rp,rm,done,r_norm,ss_mat,smat, &
               um,vecp,vecm,bp,bm,s_gm)
!
1050 format(t5,'----------------------------------------',/,&
            t7,'# target vectors:    ',i4,/,&
            t7,'# new vectors added: ',i4,/,&
            t7,'# converged vectors: ',i4,/,&
            t5,'----------------------------------------')
    return
!                      
  end subroutine caslr_eff_std
  subroutine caslr_eff_driver(verbose,n,n2,n_targ,n_max,max_iter,tol,max_dav, &
                              apbmul,ambmul,spdmul,smdmul,lrprec,eig,evec,ok)
!
!   main driver for the efficient solution to the following generalized 
!   eigenvalue problem:
!
!   / A  B \ / Y \     /  S  D \ / Y \ 
!   |      | |   | = w |       | |   |
!   \ B  A / \ Z /     \ -D -S / \ Z /
!
!   where A, B, S are symmetric matrices and D is antysimmetric.
!
!   If (w, Y, Z) are a solution, then (-w, Z, Y) is also a solution. 
!   Following J. Chem. Phys., 118, 522 (2003), we enforce this property in 
!   the iterative procedure by expanding the eigenvector as
!
!   (Y, Z) = (b+, b+) + (b-, - b-)
!
!   This routine solves the associate problem 
!
!   /  S  D \ / Y \   1 / A  B \ / Y \ 
!   |       | |   | = - |      | |   |
!   \ -D -S / \ Z /   w \ B  A / \ Z /
!
!   using the Casida matrix, which is symmetric and positive definite, as 
!   the metric. This allows us to use expansion vectors that are orthogonal
!   with respect to the dot product defined by the metric, which in turn
!   results in a Rayleigh-Ritz procedure that requires the solution of a 
!   symmetric standard eigenvalue problem
!
!   / 0  s^T \ / u+ \   1 / u+ \
!   |        | |    | = - |    |
!   \ S   0  / \ u- /   w \ u- /
!
!   which can be reduced to a half-sized eigenvalue problem
!
!   s^T s u+ = (1/w)^2 u+
!   u- = 1/w Su+
!
!   input variables
!   ===============
!
!   verbose:  logical, whether to print various information at each 
!             iteration (eigenvalues, residuals...).
!
!   n:        integer, size of the A, B, S, D matrices
!
!   n2:       integer, n2 = 2*n, size of the generalized eigenvalue
!             problem
!
!   n_targ:   integer, number of required eigenpairs.
!
!   n_max:    integer, maximum size of the search space. should be 
!             >= n_targ. note that eig and evec should be allocated 
!             n_max and (n,n_max) rather than n_targ and (n,n_targ). 
!             for better convergence, a value larger than n_targ (eg.,
!             n_targ + 10) is recommended.
!   
!   max_iter: integer, maximum allowed number of iterations.
!
!   tol:      double precision real, the convergence threshold.
!
!   maxdav:   maximum size of the expansion subspace
!
!   apbmul:   external subroutine to apply (A+B) to a vector
!
!   ambmul:   external subroutine to apply (A-B) to a vector
!
!   spdmul:   external subroutine to apply (S+D) to a vector
!
!   smdmul:   external subroutine to apply (S-D) to a vector
!
!   lrprec:   external subroutine that applies a preconditioner.
!
!   output variables:
!   =================
!
!   eig:      double precision real array of size n_max. 
!             if ok is true, contains the positive eigenvalues
!             of the generalized problem in ascending order.
!
!   evec:     double precision real array of size (n2,n_max).
!             in input, a guess for the eigenvectors.
!             if ok is true, in output the computed eigenvectors.
!
!   ok:       logical, true if caslr_eff_driver converged.
!
    use utils
    implicit none
    logical, intent(in)                          :: verbose
    integer,                       intent(in)    :: n, n2, n_targ, n_max
    integer,                       intent(in)    :: max_iter, max_dav
    real(dp),                      intent(in)    :: tol
    real(dp), dimension(n_max),    intent(inout) :: eig
    real(dp), dimension(n2,n_max), intent(inout) :: evec
    logical,                       intent(inout) :: ok
    external                                     :: apbmul, ambmul, spdmul, smdmul, &
                                                    lrprec
!
!   local variables:
!   ================
!
    integer, parameter    :: min_dav = 10
    integer               :: istat
!
!   actual expansion space size and total dimension
!
    integer               :: dim_dav, lda, lda2
!
!   number of active vectors at a given iteration, and indices to access them
!
    integer               :: n_act, ind, i_beg
!
!   current size and total dimension of the expansion space
!
    integer               :: m_dim, ldu
!
!   number of frozen (i.e. converged) vectors
!
    integer               :: n_frozen
!
    integer               :: it, i_eig
!
    real(dp)              :: sqrtn, tol_rms, tol_max
!
!   arrays to control convergence and orthogonalization
!
    logical,  allocatable :: done(:)
!
!   expansion spaces, residuals and their norms.
!
    real(dp), allocatable :: vp(:,:), vm(:,:), lvp(:,:), lvm(:,:), bvp(:,:), bvm(:,:)
    real(dp), allocatable :: rp(:,:), rm(:,:), rr(:,:), r_norm(:,:)
!
!   eigenvectors of the reduced problem and components of the ritz vectors:
!
    real(dp), allocatable :: up(:,:), um(:,:), eigp(:,:), eigm(:,:), bp(:,:), bm(:,:)
!
!   subspace matrix and eigenvalues.
!
    real(dp), allocatable :: s_red(:,:), s_copy(:,:), e_red(:)
    real(dp), allocatable :: smat(:,:)
!
!   restarting variables
!
    integer               :: n_rst
    logical               :: restart
!
!   external functions:
!   ===================
!
    real(dp)              :: dnrm2
    external              :: dcopy, dnrm2, dgemm, dsyev
!
!   compute the actual size of the expansion space, checking that
!   the input makes sense.
!   no expansion space smaller than max_dav = 10 is deemed acceptable.
!
    dim_dav = max(min_dav,max_dav)
    lda     = dim_dav*n_max
    lda2    = 2 * lda
!
!   start by allocating memory for the various lapack routines
!
    lwork = get_mem_lapack(n,n_max)
    allocate (work(lwork), tau(n_max), stat=istat)
    call check_mem(istat)
!
!   allocate memory for the expansion space, the corresponding 
!   matrix-multiplied vectors and the residual:
!
    allocate (vp(n,lda), vm(n,lda), lvp(n,lda), lvm(n,lda), bvp(n,lda), bvm(n,lda), &
              rp(n,n_max), rm(n,n_max), rr(n,n_max), stat = istat)
    call check_mem(istat)
!
!   allocate memory for convergence check
!
    allocate (done(n_max), r_norm(2,n_max), stat=istat)
    call check_mem(istat)
!
!   allocate memory for the reduced matrix and its eigenvalues:
!
    allocate (s_red(lda,lda), s_copy(lda,lda), e_red(lda2), smat(lda,lda), stat=istat)
    call check_mem(istat)
!
!   allocate memory for the plus and minus eigenvector components:
!
    allocate (up(lda,n_max), um(lda,n_max), eigp(n,n_max), eigm(n,n_max), bp(n,n_max), bm(n,n_max), stat = istat)
    call check_mem(istat)
!
!   set the tolerances and compute a useful constant to compute rms norms:
!
    sqrtn   = sqrt(real(n,dp))
    tol_rms = tol
    tol_max = 10.0_dp * tol
!
!   clean out various quantities
!
    t_diag  = zero
    t_ortho = zero
    t_mv    = zero
    t_tot   = zero
    vp      = zero
    vm      = zero
    bvp     = zero
    bvm     = zero
    lvp     = zero
    lvm     = zero
    ok      = .false.
    done    = .false.
!
    call get_time(t_tot)
!
!   move the guess into the expansion space.
!
    do i_eig = 1, n_max
      vp(:,i_eig) = evec(1:n,i_eig) + evec(n+1:n2,i_eig)
      vm(:,i_eig) = evec(1:n,i_eig) - evec(n+1:n2,i_eig)
    end do
!
!   orthogonalize the expansion space to the metric.
!
    call get_time(t1)
    call apbmul(n,n_max,vp,lvp)
    call get_time(t2)
    t_mv = t_mv + t2 - t1
    call get_time(t1)
    call b_ortho(n,n_max,vp,lvp)
    call get_time(t2)
    t_ortho = t_ortho + t2 - t1
    call get_time(t1)
    call ambmul(n,n_max,vm,lvm)
    call get_time(t2)
    t_mv = t_mv + t2 - t1
    call get_time(t1)
    call b_ortho(n,n_max,vm,lvm)
    call get_time(t2)
    t_ortho = t_ortho + t2 - t1
!
    n_act = n_max
    ind   = 1
    i_beg = 1
!
!   initialize the counter for the expansion of the subspace
!
    m_dim = 1
    ldu   = 0
!
!   initialize to false the restart
!
    restart = .false.
!
!   main loop:
!
    1030 format(t5,'Davidson-Liu iterations (tol=',d10.2,'):',/, &
                t5,'------------------------------------------------------------------',/, &
                t7,'  iter  root              eigenvalue','         rms         max ok',/, &
                t5,'------------------------------------------------------------------')
    1040 format(t9,i4,2x,i4,f24.12,2d12.4,l3)
!
    if (verbose) write(6,1030) tol
!
    n_rst   = 0
    do it = 1, max_iter
!
!     update the size of the expansion space.
!
      ldu = ldu + n_act
!
!     perform this iteration's matrix-vector multiplications:
!
      call get_time(t1)
      call spdmul(n,n_act,vp(1,i_beg),bvm(1,i_beg))
      call smdmul(n,n_act,vm(1,i_beg),bvp(1,i_beg))
      call get_time(t2)
      t_mv = t_mv + t2 - t1
!
!     update the reduced matrix 
!
      call dgemm('t','n',ldu,ldu,n,one,vm,n,bvm,n,zero,smat,lda)
!
!     save s, and assemble s^t s:
!
      s_red  = smat
      s_copy = zero
      call dgemm('t','n',ldu,ldu,ldu,one,s_red,lda,s_red,lda,zero,s_copy,lda)
!
!     diagonalize s^t s
!
      call get_time(t1)
      call dsyev('v','u',ldu,s_copy,lda,e_red,work,lwork,info)
      call get_time(t2)
      t_diag = t_diag + t2 - t1
!
!     extract the eigenvalues and compute the ritz approximation to the
!     eigenvectors 
!
      do i_eig = 1, n_max
        eig(i_eig)      = sqrt(e_red(ldu - i_eig + 1))
        up(1:ldu,i_eig) = s_copy(1:ldu,ldu - i_eig + 1)
      end do
!
!     compute the u_- eigenvectors:
!
      call dgemm('n','n',ldu,n_max,ldu,one,s_red,lda,up,lda,zero,um,lda)
      do i_eig = 1, n_max
        um(1:ldu,i_eig) = um(1:ldu,i_eig)/eig(i_eig)
      end do
!
!     assemble the symmetric and antysimmetric combinations (Y+Z) and (Y-Z)
!
      call dgemm('n','n',n,n_max,ldu,one,vp,n,up,lda,zero,eigp,n)
      call dgemm('n','n',n,n_max,ldu,one,vm,n,um,lda,zero,eigm,n)
!
!     assemble the current approximation to the eigenvectors
!
      do i_eig = 1, n_max
        evec(1:n,i_eig)    = eigp(:,i_eig) + eigm(:,i_eig)
        evec(n+1:n2,i_eig) = eigp(:,i_eig) - eigm(:,i_eig)
      end do
!
!     compute the residuals, and their rms and sup norms:
!
      call dgemm('n','n',n,n_max,ldu,one,bvp,n,um,lda,zero,rp,n)
      call dgemm('n','n',n,n_max,ldu,one,bvm,n,up,lda,zero,rm,n)
      call dgemm('n','n',n,n_max,ldu,one,lvp,n,up,lda,zero,bp,n)
      call dgemm('n','n',n,n_max,ldu,one,lvm,n,um,lda,zero,bm,n)
      
      do i_eig = 1, n_targ
!
!       if the eigenvalue is already converged, skip it.
!
        if (done(i_eig)) cycle
!
        call daxpy(n,-eig(i_eig),bp(:,i_eig),1,rp(:,i_eig),1)
        call daxpy(n,-eig(i_eig),bm(:,i_eig),1,rm(:,i_eig),1)
        r_norm(1,i_eig) = (dnrm2(n,rp(:,i_eig),1) + dnrm2(n,rm(:,i_eig),1))/(eig(i_eig)*sqrt(two)*sqrtn)
        r_norm(2,i_eig) = (maxval(abs(rp(:,i_eig))) + maxval(abs(rm(:,i_eig))))/(sqrt(two)*eig(i_eig))
      end do
!
!     check convergence. lock the first contiguous converged eigenvalues
!     by setting the logical array "done" to true.
!
      do i_eig = 1, n_targ
        if (done(i_eig)) cycle
        done(i_eig)     = r_norm(1,i_eig).lt.tol_rms .and. &
                          r_norm(2,i_eig).lt.tol_max .and. &
                          it.gt.1
        if (.not.done(i_eig)) then
          done(i_eig+1:n_max) = .false.
          exit
        end if
      end do
!
!     print some information:
!
      if (verbose) then
        do i_eig = 1, n_targ
          write(6,1040) it, i_eig, one/eig(i_eig), r_norm(:,i_eig), done(i_eig)
        end do
        write(6,*) 
      end if
!
      if (all(done(1:n_targ))) then
        ok = .true.
        do i_eig = 1, n_targ
          eig(i_eig) = one/eig(i_eig)
        end do
        exit
      end if
!
!     check whether an update is required. 
!     if not, perform a davidson restart.
!
      if (m_dim .lt. dim_dav) then
!
!       compute the preconditioned residuals using davidson's procedure
!       note that this is done with a user-supplied subroutine, that can
!       be generalized to experiment with fancy preconditioners that may
!       be more effective than the diagonal one, as in the original 
!       algorithm.
!
        m_dim = m_dim + 1
        i_beg = i_beg + n_act
        n_act = n_max
        n_frozen = 0
        do i_eig = 1, n_targ
          if (done(i_eig)) then
            n_act = n_act - 1
            n_frozen = n_frozen + 1
          else
            exit
          end if
        end do
        ind   = n_max - n_act + 1
        call lrprec(n,n_act,eig(ind),rp(1,ind),rm(1,ind),vp(1,i_beg),vm(1,i_beg))
!
!       orthogonalize the new vectors to the existing ones and then
!       orthonormalize them.
!
        call get_time(t1)
        call b_ortho_vs_x(n,ldu,n_act,vp,lvp,vp(1,i_beg))
        call get_time(t2)
        t_ortho = t_ortho + t2 - t1
        call get_time(t1)
        call apbmul(n,n_act,vp(1,i_beg),lvp(1,i_beg))
        call get_time(t2)
        t_mv = t_mv + t2 - t1
        call get_time(t1)
        call b_ortho(n,n_act,vp(1,i_beg),lvp(1,i_beg))
        call b_ortho_vs_x(n,ldu,n_act,vm,lvm,vm(1,i_beg))
        call get_time(t2)
        t_ortho = t_ortho + t2 - t1
        call get_time(t1)
        call ambmul(n,n_act,vm(1,i_beg),lvm(1,i_beg))
        call get_time(t2)
        t_mv = t_mv + t2 - t1
        call get_time(t1)
        call b_ortho(n,n_act,vm(1,i_beg),lvm(1,i_beg))
        call get_time(t2)
        t_ortho = t_ortho + t2 - t1
      else
!
        if (verbose) write(6,'(t7,a)') 'Restarting davidson.'
        restart = .true.
!
!       initialize indexes back to their starting values 
!
        ldu   = 0
        i_beg = 1
        m_dim = 1
        n_rst = 0
!        
        n_act = n_max 
        vp = zero
        vm = zero
!
!       put current eigenvectors into the first position of the 
!       expansion space
!
        do i_eig = 1, n_max
          vp(:,i_eig) = evec(1:n,i_eig) + evec(n+1:n2,i_eig)
          vm(:,i_eig) = evec(1:n,i_eig) - evec(n+1:n2,i_eig)
        end do
!
        lvp   = zero
        lvm   = zero
!
        call get_time(t1)
        call apbmul(n,n_max,vp,lvp)
        call get_time(t2)
        t_mv = t_mv + t2 - t1
        call get_time(t1)
        call b_ortho(n,n_max,vp,lvp)
        call get_time(t2)
        t_ortho = t_ortho + t2 - t1
        call get_time(t1)
        call ambmul(n,n_max,vm,lvm)
        call get_time(t2)
        t_mv = t_mv + t2 - t1
        call get_time(t1)
        call b_ortho(n,n_max,vm,lvm)
        call get_time(t2)
        t_ortho = t_ortho + t2 - t1
!
        bvp   = zero
        bvm   = zero
        s_red = zero
        smat  = zero
!        
      end if
!      
      if (verbose) write(6,1050) n_targ, n_act, n_frozen
!      
    end do
!
    call get_time(t1)
    t_tot = t1 - t_tot
    1000 format(t3,'timings for caslr_eff (cpu/wall):   ',/, &
                t3,'  matrix-vector multiplications: ',2f12.4,/, &
                t3,'  diagonalization:               ',2f12.4,/, &
                t3,'  orthogonalization:             ',2f12.4,/, &
                t3,'                                 ',24('='),/,  &
                t3,'  total:                         ',2f12.4)
    if (verbose) write(6,1000) t_mv, t_diag, t_ortho, t_tot
    deallocate(work,tau,vp,vm,lvp,lvm,bvp,bvm,rp,rm,rr,done,r_norm,s_red,s_copy,e_red,smat,up,um,eigp,eigm,bp,bm)
!
1050 format(t5,'----------------------------------------',/,&
            t7,'# target vectors:    ',i4,/,&
            t7,'# new vectors added: ',i4,/,&
            t7,'# converged vectors: ',i4,/,&
            t5,'----------------------------------------')
    return
  end subroutine caslr_eff_driver
!
  subroutine davidson_driver(verbose,n,n_targ,n_max,max_iter,tol,max_dav, &
                             shift,matvec,precnd,eig,evec,ok)
!
!   main driver for davidson-liu.
!
!   input variables:
!   ================
!
!   verbose:  logical, whether to print various information at each 
!             iteration (eigenvalues, residuals...).
!
!   n:        integer, size of the matrix to be diagonalized.
!
!   n_targ:   integer, number of required eigenpairs.
!
!   n_max:    integer, maximum size of the search space. should be 
!             >= n_targ. note that eig and evec should be allocated 
!             n_max and (n,n_max) rather than n_targ and (n,n_targ). 
!             for better convergence, a value larger than n_targ (eg.,
!             n_targ + 10) is recommended.
!   
!   max_iter: integer, maximum allowed number of iterations.
!
!   tol:      double precision real, the convergence threshold.
!
!   max_dav:  integer. maximum allowed number of iterations before a 
!             restart is forced. 
!             when n_max eigenvalues are searched, this implies a maximum
!             dimension of the expansion subspace of n_max * max_dav.
!
!   shift:    double precision real, a diagonal level shifting parameter
!
!   matvec:   external subroutine that performs the matrix-vector
!             multiplication
!
!   precnd:   external subroutine that applied a preconditioner.
!
!   output variables:
!   =================
!
!   eig:      double precision array of size n_max. if ok is true, 
!             the computed eigenvalues in asceding order.
!
!   evec:     double precision array of size (n,n_max). 
!             in input, a guess for the eigenvectors.
!             if ok is true, in output the computed eigenvectors.
!
!   ok:       logical, true if davidson converged.
!
    logical,                      intent(in)    :: verbose
    integer,                      intent(in)    :: n, n_targ, n_max
    integer,                      intent(in)    :: max_iter, max_dav
    real(dp),                     intent(in)    :: tol, shift
    real(dp), dimension(n_max),   intent(inout) :: eig
    real(dp), dimension(n,n_max), intent(inout) :: evec
    logical,                      intent(inout) :: ok
    external                                    :: matvec, precnd
!
!   local variables:
!   ================
!
    integer, parameter    :: min_dav = 10
    integer               :: istat
!
!   actual expansion space size and total dimension
!
    integer               :: dim_dav, lda
!
!   number of active vectors at a given iteration, and indices to access them
!
    integer               :: n_act, ind, i_beg
!
!   current size and total dimension of the expansion space
!
    integer               :: m_dim, ldu
!
!   number of frozen (i.e. converged) vectors
!
    integer               :: n_frozen
!
    integer               :: it, i_eig
!
    real(dp)              :: sqrtn, tol_rms, tol_max
    real(dp)              :: xx(1)
!
!   arrays to control convergence and orthogonalization
!
    logical,  allocatable :: done(:)
!
!   expansion spaces, residuals and their norms.
!
    real(dp), allocatable :: space(:,:), aspace(:,:), r(:,:), r_norm(:,:)
!
!   subspace matrix and eigenvalues.
!
    real(dp), allocatable :: a_red(:,:), a_copy(:,:), e_red(:)
!
!   restarting variables
!
    integer               :: n_rst
    logical               :: restart
!
!   external functions:
!   ===================
!
    real(dp)              :: dnrm2
    external              :: dcopy, dnrm2, dgemm, dsyev
!
!   compute the actual size of the expansion space, checking that
!   the input makes sense.
!   no expansion space smaller than max_dav = 10 is deemed acceptable.
!
    dim_dav = max(min_dav,max_dav)
    lda     = dim_dav*n_max
!
!   start by allocating memory for the various lapack routines
!
    lwork = get_mem_lapack(n,n_max)
    allocate (work(lwork), tau(n_max), stat=istat)
    call check_mem(istat)
!
!   allocate memory for the expansion space, the corresponding 
!   matrix-multiplied vectors and the residual:
!
    allocate (space(n,lda), aspace(n,lda), r(n,n_max), stat = istat)
    call check_mem(istat)
!
!   allocate memory for convergence check
!
    allocate (done(n_max), r_norm(2,n_max), stat=istat)
    call check_mem(istat)
!
!   allocate memory for the reduced matrix and its eigenvalues:
!
    allocate (a_red(lda,lda), a_copy(lda,lda), e_red(lda), stat=istat)
    call check_mem(istat)
!
!   set the tolerances and compute a useful constant to compute rms norms:
!
    sqrtn   = sqrt(real(n,dp))
    tol_rms = tol
    tol_max = 10.0_dp * tol
!
!   clean out various quantities
!
    t_diag   = zero
    t_ortho  = zero
    t_mv     = zero
    t_tot    = zero
    space   = zero
    aspace  = zero
    a_red   = zero
    ok      = .false.
    done    = .false.
!
    call get_time(t_tot)
!
!   check whether we have a guess for the eigenvectors in evec, and
!   whether it is orthonormal.
!   if evec is zero, create a random guess.
!
    call check_guess(n,n_max,evec)
!
!   move the guess into the expansion space.
!
    call dcopy(n*n_max,evec,1,space,1)
!
!   initialize the number of active vectors and the associated indices.
!
    n_act = n_max
    ind   = 1
    i_beg = 1
!
!   initialize the counter for the expansion of the subspace
!
    m_dim = 1
    ldu   = 0
!
!   initialize to false the restart
!
    restart = .false.
!
!   main loop:
!
    1030 format(t5,'Davidson-Liu iterations (tol=',d10.2,'):',/, &
                t5,'------------------------------------------------------------------',/, &
                t7,'  iter  root              eigenvalue','         rms         max ok',/, &
                t5,'------------------------------------------------------------------')
    1040 format(t9,i4,2x,i4,f24.12,2d12.4,l3)
!
    if (verbose) write(6,1030) tol
!
    n_rst   = 0
    do it = 1, max_iter
!
!     update the size of the expansion space.
!
      ldu = ldu + n_act
!
!     perform this iteration's matrix-vector multiplication:
!
      call get_time(t1)
      call matvec(n,n_act,space(1,i_beg+n_rst),aspace(1,i_beg+n_rst))
      call get_time(t2)
      t_mv = t_mv + t2 - t1
!
!     update the reduced matrix 
!
      call dgemm('t','n',ldu,n_act,n,one,space,n,aspace(1,i_beg+n_rst),n,zero,a_red(1,i_beg+n_rst),lda)
!
!     explicitly putting the first block of 
!     converged eigenvalues in the reduced matrix
!
      if (restart) then
        do i_eig = 1, n_rst
          a_red(i_eig,i_eig) = e_red(i_eig)
        end do
        restart = .false.
        n_rst   = 0
      end if
      a_copy = a_red
!
!     diagonalize the reduced matrix
!
      call get_time(t1)
      call dsyev('v','u',ldu,a_copy,lda,e_red,work,lwork,info)
      call get_time(t2)
      t_diag = t_diag + t2 - t1
!
!     extract the eigenvalues and compute the ritz approximation to the
!     eigenvectors 
!
      eig = e_red(1:n_max)
!
      call dgemm('n','n',n,n_max,ldu,one,space,n,a_copy,lda,zero,evec,n)
!
!     compute the residuals, and their rms and sup norms:
!
      call dgemm('n','n',n,n_max,ldu,one,aspace,n,a_copy,lda,zero,r,n)
!
      do i_eig = 1, n_targ
!
!       if the eigenvalue is already converged, skip it.
!
        if (done(i_eig)) cycle
!
        call daxpy(n,-eig(i_eig),evec(:,i_eig),1,r(:,i_eig),1)
        r_norm(1,i_eig) = dnrm2(n,r(:,i_eig),1)/sqrtn
        r_norm(2,i_eig) = maxval(abs(r(:,i_eig)))
      end do
!
!     check convergence. lock the first contiguous converged eigenvalues
!     by setting the logical array "done" to true.
!
      do i_eig = 1, n_targ
        if (done(i_eig)) cycle
        done(i_eig)     = r_norm(1,i_eig).lt.tol_rms .and. &
                          r_norm(2,i_eig).lt.tol_max .and. &
                          it.gt.1
        if (.not.done(i_eig)) then
          done(i_eig+1:n_max) = .false.
          exit
        end if
      end do
!
!     print some information:
!
      if (verbose) then
        do i_eig = 1, n_targ
          write(6,1040) it, i_eig, eig(i_eig) - shift, r_norm(:,i_eig), done(i_eig)
        end do
        write(6,*) 
      end if
!
      if (all(done(1:n_targ))) then
        ok = .true.
        exit
      end if
!
!     check whether an update is required. 
!     if not, perform a davidson restart.
!
      if (m_dim .lt. dim_dav) then
!
!       compute the preconditioned residuals using davidson's procedure
!       note that this is done with a user-supplied subroutine, that can
!       be generalized to experiment with fancy preconditioners that may
!       be more effective than the diagonal one, as in the original 
!       algorithm.
!
        m_dim = m_dim + 1
        i_beg = i_beg + n_act
        n_act = n_max
        n_frozen = 0
        do i_eig = 1, n_targ
          if (done(i_eig)) then
            n_act = n_act - 1
            n_frozen = n_frozen + 1
          else
            exit
          end if
        end do
        ind   = n_max - n_act + 1
        call precnd(n,n_act,-eig(ind),r(1,ind),space(1,i_beg))
!
!       orthogonalize the new vectors to the existing ones and then
!       orthonormalize them.
!
        call get_time(t1)
!        call ortho_vs_x(.false.,n,ldu,n_act,space,space(1,i_beg),xx,xx)
        call ortho_vs_x(n,ldu,n_act,space,space(1,i_beg),xx,xx)
        call get_time(t2)
        t_ortho = t_ortho + t2 - t1
      else
        if (verbose) write(6,'(t7,a)') 'Restarting davidson.'
        n_act = n_max
        space = zero
!
!       put current eigenvectors into the first position of the 
!       expansion space
!
        call dcopy(n_max*n,evec,1,space,1)
        aspace = zero
        a_red  = zero
!
!       initialize indexes back to their starting values 
!
        ldu   = 0
        i_beg = 1 
        m_dim = 1
        n_rst = 0
!
!       counting how many matvec we can skip at the next
!       iteration
!
        do i_eig = 1, n_targ
          if (done(i_eig)) then
            n_rst = n_rst + 1
          else
            exit
          end if
        end do
        restart = .true.
      end if
      if (verbose) write(6,1050) n_targ, n_act, n_frozen
!
    end do
!
    call get_time(t2)
    t_tot = t2 - t_tot
!
!   if required, print timings
!
    1000 format(t3,'timings for davidson (cpu/wall): ',/, &
                t3,'  matrix-vector multiplications: ',2f12.4,/, &
                t3,'  diagonalization:               ',2f12.4,/, &
                t3,'  orthogonalization:             ',2f12.4,/, &
                t3,'                                 ',24('='),/,  &
                t3,'  total:                         ',2f12.4)
    if (verbose) write(6,1000) t_mv, t_diag, t_ortho, t_tot
!      
!   deallocate memory
!
    deallocate (work, tau, space, aspace, r, done, r_norm, a_red, a_copy, e_red)
!
1050 format(t5,'----------------------------------------',/,&
            t7,'# target vectors:    ',i4,/,&
            t7,'# new vectors added: ',i4,/,&
            t7,'# converged vectors: ',i4,/,&
            t5,'----------------------------------------')
    return
  end subroutine davidson_driver
!
  subroutine gen_david_driver(verbose,n,n_targ,n_max,max_iter,tol,max_dav, &
                              shift,matvec,precnd,bvec,eig,evec,ok)
!
!   main driver for davidson-liu.
!
!   input variables:
!   ================
!
!   verbose:  logical, whether to print various information at each 
!             iteration (eigenvalues, residuals...).
!
!   n:        integer, size of the matrix to be diagonalized.
!
!   n_targ:   integer, number of required eigenpairs.
!
!   n_max:    integer, maximum size of the search space. should be 
!             >= n_targ. note that eig and evec should be allocated 
!             n_max and (n,n_max) rather than n_targ and (n,n_targ). 
!             for better convergence, a value larger than n_targ (eg.,
!             n_targ + 10) is recommended.
!   
!   max_iter: integer, maximum allowed number of iterations.
!
!   tol:      double precision real, the convergence threshold.
!
!   max_dav:  integer. maximum allowed number of iterations before a 
!             restart is forced. 
!             when n_max eigenvalues are searched, this implies a maximum
!             dimension of the expansion subspace of n_max * max_dav.
!
!   shift:    double precision real, a diagonal level shifting parameter
!
!   matvec:   external subroutine that performs the matrix-vector
!             multiplication
!
!   precnd:   external subroutine that applies a preconditioner.
!
!   bvec:     external subroutine that applies the metric.
!
!   output variables:
!   =================
!
!   eig:      double precision array of size n_max. if ok is true, 
!             the computed eigenvalues in asceding order.
!
!   evec:     double precision array of size (n,n_max). 
!             in input, a guess for the eigenvectors.
!             if ok is true, in output the computed eigenvectors.
!
!   ok:       logical, true if davidson converged.
!
    logical,                      intent(in)    :: verbose
    integer,                      intent(in)    :: n, n_targ, n_max
    integer,                      intent(in)    :: max_iter, max_dav
    real(dp),                     intent(in)    :: tol, shift
    real(dp), dimension(n_max),   intent(inout) :: eig
    real(dp), dimension(n,n_max), intent(inout) :: evec
    logical,                      intent(inout) :: ok
    external                                    :: matvec, precnd, bvec
!
!   local variables:
!   ================
!
    integer, parameter    :: min_dav = 10
    integer               :: istat
!
!   actual expansion space size and total dimension
!
    integer               :: dim_dav, lda
!
!   number of active vectors at a given iteration, and indices to access them
!
    integer               :: n_act, ind, i_beg
!
!   current size and total dimension of the expansion space
!
    integer               :: m_dim, ldu
!
!   number of frozen (i.e. converged) vectors
!
    integer               :: n_frozen
!
    integer               :: it, i_eig
!
    real(dp)              :: sqrtn, tol_rms, tol_max
!
!   arrays to control convergence and orthogonalization
!
    logical,  allocatable :: done(:)
!
!   expansion spaces, residuals and their norms.
!
    real(dp), allocatable :: space(:,:), aspace(:,:), bspace(:,:), r(:,:), r_norm(:,:)
!
!   subspace matrix and eigenvalues.
!
    real(dp), allocatable :: a_red(:,:), a_copy(:,:), s_red(:,:), s_copy(:,:), e_red(:)
!
!   scratch:
!
    real(dp), allocatable :: b_evec(:,:)
!
!   restarting variables
!
    integer               :: n_rst
    logical               :: restart
!
!   external functions:
!   ===================
!
    real(dp)              :: dnrm2
    external              :: dcopy, dnrm2, dgemm, dsyev
!
!   compute the actual size of the expansion space, checking that
!   the input makes sense.
!   no expansion space smaller than max_dav = 10 is deemed acceptable.
!
    dim_dav = max(min_dav,max_dav)
    lda     = dim_dav*n_max
!
!   start by allocating memory for the various lapack routines
!
    lwork = get_mem_lapack(n,n_max)
    allocate (work(lwork), tau(n_max), stat=istat)
    call check_mem(istat)
!
!   allocate memory for the expansion space, the corresponding 
!   matrix-multiplied vectors and the residual:
!
    allocate (space(n,lda), aspace(n,lda), bspace(n,lda), r(n,n_max), &
              b_evec(n,n_max), stat = istat)
    call check_mem(istat)
!
!   allocate memory for convergence check
!
    allocate (done(n_max), r_norm(2,n_max), stat=istat)
    call check_mem(istat)
!
!   allocate memory for the reduced matrix and its eigenvalues:
!
    allocate (a_red(lda,lda), a_copy(lda,lda), e_red(lda), s_red(lda,lda), &
              s_copy(lda,lda), stat=istat)
    call check_mem(istat)
!
!   set the tolerances and compute a useful constant to compute rms norms:
!
    sqrtn   = sqrt(real(n,dp))
    tol_rms = tol
    tol_max = 10.0_dp * tol
!
!   clean out various quantities
!
    t_diag   = zero
    t_ortho  = zero
    t_mv     = zero
    t_tot    = zero
    space   = zero
    aspace  = zero
    bspace  = zero
    a_red   = zero
    s_red   = zero
    ok      = .false.
    done    = .false.
!
    call get_time(t_tot)
!
!   check whether we have a guess for the eigenvectors in evec, and
!   whether it is orthonormal.
!   if evec is zero, create a random guess.
!
    call check_guess(n,n_max,evec)
!
!   move the guess into the expansion space.
!
    call dcopy(n*n_max,evec,1,space,1)
!
!   apply the b matrix and b-orthogonalize the guess:
!
    call bvec(n,n_max,space,bspace)
    call b_ortho(n,n_max,space,bspace)
!
!   initialize the number of active vectors and the associated indices.
!
    n_act = n_max
    ind   = 1
    i_beg = 1
!
!   initialize the counter for the expansion of the subspace
!
    m_dim = 1
    ldu   = 0
!
!   initialize to false the restart
!
    restart = .false.
!
!   main loop:
!
    1030 format(t5,'Generalized Davidson-Liu iterations (tol=',d10.2,'):',/, &
                t5,'------------------------------------------------------------------',/, &
                t7,'  iter  root              eigenvalue','         rms         max ok',/, &
                t5,'------------------------------------------------------------------')
    1040 format(t9,i4,2x,i4,f24.12,2d12.4,l3)
!
    if (verbose) write(6,1030) tol
!
    n_rst   = 0
    do it = 1, max_iter
!
!     update the size of the expansion space.
!
      ldu = ldu + n_act
!
!     perform this iteration's matrix-vector multiplication:
!
      call get_time(t1)
      call matvec(n,n_act,space(1,i_beg+n_rst),aspace(1,i_beg+n_rst))
!     call bvec(n,n_act,space(1,i_beg+n_rst),bspace(1,i_beg+n_rst))
      call get_time(t2)
      t_mv = t_mv + t2 - t1
!
!     update the reduced matrices 
!
      call dgemm('t','n',ldu,n_act,n,one,space,n,aspace(1,i_beg+n_rst),n,zero,a_red(1,i_beg+n_rst),lda)
!     call dgemm('t','n',ldu,n_act,n,one,space,n,bspace(1,i_beg+n_rst),n,zero,s_red(1,i_beg+n_rst),lda)
!
!     explicitly putting the first block of 
!     converged eigenvalues in the reduced matrix
!
      if (restart) then
        do i_eig = 1, n_rst
          a_red(i_eig,i_eig) = e_red(i_eig)
        end do
        restart = .false.
        n_rst   = 0
      end if
      a_copy = a_red
!     s_copy = s_red
!
!     diagonalize the reduced matrix
!
      call get_time(t1)
!     call dsygv(1,'v','u',ldu,a_copy,lda,s_copy,lda,e_red,work,lwork,info)
      call dsyev('v','u',ldu,a_copy,lda,e_red,work,lwork,info)
      call get_time(t2)
      t_diag = t_diag + t2 - t1
!
!     extract the eigenvalues and compute the ritz approximation to the
!     eigenvectors 
!
      eig = e_red(1:n_max)
!
      call dgemm('n','n',n,n_max,ldu,one,space,n,a_copy,lda,zero,evec,n)
!
!     compute the residuals, and their rms and sup norms:
!
      call dgemm('n','n',n,n_max,ldu,one,aspace,n,a_copy,lda,zero,r,n)
      call dgemm('n','n',n,n_max,ldu,one,bspace,n,a_copy,lda,zero,b_evec,n)
!
      do i_eig = 1, n_targ
!
!       if the eigenvalue is already converged, skip it.
!
        if (done(i_eig)) cycle
!
        call daxpy(n,-eig(i_eig),b_evec(:,i_eig),1,r(:,i_eig),1)
        r_norm(1,i_eig) = dnrm2(n,r(:,i_eig),1)/sqrtn
        r_norm(2,i_eig) = maxval(abs(r(:,i_eig)))
      end do
!
!     check convergence. lock the first contiguous converged eigenvalues
!     by setting the logical array "done" to true.
!
      do i_eig = 1, n_targ
        if (done(i_eig)) cycle
        done(i_eig)     = r_norm(1,i_eig).lt.tol_rms .and. &
                          r_norm(2,i_eig).lt.tol_max .and. &
                          it.gt.1
        if (.not.done(i_eig)) then
          done(i_eig+1:n_max) = .false.
          exit
        end if
      end do
!
!     print some information:
!
      if (verbose) then
        do i_eig = 1, n_targ
          write(6,1040) it, i_eig, eig(i_eig) - shift, r_norm(:,i_eig), done(i_eig)
        end do
        write(6,*) 
      end if
!
      if (all(done(1:n_targ))) then
        ok = .true.
        exit
      end if
!
!     check whether an update is required. 
!     if not, perform a davidson restart.
!
      if (m_dim .lt. dim_dav) then
!
!       compute the preconditioned residuals using davidson's procedure
!       note that this is done with a user-supplied subroutine, that can
!       be generalized to experiment with fancy preconditioners that may
!       be more effective than the diagonal one, as in the original 
!       algorithm.
!
        m_dim = m_dim + 1
        i_beg = i_beg + n_act
        n_act = n_max
        n_frozen = 0
        do i_eig = 1, n_targ
          if (done(i_eig)) then
            n_act = n_act - 1
            n_frozen = n_frozen + 1
          else
            exit
          end if
        end do
        ind   = n_max - n_act + 1
        call precnd(n,n_act,-eig(ind),r(1,ind),space(1,i_beg))
!
!       orthogonalize the new vectors to the existing ones and then
!       orthonormalize them.
!
        call get_time(t1)
        call b_ortho_vs_x(n,ldu,n_act,space,bspace,space(1,i_beg))
        call bvec(n,n_act,space(1,i_beg),bspace(1,i_beg))
        call b_ortho(n,n_act,space(1,i_beg),bspace(1,i_beg))
        call get_time(t2)
        t_ortho = t_ortho + t2 - t1
      else
        if (verbose) write(6,'(t7,a)') 'Restarting davidson.'
        n_act = n_max
        space = zero
!
!       put current eigenvectors into the first position of the 
!       expansion space
!
        call dcopy(n_max*n,evec,1,space,1)
        call dcopy(n_max*n,b_evec,1,bspace,1)
        call b_ortho(n,n_max,space,bspace)
        aspace = zero
        bspace = zero
        a_red  = zero
        s_red  = zero
!
!       initialize indexes back to their starting values 
!
        ldu   = 0
        i_beg = 1 
        m_dim = 1
        n_rst = 0
!
!       counting how many matvec we can skip at the next
!       iteration
!
        do i_eig = 1, n_targ
          if (done(i_eig)) then
            n_rst = n_rst + 1
          else
            exit
          end if
        end do
        restart = .true.
      end if
      if (verbose) write(6,1050) n_targ, n_act, n_frozen
!
    end do
!
    call get_time(t2)
    t_tot = t2 - t_tot
!
!   if required, print timings
!
    1000 format(t3,'timings for davidson (cpu/wall): ',/, &
                t3,'  matrix-vector multiplications: ',2f12.4,/, &
                t3,'  diagonalization:               ',2f12.4,/, &
                t3,'  orthogonalization:             ',2f12.4,/, &
                t3,'                                 ',24('='),/,  &
                t3,'  total:                         ',2f12.4)
    if (verbose) write(6,1000) t_mv, t_diag, t_ortho, t_tot
!      
!   deallocate memory
!
    deallocate (work, tau, space, aspace, r, done, r_norm, a_red, a_copy, e_red)
!
1050 format(t5,'----------------------------------------',/,&
            t7,'# target vectors:    ',i4,/,&
            t7,'# new vectors added: ',i4,/,&
            t7,'# converged vectors: ',i4,/,&
            t5,'----------------------------------------')
    return
  end subroutine gen_david_driver
!
! orthogonalization routines:
! ===========================
!
  subroutine ortho(n,m,u,w)
    implicit none
!
!   orthogonalize m vectors of lenght n using the QR decomposition.
!   this is done by computing U = QR and then by solving the upper 
!   triangular system U(ortho)R = U.
!
!   using this strategy allows to apply the same linear transformation
!   that orthogonalizes U to a second set of vectors that usually contain
!   the product AU, where A is some matrix that has already been applied
!   to U. 
!   this is useful when U and AU are built together without explicitly 
!   performing the matrix vector multiplication.
!
!   arguments:
!   ==========
!
    integer,                     intent(in)    :: n, m
    real(dp),  dimension(n,m),   intent(inout) :: u, w
!
!   local scratch
!   =============
!
    real(dp), allocatable :: v(:,:), work_o(:)
!
!   external functions:
!   ===================
!
    external dgeqrf, dtrsm
!
    allocate (v(n,m))
    v = u
    allocate(work_o(1))
    call dgeqrf(n,m,u,n,tau,work_o,-1,info)
    lwork = int(work_o(1))
    deallocate(work_o)
    allocate(work_o(lwork))
    call dgeqrf(n,m,u,n,tau,work_o,lwork,info)
!
    call dtrsm('r','u','n','n',n,m,one,u,n,v,n)
!
!
    u = v
!
    deallocate (v)
    return
  end subroutine ortho
!
  subroutine b_ortho(n,m,u,bu)
    implicit none
!
!   b-orthogonalize m vectors of lenght n using the cholesky decomposition
!   of the overlap matrix.
!   this is in principle not a good idea, as the u'bu matrix can be very
!   ill-conditioned, independent of how bas is b, and only works if x is
!   already orthonormal. 
!
!   arguments:
!   ==========
!
    integer,                     intent(in)    :: n, m
    real(dp),  dimension(n,m),   intent(inout) :: u, bu
!
!   local variables
!   ===============
!
    integer               :: lwork, info, i, j
    real(dp), allocatable :: metric(:,:), sigma(:), u_svd(:,:), vt_svd(:,:), &
                             temp(:,:)
    real(dp), parameter   :: tol_svd = 1.0e-5_dp
    logical,  parameter   :: use_svd = .false.
!
    real(dp), allocatable :: work(:)
!   external functions:
!   ===================
!
    external dpotrf, dtrsm, dgemm
!
    allocate (metric(m,m))
!
    call dgemm('t','n',m,m,n,one,u,n,bu,n,zero,metric,m)
!
    if (use_svd) then
!
!     debug option: use svd to b-orthonormalize, by computing
!     b**(-1/2)
!
      allocate(work(1))
      allocate (sigma(m), u_svd(m,m), vt_svd(m,m), temp(n,m))
      call dgesvd('a','a',m,m,metric,m,sigma,u_svd,m,vt_svd,m,work,-1,info)
      lwork = int(work(1)) 
      deallocate(work)
      allocate(work(lwork))
      call dgesvd('a','a',m,m,metric,m,sigma,u_svd,m,vt_svd,m,work,lwork,info)
!
!     compute sigma**(-1/2)
!
      do i = 1, m
        if (sigma(i) .gt. tol_svd) then
          sigma(i) = 1/sqrt(sigma(i))
        else
          sigma(i) = zero
        end if
      end do
!
!     compute metric ** (-1/2). first, compute sigma ** (-1/2) vt
!
      metric = zero
      do i = 1, m
        do j = 1, m
          metric(j,i) = metric(j,i) + sigma(j)*vt_svd(j,i)
        end do
      end do
!
!     now, multiply for u:
!
      vt_svd = metric
      call dgemm('n','n',m,m,m,one,u_svd,m,vt_svd,m,zero,metric,m)
!
!     metric contains s ** (-1/2), and projects out directions corresponding
!     to pathological singular values. 
!     orthogonalize u and bu:
!
      call dgemm('n','n',n,m,m,one,u,n,metric,m,zero,temp,n)
      u = temp
      call dgemm('n','n',n,m,m,one,bu,n,metric,m,zero,temp,n)
      bu = temp
!
      deallocate (sigma, u_svd, vt_svd, temp, work)
    else
!
!     compute the cholesky factorization of the metric.
!
      call dpotrf('l',m,metric,m,info)
!
!     get u * l^-T and bu * l^-T
!
      call dtrsm('r','l','t','n',n,m,one,metric,m,u,n)
      call dtrsm('r','l','t','n',n,m,one,metric,m,bu,n)
    end if
!
    deallocate (metric)
    return
  end subroutine b_ortho
!
  subroutine ortho_cd(n,m,u,growth,ok)
    implicit none
!
!   orthogonalize m vectors of lenght n using the Cholesky factorization
!   of their overlap. 
!   this is done by metric = U^t U and then by computing its cholesky 
!   decompositoin metric = L L^t. The orthogonal vectors are obtained then
!   by solving the triangular linear system
!
!     U(ortho)L^T = U
!
!   using this strategy allows to apply the same linear transformation
!   that orthogonalizes U to a second set of vectors that usually contain
!   the product AU, where A is some matrix that has already been applied
!   to U.  this is useful when U and AU are built together without explicitly 
!   performing the matrix vector multiplication.
!
!   as cholesky decomposition is not the most stable way of orthogonalizing
!   a set of vectors, the orthogonalization is refined iteratively. 
!   still, this routine can fail. 
!   a logical flag is then set to false, so that the calling program can 
!   call a more robust orthogonalization routine without aborting.
!
!   arguments:
!   ==========
!
!    logical,                   intent(in)    :: do_other
    integer,                   intent(in)    :: n, m
!    real(dp),  dimension(n,m), intent(inout) :: u, w
    real(dp),  dimension(n,m), intent(inout) :: u
    real(dp),                  intent(inout) :: growth
    logical,                   intent(inout) :: ok
!
!   local variables
!   ===============
!
    integer               :: it, it_micro
!    real(dp)              :: metric_norm, dnrm2, alpha, unorm, shift
!    logical               :: macro_done, micro_done
    real(dp)              :: error, dnrm2, alpha, unorm, shift
    real(dp)              :: rcond, l_norm, linv_norm
    logical               :: macro_done, micro_done
    real(dp), parameter   :: tol_ortho_cd = two * epsilon(one)
    integer, parameter    :: maxit = 10
!
!   local scratch
!   =============
!
    real(dp), allocatable :: metric(:,:), msave(:,:)
!
!   external functions:
!   ===================
!
    external              :: dgemm, dgeqrf, dtrsm, dnrm2, dtrcon
!
!   get memory for the metric.
!
    allocate (metric(m,m), msave(m,m))
    metric = zero
    macro_done = .false.
!
!   assemble the metric
!
    it = 0
    growth = one
    do while(.not. macro_done)
      it = it + 1
      if (it .gt. maxit) then
!
!       ortho_cd failed. return with an error message
!
        ok = .false.
        write(6,100) ' maximum number of iterations reached.'
        return
      end if
      call dgemm('t','n',m,m,n,one,u,n,u,n,zero,metric,m)
      msave = metric
!
!   compute the cholesky factorization of the metric.
!
      call dpotrf('l',m,metric,m,info)
!
!     if dpotrf failed, try a second time, after level-shifting the diagonal of the metric.
!
      if (info.ne.0) then
!
        alpha      = 100.0_dp
        unorm      = dnrm2(n*m,u,1)
        it_micro   = 0
        micro_done = .false.
!
!       add larger and larger shifts to the diagonal until dpotrf manages to factorize it.
!
        do while (.not. micro_done)
          it_micro = it_micro + 1
          if (it_micro.gt.maxit) then
!
!           something went very wrong. return with an error status, the orthogonalization
!           will be carried out using a different algorithm.
!
            ok = .false.
            write(6,100) ' maximum number of iterations for factorization reached.'
            stop
            return
          end if
!
          shift = max(epsilon(one)*alpha*unorm,tol_ortho)
          metric = msave
          call diag_shift(m,shift,metric)
          call dpotrf('l',m,metric,m,info)
          alpha = alpha * 10.0_dp
          micro_done = info.eq.0
        end do
!
      end if
!
!     we assume that the error on the orthogonality is of order k(l)^2 * eps,
!     where eps is the machine precision.
!     the condition number k(l) is estimated by computing
!
!     k(l) ||l|| ||l^-1||,
!
!     where the norm used is the following (see norm_estimate):
!
!     || A || = || D + O || <= || D ||_inf + || O ||_2
!
!     compute l^-1, using msave to store the inverse cholesky factor
!
!     orthogonalize u by computing a solution to u(ortho) l^t = u
!     if required, apply the same transformation to w.
!
      msave = metric
      call dtrtri('l','n',m,msave,m,info)
!
!     compute the norm of l, l^-1 and the condition number:
!
      l_norm    = norm_est(m,metric)
      linv_norm = norm_est(m,msave)
      rcond     = l_norm * linv_norm
!
!     in each iteration of ortho_cd, we apply l^-t to u, which introduces
!     a numerical error of order ||l^-1||. 
!     this error is saved in growth and used in ortho_vs_x to check how much
!     ortho_cd spoiled the previously computed orthogonality to x. 
!
      growth    = growth * linv_norm
!
!     orthogonalize u by applying l^(-t)
!    
      call dtrmm('r','l','t','n',n,m,one,msave,m,u,n)
!
!     check the error:
!
      error      = epsilon(one) * rcond*rcond
      macro_done = error .lt. tol_ortho_cd
    end do
!
    100 format(t3,'ortho_cd failed with the following error:',a)
!
    ok = .true.
!
    deallocate (metric)

!      call dtrsm('r','l','t','n',n,m,one,metric,m,u,n)
!      if (do_other) call dtrsm('r','l','t','n',n,m,one,metric,m,w,n)
!
!     check that the vectors in v are really orthogonal
!
!      call dgemm('t','n',m,m,n,one,u,n,u,n,zero,metric,m)
!      metric_norm = abs(dnrm2(m*m,metric,1) - sqrt(dble(m)))
!      macro_done = metric_norm .lt. tol_ortho
!    end do
!
!    100 format(t3,'ortho_cd failed with the following error:',a)
!
!    ok = .true.
!
!    deallocate (metric)
    return
  end subroutine ortho_cd
!
  real(dp) function norm_est(m,a)
!
!   compute a cheap estimate of the norm of a lower triangular matrix.
!   let a = d + o, where d = diag(a). as
!
!   || a || <= || d || + || o ||
!
!   we compute || d || as max_i |d(i)| and || o || as its frobenius norm.
!   this is tight enough, and goes to 1 when a approaches the identity.
!
    implicit none
    integer,                  intent(in) :: m
    real(dp), dimension(m,m), intent(in) :: a
!
    integer  :: i, j
    real(dp) :: diag_norm, od_norm
!
    diag_norm = zero
    do i = 1, m
      diag_norm = max(diag_norm,abs(a(i,i)))
    end do
!
    od_norm = zero
    do i = 1, m
      do j = 1, i - 1
        od_norm = od_norm + a(i,j)**2
      end do
    end do
    od_norm = sqrt(od_norm)
!
    norm_est = diag_norm + od_norm
    return
  end function norm_est
!
  subroutine ortho_vs_x(n,m,k,x,u,ax,au,dropping)
    implicit none
!
!   given two sets x(n,m) and u(n,k) of vectors, where x 
!   is assumed to be orthogonal, orthogonalize u against x.
!
!   if required, orthogonalize au to ax using the same linear
!   transformation, where ax and au are the results of the 
!   application of a matrix a to both x and u.
!
!   furthermore, orthonormalize u and, if required, apply the
!   same transformation to au.
!
!   this routine performs the u vs x orthogonalization and the
!   subsequent orthonormalization of u iteratively, until the
!   overlap between x and the orthogonalized u is smaller than
!   a (tight) threshold. 
!
!   arguments:
!   ==========
!
    integer,                   intent(in)    :: n, m, k
    real(dp),  dimension(n,m), intent(in)    :: x, ax
    real(dp),  dimension(n,k), intent(inout) :: u, au
    logical,   optional,       intent(in)    :: dropping
!
!   local variables:
!   ================
!
    logical                :: done, ok
    integer                :: it 
    real(dp)               :: xu_norm, growth
    real(dp),  allocatable :: xu(:,:)
!
!   external functions:
!   ===================
!
    intrinsic              :: random_number
    real(dp)               :: dnrm2
    external               :: dnrm2, dgemm
!   
    integer, parameter     :: maxit = 10
    logical, parameter     :: useqr = .false.
!
!   allocate space for the overlap between x and u.
!
    ok = .false.
    allocate (xu(m,k))
    done = .false.
    it   = 0
!
!   start with an initial orthogonalization to improve conditioning.
!
    if(.not.present(dropping)) then
      if (.not. useqr) call ortho_cd(n,k,u,growth,ok)
      if (.not. ok .and. useqr) call ortho(n,k,u,au)
    end if
!
!   iteratively orthogonalize u against x, and then orthonormalize u.
!
    do while (.not. done)
      it = it + 1
!
!     u = u - x (x^t u)
!
      call dgemm('t','n',m,k,n,one,x,n,u,n,zero,xu,m)
      call dgemm('n','n',n,k,m,-one,x,n,xu,m,one,u,n)
!
      if (present(dropping)) then
              
      endif
!
!     now, orthonormalize u.
!     
      if (.not. useqr) call ortho_cd(n,k,u,growth,ok)
      if (.not. ok .or. useqr) call ortho(n,k,u,au)
!
!     the orthogonalization has introduced an error that makes the new
!     vector no longer fully orthogonal to x. assuming that u was
!     orthogonal to x to machine precision before, we estimate the
!     error with growth * eps, where growth is the product of the norms
!     of all the linear transformations applied to u.
!     if ortho_cd has failed, we just compute the overlap and its norm.
!
!
      if (.not. ok .or. useqr) then
        call dgemm('t','n',m,k,n,one,x,n,u,n,zero,xu,m)
        xu_norm = dnrm2(m*k,xu,1)
      else
        xu_norm = growth * epsilon(one)
      end if
      done    = xu_norm.lt.tol_ortho
!
!     if things went really wrong, abort.
!
      if (it.gt.maxit) stop ' catastrophic failure of ortho_vs_x'
    end do
!
    deallocate(xu)
!
    return
  end subroutine ortho_vs_x
!
  subroutine b_ortho_vs_x(n,m,k,x,bx,u,dropping,tol_o)
    implicit none
!
!   given two sets x(n,m) and u(n,k) of vectors, where x 
!   is assumed to be orthogonal, b-orthogonalize u against x.
!   furthermore, orthonormalize u.
!
!   this routine performs the u vs x orthogonalization and the
!   subsequent orthonormalization of u iteratively, until the
!   overlap between x and the orthogonalized u is smaller than
!   a (tight) threshold. 
!
!   arguments:
!   ==========
!
    integer,                   intent(in)    :: n, m
    integer,                   intent(inout) :: k 
    real(dp), dimension(n,m),  intent(in)    :: x, bx
    real(dp), dimension(n,k),  intent(inout) :: u
    logical,  optional,        intent(in)    :: dropping
    real(dp), optional,        intent(in)    :: tol_o
!
!   local variables:
!   ================
!
    logical                :: done, ok
    integer                :: it, i, j 
    real(dp)               :: xu_norm, growth, xx(1), u_norm
    real(dp),  allocatable :: xu(:,:)
!
!   external functions:
!   ===================
!
    intrinsic              :: random_number
    real(dp)               :: dnrm2
    external               :: dnrm2, dgemm
!   
    integer, parameter     :: maxit = 10
    logical, parameter     :: useqr = .false.
!
!   allocate space for the overlap between x and u.
!
    ok = .false.
    allocate (xu(m,k))
    done = .false.
    it = 0
!
!   start with an initial orthogonalization to improve conditioning.
!
    if (.not.present(dropping)) then
      if (.not. useqr) call ortho_cd(n,k,u,growth,ok)
      if (.not. ok .or. useqr) call ortho(n,k,u,xx)
    end if
!
!   iteratively orthogonalize u against x, and then orthonormalize u.
!
    do while (.not. done)
      it = it + 1
!
!     u = u - x (bx^t u)
!
      call dgemm('t','n',m,k,n,one,bx,n,u,n,zero,xu,m)
      call dgemm('n','n',n,k,m,-one,x,n,xu,m,one,u,n)
!
!     perform the vector drop
!     
      if (present(dropping)) then
        call vector_drop(u,n,k,tol_o)  
      end if
!
!     now, orthonormalize u.
!
      if (.not. useqr) call ortho_cd(n,k,u,growth,ok)
      if (.not. ok .or. useqr) call ortho(n,k,u,xx)
!
!     compute the overlap between the orthonormalized u and x and decide
!     whether the orthogonalization procedure converged.
!
!     note that, if we use ortho_cd, we estimate the norm of the overlap
!     using the growth factor returned in growth.
!     see ortho_vs_x for more information.
!
     if (.not. ok .or. useqr) then
       call dgemm('t','n',m,k,n,one,bx,n,u,n,zero,xu,m)
       xu_norm = dnrm2(m*k,xu,1)
     else
       xu_norm = growth * epsilon(one)
     end if  
       done    = xu_norm.lt.tol_ortho
!
!     if things went really wrong, abort.
!
      if (it.gt.maxit) stop ' catastrophic failure of b_ortho_vs_x'
    end do
!
    deallocate(xu)
!
    return
  end subroutine b_ortho_vs_x
!
! utilities
! =========
!
  subroutine diag_shift(n,shift,a)
    implicit none
!  
!   add shift to the diagonal elements of the matric a
!  
    integer,                  intent(in)    :: n
    real(dp),                 intent(in)    :: shift
    real(dp), dimension(n,n), intent(inout) :: a
!  
    integer :: i
!  
    do i = 1, n
      a(i,i) = a(i,i) + shift
    end do
!  
    return
  end subroutine diag_shift
!
  subroutine get_coeffs(len_a,len_u,n_max,n_act,a_red,u_x,u_p)
    implicit none
!
!   given the eigenvetors of the reduced matrix in a_red, extract
!   the expansion coefficients for x_new (u_x) and assemble the 
!   ones for p_new in u_p.
!
!   the coefficients u_p are computed as the difference between the
!   coefficients for x_new and x_old, and only the columns associated
!   with active eigenvectors are considered. 
!   u_p is then orthogonalized to u_x: this not only guarantees that
!   the p_new vectors will be orthogonal to x_new, but also allows one
!   to reuse the ax, aw, and ap vectors to compute ap_new, without
!   loosing numerical precision.
!
    integer,                          intent(in)    :: len_a, len_u, n_max, n_act
    real(dp), dimension(len_a,len_a), intent(in)    :: a_red
    real(dp), dimension(len_u,n_max), intent(inout) :: u_x
    real(dp), dimension(len_u,n_act), intent(inout) :: u_p
!
    integer               :: ind_x, off_x, i_eig
    real(dp)              :: xx(1)
!
    off_x = n_max - n_act
    ind_x = off_x + 1
!
    u_x(1:len_u,1:n_max) = a_red(1:len_u,1:n_max)
!
!   u_p = u_x for the active vectors only
!
    u_p = u_x(:,ind_x:n_max)
!
!   remove the coefficients for x from u_p
!
    do i_eig = 1, n_act
      u_p(off_x + i_eig,i_eig) = u_p(off_x + i_eig,i_eig) - one
    end do
!
!   orthogonalize:
!
!    call ortho_vs_x(.false.,len_u,n_max,n_act,u_x,u_p,xx,xx)
    call ortho_vs_x(len_u,n_max,n_act,u_x,u_p,xx,xx)
!
!   all done.
!
    return
!
  end subroutine get_coeffs
!
  subroutine check_guess(n,m,evec)
    implicit none
    integer,                  intent(in)    :: n, m
    real(dp), dimension(n,m), intent(inout) :: evec
!
    integer               :: i, j, istat
!    real(dp)              :: fac, diag_norm, out_norm, xx(1)
    real(dp)              :: fac, diag_norm, out_norm, growth, xx(1)
    logical               :: ok
!
    real(dp), allocatable :: overlap(:,:)
    real(dp)              :: dnrm2
    external              :: dnrm2
!
!   check whether evec is zero.
!
    fac = dnrm2(n*m,evec,1)
    if (fac.eq.zero) then
!
!     no luck. make a random guess, then orthonormalize it.
!
      call random_number(evec)
!      call ortho(.false.,n,m,evec,xx)
      call ortho_cd(n,m,evec,growth,ok)
    else
!
!     compute the overlap and check that the vectors are orthonormal.
!
      allocate (overlap(m,m), stat=istat)
      call check_mem(istat)
      call dgemm('t','n',m,m,n,one,evec,n,evec,n,zero,overlap,m)
      diag_norm = zero
      out_norm  = zero
      do i = 1, m
        diag_norm = diag_norm + overlap(i,i)**2
        do j = 1, i-1
          out_norm = out_norm + overlap(j,i)**2
        end do
      end do
!
      diag_norm = diag_norm/real(m,dp)
!
      if (diag_norm .ne. one .or. out_norm.ne.zero) then
!
!       orthogonalize the guess:
!
!        call ortho(.false.,n,m,evec,xx)
        call ortho_cd(n,m,evec,growth,ok)
      end if
!
      deallocate (overlap, stat = istat)
      call check_mem(istat)
    end if
!
    return
  end subroutine check_guess
!
!
  subroutine check_mem(istat)
    implicit none
!
!   silly routine that checks that istat is zero after memory allocation
!   and aborts otherwise printing an error message.
!
    integer, intent(in) :: istat
!
    1000 format(t3,'memory allocation failed. istat = ',i8)
    if (istat.ne.0) then
      write(6,1000) istat
      stop
    end if
    return
  end subroutine check_mem
!
  integer function get_mem_lapack(n,n_max)
    integer, intent(in)    :: n, n_max
!
    integer           :: lwork1, lwork2, len_rr, len_qr, nb
!fl
    integer           :: lwork3
    integer, external :: ilaenv
!
!   maximum size of the rayleigh-ritz matrix:
!
    len_rr = 3*n_max
!
!   maximum size of the space to orthonormalize:
!
    len_qr = 6*n_max
!
!   use lapack query routines to compute the optimal memory required
!   for diagonalization and QR decomposition.
!
    nb     = ilaenv( 1, 'DSYTRD', 'l', len_rr, -1, -1, -1 )
    lwork1 = len_rr * nb
!
    nb     = ilaenv( 1, 'DGEQRF', 'l', n, len_qr, -1, -1 )
    lwork2 = len_qr*nb
!
    nb     = ilaenv( 1, 'DSYTRD', 'l', len_rr, -1, -1, -1)
    lwork3 = len_rr * nb
!
    get_mem_lapack = max(lwork1,lwork2,lwork3)
    return
  end function get_mem_lapack
! 
  subroutine prtmat(n,lda,mat,iform)
!
! prints a symmetric matrix
!   
    integer,                      intent(in) :: n, lda, iform
    real(dp), dimension(lda,lda), intent(in) :: mat
!
    integer :: i
!
    100 format(t3,20f12.6)
    200 format(t3,20d12.4)
!
    do i = 1, n
      if (iform.eq.1) write(6,100) mat(i,1:n)
      if (iform.eq.2) write(6,200) mat(i,1:n)
    end do
    return
  end subroutine prtmat
!
  subroutine prt_generic(n,m,nline,ncol,mat,iform,mat_name)
!
!  prints nxm matrix derived from a generic nlinexncol matrix
!
   implicit none
   integer, intent(in)                          :: n, m, nline, ncol, iform
   real(dp), dimension(nline, ncol), intent(in) :: mat
   character(len=*), intent(in)                 :: mat_name 
!
   integer :: i, j
   character(len=80) :: fmt
!
! Print the name of the matrix
!
    print *, "Matrix: ", trim(mat_name)
    print *, ""
!
!  Define the formats for the output
!
   if (iform .eq. 1) then
     fmt = '(t3, ' // trim(adjustl(itoa(m))) // 'f12.6)'
   else if (iform .eq. 2) then
     fmt = '(t3, ' // trim(adjustl(itoa(m))) // 'd12.4)'
   end if
!
!  Loop over the rows of the matrix and print each row
!
   do i = 1, n
     write(6, fmt) (mat(i, j), j = 1, m)
   end do
!
   return
  end subroutine prt_generic

! Helper function to convert integer to string
  function itoa(i)
   implicit none
   integer, intent(in) :: i
   character(len=20) :: itoa
   write(itoa, '(I0)') i
  end function itoa
!  
  subroutine get_time(t)
    real(dp), dimension(2), intent(inout) :: t
!
!$  real(dp) :: omp_get_wtime
!$  external :: omp_get_wtime
!
!   get cpu and (if openmp is available) wall time.
!
    t = zero
    call cpu_time(t(1))
!$  t(2) = omp_get_wtime()
!
    return
  end subroutine get_time
!
 subroutine vector_drop(u,n,k,tol_o)
   implicit none
!
!  checks a vector norm and, based on a set tolerance, keeps or drops
!  the said vector. A logical is defined to keep track of the dropped
!  vector's positions.

!  imput variables:
!  ================
!
   integer,  intent(inout)                  :: k
   integer,  intent(in)                     :: n
   real(dp), dimension(n,k), intent(inout)  :: u
   real(dp), optional, intent(in)           :: tol_o
!
!  local variables:
!  ================
!
   real(dp)                                :: u_norm, dnrm2, tol_drop
   integer                                 :: i, j
   integer                                 :: n_drop
   logical                                 :: restart
!
!  external functions:
!  ===================
   external                                :: dnrm2 
!
!  set the tolerance to a certain value and
!  initialize the integer that keeps track
!  of the number of vectors removed:
!
   tol_drop = min(tol_o*1.0e-2_dp,1.0e-8_dp)
   n_drop   = 0
! 
!  set the restart logical to false. If true,
!  the routine will restart with a different
!  tolerance.
!
   restart = .false.
!
!  check vector norms and decide whether to drop 
!  them or not:
! 
   u_norm = 0.0_dp
   do
     do j = 1, k
       if (restart) then
         restart = .false.
         exit
       end if
       u_norm = dnrm2(n,u(:,j-n_drop),1) 
       if (u_norm .lt. tol_drop .and.(k-n_drop).gt.1) then
         u(:,j-n_drop) = u(:,k-n_drop)
         n_drop = n_drop + 1
       else if (u_norm .lt. tol_drop .and.(k-n_drop).eq.1) then
         restart  = .true.
         exit
       end if   
     end do
     if(.not.restart) exit
   end do
   k = k - n_drop
   print*, "n_drop=", n_drop
!
  end subroutine vector_drop  
end module diaglib
