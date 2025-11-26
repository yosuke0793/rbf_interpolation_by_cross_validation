module util_math
implicit none
contains
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine invmat(N,A)
  !----------------------------------------------------------------------
  use mpival
  ! --- global
  implicit none
  integer,intent(in)::N
  double precision,intent(inout)::A(N,N)
  ! --- local
  integer::i,info,lwork,ipiv(1:N)
  double precision::tw(1:1)
  double precision,allocatable::work(:)
  !======================================================================
  ! --- LU DECOMPOSITION
  info=0
  call dgetrf(N,N,A,N,ipiv,info)
  if(info.ne.0)then
  !WRITE(*,*)"ERROR AT INVERSE MATRIX DGETRF, INFO",INFO
    A = 0.d0
    do i=1,N
      A(i,i) = 1.d0
    enddo
    return
  endif
  !
  ! --- GET SIZE OF WORK
  call dgetri(N,A,N,ipiv,tw,-1,info)
  lwork=nint(tw(1))
  allocate(work(1:lwork+1))
  work=0.d0
  !
  ! --- INVERSE
  call dgetri(N,A,N,ipiv,work,lwork,info)
  if(info.ne.0)then
  !WRITE(*,*)"ERROR AT INVERSE MATRIX DGETRI, INFO",INFO
    A = 0.d0
    do i=1,N
      A(i,i) = 1.d0
    enddo
    return
  endif
  !
  !======================================================================
  return
  !----------------------------------------------------------------------
  end subroutine invmat
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine pardisoT(neq, nindx, AA, indx, iptr, BB, xx)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! solve Ax=B for x
  implicit none
  integer, intent(in) :: neq, nindx, indx(nindx), iptr(neq+1)
  double precision, intent(inout) :: AA(nindx), BB(neq), xx(neq)
  !
  integer :: iparm(128), pt(128)
  integer :: maxfct, mnum, mtype, phase, nrhs, error, msglvl
  integer :: idum
  !
  double precision, allocatable :: yy(:)
  !======================================================================
  allocate(yy(neq)); yy=0
  !
  mtype = 11 
  call pardisoinit(pt,mtype,iparm)
  iparm(1) = 0
  phase = 13
  !
  nrhs = 1
  msglvl = 0
  maxfct = 1
  mnum = 1
  !
  call pardiso(pt, maxfct, mnum, mtype, phase, neq, AA, iptr, indx, &
              idum, nrhs, iparm, msglvl, BB, yy, error)
  !
  phase = -1
  call pardiso(pt, maxfct, mnum, mtype, phase, neq, AA, iptr, indx, &
              idum, nrhs, iparm, msglvl, BB, yy, error)
  !
  xx(:) = yy(:)
  deallocate(yy)
  !
  !----------------------------------------------------------------------
  return
  end subroutine pardisoT
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
end module util_math