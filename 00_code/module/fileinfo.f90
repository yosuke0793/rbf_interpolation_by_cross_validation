!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module fileinfo
!================================================================
implicit none
! --- FILENAME
character(99) :: filename
! --- INPUTS
integer :: inp_condt , inp_origin
integer :: inp_pod   , inp_expln , inp_objct, inp_restat
! --- OUTPUTS
integer :: out_reprod, out_hist_err, out_best, out_log, out_first
integer :: out_stplog, out_restat,   out_rbf,  out_trnval
! --- INTERMIDIATE
integer :: intermid
integer :: random_table
!================================================================
end module fileinfo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
