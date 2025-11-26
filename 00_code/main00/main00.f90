!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 差分進化法による最適化プログラム
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program differential_evolution
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! --- 必要なモジュールの読み込み
use fileinfo      ! ファイル入出力関連の情報
use commonval     ! 共通の変数・定数
use mpival        ! MPI並列計算用の値
use mpi           ! MPIライブラリ（並列計算用）
use util_common   ! 共通ユーティリティ関数
use util_de       ! 差分進化法のユーティリティ
use util_rbf      ! 動径基底関数（Radial Basis Function）のユーティリティ
! -------------------------------------------------------------------------
implicit none

! --- フラグ変数群
! 各種処理の状態を管理するフラグ
integer          :: flag_post      ! 後処理フラグ
integer          :: flag_strategy  ! DE戦略の選択フラグ（DE/rand/1, DE/best/1等）
integer          :: flag_update    ! パラメータ更新フラグ
integer          :: flag_finish    ! 計算終了フラグ

! --- 最適化イテレーション
integer          :: nopt           ! 最適化の総イテレーション数
integer          :: iopt           ! 現在のイテレーション番号

! --- データ関連
integer          :: ndata          ! データ点の総数
integer          :: nexp           ! 実験データの数
integer          :: nobj           ! 目的関数の数（多目的最適化の場合）

! --- エージェント（個体群）関連
! DEでは複数の解候補（エージェント）を持つ集団ベース手法
integer          :: npara          ! 最適化パラメータの総数
integer          :: np_fix         ! 固定パラメータの数
integer          :: ip_fix         ! 固定パラメータのインデックス
integer          :: np_target      ! ターゲットパラメータの数
integer          :: ip_target      ! ターゲットパラメータのインデックス
integer          :: nagent         ! エージェントの総数
integer          :: iagent         ! 現在のエージェント番号
integer          :: nbest          ! ベストエージェントのインデックス

! --- エージェントの探索パラメータ
double precision :: mf             ! 突然変異係数
double precision :: sf             ! 交叉率
double precision :: mcr            ! 突然変異確率
double precision :: scr            ! 交叉確率

! --- エージェントおよびパラメータの配列
double precision, allocatable :: f1(:), f2(:), cr(:)          ! 探索パラメータ
double precision, allocatable :: agent(:,:), para(:)         ! エージェントとそのパラメータ
double precision, allocatable :: err_agent(:), agent_string(:) ! エージェントの誤差と文字列表現
double precision, allocatable :: agent_best(:)               ! ベストエージェント
double precision, allocatable :: recv_buff(:)                ! 受信バッファ

! --- 計算時間のカウンタ
integer          :: count1, count2, count3,  count_rate
double precision :: time

! --- 収束関連
double precision :: tole            ! 収束許容誤差
double precision :: dev_percent     ! 偏差率
double precision :: err_gbest       ! グローバルベストエージェントの誤差
double precision :: derr_best       ! ベストエージェントの誤差変化量
double precision :: error           ! 誤差
double precision :: dev_error       ! 誤差の偏差
double precision :: mean_error      ! 誤差の平均値

! --- 元データ
double precision, allocatable :: original(:,:)
double precision, allocatable :: original_s(:)
! =========================================================================
!
! --- INITIALIZE MPI
call MPI_INIT(ierr) ! MPIの初期化
call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr) ! プロセス数の取得
call MPI_COMM_RANK(MPI_COMM_WORLD,irank,ierr)  ! プロセスランクの取得
!
! --- Set files
call setfile ! ファイル設定の初期化
!
! --- Read DE conditions
call readin(flag_strategy, nopt, nagent, ncv, ndata, nexp, nobj, &
            npara, np_target, np_fix, mf, sf, mcr, scr, tole  ) ! DE条件の読み込み
!
! --- Initialization 
allocate( f1(nagent),f2(nagent), cr(nagent) )
allocate( agent(nagent,np_target) ) 
allocate( agent_best(np_target), err_agent(nagent)    )
allocate( para(npara), agent_string(nagent*np_target) )
allocate( original(nexp+nobj,ndata)   )
allocate( original_s((nexp+nobj)*ndata)   )
flag_finish = 0
flag_update = 0
flag_post   = 0
err_gbest   = 1.D0
nstart = 1+irank *nagent/nprocs
nend   =(1+irank)*nagent/nprocs
!
! --- Read data
call read_data(ndata, nexp, nobj, original,original_s) ! データの読み込み
!
! --- Create agents
call createagent(nagent, np_target, mf, mcr, f1, f2, cr, agent, agent_string )! エージェントの生成
call bcast_agents(nagent, np_target, ndata, nexp, nobj, &
                  f1, f2, cr, agent, agent_string, &
                  original, original_s ) ! エージェントのブロードキャスト
!
! --- Set clock
call stplog(0, 0, 0.d0, 0.d0, 1.d0, count1, count2, count3, count_rate, time)
!
! --- Parameter fitting
do iopt = 1,nopt
  call stplog(1, 0, 0.d0, 0.d0, 1.d0, count1, count2, count3, count_rate, time)
  !
  do iagent = nstart,nend
    ! --- SET EXPLORATION PARAMETERS
    do ip_target = 1,np_target
      para(nint(maxpara(1,ip_target))) = agent(iagent,ip_target)
    enddo
    !
    ! --- SET FIXED PARAMETERS
    do ip_fix = 1,np_fix
      para(nint(parafix(1,ip_fix))) = parafix(2,ip_fix)
    enddo
    !
    ! --- ERROR ESTIMATION
    call cv_rbf(flag_post,flag_update,ncv, ndata, nexp, nobj, npara,para,original,error)
    err_agent(iagent) = error
    !
    ! --- Avoid NaN & overflow
    if( err_agent(iagent).ne.err_agent(iagent) .or. err_agent(iagent).gt.1.0D0/zero )then
      err_agent(iagent)= 1.0D0/zero
    endif
    !
  enddo
  !
  ! --- Gather error to HEAD process
  allocate( recv_buff(nagent) )
  recv_buff = 0.d0
  call MPI_GATHER( err_agent(nstart), nagent/nprocs, MPI_DOUBLE_PRECISION, &
                   recv_buff(1),      nagent/nprocs, MPI_DOUBLE_PRECISION, &
                   0, MPI_COMM_WORLD, ierr )
  if( irank.eq.0 )err_agent = recv_buff
  deallocate( recv_buff )
  !
  ! --- Find best
  call findbest( iopt, npara, nagent, np_target, flag_update, nbest,         &
                 flag_finish, tole, err_agent, agent, derr_best, mean_error, &
                 dev_error, agent_best , err_gbest, dev_percent, para        )
  !
  ! --- Output best results
  if( irank.eq.0 )then
    flag_post = 1
    call cv_rbf(flag_post,flag_update,ncv, ndata, nexp, nobj, npara,para,original,error)
    flag_post = 0
  endif
  !
  ! --- Update agents
  call updateagent(flag_strategy, nagent, np_target, mf, sf, mcr, scr,   &
                   err_agent, agent_best, f1, f2, cr, agent, agent_string)
  !
  ! --- Scatter agents from head process to others
  call bcast_agents(nagent, np_target, ndata, nexp, nobj, f1, f2, cr, agent, agent_string, original, original_s )
  !
  ! --- Record COMPUTATION TIME
  call stplog(2, iopt, err_gbest, dev_error, dev_percent, &
              count1, count2, count3, count_rate, time  )
  !
  ! --- Exit update loop
  call MPI_BCAST(flag_finish, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
  if( flag_finish.eq.1 )exit
  !
enddo
call stplog(3, iopt, err_gbest, dev_error, dev_percent, &
            count1, count2, count3, count_rate, time  )
call stop00
! =========================================================================
end program differential_evolution
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
