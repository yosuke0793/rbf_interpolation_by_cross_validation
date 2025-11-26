module util_de
use util_common, only: stop00
implicit none
contains
  ! ==========================================================================
  subroutine mutation(flag_strategy,nagent,np_target,mf,sf,mcr,scr, &
                      agent_best,f1,f2,cr,agent)
  ! サブルーチンmutation: 差分進化アルゴリズムに基づいてエージェントの突然変異を実行する
  ! 引数:
  !   flag_strategy : 整数, 使用する差分進化戦略のフラグ
  !   nagent        : 整数, エージェントの数
  !   np_target     : 整数, 目的変数の数
  !   mf, sf        : 実数, 突然変異係数の平均と標準偏差
  !   mcr, scr      : 実数, 交差率の平均と標準偏差
  !   agent_best    : 実数配列(1:np_target), 最良エージェントのパラメータ
  !   f1, f2        : 実数配列(1:nagent), 各エージェントの突然変異係数
  !   cr            : 実数配列(1:nagent), 各エージェントの交差率
  !   agent         : 実数配列(1:nagent,1:np_target), エージェントのパラメータ
  ! 動作:
  !   指定された差分進化戦略に基づいて、エージェントの突然変異を実行し、エージェント配列を更新する。
  ! ==========================================================================
  implicit none
  integer, intent(in)             :: flag_strategy, nagent, np_target
  double precision, intent(in)    :: mf,sf,mcr,scr 
  double precision, intent(in)    :: agent_best(np_target)
  double precision, intent(inout) :: f1(nagent), f2(nagent), cr(nagent)
  double precision, intent(inout) :: agent(nagent,np_target)
  integer                         :: ipara      
  integer                         :: iagent     , nagent_half
  integer                         :: seedsize   , n_rand     , clock
  integer                         :: np_mutation, n_inherit  , index_mutant(3)
  integer         , allocatable   :: seed(:)
  double precision                :: r1, r2, r3, rand1, rand2, rand3
  double precision                :: rand_cr, rand_mutant(3)
  double precision, allocatable   :: mutant(:)
  double precision, parameter     :: pi = 4.d0*atan(1.d0)
  !
  ! --- ALLOCATE
  allocate( mutant(np_target) )
  !
  ! --- SET SEED FOR RANDOM NUMBER 
  call random_seed(size=seedsize)
  allocate(seed(seedsize))
  call system_clock(count=clock)
  seed(:) = clock
  !
  ! --- HALF OF NUMBER OF AGENTS
  nagent_half = int(nagent/2)
  !
  ! --- SET MUTANTS
  do iagent = 1,nagent_half
    ! --- INITIALIZATION
    agent(iagent+nagent_half,:) = agent(iagent,:)
    !
    ! --- SELECT 3 AGENT RANDOMLY
    do
      seed(:) = seed(:) + 1
      call random_seed(put=seed(:))
      call random_number(rand_mutant)
      index_mutant(:) = int(rand_mutant(:)*nagent_half)
      do ipara=1,3
        if( index_mutant(ipara).eq.0 )index_mutant(ipara) = 1
      enddo
      if( index_mutant(1).ne.index_mutant(2) .and. &
          index_mutant(2).ne.index_mutant(3) .and. &
          index_mutant(3).ne.index_mutant(1) )then
        exit
      endif
    enddo
    !
    ! --- GENERATE SCALE FACTOR AND CROSSOVER RATE (jDE)
    call random_number(r1)
    call random_number(r2)
    call random_number(r3)
    rand1 = sqrt(-2.d0*log(r1))*dcos(2.d0*pi*r2)
    rand2 = sqrt(-2.d0*log(r2))*dsin(2.d0*pi*r3)
    rand3 = sqrt(-2.d0*log(r3))*dcos(2.d0*pi*r1)
    f1(iagent) = mf  +  sf*rand1 
    f2(iagent) = mf  +  sf*rand2
    cr(iagent) = mcr + scr*rand3
    f1(iagent+nagent_half) = f1(iagent)
    f2(iagent+nagent_half) = f2(iagent)
    cr(iagent+nagent_half) = cr(iagent)
    !
    ! --- GENETIC STRATEGY
    select case(flag_strategy)
    case(1) ! DE/RAND/1 
      mutant(:) = agent(index_mutant(1),:)                                 &
                + f1(iagent)*( agent(index_mutant(2),:)-agent(index_mutant(3),:) )
    case(2) ! DE/BEST/1
      mutant(:) = agent_best(:)                                            &
                + f1(iagent)*( agent(index_mutant(2),:)-agent(index_mutant(3),:) )
    case(3) ! DE/CURRENT-TO/1
      mutant(:) = agent(iagent,:)                                          &
                + f2(iagent)*( agent(index_mutant(2),:)-agent(index_mutant(3),:) )
    case(4) ! DE/CURRENT-TO-BEST/1
      mutant(:) = agent(iagent,:)                                          &
                + f1(iagent)*( agent_best(:)-agent(iagent,:) )            &
                + f2(iagent)*( agent(index_mutant(2),:)-agent(index_mutant(3),:) )
    case(5) ! DE/RAND-TO-BEST
      mutant(:) = agent(index_mutant(1),:)                                          &
                + f1(iagent)*( agent_best(:)-agent(index_mutant(1),:) )            &
                + f2(iagent)*( agent(index_mutant(2),:)-agent(index_mutant(3),:) )
    case default
      write(*,*)'FLAG_DEMETHOD is not correct! @main00/mutaiton.f90'
      call stop00
    end select
    !
    ! --- CREATE TRIAL AGENTS BY CROSS OVER
    seed(:) = seed(:) + 1
    call random_seed(put=seed(:))
    call random_number(rand_cr)
    n_rand = nint(rand_cr*np_target)
    if( n_rand.eq.0 )n_rand=1
    n_inherit = 1
    do ipara = 1,np_target
      seed(:) = seed(:) + 1
      call random_seed(put=seed(:))
      call random_number(rand_cr)
      if( rand_cr.le.cr(iagent) )then
        n_inherit = n_inherit + 1
      else
        exit
      endif
    enddo
    !
    ! --- SET LATTER HALF OF AGENTS
    do ipara = n_rand,n_rand+n_inherit-1
      np_mutation = mod(ipara,np_target)
      if(np_mutation.eq.0) np_mutation = np_target
      agent(iagent+nagent_half,np_mutation) = mutant(np_mutation)
    enddo
    !
  enddo
  ! --- DEALLOCATE
  deallocate( mutant,seed )
  !
  return
  end subroutine mutation
  ! ==========================================================================

  ! ==========================================================================
  subroutine updateagent(flag_strategy, nagent, np_target, mf, sf, mcr, scr, &
                        err_agent, agent_best, f1, f2, cr, agent, agent_string)
  ! サブルーチンupdateagent: 差分進化アルゴリズムに基づいてエージェントを更新する
  ! 引数:
  !   flag_strategy : 整数, 使用する差分進化戦略のフラグ
  !   nagent        : 整数, エージェントの数
  !   np_target     : 整数, 目的変数の数
  !   mf, sf        : 実数, 突然変異係数の平均と標準偏差
  !   mcr, scr      : 実数, 交差率の平均と標準偏差
  !   err_agent     : 実数配列(1:nagent), 各エージェントの誤差
  !   agent_best    : 実数配列(1:np_target), 最良エージェントのパラメータ
  !   f1, f2        : 実数配列(1:nagent), 各エージェントの突然変異係数
  !   cr            : 実数配列(1:nagent), 各エージェントの交差率
  !   agent         : 実数配列(1:nagent,1:np_target), エージェントのパラメータ
  !   agent_string  : 実数配列(1:nagent*np_target), MPI通信用のエージェントパラメータ
  ! 動作:
  !   指定された差分進化戦略に基づいて、エージェントを更新し、エージェント配列とエージェント文字列を更新する。
  ! ==========================================================================
  use commonval, only: minpara, maxpara, flag_intgr
  use mpival, only: irank, nprocs
  implicit none
  integer,          intent(in)    :: flag_strategy
  integer,          intent(in)    :: nagent, np_target
  double precision, intent(in)    :: mf, sf, mcr, scr
  double precision, intent(in)    :: agent_best(np_target)
  double precision, intent(inout) :: f1(nagent), f2(nagent), cr(nagent)
  double precision, intent(inout) :: agent(nagent,np_target), err_agent(nagent)
  double precision, intent(inout) :: agent_string(nagent*np_target)
  integer                         :: ipara, ip_target
  integer                         :: iagent, iagent_s, nagent_half
  !
  if( irank.ne.0 )then
    return
  endif
  !
  ! --- Selection: Parent or child
  nagent_half = int(nagent/2)
  do iagent = 1,nagent_half
    if( err_agent(iagent).gt.err_agent(iagent+nagent_half) )then
      agent(iagent,:) = agent(iagent+nagent_half,:)
    endif
  enddo
  !
  ! --- MUTATION
  call mutation(flag_strategy,nagent,np_target,mf,sf,mcr,scr,agent_best,f1,f2,cr,agent)
  !
  ! --- CONSTRAINT CONDITION
  do iagent = 1,nagent
    do ipara = 1,np_target
      if( agent(iagent,ipara).lt.minpara(2,ipara) )then
        agent(iagent,ipara) = minpara(2,ipara)
      elseif( agent(iagent,ipara).gt.maxpara(2,ipara) )then
        agent(iagent,ipara) = maxpara(2,ipara)
      endif
      if( flag_intgr(ipara).eq.1 )then
        agent(iagent,ipara)=int(agent(iagent,ipara))
      endif
    enddo
  enddo
  !
  ! --- MAKE AGENT STRING
  if( nprocs.gt.1 )then
    do iagent=1,nagent
      do ip_target = 1,np_target
        iagent_s = ip_target + (iagent-1)*np_target
        agent_string(iagent_s) = agent(iagent,ip_target)
      enddo
    enddo
  endif
  !
  return
  end subroutine updateagent
  ! ==========================================================================

  ! ==========================================================================
  subroutine createagent(nagent, np_target, mf, mcr, f1, f2, cr, agent, agent_string )
  ! サブルーチンcreateagent: ラテンハイパーキューブサンプリングに基づいてエージェントを初期化する
  ! 引数:
  !   nagent       : 整数, エージェントの数
  !   np_target    : 整数, 目的変数の数
  !   mf           : 実数, 突然変異係数の平均
  !   mcr          : 実数, 交差率の平均
  !   f1, f2       : 実数配列(1:nagent), 各エージェントの突然変異係数
  !   cr           : 実数配列(1:nagent), 各エージェントの交差率
  !   agent        : 実数配列(1:nagent,1:np_target), エージェントのパラメータ
  !   agent_string : 実数配列(1:nagent*np_target), MPI通信用のエージェントパラメータ
  ! 動作:
  !   ラテンハイパーキューブサンプリングに基づいてエージェントを初期化し、エージェント配列とエージェント文字列を設定する。
  ! ==========================================================================
  use commonval, only: minpara, maxpara, flag_intgr
  use fileinfo
  use mpival
  implicit none
  ! --- Global
  integer, intent(in)    :: nagent,np_target
  double precision, intent(in)  :: mf,mcr
  double precision, intent(out) :: f1(nagent), f2(nagent), cr(nagent)
  double precision, intent(out) :: agent(nagent,np_target)
  double precision, intent(out) :: agent_string(nagent*np_target)
  ! --- Local
  integer                      :: ipara, ip_target
  integer                      :: iagent, iagent_s
  integer                      :: nseed, id
  integer                      :: seedsize, clock
  integer, allocatable         :: seed(:), ids(:,:)
  double precision             :: rand
  !
  if( irank.ne.0 )then
    return
  endif
  !
  ! --- Latin Hypercube Sampling
  nseed = 0.d0
  !
  ! --- SET SEED FOR RANDOM NUMBER 
  call random_seed(size=seedsize)
  allocate(seed(seedsize))
  call system_clock(count=clock)
  seed(:) = clock
  call random_seed(put=seed(:))
  !
  ! --- MAKE IDS ARRAY
  allocate( ids(nagent,np_target) )
  do iagent = 1,nagent
    ids(iagent,:) = iagent
  enddo
  !
  ! --- MAKE UNIFORM RANDOM NUMBER
  agent =0.d0
  do iagent = 1,nagent
    do ipara=1,np_target
      do 
        call random_number(rand)
        id = int(rand*nagent)+1
        if( ids(id,ipara).ne.0 )then
          ids(id,ipara)=0
          exit
        endif
      enddo
      call random_number(rand)
      agent(iagent,ipara)=(rand+id-1)/dble(nagent)
    enddo
  enddo
  !
  ! --- MAKE PARAMETERS
  do iagent = 1,nagent
    agent(iagent,:) = minpara(2,:) + ( maxpara(2,:)-minpara(2,:) )*agent(iagent,:)
    do ipara = 1,np_target
      if( agent(iagent,ipara).lt.minpara(2,ipara) )then
        agent(iagent,ipara)=minpara(2,ipara)
      elseif( agent(iagent,ipara).gt.maxpara(2,ipara) )then
        agent(iagent,ipara)=maxpara(2,ipara)
      endif
      if( flag_intgr(ipara).eq.1 )then
        agent(iagent,ipara)=int(agent(iagent,ipara))
      endif
    enddo
    ! --- Coefficient for mutation and cross-over
    f1(iagent) = mf
    f2(iagent) = mf
    cr(iagent) = mcr
  enddo
  !
  ! --- Make agent array for MPI
  if( nprocs.gt.1 )then
    do iagent=1,nagent
      do ip_target = 1,np_target
        iagent_s = ip_target + (iagent-1)*np_target
        agent_string(iagent_s) = agent(iagent,ip_target)
      enddo
    enddo
  endif
  !
  return
  end subroutine createagent
  ! ==========================================================================

  ! ==========================================================================
  subroutine findbest( iopt       , npara    , nagent     , np_target, & 
                      flag_update, nbest    , flag_finish,            &
                      tol_global , err_agent, agent      , derr_best, &
                      mean_error , dev_error, agent_best , err_gbest, &
                      dev_percent, para                               )
  ! サブルーチンfindbest: 最良エージェントを特定し、収束を評価する
  ! 引数:
  !   iopt         : 整数, 最適化の反復回数
  !   npara        : 整数, パラメータの数
  !   nagent       : 整数, エージェントの数
  !   np_target    : 整数, 目的変数の数
  !   flag_update  : 整数, 最良エージェントが更新されたかどうかのフラグ
  !   nbest        : 整数, 最良エージェントのインデックス
  !   flag_finish  : 整数, 収束フラグ
  !   tol_global   : 実数, 収束の許容誤差
  !   err_agent    : 実数配列(1:nagent), 各エージェントの誤差
  !   agent        : 実数配列(1:nagent,1:np_target), エージェントのパラメータ
  !   derr_best    : 実数, 最良エージェントの誤差の変化量
  !   mean_error   : 実数, エージェントの誤差の平均
  !   dev_error    : 実数, エージェントの誤差の標準偏差
  !   agent_best   : 実数配列(1:np_target), 最良エージェントのパラメータ
  !   err_gbest    : 実数, 最良エージェントの誤差
  !   dev_percent  : 実数, エージェントのパラメータの相対的な標準偏差
  !   para         : 実数配列(1:npara), 最良パラメータ
  ! 動作:
  !   最良エージェントを特定し、収束を評価し、必要に応じて最良パラメータを出力する。 
  use commonval, only: minpara, maxpara, flag_expnd
  use fileinfo, only: out_best, out_hist_err
  use mpival
  implicit none
  integer         , intent(in)    :: iopt
  integer         , intent(in)    :: npara, nagent, np_target
  integer         , intent(inout) :: flag_update, nbest
  integer         , intent(inout) :: flag_finish
  double precision, intent(in)    :: tol_global
  double precision, intent(in)    :: err_agent(nagent), agent(nagent,np_target)
  double precision, intent(inout) :: derr_best, mean_error, dev_error
  double precision, intent(inout) :: agent_best(np_target)
  double precision, intent(inout) :: err_gbest, dev_percent, para(npara)
  integer                         :: ipara, ip_target, iagent
  double precision                :: gbest_old
  double precision                :: dev_percent_0
  double precision                :: dev_agent(np_target), dev_agent_0(np_target)
  double precision                :: mean_agent(np_target)
  double precision, parameter     :: expnd = 1.0D-01
  !
  if( irank.ne.0 )then
    return
  else
    flag_update = 0
    gbest_old = err_gbest
    do iagent = 1,nagent
      if( err_agent(iagent).lt.err_gbest .or. iagent*iopt.eq.1 )then
        flag_update = 1
        nbest  = iagent
        err_gbest = err_agent(iagent)
        do ip_target = 1,np_target
          agent_best(ip_target) = agent(iagent,ip_target)
        enddo
        ! --- Check if best parameter reach limit of exploration range and expand 
        do ip_target = 1,np_target
          if( agent_best(ip_target).eq.minpara(2,ip_target) .and. flag_expnd(ip_target).eq.1 )then
            minpara(2,ip_target) = minpara(2,ip_target) - expnd*dabs(minpara(2,ip_target))
          elseif( agent_best(ip_target).eq.maxpara(2,ip_target) .and. flag_expnd(ip_target).eq.1 )then
            maxpara(2,ip_target) = maxpara(2,ip_target) + expnd*dabs(maxpara(2,ip_target))
          endif
        enddo
      endif
    enddo
    derr_best = gbest_old - err_gbest
    !
    ! --- Means
    mean_agent = 0.d0
    mean_error = 0.d0
    do iagent = 1,nagent
      mean_error = mean_error + err_agent(iagent)/nagent
      do ip_target=1,np_target
        mean_agent(ip_target) = mean_agent(ip_target) + agent(iagent,ip_target)/nagent
      enddo
    enddo
    !
    ! --- Deviations
    dev_error     = 0.d0
    dev_agent     = 0.d0
    dev_agent_0   = 0.d0
    dev_percent   = 0.d0
    dev_percent_0 = 0.d0
    do iagent = 1,nagent
      dev_error = dev_error + ( err_agent(iagent)-mean_error )**2.d0
      do ip_target=1,np_target
        dev_agent(ip_target) = dev_agent(ip_target) + ( agent(iagent,ip_target)-mean_agent(ip_target) )**2.d0
      enddo
    enddo
    dev_error = dsqrt(dev_error/nagent)
    do ip_target=1,np_target
      dev_agent(ip_target)   = dev_agent(ip_target)/nagent
      dev_agent_0(ip_target) = ( (maxpara(2,ip_target)-minpara(2,ip_target))/2.d0 )**2.d0
      dev_percent   = dev_percent   + dev_agent(ip_target)**2.d0
      dev_percent_0 = dev_percent_0 + dev_agent_0(ip_target)**2.d0
    enddo
    dev_percent = dsqrt(dev_percent)/dsqrt(dev_percent_0)
    !
    ! --- CONVERGENCE JUDGEMENT
    if( dev_percent.lt.tol_global )then
      flag_finish = 1
    else
      flag_finish = 0
    endif
    !
    ! --- Output best parameters
    if( flag_update.eq.1 )then
      write(out_hist_err,'(2i10,99ES15.5)')iopt,nbest,err_gbest,derr_best,&
                                           ( agent(nbest,ip_target),ip_target=1,np_target)
      do ip_target = 1,np_target
        para(nint(maxpara(1,ip_target))) = agent_best(ip_target)
      enddo
      open(out_best,file='./results/best.txt',status='replace')
      write(out_best,'(a15)')'Best error'
      write(out_best,'(99E15.5)')err_gbest
      write(out_best,'(a15)')'Parameters'
      do ipara = 1,npara
        write(out_best,'(99ES15.5)')para(ipara)
      enddo
      close(out_best)
    endif
    !
  endif
  !
  return
  end subroutine findbest 
  ! ==========================================================================
  
  ! ==========================================================================
  subroutine bcast_agents(nagent, np_target,  ndata, nexp, nobj, f1, f2, cr, &
                          agent,  agent_string, original, original_s )
  ! サブルーチンbcast_agents: MPIを使用してエージェントとオリジナルデータをブロードキャストする
  ! 引数:
  !   nagent       : 整数, エージェントの数
  !   np_target    : 整数, 目的変数の数
  !   ndata        : 整数, データの数
  !   nexp         : 整数, 実験データの数
  !   nobj         : 整数, 目的関数の数
  !   f1, f2       : 実数配列(1:nagent), 各エージェントの突然変異係数
  !   cr           : 実数配列(1:nagent), 各エージェントの交差率
  !   agent        : 実数配列(1:nagent,1:np_target), エージェントのパラメータ
  !   agent_string : 実数配列(1:nagent*np_target), MPI通信用のエージェントパラメータ
  !   original     : 実数配列(1:nexp+nobj,1:ndata), オリジナルデータ
  !   original_s   : 実数配列(1:(nexp+nobj)*ndata), MPI通信用のオリジナルデータ
  ! 動作:
  !   MPIを使用してエージェントとオリジナルデータをブロードキャストする。
  ! ==========================================================================
  integer         :: ierr
  use mpi
  use commonval
  use mpival
  implicit none
  integer         , intent(in)    :: nagent,np_target
  integer         , intent(in)    :: ndata, nexp, nobj
  double precision, intent(inout) :: f1(nagent), f2(nagent), cr(nagent)
  double precision, intent(inout) :: agent(nagent,np_target)
  double precision, intent(inout) :: agent_string(nagent*np_target)
  double precision, intent(inout) :: original(nexp+nobj,ndata)
  double precision, intent(inout) :: original_s((nexp+nobj)*ndata)
  integer                      :: i, idata, ip_target, jdata
  integer                      :: iagent, iagent_s
  ! 
  if( nprocs.gt.1 )then
    jdata = 0
    do idata = 1,ndata
      do i = 1,nexp+nobj
        jdata = jdata+1
        original_s(jdata) = original(i,idata)
      enddo
    enddo
    !
    call MPI_BCAST( original_s(1),   (nexp+nobj)*ndata, MPI_DOUBLE_PRECISION, &
                    0              , MPI_COMM_WORLD  , ierr )
    call MPI_BCAST( agent_string(1), nagent*np_target, MPI_DOUBLE_PRECISION, &
                    0              , MPI_COMM_WORLD  , ierr )
    call MPI_BCAST( f1(1)          , nagent          , MPI_DOUBLE_PRECISION, &
                    0              , MPI_COMM_WORLD  , ierr )
    call MPI_BCAST( f2(1)          , nagent          , MPI_DOUBLE_PRECISION, &
                    0              , MPI_COMM_WORLD  , ierr )
    call MPI_BCAST( cr(1)          , nagent          , MPI_DOUBLE_PRECISION, &
                    0              , MPI_COMM_WORLD  , ierr )
    do iagent=1,nagent
      do ip_target = 1,np_target
        iagent_s = ip_target + (iagent-1)*np_target
        agent(iagent,ip_target) = agent_string(iagent_s) 
      enddo
    enddo
    !
    jdata = 0
    do idata = 1,ndata
      do i = 1,nexp+nobj
        jdata = jdata+1
        original(i,idata) = original_s(jdata)
      enddo
    enddo
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  endif
  !
  return
  end subroutine bcast_agents
  ! ==========================================================================

end module util_de