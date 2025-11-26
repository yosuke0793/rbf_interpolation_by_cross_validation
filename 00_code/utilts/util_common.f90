module util_common
implicit none
contains

  ! ==========================================================================
  subroutine outlog( iopt, nagent, np_target, f1, f2, cr, agent, err_agent )
  ! サブルーチンoutlog: エージェントの情報をログファイルに出力する
  ! 
  !　引数:
  !   iopt        : 整数, 現在の最適化ステップ
  !   nagent      : 整数, エージェントの総数
  !   np_target   : 整数, 最適化対象パラメータ数
  !   f1          : 実数配列(nagent), 各エージェントの突然変異係数1
  !   f2          : 実数配列(nagent), 各エージェントの突然変異係数2
  !   cr          : 実数配列(nagent), 各エージェントの交差率
  !   agent       : 実数2次元配列(nagent,np_target), 各エージェントのパラメータ値
  !   err_agent   : 実数配列(nagent), 各エージェントの誤差値
  ! 使用モジュール:
  !   fileinfo    : ファイル入出力に関する情報を管理
  !   mpival      : MPIに関する情報を管理
  ! 動作:
  !   エージェントごとの情報をログファイルに出力する。MPIのランク0プロセスのみが出力を担当する。
  ! ==========================================================================
  use mpi
  use fileinfo, only: out_log
  use mpival, only: irank
  implicit none
  integer         , intent(in) :: iopt, nagent, np_target
  double precision, intent(in) :: f1(nagent), f2(nagent), cr(nagent)
  double precision, intent(in) :: agent(nagent,np_target), err_agent(nagent)
  integer                      :: ip_target, iagent
  !
  if( irank.ne.0 )then
    return
  else
    do iagent=1,nagent
      write(out_log,'(2i10,99es15.5)')iopt, iagent, err_agent(iagent),&
                                      f1(iagent), f2(iagent), cr(iagent), &
                                      ( agent(iagent,ip_target), ip_target=1,np_target )
    enddo
  endif
  !
  return
  end subroutine outlog
  ! ==========================================================================

  ! ==========================================================================
  subroutine read_data(ndata, nexp, nobj, original, original_s)
  ! サブルーチンread_data: データセットファイルから実験データを読み込む
  !　引数:
  !   ndata       : 整数, データ点の総数
  !   nexp        : 整数, 説明変数の数
  !   nobj        : 整数, 目的変数の数
  !   original    : 実数2次元配列(nexp+nobj,ndata), 読み込んだ実験データを格納
  !   original_s  : 実数1次元配列((nexp+nobj)*ndata), 読み込んだ実験データを1次元配列に変換して格納
  ! 使用モジュール:
  !   fileinfo    : ファイル入出力に関する情報を管理
  ! 動作:
  !   指定されたデータセットファイルから実験データを読み込み、2次元配列と1次元配列の両方に格納する。
  ! ==========================================================================
  use fileinfo, only: inp_origin
  implicit none
  integer, intent(in) :: ndata, nexp, nobj
  double precision, intent(out) :: original(nexp+nobj,ndata)
  double precision, intent(out) :: original_s((nexp+nobj)*ndata)
  integer :: i, idata, jdata
  original = 0.d0
  open(inp_origin,file='./inputs/dataset.txt',status='old')
  do idata = 1,ndata
     read(inp_origin,*)( original(i,idata), i=1,nexp+nobj )
  enddo
  close(inp_origin)
  !
  jdata = 0
  do idata = 1,ndata
    do i = 1,nexp+nobj
      jdata = jdata+1
      original_s(jdata) = original(i,idata)
    enddo
  enddo
  return
  end subroutine read_data
  ! ==========================================================================

  ! ==========================================================================
  subroutine heapsort2(n,array,turn)
  ! サブルーチンheapsort2: ヒープソートアルゴリズムを使用して配列を降順にソートし、元のインデックスを追跡する
  !　引数:
  !   n           : 整数, ソートする要素の数
  !   array       : 実数配列(1:n), ソート対象の配列
  !   turn        : 整数配列(1:n), 元のインデックスを格納する配列
  ! 動作:
  !   ヒープソートアルゴリズムを使用して、指定された配列を降順にソートし、元のインデックスを追跡する。
  ! ==========================================================================
  implicit none
  integer,intent(in)::n
  integer,intent(out)::turn(1:n)
  double precision,intent(inout)::array(1:n)
  integer::i,k,j,l,m
  double precision::t
  if(n.le.0)then
     write(6,*)"Error, at heapsort"; stop
  endif
  if(n.eq.1)return
  do i=1,n
     turn(i)=i
  enddo
  l=n/2+1
  k=n
  do while(k.ne.1)
     if(l.gt.1)then
        l=l-1
        t=array(l)
        m=turn(l)
     else
        t=array(k)
        m=turn(k)
        array(k)=array(1)
        turn(k)=turn(1)
        k=k-1
        if(k.eq.1) then
           array(1)=t
           turn(1)=m
           exit
        endif
     endif
     i=l
     j=l+l
     do while(j.le.k)
        if(j.lt.k)then
           if(array(j).lt.array(j+1))j=j+1
        endif
        if (t.lt.array(j))then
           array(i)=array(j)
           turn(i)=turn(j)
           i=j
           j=j+j
        else
           j=k+1
        endif
     enddo
     array(i)=t
     turn(i)=m
  enddo
  return
  end subroutine heapsort2
  ! ==========================================================================

  ! ==========================================================================
  subroutine readin( flag_strategy, nopt,  nagent, ncv, ndata, nexp, nobj, &
                     npara, np_target, np_fix, mf, sf, mcr, scr, tole)
  ! サブルーチンreadin: 条件ファイルから最適化の設定を読み込む
  !　引数:
  !   flag_strategy : 整数, 最適化戦略のフラグ
  !   nopt          : 整数, 最適化ステップ数
  !   nagent        : 整数, エージェントの総数
  !   ncv           : 整数, クロスバリデーションの分割数
  !   ndata         : 整数, データ点の総数
  !   nexp          : 整数, 説明変数の数
  !   nobj          : 整数, 目的変数の数
  !   npara         : 整数, 最適化パラメータの総数
  !   np_target     : 整数, 最適化対象パラメータ数
  !   np_fix        : 整数, 固定パラメータ数
  !   mf            : 実数, 突然変異係数1の平均値
  !   sf            : 実数, 突然変異係数1の標準偏差
  !   mcr           : 実数, 交差率の平均値
  !   scr           : 実数, 交差率の標準偏差
  !   tole          : 実数, 収束判定の許容誤差
  ! 使用モジュール:
  !   commonval    : 共通変数を管理
  !   fileinfo     : ファイル入出力に関する情報を管理
  !   mpival       : MPIに関する情報を管理
  ! 動作:
  !   条件ファイルから最適化の各種設定を読み込み、対応する変数に格納する。
  ! ==========================================================================
  use commonval
  use fileinfo, only: inp_condt
  use mpival, only: nprocs
  ! --- Arguments
  implicit none
  integer, intent(out) :: flag_strategy
  integer, intent(out) :: nopt,  nagent, ncv
  integer, intent(out) :: ndata, nexp, nobj
  integer, intent(out) :: npara, np_target, np_fix
  double precision, intent(out) :: mf, sf, mcr, scr, tole
  ! --- Local variables
  integer          :: i, flag_fixed, nagent_half
  integer          :: ipara, ipara_fix, ipara_target, intbuf(2)
  double precision :: buf(2)
  double precision, parameter :: pi=4.0*atan(1.d0)
  character(len=7) :: tag
  ! -------------------------------------------------------------------------
  ! --- OPEN CONDITION FILE
  open(inp_condt,file='./inputs/condition_training.txt',action='read',status='old')
  ! --- READ SEVERAL CONDITIONS FROM de_condt.txt
  do 
    read(inp_condt,'(a7)')tag
    select case(tag)
    case('/CONDT/')
      read(inp_condt,*)flag_strategy, tole
      read(inp_condt,*)nopt,          nagent_half, ncv
      read(inp_condt,*)mf,            sf,          mcr,         scr
      nagent = 2 * nagent_half
    case('/NDATA/')
      read(inp_condt,*)nexp , nobj
      read(inp_condt,*)ndata
    case('/PARAM/')
      read(inp_condt,*)npara, np_fix
      np_target = npara - np_fix
      allocate( parafix(2,np_fix)      )
      allocate( flag_expnd(np_target), flag_intgr(np_target) )
      allocate( maxpara(2,np_target),  minpara(2,np_target)  )
      ipara=0; ipara_fix=0; ipara_target=0
      do i = 1, npara
        ipara = ipara + 1
        read(inp_condt,*) flag_fixed, intbuf(1), intbuf(2), buf(1), buf(2)
        select case(flag_fixed)
        case(1)
          ipara_fix = ipara_fix + 1
          parafix(1,ipara_fix)  = ipara
          parafix(2,ipara_fix)  = buf(1)
        case(0)
          ipara_target = ipara_target + 1
          maxpara(1,ipara_target)  = ipara
          maxpara(2,ipara_target)  = buf(1)
          minpara(1,ipara_target)  = ipara
          minpara(2,ipara_target)  = buf(2)
          flag_expnd(ipara_target) = intbuf(1)
          flag_intgr(ipara_target) = intbuf(2)
          if( maxpara(2,ipara_target).lt.minpara(2,ipara_target) )then
            write(*,*)'ERROR! Maximum value is set less than minimumvalue!',ipara
            call stop00
          endif
          if( flag_intgr(ipara_target).eq.1 )then
            maxpara(2,ipara_target) = buf(1) + cos(pi/1.8D+03)
          endif
        case default
          write(*,*)'ERROR! flag_fixed have to be 0 or 1!'
          call stop00
        end select
      enddo
    case('/ENDOF/')
      exit
    case default
      write(*,*)'ERROR! Unknown tag',tag,' was found in de_condt.txt!'    
      call stop00
    end select
  enddo
  close(inp_condt)
  !
  if( mod(nagent,nprocs).ne.0 )then
    write(*,*)'NAGENT cannot be devided by NPROCS!'
    call stop00  
  endif
  !
  return
  end subroutine readin
  ! ==========================================================================

  ! ==========================================================================
  subroutine setfile
  ! サブルーチンsetfile: 出力ファイルを設定し、必要に応じて初期化する
  ! 使用モジュール:
  !   fileinfo    : ファイル入出力に関する情報を管理
  !   mpival      : MPIに関する情報を管理
  ! 動作:
  !   出力ファイルのユニット番号を設定し、ランク0プロセスがログファイルと履歴ファイルを初期化する。
  ! ===========================================================================
  use fileinfo
  use mpival, only: irank
  implicit none
  ! --- Input file setting 
  inp_condt  = 11
  inp_origin = 12
  inp_pod    = 13
  inp_expln  = 14
  ! --- Output file settting
  out_reprod   = 21
  out_hist_err = 22
  out_best     = 23
  out_log      = 24
  out_rbf      = 25
  out_stplog   = 26
  out_first    = 27
  out_trnval   = 28
  intermid     = 30
  random_table = 31
  if( irank.eq.0 )then
    ! --- Set files
    open(out_log,      file='./results/log.txt',     status='replace')
    open(out_stplog,   file='./results/stplog.txt',  status='replace')
    open(out_hist_err, file='./results/history.txt', status='replace')
    ! --- Write header line
    write(out_hist_err,'(2a10,5a15)')'Step','ID_agent','Loss','Diff.Error','Parameters'
  endif
  open(out_trnval,  file='./results/trnval.txt',  status='replace')
  return
  end subroutine setfile
  ! ==========================================================================

  ! ==========================================================================
  subroutine stop00
  ! サブルーチンstop00: MPIを終了し、プログラムを停止する
  ! 使用モジュール:
  !   mpival      : MPIに関する情報を管理
  !   mpi         : MPIライブラリ
  ! 動作:
  !   MPIを終了し、プログラムの実行を停止する。
  ! ==========================================================================
  use mpival, only: ierr
  use mpi
  implicit none
  call MPI_FINALIZE(ierr)
  stop
  end subroutine stop00
  ! ==========================================================================

  ! ==========================================================================
  subroutine stplog(flag_stplog, iopt,   err_gbest, dev_error,  dev_percent, &
                    count1,      count2, count3,    count_rate, time         )
  ! サブルーチンstplog: 最適化ステップごとのログを出力する
  !　引数:
  !   flag_stplog : 整数, ログ出力のフラグ
  !   iopt         : 整数, 現在の最適化ステップ
  !   err_gbest    : 実数, 現 在の最良誤差値
  !   dev_error    : 実数, 現在の誤差の偏差
  !   dev_percent  : 実数, 現在の誤差の偏差率
  !   count1      : 整数, システムクロックの開始カウント
  !   count2      : 整数, システムクロックの終了カウント
  !   count3      : 整数, システムクロックの最大カウント
  !   count_rate  : 整数, システムクロックのレート
  !   time         : 実数, 経過時間（秒）
  ! 使用モジュール:
  !   mpi         : MPIライブラリ
  !   fileinfo    : ファイル入出力に関する情報を管理
  !   mpival      : MPIに関する情報を管理
  ! 動作:
  !   最適化ステップごとにログを出力する。MPIのランク0プロセスのみが出力を担当する。
  ! ==========================================================================
  use mpi
  use fileinfo, only: out_stplog
  use mpival, only: irank
  implicit none
  integer         , intent(in)    :: flag_stplog, iopt
  integer         , intent(inout) :: count1     , count2    , count3, count_rate
  double precision, intent(inout) :: time
  double precision, intent(in)    :: err_gbest  , dev_error , dev_percent
  integer :: count_fin, hour
  double precision :: dtime
  !
  if( irank.ne.0 )then
    return
  else
    select case(flag_stplog)
    case(0)
      write(out_stplog,'(a7,4a15,3a5)') &
                        'Step','Loss','Dev.Err','Dev.Para','Time [s]','H','M','S'
      time = 0.d0
    case(1)
      call system_clock(count1, count_rate, count3)
    case(2)
      call system_clock(count2)
      if( count2.lt.count1 )then
        count2 = count2 + count3
      endif
      dtime = dble( count2 - count1 )/dble(count_rate)
      time  = time + dtime
      if( time.ge.3600 )then
        hour = int(time/3600)
      else
        hour = 0
      endif
      write(out_stplog,'(i7,4ES15.5,3i5)')&
                          iopt,err_gbest,dev_error,dev_percent,time,&
                          hour,int(mod(int(time),3600)/60),int(mod(mod(int(time),3600),60))
    case(3)
      call system_clock(count_fin)
      if( count_fin.lt.count2 )then
        count_fin = count_fin + count3
      endif
      dtime = int( count_fin - count2 )/count_rate
      time  = time + dtime
      if( time.ge.3600 )then
        hour = int(time/3600)
      else
        hour = 0
      endif
      write(out_stplog,'(a,i5,a5,i5,a4,i5,a4)')'Computation finished in',       &
                                              hour,' hour', &
                                              int(mod(int(time),3600)/60) ,' min' , &
                                              int(mod(mod(int(time),3600),60)) ,' sec'
    end select
  endif
  return
  end subroutine stplog 
  ! ==========================================================================

end module util_common
