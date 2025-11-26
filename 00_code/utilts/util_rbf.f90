module util_rbf
use util_math
use util_common, only: stop00
implicit none
contains

  ! =========================================================================
  subroutine RadialBasisFunction_select(ftyp,enorm,gaus,rbf)
  ! ==============================================================================
  ! サブルーチン: RBF関数の計算
  ! 
  ! 目的:
  !   指定された関数タイプに応じて、放射基底関数(Radial Basis Function, RBF)
  !   の値を計算します。
  !
  ! 引数:
  !   ftyp  (in)  : 関数タイプ (1=Multiquadric, 
  !                           2=Linear,
  !                           3=Inverse multiquadric,
  !                           4=Negative linear with Macaulay brackets,
  !                           5=Gaussian,
  !                           6=Logarithmic)
  !   enorm (in)  : ユークリッドノルム（距離）
  !   gaus  (in)  : 形状パラメータ（スケーリングファクター）
  !   rbf   (out) : 計算されたRBF値
  !
  ! 関数タイプの詳細:
  !   1. Multiquadric: sqrt((r/ε)^2 + 1)
  !      - 滑らかで全体的な近似に適している
  !   2. Linear: r/ε
  !      - 最もシンプルな線形関数
  !   3. Inverse multiquadric: 1/sqrt((r/ε)^2 + 1)
  !      - コンパクトな影響範囲を持つ
  !   4. Negative linear with Macaulay brackets: 1 - r/ε
  !      - 有限サポートを持つ関数
  !   5. Gaussian: exp(-π*(r/ε)^2)
  !      - 非常に滑らかで、局所的な近似に適している
  !   6. Logarithmic: log( -2/ε*r + exp(1) - 1 )
  !      - 特殊な形状を持つ関数
  !
  ! 注意:
  !   - enorm: データ点間の距離（正の値）
  !   - gaus: 形状パラメータ。大きいほど関数は平坦になる
  ! ==============================================================================
  implicit none
  integer,          intent(in)  :: ftyp
  double precision, intent(in)  :: enorm, gaus
  double precision, intent(out) :: rbf
  double precision              :: coef, func
  double precision, parameter   :: pi = 4.d0*atan(1.d0)
  ! ---
  rbf  = 0.d0
  coef = 0.d0
  func = 0.d0
  !
  select case(ftyp)
  case(1) ! Multiquadric
    rbf = dsqrt( (enorm/gaus)**2.d0 + 1.d0 )
    !
  case(2) ! Linear
    rbf = enorm/gaus
    !
  case(3) ! Inverse multiquadric
    rbf = 1.d0/dsqrt( (enorm/gaus)**2.d0 + 1.d0 )
    !
  case(4) ! Negative linear with Macaulay brackets
    rbf  = 1.d0 - enorm/gaus
    if( rbf.lt.0.d0 )  rbf = 0.d0
    !
  case(5) ! Gaussian
    rbf = dexp( -pi*(enorm/gaus)**2.d0 ) 
    !
  case(6) ! --- Logarithmic
    func = -2.d0/gaus*enorm+dexp(1.d0)-1.d0
    if( func.le.0.d0 ) func = 0.d0
    rbf  = dlog( func+1.d0 ) 
    if( rbf.le.0.d0  ) rbf = 0.d0
    !
  case default
    write(*,*)'RBF TYPE is not set appropriately! @subroutine RadialBasisFunction_select in utilts/utilmt.f90'
    write(*,*)' Current RBF TYPE = ', ftyp
    write(*,*)' Please set 1:Gaussian, 2:Linear, 3:Logarithmic'
    call stop00
  end select
    !
  ! ---
  return
  end subroutine RadialBasisFunction_select
  ! =========================================================================

  ! =========================================================================
  subroutine weight_rbf_ridge( flag_full,  ndata, nexp, nobj, &
                               ftyp, gaus, smth, expln, objct, wght )
  ! ==============================================================================
  ! サブルーチン: リッジ回帰によるRBF重み係数の計算
  ! 
  ! 目的:
  !   Radial Basis Function (RBF) 補間のための重み係数を、リッジ回帰を用いて
  !   計算します。正則化されたカーネル行列を構築し、その逆行列を求めることで
  !   重み係数を決定します。
  !
  ! 引数:
  !   flag_full (in)  : 行列解法の選択フラグ (0=疎行列ソルバー, 1=密行列逆行列)
  !   ndata     (in)  : データ点数
  !   nexp      (in)  : 説明変数の次元数
  !   nobj      (in)  : 目的変数の次元数
  !   ftyp      (in)  : RBF関数タイプ (1-6)
  !   gaus      (in)  : RBF形状パラメータ (ε)
  !   smth      (in)  : リッジ回帰の正則化パラメータ (λ)
  !   expln     (in)  : 説明変数配列 (nexp×ndata)
  !   objct     (in)  : 目的変数配列 (nobj×ndata)
  !   wght      (out) : 計算された重み係数行列 (nobj×ndata)
  !
  ! アルゴリズム:
  !   1. カーネル行列 K の計算: K_ij = RBF(||x_i - x_j||)
  !   2. 正規化: 各行を行和で除算
  !   3. 共分散行列の構築: C = K^T K + λI
  !   4. 重み係数の計算: W = C^(-1) K^T Y
  !      ここで、Y は目的変数行列
  !
  ! 注意:
  !   - flag_full=0 の場合、疎行列形式(CSR)とPardisoソルバーを使用
  !   - flag_full=1 の場合、密行列の直接逆行列計算を使用
  ! ==============================================================================
  use commonval, only : zero 
  implicit none
  integer,          intent(in)    :: flag_full
  integer,          intent(in)    :: ndata, nexp, nobj, ftyp
  double precision, intent(in)    :: gaus, smth
  double precision, intent(in)    :: expln(nexp,ndata), objct(nobj,ndata)
  double precision, intent(out)   :: wght(nobj,ndata)
  integer                         :: id, jd, kd, io
  integer                         :: nindx,       iindx
  integer,          allocatable   :: iptr(:),     indx(:)
  double precision                :: sum_kernel,  enorm
  double precision                :: rbf
  double precision, allocatable   :: xx(:),       bb(:)
  integer                         :: locater(ndata,ndata)
  double precision                :: kernel(ndata,ndata)
  double precision                :: tkernel(ndata,ndata)
  double precision                :: covariant(ndata,ndata), kobj(nobj,ndata)
  double precision, allocatable   :: csr_cov(:)
  integer :: count1, rate_count
  !
  call system_clock(count1, rate_count)
  ! --- Kernel matrix
  kernel = 0.d0
  tkernel = 0.d0
  locater = 0
  do id = 1,ndata
    sum_kernel = 0.d0
    do jd = 1,ndata
      ! --- Euclidean norm
      enorm = 0.d0
      do io = 1,nexp
        enorm = enorm + ( expln(io,id)-expln(io,jd) )**2.d0
      enddo
      enorm = dsqrt( enorm )
      ! --- Components of regularized kernel matrix
      call RadialBasisFunction_select(ftyp,enorm,gaus,rbf)
      kernel(id,jd) = rbf
      sum_kernel = sum_kernel + rbf
      if( rbf.gt.zero )locater(id,jd)=1
    enddo
    kernel(id,:) = kernel(id,:)/sum_kernel
    tkernel(:,id) = kernel(id,:)
  enddo
  !
  ! --- Covariant
  locater = 0
  covariant = 0.d0
  nindx   = 0
  iindx   = 0
  do id = 1,ndata
    covariant(id,id) = 1.d0
  enddo
  CALL DGEMM('N','N',ndata,ndata,ndata,1.d0,tkernel,ndata,kernel,ndata,smth,covariant,ndata)
  do jd=1,ndata
    do id=1,ndata
      if( covariant(id,jd).gt.zero )then
        iindx = iindx+1
        locater(id,jd) = 1
      endif
    enddo
  enddo
  !
  nindx = iindx
  if( flag_full.eq.0 )then
    allocate( iptr(ndata+1), indx(nindx), csr_cov(nindx) )
    indx    = 0
    iptr    = 0
    iptr(1) = 1 ! One-based indexing
    csr_cov = 0.d0
    iindx   = 0
    do id=1,ndata
      do jd=1,ndata
        if( covariant(id,jd).gt.zero )then
          iindx = iindx+1
          csr_cov(iindx) = covariant(id,jd) 
          indx(iindx)    = jd ! One-based indexing
        endif
      enddo
      iptr(id+1)=iindx+1 ! One-based indexing
    enddo
  endif
  !
  ! --- Convert objects
  kobj=0.d0
  CALL DGEMM('N','N',nobj,ndata,ndata,1.d0,objct,nobj,kernel,ndata,0.d0,kobj,nobj)
  !
  ! --- Make weight
  select case(flag_full)
    case(1)
    ! --- Inverse of regularized kernel matrix
    call invmat( ndata,covariant )
    ! --- Weight coefficient matrix
    wght = 0.d0
    do kd=1,ndata
      do id=1,ndata
        do io=1,nobj
          wght(io,id) = wght(io,id) + covariant(id,kd)*kobj(io,kd)
        enddo
      enddo
    enddo
    !
  case(0)
    ! --- Weight coefficient matrix
    allocate( xx(ndata), bb(ndata) )
    do io=1,nobj
      xx(:) = 0.d0
      bb(:) = kobj(io,:)
      call pardisoT(ndata, nindx, csr_cov(1:nindx), indx(1:nindx), iptr, bb, xx)
      wght(io,:) = xx(:)
    enddo
    deallocate( xx, bb )
    !
  case default
    write(*,*)'FlAG_FULL is not set appropriately! @subroutine weight_rbf in utilts/utilmt.f90'
    call stop00
  end select
  !
  return
  end subroutine weight_rbf_ridge
  ! =========================================================================
  
  ! =========================================================================
  subroutine regres_rbf( flag_post, flag_update, nexp, nobj, ntrn, ntst, ftyp, gaus, & 
                         expln_lrn, expln_tst, objct_tst, wght, error)
  ! サブルーチン: RBF回帰による予測
  ! 
  ! 目的:
  !   学習済みのRBFモデル（重み係数）を用いて、テストデータに対する予測を行い、
  !   予測誤差を計算します。
  !
  ! 引数:
  !   flag_post   (in)  : 後処理フラグ (0=通常, 1=後処理モード)
  !   flag_update (in)  : 出力更新フラグ (0=出力なし, 1=結果ファイル出力)
  !   nexp        (in)  : 説明変数の次元数
  !   nobj        (in)  : 目的変数の次元数
  !   ntrn        (in)  : 学習データ点数
  !   ntst        (in)  : テストデータ点数
  !   ftyp        (in)  : RBF関数タイプ (1-6)
  !   gaus        (in)  : RBF形状パラメータ (ε)
  !   expln_lrn   (in)  : 学習データの説明変数配列 (nexp×ntrn)
  !   expln_tst   (in)  : テストデータの説明変数配列 (nexp×ntst)
  !   objct_tst   (in)  : テストデータの目的変数配列 (nobj×ntst)
  !   wght        (in)  : 学習済み重み係数行列 (nobj×ntrn)
  !   error       (out) : 予測誤差（RRMSE: Root Relative Mean Squared Error）
  !
  ! アルゴリズム:
  !   1. 各テスト点に対して、学習データ点とのRBFカーネル値を計算
  !   2. カーネル値を正規化
  !   3. 重み係数との内積により目的変数を予測
  !   4. 真値との誤差を計算し、RRMSEを求める
  !
  ! 出力ファイル:
  !   - results/reprod.txt: 説明変数、真値、予測値
  !   - results/errors.txt: 説明変数、予測誤差
  !
  ! 注意:
  !   - flag_update=1 の場合のみ結果ファイルを出力
  ! ==============================================================================
  use commonval, only : zero 
  use fileinfo,  only : filename
  implicit none
  integer,          intent(in)  :: flag_post, flag_update, nexp, nobj, ntrn, ntst, ftyp
  double precision, intent(in)  :: gaus
  double precision, intent(in)  :: expln_lrn(nexp,ntrn), expln_tst(nexp,ntst)
  double precision, intent(in)  :: objct_tst(nobj,ntst)
  double precision, intent(in)  :: wght(nobj,ntrn)
  double precision, intent(out) :: error 
  integer                       :: it, jl, k, nindx
  double precision              :: enorm, rbf, sum_rbf
  double precision              :: objct_rbf(nobj,ntst), kernel_ij(ntrn)
  double precision              :: sqnorm_rrmse, rrmse(nobj)
  double precision              :: error_step(nobj), error_sum(nobj)
  character(99) :: outfmt
  ! 
  outfmt='(999es15.5)'
  ! --- Set file
  if( flag_update.eq.1 )then
    write(filename,'(a18)')'results/reprod.txt'
    open(91,file=trim(filename),status='replace')
    write(filename,'(a18)')'results/errors.txt'
    open(92,file=trim(filename),status='replace')
  endif
  !
  ! --- Step loop
  error_sum = 0.d0
  objct_rbf = 0.d0
  do it = 1,ntst
    nindx     = 0
    sum_rbf   = 0.d0
    ! --- Kernel values between test and training data
    do jl = 1,ntrn
      enorm = 0.d0
      do k = 1,nexp
        enorm = enorm + ( expln_tst(k,it)-expln_lrn(k,jl) )**2.d0
      enddo
      enorm = dsqrt( enorm )
      call RadialBasisFunction_select(ftyp,enorm,gaus,rbf)
      kernel_ij(jl) = rbf
      sum_rbf = sum_rbf + rbf
      if( kernel_ij(jl).gt.zero )nindx=nindx+1
      ! --- Objective variables
      do k = 1,nobj
        objct_rbf(k,it) = objct_rbf(k,it) + wght(k,jl)*rbf
      enddo
    enddo
    if( sum_rbf.lt.zero )sum_rbf=1.d0
    !
    ! --- Normalization
    objct_rbf(:,it) = objct_rbf(:,it)/sum_rbf
    if( flag_post.eq.1 )then
      do jl=1,ntrn
        kernel_ij(jl)=kernel_ij(jl)/sum_rbf
        if( kernel_ij(jl).lt.zero )then
          kernel_ij(jl)=0.d0
        endif
      enddo
    endif
    !
    ! --- Avoid quite small values
    do k=1,nobj
      if( dabs(objct_rbf(k,it)).lt.zero )objct_rbf(k,it)=0.d0
    enddo
    !
    ! --- Error for each component of objective variables
    do k = 1,nobj
      error_step(k) = objct_rbf(k,it) - objct_tst(k,it)
      error_sum(k) = error_sum(k) + ( error_step(k) )**2.d0
    enddo
    if( flag_update.eq.1 )then
      write(92,trim(outfmt))( expln_tst(k,it), k=1,nexp ), &
                            ( error_step(k),   k=1,nobj )
    endif
    !
    ! --- Avoid overflow
    do k = 1,nobj
      if( dabs(error_sum(k)).gt.1.d0/zero )then
        error = 1.d0/zero
        return
      endif
    enddo
    !
  enddo
  !
  ! --- Output results
  if( flag_update.eq.1 )then
    do it=1,ntst
      write(91,trim(outfmt))( expln_tst(k,it),  k=1,nexp ), &
                                ( objct_tst(k,it),  k=1,nobj ), &
                                ( objct_rbf(k,it),  k=1,nobj )
    enddo
    close(91)
    close(92)
  endif
  !
  ! --- Error in one cv-case
  sqnorm_rrmse = 0.d0
  do k=1,nobj
    rrmse(k)     = dsqrt(error_sum(k)/ntst)
    sqnorm_rrmse = sqnorm_rrmse + rrmse(k)**2.d0
  enddo
  error = dsqrt(sqnorm_rrmse)
  !
  return
  end subroutine regres_rbf
  ! =========================================================================

  ! =========================================================================
  subroutine output_rbf( ntrn, nexp, nobj, ftyp, gaus, smth, expln, wght ) 
  ! サブルーチン: RBFモデル情報の出力
  ! 
  ! 目的:
  !   学習済みのRBFモデル情報（説明変数、重み係数、パラメータ）を
  !   バイナリおよびテキスト形式で出力します。
  ! 引数:
  !   ntrn   (in)  : 学習データ点数
  !   nexp   (in)  : 説明変数の次元数
  !   nobj   (in)  : 目的変数の次元数
  !   ftyp   (in)  : RBF関数タイプ (1-6)
  !   gaus   (in)  : RBF形状パラメータ (ε)
  !   smth   (in)  : リッジ回帰の正則化パラメータ (λ)
  !   expln  (in)  : 学習データの説明変数配列 (nexp×ntrn)
  !   wght   (in)  : 学習済み重み係数行列 (nobj×ntrn)
  ! 出力ファイル:
  !   - results/rbf.bin: バイナリ形式のRBFモデル情報
  !   - results/rbf.txt: テキスト形式のRBFモデル情報
  ! ==============================================================================
  use fileinfo,  only : filename, out_rbf
  implicit none
  integer, intent(in)          :: ntrn, nexp, nobj, ftyp
  double precision, intent(in) :: gaus, smth
  double precision, intent(in) :: expln(nexp,ntrn), wght(nobj,ntrn)
  integer :: k,id
  ! ---
  !
  ! --- 'rbf.bin
  write(filename,'(a15)')'results/rbf.bin'
  open(out_rbf,file=trim(filename),form='unformatted',status='replace')
  write(out_rbf)ftyp, 1
  write(out_rbf)ntrn
  write(out_rbf)nexp, nobj
  write(out_rbf)gaus
  write(out_rbf)smth
  do id = 1,ntrn
    write(out_rbf)( expln(k,id), k=1,nexp )
  enddo
  do k=1,nobj
    write(out_rbf)( wght(k,id), id=1,ntrn )
  enddo
  close(out_rbf)
  !
  ! --- 'rbf.txt'
  write(filename,'(a15)')'results/rbf.txt'
  open(out_rbf,file=trim(filename),status='replace')
  write(out_rbf,'(2i15)')ftyp, 1
  write(out_rbf,'(2i15)')ntrn
  write(out_rbf,'(2i15)')nexp, nobj
  write(out_rbf,'(3ES15.5)')gaus
  write(out_rbf,'(1ES15.5)')smth
  do id = 1,ntrn
    write(out_rbf,'(99ES15.5)')( expln(k,id), k=1,nexp ) 
  enddo
  do k=1,nobj
    write(out_rbf,'(99999ES15.5)')( wght(k,id), id=1,ntrn )
  enddo
  close(out_rbf)
  !
  return
  end subroutine output_rbf
  ! =========================================================================

  ! =========================================================================
  subroutine rbfsurrogate(flag_post,flag_update,ndata,nexp,nobj,npara,ntrn,para,origin,error)
  ! サブルーチン: RBFサロゲートモデルの構築と評価
  ! 
  ! 目的:
  !   Radial Basis Function (RBF) サロゲートモデルを構築し、与えられたデータに対して
  !   予測誤差を評価します。リッジ回帰を用いて重み係数を計算し、RBF補間を通じて予測を行います。
  ! 引数:
  !   flag_post   (in)  : 後処理フラグ (0=通常, 1=後処理モード)
  !   flag_update (in)  : 出力更新フラグ (0=出力なし, 1=結果ファイル出力)
  !   ndata       (in)  : 全データ点数
  !   nexp        (in)  : 説明変数の次元数
  !   nobj        (in)  : 目的変数の次元数
  !   npara       (in)  : パラメータ数
  !   ntrn        (in)  : 学習データ点数
  !   para        (in)  : RBFモデルパラメータ配列 (npara)
  !   origin      (in)  : 元データ配列 (nexp+nobj, ndata)
  !   error       (out) : 予測誤差（RRMSE: Root Relative Mean Squared Error）
  ! ==============================================================================
  implicit none
  integer         , intent(in)  :: flag_post, flag_update
  integer         , intent(in)  :: ndata, nexp, nobj, npara, ntrn
  double precision, intent(in)  :: para(npara)
  double precision, intent(in)  :: origin(nexp+nobj,ndata)
  double precision, intent(out) :: error

  integer :: flag_exit
  integer                       :: ftyp, flag_full
  double precision              :: gaus, smth
  double precision, allocatable :: wght(:,:)
  double precision, allocatable :: expln_ttl(:,:), expln_trn(:,:)
  double precision, allocatable :: objct_trn(:,:), objct_ttl(:,:)

  ! --- Fitting parameters
  ftyp = int(para(1))
  gaus = 10.d0**para(2)
  smth = 10.d0**para(3)
  flag_full = 0
  flag_exit = 0
  !
  ! --- Weight coefficient
  allocate( wght(nobj,ntrn) )
  allocate( expln_trn(nexp,ntrn),  objct_trn(nobj,ntrn)  )
  allocate( expln_ttl(nexp,ndata), objct_ttl(nobj,ndata) )
  expln_trn = 0.d0; objct_trn=0.d0
  expln_ttl = 0.d0; objct_ttl=0.d0
  !
  expln_trn = origin(1:nexp,1:ntrn)
  objct_trn = origin(nexp+1:nexp+nobj,1:ntrn)
  expln_ttl = origin(1:nexp,:)
  objct_ttl = origin(nexp+1:nexp+nobj,:)
  !
  ! --- Weight coefficients
  call weight_rbf_ridge( flag_full, ntrn, nexp, nobj, ftyp, gaus, smth, &
                         expln_trn, objct_trn, wght)
  !
  ! --- Regression by RBF interpolation
  call regres_rbf( flag_post, flag_update, nexp, nobj, ntrn, ndata, ftyp, gaus, &
                  expln_trn, expln_ttl, objct_ttl, wght, error )
  !
  ! --- Output model information for Post-process
  if( flag_post*flag_update.eq.1 )then
    call output_rbf( ntrn, nexp, nobj, ftyp, gaus, smth, expln_trn, wght )
  endif
  !
  deallocate( wght )
  deallocate( expln_trn, objct_trn )
  deallocate( expln_ttl, objct_ttl )
  return
  end subroutine
  ! =========================================================================

  ! =========================================================================
  subroutine cv_rbf(flag_post,flag_update,ncv,ndata,nexp,nobj,npara,para,origin,loss)
  ! サブルーチン: RBFサロゲートモデルの交差検証
  ! 
  ! 目的:
  !   Radial Basis Function (RBF) サロゲートモデルの交差検証を実施し、
  !   モデルの汎化性能を評価します。
  ! 引数:
  !   flag_post   (in)  : 後処理フラグ (0=通常, 1=後処理モード)
  !   flag_update (in)  : 出力更新フラグ (0=出力なし, 1=結果ファイル出力)
  !   ncv         (in)  : 交差検証の分割数
  !   ndata       (in)  : 全データ点数
  !   nexp        (in)  : 説明変数の次元数
  !   nobj        (in)  : 目的変数の次元数
  !   npara       (in)  : パラメータ数
  !   para        (in)  : RBFモデルパラメータ配列 (npara)
  !   origin      (in)  : 元データ配列 (nexp+nobj, ndata)
  !   loss        (out) : 交差検証による損失関数値（RRMSE）
  ! ==============================================================================
  use mpival
  use mpi
  ! === Global variables
  implicit none
  integer, intent(in)           :: flag_post, flag_update
  integer, intent(in)           :: ncv, ndata, nexp, nobj, npara
  double precision, intent(in)  :: para(npara), origin(nexp+nobj,ndata)
  double precision, intent(out) :: loss
  ! === Local variables
  integer :: i,k,id,icv,ibatch
  !
  ! --- Cross-validation
  integer                       :: batchsize, nbatch, ntrn, ntst
  integer                       :: nstart, nend
  double precision              :: error
  double precision, allocatable :: batches(:,:,:), alldata(:,:)
  double precision, allocatable :: error_cvcase(:)
  !==========================================================================
  !
  ! --- Set batch size
  nbatch = ncv
  batchsize = ndata/nbatch
  ntrn = max(1,(nbatch-1))*batchsize
  ntst = batchsize
  !
  ! --- Set explanatory and objective
  if( flag_post.eq.0 )then
    ! --- Allocate arrays for cross-validation
    allocate( batches(nexp+nobj,batchsize,nbatch) )
    allocate( alldata(nexp+nobj,batchsize*nbatch) )
    allocate( error_cvcase(ncv) )
    error_cvcase = 0.d0
    !
    ! --- Make bacthes from original data
    do ibatch = 1,nbatch
      if ( ibatch.lt.nbatch ) then
        nstart = batchsize*( ibatch-1 ) + 1
        nend   = batchsize*ibatch
      else
        nstart = ndata - batchsize + 1
        nend   = ndata
      endif
      k = 0
      do id = nstart,nend
        k = k + 1
        batches(:,k,ibatch) = origin(:,id)
      enddo
    enddo
    !
    ! --- Cross-validation loop
    do icv = 1,ncv
      ! --- Make alldata for cv-case
      id = 0
      do ibatch = 1,nbatch
        if( ibatch.eq.icv )then
          nstart = batchsize*(nbatch - 1) + 1
          nend   = batchsize*nbatch
          alldata(:,nstart:nend) = batches(:,:,ibatch)
        else
          do i = 1,batchsize
            id = id + 1
            alldata(:,id) = batches(:,i,ibatch)
          enddo
        endif
      enddo
      !
      ! --- RBF regression for each cv-case
      call rbfsurrogate( flag_post, flag_update, ndata, nexp, nobj, npara, ntrn, para, alldata, error )
      error_cvcase(icv) = error
      !
    enddo
    !
    ! --- Loss function for training data
    loss = 0.d0
    do icv = 1,ncv
      loss = loss + error_cvcase(icv)**2.d0
    enddo
    loss = dsqrt( loss/ncv )
    !
    ! --- Deallocation
    deallocate( batches )
    deallocate( alldata )
    deallocate( error_cvcase )
    !
  elseif( flag_post.eq.1 )then
    call rbfsurrogate( flag_post, flag_update, ndata, nexp, nobj, npara, ndata, para, origin, loss )
  endif
  !
  return
  end subroutine
  ! =========================================================================

end module util_rbf