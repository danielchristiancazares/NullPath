module geokerr_ref_api
  use, intrinsic :: iso_c_binding
  implicit none
contains

  pure real(c_double) function delta_fn(r, a) result(delta)
    use, intrinsic :: iso_c_binding
    implicit none
    real(c_double), value :: r, a
    delta = r*r - 2.0_c_double*r + a*a
  end function delta_fn

  pure real(c_double) function R_of_r(r, a, Lz, Q2) result(Rval)
    use, intrinsic :: iso_c_binding
    implicit none
    real(c_double), value :: r, a, Lz, Q2
    real(c_double) :: Delta, term1, term2
    ! E = 1 for photons
    Delta = delta_fn(r, a)
    term1 = (r*r + a*a) - a*Lz
    term1 = term1*term1
    term2 = Delta * (Q2 + (Lz - a)*(Lz - a))
    Rval = term1 - term2
  end function R_of_r

  pure real(c_double) function Theta_mu(mu, a, Lz, Q2) result(Th)
    use, intrinsic :: iso_c_binding
    implicit none
    real(c_double), value :: mu, a, Lz, Q2
    real(c_double) :: one_minus_mu2
    ! In Mino time: (dmu/dλ)^2 = (1-μ^2)Q - μ^2 Lz^2 (E=1)
    one_minus_mu2 = max(0.0_c_double, 1.0_c_double - mu*mu)
    Th = one_minus_mu2 * Q2 - (mu*mu) * (Lz*Lz)
  end function Theta_mu

  pure real(c_double) function dphi_dlam(r, mu, a, Lz) result(dphi)
    use, intrinsic :: iso_c_binding
    implicit none
    real(c_double), value :: r, mu, a, Lz
    real(c_double) :: Delta, sin2
    ! E = 1, unwrapped phi
    Delta = max(1.0e-30_c_double, delta_fn(r, a))
    sin2 = max(1.0e-30_c_double, 1.0_c_double - mu*mu)
    dphi = a * ( ( (r*r + a*a) - a*Lz ) / Delta ) + Lz / sin2 - a
  end function dphi_dlam

  pure real(c_double) function dt_dlam(r, mu, a, Lz) result(dt)
    use, intrinsic :: iso_c_binding
    implicit none
    real(c_double), value :: r, mu, a, Lz
    real(c_double) :: Delta
    ! E = 1
    Delta = max(1.0e-30_c_double, delta_fn(r, a))
    dt = ((r*r + a*a) * ( (r*r + a*a) - a*Lz ) ) / Delta + a * ( Lz - a * (1.0_c_double - mu*mu) )
  end function dt_dlam

  subroutine geokerr_eval_lambda(u0, mu0, uf, a, Lz, Q2, t0, &
                                 lam0, dlam, S,                &
                                 u, mu, t, phi, affine,        &
                                 tpmi, tpri, info) bind(C,name="geokerr_eval_lambda")
    use, intrinsic :: iso_c_binding
    implicit none
    ! inputs (by value for C ABI)
    real(c_double), value :: u0, mu0, uf, a, Lz, Q2, t0
    real(c_double), value :: lam0, dlam
    integer(c_int),  value :: S
    ! outputs (C passes pointers) — assumed-size for C interoperability
    real(c_double)          :: u(*), mu(*), t(*), phi(*), affine(*)
    integer(c_int)          :: tpmi, tpri
    integer(c_int)          :: info  ! 0 ok; >0 code+reason
    ! Analytic sampler via GEOKERR + resample to uniform lambda grid
    interface
      subroutine GEOKERR(U0,UF,UOUT,MU0,MUF,A,L,Q2,ALPHA,BETA,TPM,TPR,SU,SM,NUP,OFFSET,PHIT,USEGEOR,MUFILL,NCASE,KEXT,NPTS,&
                         UFI,MUFI,DTI,DPHI,TPMIv,TPRIv,LAMBDAI)
        logical PHIT, USEGEOR, MUFILL
        integer NUP,NCASE,KEXT,NPTS,TPM,TPR
        double precision U0,UF,UOUT,MU0,MUF,A,L,Q2,ALPHA,BETA,SU,SM,OFFSET
        double precision UFI(*), MUFI(*), DTI(*), DPHI(*), LAMBDAI(*)
        integer TPMIv(*), TPRIv(*)
      end subroutine GEOKERR
    end interface

    integer :: i, j, NUP, NCASE, KEXT, NPTS
    integer :: TPM, TPR
    real(c_double) :: UOUT, MUF, SU, SM, OFFSET
    real(c_double) :: ALPHA, BETA
    logical :: PHIT, USEGEOR, MUFILL

  real(c_double), allocatable :: UFI(:), MUFI(:), DTI(:), DPHI(:), LAMI(:)
  integer, allocatable :: TPMI_arr(:), TPRI_arr(:)
    real(c_double) :: lamk, lam_base, lam_prev
    integer :: tp_idx_mu, tp_idx_r
  real(c_double) :: x0, x1
  real(c_double) :: w

    info = 0_c_int
    if (S <= 0_c_int) then
      info = 1_c_int
      return
    end if

  NUP = S
    UOUT = u0
    MUF = mu0
    ALPHA = 0.0_c_double
    BETA  = 0.0_c_double
    TPM = 0
    TPR = 1
    SU = merge(1.0_c_double, -1.0_c_double, (uf - u0) >= 0.0_c_double)
    SM = merge(1.0_c_double, -1.0_c_double, mu0 <= 0.0_c_double)
    OFFSET = 0.5_c_double
    PHIT = .true.
    USEGEOR = .false.
    MUFILL = .true.
    NCASE = 0
    KEXT = 0
    NPTS = 0

  allocate(UFI(NUP), MUFI(NUP), DTI(NUP), DPHI(NUP), LAMI(NUP))
  allocate(TPMI_arr(NUP), TPRI_arr(NUP))

  call GEOKERR(u0, uf, UOUT, mu0, MUF, a, Lz, Q2, ALPHA, BETA, TPM, TPR, SU, SM, NUP, OFFSET, &
         PHIT, USEGEOR, MUFILL, NCASE, KEXT, NPTS, UFI, MUFI, DTI, DPHI, TPMI_arr, TPRI_arr, LAMI)

    ! Normalize lambda to start at zero
    lam_base = LAMI(1)
    lam_prev = lam0
    tp_idx_mu = -1
    tp_idx_r  = -1

    ! Resample to uniform lambda grid using linear interpolation
    j = 1
    do i = 1, S
      lamk = lam0 + dlam * real(i-1, kind=c_double)
      affine(i) = lamk

      ! Advance j until LAMI(j) <= lamk <= LAMI(j+1)
      do while (j < NUP .and. (LAMI(j+1)-lam_base) < lamk - 1.0e-30_c_double)
        j = j + 1
      end do

      if (lamk <= (LAMI(1) - lam_base)) then
  u(i) = UFI(1)
  mu(i) = MUFI(1)
        t(i) = t0 + DTI(1)
        phi(i) = DPHI(1)
      else if (lamk >= (LAMI(NUP) - lam_base)) then
        u(i) = UFI(NUP)
        mu(i) = MUFI(NUP)
        t(i) = t0 + DTI(NUP)
        phi(i) = DPHI(NUP)
      else
  ! Linear interpolation between j and j+1
        x0 = LAMI(j)   - lam_base
        x1 = LAMI(j+1) - lam_base
        if (x1 <= x0 + 1.0e-30_c_double) then
          w = 0.0_c_double
        else
          w = (lamk - x0) / (x1 - x0)
        end if
        u(i)  = (1.0_c_double - w) * UFI(j)  + w * UFI(j+1)
        mu(i) = (1.0_c_double - w) * MUFI(j) + w * MUFI(j+1)
        t(i)  = t0 + ((1.0_c_double - w) * DTI(j)  + w * DTI(j+1))
        phi(i)= (1.0_c_double - w) * DPHI(j) + w * DPHI(j+1)
      end if

      ! Force exact seeding at the first uniform sample to the provided inputs
      if (i == 1) then
        u(i)  = u0
        mu(i) = mu0
      end if
    end do

  ! Map first turning point indices to nearest lambda grid indices (0-based)
    tp_idx_mu = -1
    tp_idx_r  = -1
    do j = 1, NUP
  if (tp_idx_mu < 0 .and. TPMI_arr(j) > 0) then
        ! Find nearest i where affine(i) >= (LAMI(j)-lam_base)
        do i = 1, S
          if (affine(i) >= (LAMI(j) - lam_base)) then
            tp_idx_mu = i - 1
            exit
          end if
        end do
      end if
  if (tp_idx_r < 0 .and. TPRI_arr(j) > 0) then
        do i = 1, S
          if (affine(i) >= (LAMI(j) - lam_base)) then
            tp_idx_r = i - 1
            exit
          end if
        end do
      end if
      if (tp_idx_mu >= 0 .and. tp_idx_r >= 0) exit
    end do

  tpmi = tp_idx_mu
  tpri = tp_idx_r

  deallocate(UFI, MUFI, DTI, DPHI, LAMI, TPMI_arr, TPRI_arr)
    info = 0_c_int
  end subroutine geokerr_eval_lambda

  subroutine geokerr_estimate_periods(u0, mu0, a, Lz, Q2, Tr, Ttheta, info) bind(C,name="geokerr_estimate_periods")
    use, intrinsic :: iso_c_binding
    implicit none
    real(c_double), value :: u0, mu0, a, Lz, Q2
    real(c_double) :: Tr, Ttheta
    integer(c_int) :: info

    ! Analytic estimation via GEOKERR sampling: measure lambda between successive turning points
    interface
      subroutine GEOKERR(U0,UF,UOUT,MU0,MUF,A,L,Q2,ALPHA,BETA,TPM,TPR,SU,SM,NUP,OFFSET,PHIT,USEGEOR,MUFILL,NCASE,KEXT,NPTS,&
                         UFI,MUFI,DTI,DPHI,TPMIv,TPRIv,LAMBDAI)
        logical PHIT, USEGEOR, MUFILL
        integer NUP,NCASE,KEXT,NPTS,TPM,TPR
        double precision U0,UF,UOUT,MU0,MUF,A,L,Q2,ALPHA,BETA,SU,SM,OFFSET
        double precision UFI(*), MUFI(*), DTI(*), DPHI(*), LAMBDAI(*)
        integer TPMIv(*), TPRIv(*)
      end subroutine GEOKERR
    end interface

    integer :: NUP, NCASE, KEXT, NPTS
    integer :: TPM, TPR
    real(c_double) :: UOUT, MUF, SU, SM, OFFSET
    real(c_double) :: ALPHA, BETA
    logical :: PHIT, USEGEOR, MUFILL
    real(c_double), allocatable :: UFI(:), MUFI(:), DTI(:), DPHI(:), LAMI(:)
    integer, allocatable :: TPMI_arr(:), TPRI_arr(:)
    real(c_double) :: lam_base
    integer :: j
    integer :: prev_mu, prev_r
    integer :: idx_mu1, idx_mu2, idx_r1, idx_r2

    ! Choose a reasonably fine sampling to reliably hit both turning points
    NUP   = 4096
    UOUT  = u0
    MUF   = mu0
    ALPHA = 0.0_c_double
    BETA  = 0.0_c_double
    TPM   = 0
    TPR   = 1
    SU    = 1.0_c_double
    SM    = merge(1.0_c_double, -1.0_c_double, mu0 <= 0.0_c_double)
    OFFSET= 0.5_c_double
    PHIT  = .true.
    USEGEOR = .false.
    MUFILL  = .true.
    NCASE = 0
    KEXT  = 0
    NPTS  = 0

    allocate(UFI(NUP), MUFI(NUP), DTI(NUP), DPHI(NUP), LAMI(NUP))
    allocate(TPMI_arr(NUP), TPRI_arr(NUP))

    call GEOKERR(u0, u0, UOUT, mu0, MUF, a, Lz, Q2, ALPHA, BETA, TPM, TPR, SU, SM, NUP, OFFSET, &
                 PHIT, USEGEOR, MUFILL, NCASE, KEXT, NPTS, UFI, MUFI, DTI, DPHI, TPMI_arr, TPRI_arr, LAMI)

    lam_base = LAMI(1)
    idx_mu1 = -1; idx_mu2 = -1
    idx_r1  = -1; idx_r2  = -1
    prev_mu = 0
    prev_r  = 0
    do j = 1, NUP
      if (idx_mu1 < 0 .and. prev_mu < 1 .and. TPMI_arr(j) >= 1) then
        idx_mu1 = j
      end if
      if (idx_mu2 < 0 .and. prev_mu < 2 .and. TPMI_arr(j) >= 2) then
        idx_mu2 = j
      end if
      if (idx_r1 < 0 .and. prev_r < 1 .and. TPRI_arr(j) >= 1) then
        idx_r1 = j
      end if
      if (idx_r2 < 0 .and. prev_r < 2 .and. TPRI_arr(j) >= 2) then
        idx_r2 = j
      end if
      prev_mu = max(prev_mu, TPMI_arr(j))
      prev_r  = max(prev_r,  TPRI_arr(j))
      if (idx_mu2 >= 0 .and. idx_r2 >= 0) exit
    end do

    if (idx_r1 >= 0 .and. idx_r2 >= 0) then
      Tr = 2.0_c_double * max(1.0e-18_c_double, (LAMI(idx_r2)-lam_base) - (LAMI(idx_r1)-lam_base))
    else
      Tr = -1.0_c_double
    end if

    if (idx_mu1 >= 0 .and. idx_mu2 >= 0) then
      Ttheta = 2.0_c_double * max(1.0e-18_c_double, (LAMI(idx_mu2)-lam_base) - (LAMI(idx_mu1)-lam_base))
    else
      Ttheta = -1.0_c_double
    end if

    deallocate(UFI, MUFI, DTI, DPHI, LAMI, TPMI_arr, TPRI_arr)
    info = 0_c_int
  end subroutine geokerr_estimate_periods

  real(c_double) function geokerr_phi_omega(u0, mu0, a, Lz, Q2) bind(C,name="geokerr_phi_omega")
    use, intrinsic :: iso_c_binding
    implicit none
    real(c_double), value :: u0, mu0, a, Lz, Q2
    real(c_double) :: r
    r = max(1.0e-12_c_double, 1.0_c_double / max(1.0e-12_c_double, u0))
    geokerr_phi_omega = dphi_dlam(r, mu0, a, Lz)
  end function geokerr_phi_omega

end module geokerr_ref_api
