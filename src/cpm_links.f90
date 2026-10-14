! link functions and shared numerics for every routine in the package

module cpm_links
use, intrinsic :: ISO_FORTRAN_ENV, only: dp => real64, int32
implicit none

real(dp), parameter :: INV_SQRT2   = 0.70710678118654752440_dp
real(dp), parameter :: INV_SQRT2PI = 0.39894228040143267794_dp
real(dp), parameter :: INV_PI      = 0.31830988618379067154_dp
real(dp), parameter :: PIVOT_EPS   = 1.0e-14_dp
integer(int32), parameter :: MAX_HALVING = 30

contains

! interval probability of each observed category, with the link density and its
! derivative at both ends of the interval
subroutine compute_probs(n, k, alpha, lp, y, link, d, pdf1, pdf2, dpdf1, dpdf2)
    integer(int32), intent(in) :: n, k, link, y(n)
    real(dp), intent(in)  :: alpha(k), lp(n)
    real(dp), intent(out) :: d(n), pdf1(n), pdf2(n), dpdf1(n), dpdf2(n)
    integer(int32) :: i, j
    real(dp) :: eta1, eta2, cdf1, cdf2, e1, e2, t1, t2

    select case (link)

    case (1)  ! logistic
        do i = 1, n
            j = y(i)
            if (j == 0) then
                eta1 = alpha(1) + lp(i)
                e1 = exp(-eta1); cdf1 = 1.0_dp / (1.0_dp + e1)
                d(i) = 1.0_dp - cdf1; pdf1(i) = cdf1 * d(i)
                dpdf1(i) = pdf1(i) * (1.0_dp - 2.0_dp * cdf1)
                pdf2(i) = 0.0_dp; dpdf2(i) = 0.0_dp
            else if (j == k) then
                eta1 = alpha(k) + lp(i)
                e1 = exp(-eta1); cdf1 = 1.0_dp / (1.0_dp + e1)
                d(i) = cdf1; pdf1(i) = cdf1 * (1.0_dp - cdf1)
                dpdf1(i) = pdf1(i) * (1.0_dp - 2.0_dp * cdf1)
                pdf2(i) = 0.0_dp; dpdf2(i) = 0.0_dp
            else
                eta1 = alpha(j) + lp(i); eta2 = alpha(j+1) + lp(i)
                e1 = exp(-eta1); e2 = exp(-eta2)
                cdf1 = 1.0_dp / (1.0_dp + e1); cdf2 = 1.0_dp / (1.0_dp + e2)
                d(i) = cdf1 - cdf2
                pdf1(i) = cdf1 * (1.0_dp - cdf1); pdf2(i) = cdf2 * (1.0_dp - cdf2)
                dpdf1(i) = pdf1(i) * (1.0_dp - 2.0_dp * cdf1)
                dpdf2(i) = pdf2(i) * (1.0_dp - 2.0_dp * cdf2)
            end if
        end do

    case (2)  ! probit
        do i = 1, n
            j = y(i)
            if (j == 0) then
                eta1 = alpha(1) + lp(i)
                cdf1 = 0.5_dp * (1.0_dp + erf(eta1 * INV_SQRT2))
                d(i) = 1.0_dp - cdf1
                pdf1(i) = INV_SQRT2PI * exp(-0.5_dp * eta1 * eta1)
                dpdf1(i) = -pdf1(i) * eta1
                pdf2(i) = 0.0_dp; dpdf2(i) = 0.0_dp
            else if (j == k) then
                eta1 = alpha(k) + lp(i)
                cdf1 = 0.5_dp * (1.0_dp + erf(eta1 * INV_SQRT2))
                d(i) = cdf1
                pdf1(i) = INV_SQRT2PI * exp(-0.5_dp * eta1 * eta1)
                dpdf1(i) = -pdf1(i) * eta1
                pdf2(i) = 0.0_dp; dpdf2(i) = 0.0_dp
            else
                eta1 = alpha(j) + lp(i); eta2 = alpha(j+1) + lp(i)
                cdf1 = 0.5_dp * (1.0_dp + erf(eta1 * INV_SQRT2))
                cdf2 = 0.5_dp * (1.0_dp + erf(eta2 * INV_SQRT2))
                d(i) = cdf1 - cdf2
                pdf1(i) = INV_SQRT2PI * exp(-0.5_dp * eta1 * eta1)
                pdf2(i) = INV_SQRT2PI * exp(-0.5_dp * eta2 * eta2)
                dpdf1(i) = -pdf1(i) * eta1; dpdf2(i) = -pdf2(i) * eta2
            end if
        end do

    case (3)  ! log-log
        do i = 1, n
            j = y(i)
            if (j == 0) then
                eta1 = alpha(1) + lp(i)
                e1 = exp(-eta1); cdf1 = exp(-e1)
                d(i) = 1.0_dp - cdf1; pdf1(i) = e1 * cdf1
                dpdf1(i) = pdf1(i) * (-1.0_dp + e1)
                pdf2(i) = 0.0_dp; dpdf2(i) = 0.0_dp
            else if (j == k) then
                eta1 = alpha(k) + lp(i)
                e1 = exp(-eta1); cdf1 = exp(-e1)
                d(i) = cdf1; pdf1(i) = e1 * cdf1
                dpdf1(i) = pdf1(i) * (-1.0_dp + e1)
                pdf2(i) = 0.0_dp; dpdf2(i) = 0.0_dp
            else
                eta1 = alpha(j) + lp(i); eta2 = alpha(j+1) + lp(i)
                e1 = exp(-eta1); e2 = exp(-eta2)
                cdf1 = exp(-e1); cdf2 = exp(-e2)
                d(i) = cdf1 - cdf2
                pdf1(i) = e1 * cdf1; pdf2(i) = e2 * cdf2
                dpdf1(i) = pdf1(i) * (-1.0_dp + e1); dpdf2(i) = pdf2(i) * (-1.0_dp + e2)
            end if
        end do

    case (4)  ! complementary log-log
        do i = 1, n
            j = y(i)
            if (j == 0) then
                eta1 = alpha(1) + lp(i)
                e1 = exp(eta1); cdf1 = 1.0_dp - exp(-e1)
                d(i) = 1.0_dp - cdf1; pdf1(i) = e1 * exp(-e1)
                dpdf1(i) = pdf1(i) * (1.0_dp - e1)
                pdf2(i) = 0.0_dp; dpdf2(i) = 0.0_dp
            else if (j == k) then
                eta1 = alpha(k) + lp(i)
                e1 = exp(eta1); cdf1 = 1.0_dp - exp(-e1)
                d(i) = cdf1; pdf1(i) = e1 * exp(-e1)
                dpdf1(i) = pdf1(i) * (1.0_dp - e1)
                pdf2(i) = 0.0_dp; dpdf2(i) = 0.0_dp
            else
                eta1 = alpha(j) + lp(i); eta2 = alpha(j+1) + lp(i)
                e1 = exp(eta1); e2 = exp(eta2)
                cdf1 = 1.0_dp - exp(-e1); cdf2 = 1.0_dp - exp(-e2)
                d(i) = cdf1 - cdf2
                pdf1(i) = e1 * exp(-e1); pdf2(i) = e2 * exp(-e2)
                dpdf1(i) = pdf1(i) * (1.0_dp - e1); dpdf2(i) = pdf2(i) * (1.0_dp - e2)
            end if
        end do

    case (5)  ! cauchit
        do i = 1, n
            j = y(i)
            if (j == 0) then
                eta1 = alpha(1) + lp(i)
                t1 = 1.0_dp + eta1 * eta1; cdf1 = INV_PI * atan(eta1) + 0.5_dp
                d(i) = 1.0_dp - cdf1; pdf1(i) = INV_PI / t1
                dpdf1(i) = -2.0_dp * eta1 * INV_PI / (t1 * t1)
                pdf2(i) = 0.0_dp; dpdf2(i) = 0.0_dp
            else if (j == k) then
                eta1 = alpha(k) + lp(i)
                t1 = 1.0_dp + eta1 * eta1; cdf1 = INV_PI * atan(eta1) + 0.5_dp
                d(i) = cdf1; pdf1(i) = INV_PI / t1
                dpdf1(i) = -2.0_dp * eta1 * INV_PI / (t1 * t1)
                pdf2(i) = 0.0_dp; dpdf2(i) = 0.0_dp
            else
                eta1 = alpha(j) + lp(i); eta2 = alpha(j+1) + lp(i)
                t1 = 1.0_dp + eta1 * eta1; t2 = 1.0_dp + eta2 * eta2
                cdf1 = INV_PI * atan(eta1) + 0.5_dp; cdf2 = INV_PI * atan(eta2) + 0.5_dp
                d(i) = cdf1 - cdf2
                pdf1(i) = INV_PI / t1; pdf2(i) = INV_PI / t2
                dpdf1(i) = -2.0_dp * eta1 * INV_PI / (t1 * t1)
                dpdf2(i) = -2.0_dp * eta2 * INV_PI / (t2 * t2)
            end if
        end do

    end select
end subroutine compute_probs

! first and second derivatives of the negative log-likelihood in the linear predictor
subroutine score_weight(n, k, y, d, pdf1, pdf2, dpdf1, dpdf2, s, w)
    integer(int32), intent(in) :: n, k, y(n)
    real(dp), intent(in)  :: d(n), pdf1(n), pdf2(n), dpdf1(n), dpdf2(n)
    real(dp), intent(out) :: s(n), w(n)
    integer(int32) :: i, j
    real(dp) :: dd
    do i = 1, n
        j = y(i); dd = d(i)
        if (j == 0) then
            s(i) = pdf1(i) / dd
            w(i) = (pdf1(i) * pdf1(i) + dd * dpdf1(i)) / (dd * dd)
        else if (j == k) then
            s(i) = -pdf1(i) / dd
            w(i) = (pdf1(i) * pdf1(i) - dd * dpdf1(i)) / (dd * dd)
        else
            s(i) = -(pdf1(i) - pdf2(i)) / dd
            w(i) = ((pdf1(i) - pdf2(i)) ** 2 - dd * (dpdf1(i) - dpdf2(i))) / (dd * dd)
        end if
    end do
end subroutine score_weight

subroutine solve_tridiag(nn, diag, offdiag, rhs, solution, status)
    integer(int32), intent(in)  :: nn
    real(dp), intent(in)  :: diag(nn), offdiag(max(nn-1,1)), rhs(nn)
    real(dp), intent(out) :: solution(nn)
    integer(int32), intent(out) :: status
    real(dp) :: wd(nn), wr(nn), m
    integer(int32) :: i
    status = 0; wd(1) = diag(1); wr(1) = rhs(1)
    if (abs(wd(1)) < PIVOT_EPS) then; status = 1; return; end if
    do i = 2, nn
        m = offdiag(i-1) / wd(i-1)
        wd(i) = diag(i) - m * offdiag(i-1)
        if (abs(wd(i)) < PIVOT_EPS) then; status = i; return; end if
        wr(i) = rhs(i) - m * wr(i-1)
    end do
    solution(nn) = wr(nn) / wd(nn)
    do i = nn-1, 1, -1
        solution(i) = (wr(i) - offdiag(i) * solution(i+1)) / wd(i)
    end do
end subroutine solve_tridiag

! P(Y <= y_j | x) = 1 - F(alpha_j + x'beta), the model in the rms orientation
pure function link_cdf(eta, link) result(c)
    real(dp), intent(in) :: eta
    integer(int32), intent(in) :: link
    real(dp) :: c
    select case (link)
    case (1)
        c = 1.0_dp - 1.0_dp / (1.0_dp + exp(-eta))
    case (2)
        c = 1.0_dp - 0.5_dp * (1.0_dp + erf(eta * INV_SQRT2))
    case (3)
        c = 1.0_dp - exp(-exp(-eta))
    case (4)
        c = exp(-exp(eta))
    case (5)
        c = 0.5_dp - atan(eta) * INV_PI
    case default
        c = 0.0_dp
    end select
end function link_cdf

end module cpm_links
