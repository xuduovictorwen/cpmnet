subroutine ormll_beta_ccd_enet(n, k, p, x, y, wt, link, alpha, beta, &
                          lambda, alpha_en, logL, u, hb, salloc)

use, intrinsic :: ISO_FORTRAN_ENV, only: dp => real64, int32
implicit none

real(dp), parameter :: PI = 3.14159265358979323846_dp
real(dp), parameter :: INV_SQRT2 = 0.70710678118654752440_dp
real(dp), parameter :: INV_SQRT2PI = 0.39894228040143267794_dp
real(dp), parameter :: INV_PI = 0.31830988618379067154_dp

! lambda_max only

integer(int32), intent(in)     :: n, y(n), k, p, link
real(dp),       intent(in)     :: x(n, p), wt(n), alpha(k), lambda, alpha_en
real(dp),       intent(out)    :: logL, u(p), hb(p)
real(dp),       intent(inout)  :: beta(p)
integer(int32), intent(out)    :: salloc

integer(int32)  :: i, j, l
real(dp)        :: z, s, beta_old, lambda1, lambda2, z_j, inv_n, inv_d, inv_d2
real(dp), allocatable :: lp(:), d(:), pdf1(:), pdf2(:), dpdf1(:), dpdf2(:)

lambda1 = lambda * alpha_en
lambda2 = lambda * (1.0_dp - alpha_en)
inv_n = 1.0_dp / real(n, dp)

allocate(lp(n), d(n), pdf1(n), pdf2(n), dpdf1(n), dpdf2(n), stat=salloc)
if (salloc /= 0) return


lp = matmul(x, beta)


call compute_probs_fast(n, k, alpha, lp, y, link, d, pdf1, pdf2, dpdf1, dpdf2)

! all derivatives before any update
do l = 1, p
    u(l) = 0.0_dp
    hb(l) = 0.0_dp
    
    do i = 1, n
        j = y(i)
        inv_d = wt(i) / d(i)
        inv_d2 = inv_d / d(i)
        
        if (j == 0) then
            u(l) = u(l) + inv_d * pdf1(i) * x(i, l)
            z = pdf1(i) * pdf1(i) + d(i) * dpdf1(i)
        else if (j == k) then
            u(l) = u(l) - inv_d * pdf1(i) * x(i, l)
            z = pdf1(i) * pdf1(i) - d(i) * dpdf1(i)
        else
            u(l) = u(l) - inv_d * x(i, l) * (pdf1(i) - pdf2(i))
            z = (pdf1(i) - pdf2(i)) ** 2 - d(i) * (dpdf1(i) - dpdf2(i))
        end if
        
        hb(l) = hb(l) + inv_d2 * z * x(i, l) * x(i, l)
    end do
    
    u(l) = u(l) * inv_n
    hb(l) = hb(l) * inv_n
end do


do l = 1, p
    beta_old = beta(l)
    
    z_j = hb(l) * beta_old - u(l)
    
    if (abs(z_j) <= lambda1) then
        beta(l) = 0.0_dp
    else
        beta(l) = sign(1.0_dp, z_j) * (abs(z_j) - lambda1) / (hb(l) + lambda2)
    end if
    
    s = beta(l) - beta_old
    if (s /= 0.0_dp) then
        do i = 1, n
            lp(i) = lp(i) + x(i, l) * s
        end do
        call compute_probs_fast(n, k, alpha, lp, y, link, d, pdf1, pdf2, dpdf1, dpdf2)
    end if
end do


logL = 0.0_dp
do i = 1, n
    logL = logL + wt(i) * log(d(i))
end do
logL = -2.0_dp * logL

deallocate(lp, d, pdf1, pdf2, dpdf1, dpdf2)
return

contains

subroutine compute_probs_fast(nn, kk, alpha, lp, yy, link_type, dd, pdf1, pdf2, dpdf1, dpdf2)
    integer(int32), intent(in) :: nn, kk, link_type
    integer(int32), intent(in) :: yy(nn)
    real(dp), intent(in) :: alpha(kk), lp(nn)
    real(dp), intent(out) :: dd(nn), pdf1(nn), pdf2(nn), dpdf1(nn), dpdf2(nn)
    
    integer(int32) :: ii, jj
    real(dp) :: eta1, eta2, cdf1, cdf2, e1, e2, t1, t2
    
    select case(link_type)
    
    case(1)  ! Logistic
        do ii = 1, nn
            jj = yy(ii)
            if (jj == 0) then
                eta1 = alpha(1) + lp(ii)
                e1 = exp(-eta1)
                cdf1 = 1.0_dp / (1.0_dp + e1)
                dd(ii) = 1.0_dp - cdf1
                pdf1(ii) = cdf1 * dd(ii)
                dpdf1(ii) = pdf1(ii) * (1.0_dp - 2.0_dp * cdf1)
                pdf2(ii) = 0.0_dp
                dpdf2(ii) = 0.0_dp
            else if (jj == kk) then
                eta1 = alpha(kk) + lp(ii)
                e1 = exp(-eta1)
                cdf1 = 1.0_dp / (1.0_dp + e1)
                dd(ii) = cdf1
                pdf1(ii) = cdf1 * (1.0_dp - cdf1)
                dpdf1(ii) = pdf1(ii) * (1.0_dp - 2.0_dp * cdf1)
                pdf2(ii) = 0.0_dp
                dpdf2(ii) = 0.0_dp
            else
                eta1 = alpha(jj)   + lp(ii)
                eta2 = alpha(jj+1) + lp(ii)
                e1 = exp(-eta1)
                e2 = exp(-eta2)
                cdf1 = 1.0_dp / (1.0_dp + e1)
                cdf2 = 1.0_dp / (1.0_dp + e2)
                dd(ii) = cdf1 - cdf2
                pdf1(ii) = cdf1 * (1.0_dp - cdf1)
                pdf2(ii) = cdf2 * (1.0_dp - cdf2)
                dpdf1(ii) = pdf1(ii) * (1.0_dp - 2.0_dp * cdf1)
                dpdf2(ii) = pdf2(ii) * (1.0_dp - 2.0_dp * cdf2)
            end if
        end do
        
    case(2)  ! Probit
        do ii = 1, nn
            jj = yy(ii)
            if (jj == 0) then
                eta1 = alpha(1) + lp(ii)
                cdf1 = 0.5_dp * (1.0_dp + erf(eta1 * INV_SQRT2))
                dd(ii) = 1.0_dp - cdf1
                pdf1(ii) = INV_SQRT2PI * exp(-0.5_dp * eta1 * eta1)
                dpdf1(ii) = -pdf1(ii) * eta1
                pdf2(ii) = 0.0_dp
                dpdf2(ii) = 0.0_dp
            else if (jj == kk) then
                eta1 = alpha(kk) + lp(ii)
                cdf1 = 0.5_dp * (1.0_dp + erf(eta1 * INV_SQRT2))
                dd(ii) = cdf1
                pdf1(ii) = INV_SQRT2PI * exp(-0.5_dp * eta1 * eta1)
                dpdf1(ii) = -pdf1(ii) * eta1
                pdf2(ii) = 0.0_dp
                dpdf2(ii) = 0.0_dp
            else
                eta1 = alpha(jj)   + lp(ii)
                eta2 = alpha(jj+1) + lp(ii)
                cdf1 = 0.5_dp * (1.0_dp + erf(eta1 * INV_SQRT2))
                cdf2 = 0.5_dp * (1.0_dp + erf(eta2 * INV_SQRT2))
                dd(ii) = cdf1 - cdf2
                pdf1(ii) = INV_SQRT2PI * exp(-0.5_dp * eta1 * eta1)
                pdf2(ii) = INV_SQRT2PI * exp(-0.5_dp * eta2 * eta2)
                dpdf1(ii) = -pdf1(ii) * eta1
                dpdf2(ii) = -pdf2(ii) * eta2
            end if
        end do
        
    case(3)  ! Log-log
        do ii = 1, nn
            jj = yy(ii)
            if (jj == 0) then
                eta1 = alpha(1) + lp(ii)
                e1 = exp(-eta1)
                cdf1 = exp(-e1)
                dd(ii) = 1.0_dp - cdf1
                pdf1(ii) = e1 * cdf1
                dpdf1(ii) = pdf1(ii) * (-1.0_dp + e1)
                pdf2(ii) = 0.0_dp
                dpdf2(ii) = 0.0_dp
            else if (jj == kk) then
                eta1 = alpha(kk) + lp(ii)
                e1 = exp(-eta1)
                cdf1 = exp(-e1)
                dd(ii) = cdf1
                pdf1(ii) = e1 * cdf1
                dpdf1(ii) = pdf1(ii) * (-1.0_dp + e1)
                pdf2(ii) = 0.0_dp
                dpdf2(ii) = 0.0_dp
            else
                eta1 = alpha(jj)   + lp(ii)
                eta2 = alpha(jj+1) + lp(ii)
                e1 = exp(-eta1)
                e2 = exp(-eta2)
                cdf1 = exp(-e1)
                cdf2 = exp(-e2)
                dd(ii) = cdf1 - cdf2
                pdf1(ii) = e1 * cdf1
                pdf2(ii) = e2 * cdf2
                dpdf1(ii) = pdf1(ii) * (-1.0_dp + e1)
                dpdf2(ii) = pdf2(ii) * (-1.0_dp + e2)
            end if
        end do
        
    case(4)  ! Complementary log-log
        do ii = 1, nn
            jj = yy(ii)
            if (jj == 0) then
                eta1 = alpha(1) + lp(ii)
                e1 = exp(eta1)
                cdf1 = 1.0_dp - exp(-e1)
                dd(ii) = 1.0_dp - cdf1
                pdf1(ii) = e1 * exp(-e1)
                dpdf1(ii) = pdf1(ii) * (1.0_dp - e1)
                pdf2(ii) = 0.0_dp
                dpdf2(ii) = 0.0_dp
            else if (jj == kk) then
                eta1 = alpha(kk) + lp(ii)
                e1 = exp(eta1)
                cdf1 = 1.0_dp - exp(-e1)
                dd(ii) = cdf1
                pdf1(ii) = e1 * exp(-e1)
                dpdf1(ii) = pdf1(ii) * (1.0_dp - e1)
                pdf2(ii) = 0.0_dp
                dpdf2(ii) = 0.0_dp
            else
                eta1 = alpha(jj)   + lp(ii)
                eta2 = alpha(jj+1) + lp(ii)
                e1 = exp(eta1)
                e2 = exp(eta2)
                cdf1 = 1.0_dp - exp(-e1)
                cdf2 = 1.0_dp - exp(-e2)
                dd(ii) = cdf1 - cdf2
                pdf1(ii) = e1 * exp(-e1)
                pdf2(ii) = e2 * exp(-e2)
                dpdf1(ii) = pdf1(ii) * (1.0_dp - e1)
                dpdf2(ii) = pdf2(ii) * (1.0_dp - e2)
            end if
        end do
        
    case(5)  ! Cauchit
        do ii = 1, nn
            jj = yy(ii)
            if (jj == 0) then
                eta1 = alpha(1) + lp(ii)
                t1 = 1.0_dp + eta1 * eta1
                cdf1 = INV_PI * atan(eta1) + 0.5_dp
                dd(ii) = 1.0_dp - cdf1
                pdf1(ii) = INV_PI / t1
                dpdf1(ii) = -2.0_dp * eta1 * INV_PI / (t1 * t1)
                pdf2(ii) = 0.0_dp
                dpdf2(ii) = 0.0_dp
            else if (jj == kk) then
                eta1 = alpha(kk) + lp(ii)
                t1 = 1.0_dp + eta1 * eta1
                cdf1 = INV_PI * atan(eta1) + 0.5_dp
                dd(ii) = cdf1
                pdf1(ii) = INV_PI / t1
                dpdf1(ii) = -2.0_dp * eta1 * INV_PI / (t1 * t1)
                pdf2(ii) = 0.0_dp
                dpdf2(ii) = 0.0_dp
            else
                eta1 = alpha(jj)   + lp(ii)
                eta2 = alpha(jj+1) + lp(ii)
                t1 = 1.0_dp + eta1 * eta1
                t2 = 1.0_dp + eta2 * eta2
                cdf1 = INV_PI * atan(eta1) + 0.5_dp
                cdf2 = INV_PI * atan(eta2) + 0.5_dp
                dd(ii) = cdf1 - cdf2
                pdf1(ii) = INV_PI / t1
                pdf2(ii) = INV_PI / t2
                dpdf1(ii) = -2.0_dp * eta1 * INV_PI / (t1 * t1)
                dpdf2(ii) = -2.0_dp * eta2 * INV_PI / (t2 * t2)
            end if
        end do
        
    end select
    
end subroutine compute_probs_fast

end subroutine ormll_beta_ccd_enet