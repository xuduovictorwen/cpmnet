! Compile: gfortran -shared -fPIC -O2 -o ormll_path_enet.so ormll_path_enet.f90
! Note: -ffast-math and -march=native are intentionally omitted for
! portability and IEEE correctness. The numerical guards below (Newton
! line search, pivot checks) rely on standard IEEE semantics.

subroutine ormll_path_enet(n, k, p, x, y, wt, link, &
                           nlambda, lambda_path, alpha_en, &
                           penalty_factor , &
                           alpha_init, beta_init, &
                           max_iter_outer, tol, &
                           alpha_mat, beta_mat, logL_vec, &
                           converged_vec, niter_vec, salloc)

use, intrinsic :: ISO_FORTRAN_ENV, only: dp => real64, int32
implicit none

! Precomputed constants
real(dp), parameter :: PI = 3.14159265358979323846_dp
real(dp), parameter :: SQRT2 = 1.41421356237309504880_dp
real(dp), parameter :: INV_SQRT2 = 0.70710678118654752440_dp
real(dp), parameter :: INV_SQRT2PI = 0.39894228040143267794_dp
real(dp), parameter :: INV_PI = 0.31830988618379067154_dp

! Numerical safety constants
! The intercept update is a damped Newton step (backtracking line search,
! MAX_HALVING below). Each trial step is accepted only if every interval
! probability d(i) stays > 0 AND the deviance does not increase, so the
! fitted intercepts remain strictly ordered and d(i) is positive by
! construction. No probability floor is therefore needed (or used): d(i)
! enters the gradient, Hessian, and log-likelihood directly. This is what
! lets the heavy-tailed cauchit link converge -- an undamped Newton step
! overshoots there (the Cauchy log-likelihood is not concave) and used to
! drive the intercepts to +/-1e4, which a probability floor only masked.
! The same d(i) > 0 guard also makes a separate eta clamp unnecessary: for
! the Gumbel links (loglog, cloglog) the interval probability underflows to
! 0 by |eta| ~ 6.6, so any step heading toward the inner-exp() overflow
! threshold (eta ~ 709.8) is rejected long before eta could get there.
! PIVOT_EPS: minimum absolute pivot in the tridiagonal Thomas solver.
!          Below this we flag the Hessian as effectively singular and
!          return a nonzero salloc.
! MAX_HALVING: maximum backtracking steps in the intercept line search.
integer(int32), parameter :: MAX_HALVING = 30
real(dp), parameter :: PIVOT_EPS = 1.0e-14_dp

integer(int32), intent(in) :: n, k, p, nlambda, link, max_iter_outer
integer(int32), intent(in) :: y(n)
real(dp), intent(in) :: x(n, p), wt(n), lambda_path(nlambda), alpha_en, tol
real(dp), intent(in) :: penalty_factor(p)
real(dp), intent(in) :: alpha_init(k), beta_init(p)

real(dp), intent(out) :: alpha_mat(k, nlambda), beta_mat(p, nlambda), logL_vec(nlambda)
logical, intent(out) :: converged_vec(nlambda)
integer(int32), intent(out) :: niter_vec(nlambda), salloc

integer(int32) :: i_lambda, iter, i, j, l, tri_status, ls
logical :: accepted, dpos
real(dp) :: lambda, lambda1, lambda2, max_change, beta_old, s, w, z, z_j
real(dp) :: inv_n, inv_d, inv_d2, dev_old, dev_new, step_t

real(dp), allocatable :: alpha_curr(:), beta_curr(:), alpha_prev(:), beta_prev(:)
real(dp), allocatable :: lp(:), d(:), pdf1(:), pdf2(:), dpdf1(:), dpdf2(:)
real(dp), allocatable :: grad_alpha(:), hess_diag(:), hess_offdiag(:), alpha_step(:)
real(dp), allocatable :: eta1(:), eta2(:)  ! Precomputed eta values

allocate(alpha_curr(k), beta_curr(p), alpha_prev(k), beta_prev(p), &
         lp(n), d(n), pdf1(n), pdf2(n), dpdf1(n), dpdf2(n), &
         grad_alpha(k), hess_diag(k), hess_offdiag(k-1), alpha_step(k), &
         eta1(n), eta2(n), stat=salloc)
if (salloc /= 0) return

alpha_curr = alpha_init
beta_curr = beta_init
inv_n = 1.0_dp / real(n, dp)

do i_lambda = 1, nlambda

    lambda = lambda_path(i_lambda)
    lambda1 = lambda * alpha_en
    lambda2 = lambda * (1.0_dp - alpha_en)

    ! Compute initial linear predictor
    lp = 0.0_dp
    do l = 1, p
        if (beta_curr(l) /= 0.0_dp) then
            do i = 1, n
                lp(i) = lp(i) + x(i, l) * beta_curr(l)
            end do
        end if
    end do

    converged_vec(i_lambda) = .false.

    do iter = 1, max_iter_outer

        alpha_prev = alpha_curr
        beta_prev = beta_curr

        ! Compute probabilities using optimized routine
        call compute_probs_fast(n, k, alpha_curr, lp, y, link, d, pdf1, pdf2, dpdf1, dpdf2, eta1, eta2)

        ! Compute gradient for alpha (vectorized where possible)
        grad_alpha = 0.0_dp
        do i = 1, n
            j = y(i)
            inv_d = wt(i) / d(i)
            if (j == 0) then
                grad_alpha(1) = grad_alpha(1) + inv_d * pdf1(i)
            else if (j == k) then
                grad_alpha(k) = grad_alpha(k) - inv_d * pdf1(i)
            else
                grad_alpha(j) = grad_alpha(j) - inv_d * pdf1(i)
                grad_alpha(j+1) = grad_alpha(j+1) + inv_d * pdf2(i)
            end if
        end do

        ! Compute Hessian for alpha
        hess_diag = 0.0_dp
        hess_offdiag = 0.0_dp
        do i = 1, n
            j = y(i)
            inv_d2 = wt(i) / (d(i) * d(i))
            if (j == 0) then
                z = pdf1(i) * pdf1(i) + d(i) * dpdf1(i)
                hess_diag(1) = hess_diag(1) + inv_d2 * z
            else if (j == k) then
                z = pdf1(i) * pdf1(i) - d(i) * dpdf1(i)
                hess_diag(k) = hess_diag(k) + inv_d2 * z
            else
                hess_diag(j) = hess_diag(j) + inv_d2 * (pdf1(i) * pdf1(i) - d(i) * dpdf1(i))
                hess_diag(j+1) = hess_diag(j+1) + inv_d2 * (pdf2(i) * pdf2(i) + d(i) * dpdf2(i))
                hess_offdiag(j) = hess_offdiag(j) - inv_d2 * pdf1(i) * pdf2(i)
            end if
        end do

        call solve_tridiag(k, hess_diag, hess_offdiag, grad_alpha, alpha_step, tri_status)
        if (tri_status /= 0) then
            salloc = 1000 + tri_status  ! singular Hessian, surface to R
            deallocate(alpha_curr, beta_curr, alpha_prev, beta_prev, &
                       lp, d, pdf1, pdf2, dpdf1, dpdf2, &
                       grad_alpha, hess_diag, hess_offdiag, alpha_step, eta1, eta2)
            return
        end if

        ! Deviance at alpha_prev. d still holds the alpha_prev probabilities
        ! computed above; for a strictly ordered alpha every d(i) > 0. Guard
        ! against a non-ordered start so the line search can still recover.
        dev_old = 0.0_dp
        do i = 1, n
            if (d(i) <= 0.0_dp) then
                dev_old = huge(1.0_dp)
                exit
            end if
            dev_old = dev_old - wt(i) * log(d(i))
        end do

        ! Damped Newton: backtrack along the full step until every interval
        ! probability is positive AND the deviance does not increase. This
        ! keeps the intercepts strictly ordered (so no probability floor is
        ! needed) and makes the update globally convergent for every link,
        ! including the non-log-concave cauchit, whose undamped step diverges.
        step_t = 1.0_dp
        accepted = .false.
        do ls = 1, MAX_HALVING
            alpha_curr = alpha_prev - step_t * alpha_step
            call compute_probs_fast(n, k, alpha_curr, lp, y, link, d, pdf1, pdf2, dpdf1, dpdf2, eta1, eta2)
            dpos = .true.
            do i = 1, n
                if (d(i) <= 0.0_dp) then
                    dpos = .false.
                    exit
                end if
            end do
            if (dpos) then
                dev_new = 0.0_dp
                do i = 1, n
                    dev_new = dev_new - wt(i) * log(d(i))
                end do
                if (dev_new <= dev_old) then
                    accepted = .true.
                    exit
                end if
            end if
            step_t = 0.5_dp * step_t
        end do

        if (.not. accepted) then
            ! No improving step along the Newton direction (indefinite Hessian
            ! far from the optimum); leave the intercepts unchanged and let the
            ! beta sweep / next outer iteration make progress.
            alpha_curr = alpha_prev
            call compute_probs_fast(n, k, alpha_curr, lp, y, link, d, pdf1, pdf2, dpdf1, dpdf2, eta1, eta2)
        end if

        ! Coordinate descent on beta
        do l = 1, p
            beta_old = beta_curr(l)

            ! Compute gradient and Hessian for beta_l in single pass
            z_j = 0.0_dp
            w = 0.0_dp
            do i = 1, n
                j = y(i)
                inv_d = 1.0_dp / d(i)
                inv_d2 = inv_d * inv_d

                if (j == 0) then
                    z_j = z_j + wt(i) * pdf1(i) * x(i, l) * inv_d
                    z = pdf1(i) * pdf1(i) + d(i) * dpdf1(i)
                else if (j == k) then
                    z_j = z_j - wt(i) * pdf1(i) * x(i, l) * inv_d
                    z = pdf1(i) * pdf1(i) - d(i) * dpdf1(i)
                else
                    z_j = z_j - wt(i) * x(i, l) * (pdf1(i) - pdf2(i)) * inv_d
                    z = (pdf1(i) - pdf2(i)) * (pdf1(i) - pdf2(i)) - d(i) * (dpdf1(i) - dpdf2(i))
                end if
                w = w + wt(i) * inv_d2 * z * x(i, l) * x(i, l)
            end do
            z_j = z_j * inv_n
            w = w * inv_n

            ! Elastic net update
            z = w * beta_old - z_j
            if (abs(z) <= lambda1 * penalty_factor(l)) then
                beta_curr(l) = 0.0_dp
            else
                beta_curr(l) = sign(1.0_dp, z) * (abs(z) - lambda1 * penalty_factor(l)) / (w + lambda2 * penalty_factor(l))
            end if

            ! Update linear predictor and recompute probabilities if beta changed
            s = beta_curr(l) - beta_old
            if (s /= 0.0_dp) then
                do i = 1, n
                    lp(i) = lp(i) + x(i, l) * s
                end do
                call compute_probs_fast(n, k, alpha_curr, lp, y, link, d, pdf1, pdf2, dpdf1, dpdf2, eta1, eta2)
            end if
        end do

        ! Check convergence
        max_change = 0.0_dp
        do j = 1, k
            if (abs(alpha_curr(j) - alpha_prev(j)) > max_change) then
                max_change = abs(alpha_curr(j) - alpha_prev(j))
            end if
        end do
        do l = 1, p
            if (abs(beta_curr(l) - beta_prev(l)) > max_change) then
                max_change = abs(beta_curr(l) - beta_prev(l))
            end if
        end do

        if (max_change < tol) then
            converged_vec(i_lambda) = .true.
            niter_vec(i_lambda) = iter
            exit
        end if

        if (iter == max_iter_outer) then
            niter_vec(i_lambda) = max_iter_outer
        end if

    end do

    ! Final log-likelihood
    call compute_probs_fast(n, k, alpha_curr, lp, y, link, d, pdf1, pdf2, dpdf1, dpdf2, eta1, eta2)
    logL_vec(i_lambda) = 0.0_dp
    do i = 1, n
        logL_vec(i_lambda) = logL_vec(i_lambda) + wt(i) * log(d(i))
    end do
    logL_vec(i_lambda) = -2.0_dp * logL_vec(i_lambda)

    alpha_mat(:, i_lambda) = alpha_curr
    beta_mat(:, i_lambda) = beta_curr

end do

deallocate(alpha_curr, beta_curr, alpha_prev, beta_prev, &
           lp, d, pdf1, pdf2, dpdf1, dpdf2, &
           grad_alpha, hess_diag, hess_offdiag, alpha_step, eta1, eta2)

return

contains

! Optimized probability computation - computes CDF, PDF, dPDF in single pass
! Avoids redundant exp() calls by reusing intermediate values
subroutine compute_probs_fast(nn, kk, alpha, lp, yy, link_type, dd, pdf1, pdf2, dpdf1, dpdf2, eta1, eta2)
    integer(int32), intent(in) :: nn, kk, link_type
    integer(int32), intent(in) :: yy(nn)
    real(dp), intent(in) :: alpha(kk), lp(nn)
    real(dp), intent(out) :: dd(nn), pdf1(nn), pdf2(nn), dpdf1(nn), dpdf2(nn)
    real(dp), intent(out) :: eta1(nn), eta2(nn)  ! Workspace for eta values

    integer(int32) :: ii, jj
    real(dp) :: cdf1, cdf2, e1, e2, t1, t2

    ! Precompute eta values to improve cache locality. No clamp: the d(i) > 0
    ! line-search guard rejects any step that would drive |eta| toward the
    ! inner-exp() overflow threshold (see header).
    do ii = 1, nn
        jj = yy(ii)
        if (jj == 0) then
            eta1(ii) = alpha(1) + lp(ii)
        else if (jj == kk) then
            eta1(ii) = alpha(kk) + lp(ii)
        else
            eta1(ii) = alpha(jj)   + lp(ii)
            eta2(ii) = alpha(jj+1) + lp(ii)
        end if
    end do

    ! Compute CDF, PDF, dPDF based on link function
    select case(link_type)

    case(1)  ! Logistic
        do ii = 1, nn
            jj = yy(ii)
            if (jj == 0) then
                e1 = exp(-eta1(ii))
                cdf1 = 1.0_dp / (1.0_dp + e1)
                dd(ii) = 1.0_dp - cdf1
                pdf1(ii) = cdf1 * dd(ii)
                dpdf1(ii) = pdf1(ii) * (1.0_dp - 2.0_dp * cdf1)
                pdf2(ii) = 0.0_dp
                dpdf2(ii) = 0.0_dp
            else if (jj == kk) then
                e1 = exp(-eta1(ii))
                cdf1 = 1.0_dp / (1.0_dp + e1)
                dd(ii) = cdf1
                pdf1(ii) = cdf1 * (1.0_dp - cdf1)
                dpdf1(ii) = pdf1(ii) * (1.0_dp - 2.0_dp * cdf1)
                pdf2(ii) = 0.0_dp
                dpdf2(ii) = 0.0_dp
            else
                e1 = exp(-eta1(ii))
                e2 = exp(-eta2(ii))
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
                cdf1 = 0.5_dp * (1.0_dp + erf(eta1(ii) * INV_SQRT2))
                dd(ii) = 1.0_dp - cdf1
                pdf1(ii) = INV_SQRT2PI * exp(-0.5_dp * eta1(ii) * eta1(ii))
                dpdf1(ii) = -pdf1(ii) * eta1(ii)
                pdf2(ii) = 0.0_dp
                dpdf2(ii) = 0.0_dp
            else if (jj == kk) then
                cdf1 = 0.5_dp * (1.0_dp + erf(eta1(ii) * INV_SQRT2))
                dd(ii) = cdf1
                pdf1(ii) = INV_SQRT2PI * exp(-0.5_dp * eta1(ii) * eta1(ii))
                dpdf1(ii) = -pdf1(ii) * eta1(ii)
                pdf2(ii) = 0.0_dp
                dpdf2(ii) = 0.0_dp
            else
                cdf1 = 0.5_dp * (1.0_dp + erf(eta1(ii) * INV_SQRT2))
                cdf2 = 0.5_dp * (1.0_dp + erf(eta2(ii) * INV_SQRT2))
                dd(ii) = cdf1 - cdf2
                pdf1(ii) = INV_SQRT2PI * exp(-0.5_dp * eta1(ii) * eta1(ii))
                pdf2(ii) = INV_SQRT2PI * exp(-0.5_dp * eta2(ii) * eta2(ii))
                dpdf1(ii) = -pdf1(ii) * eta1(ii)
                dpdf2(ii) = -pdf2(ii) * eta2(ii)
            end if
        end do

    case(3)  ! Log-log
        do ii = 1, nn
            jj = yy(ii)
            if (jj == 0) then
                e1 = exp(-eta1(ii))
                cdf1 = exp(-e1)
                dd(ii) = 1.0_dp - cdf1
                pdf1(ii) = e1 * cdf1
                dpdf1(ii) = pdf1(ii) * (-1.0_dp + e1)
                pdf2(ii) = 0.0_dp
                dpdf2(ii) = 0.0_dp
            else if (jj == kk) then
                e1 = exp(-eta1(ii))
                cdf1 = exp(-e1)
                dd(ii) = cdf1
                pdf1(ii) = e1 * cdf1
                dpdf1(ii) = pdf1(ii) * (-1.0_dp + e1)
                pdf2(ii) = 0.0_dp
                dpdf2(ii) = 0.0_dp
            else
                e1 = exp(-eta1(ii))
                e2 = exp(-eta2(ii))
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
                e1 = exp(eta1(ii))
                cdf1 = 1.0_dp - exp(-e1)
                dd(ii) = 1.0_dp - cdf1
                pdf1(ii) = e1 * exp(-e1)
                dpdf1(ii) = pdf1(ii) * (1.0_dp - e1)
                pdf2(ii) = 0.0_dp
                dpdf2(ii) = 0.0_dp
            else if (jj == kk) then
                e1 = exp(eta1(ii))
                cdf1 = 1.0_dp - exp(-e1)
                dd(ii) = cdf1
                pdf1(ii) = e1 * exp(-e1)
                dpdf1(ii) = pdf1(ii) * (1.0_dp - e1)
                pdf2(ii) = 0.0_dp
                dpdf2(ii) = 0.0_dp
            else
                e1 = exp(eta1(ii))
                e2 = exp(eta2(ii))
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
                t1 = 1.0_dp + eta1(ii) * eta1(ii)
                cdf1 = INV_PI * atan(eta1(ii)) + 0.5_dp
                dd(ii) = 1.0_dp - cdf1
                pdf1(ii) = INV_PI / t1
                dpdf1(ii) = -2.0_dp * eta1(ii) * INV_PI / (t1 * t1)
                pdf2(ii) = 0.0_dp
                dpdf2(ii) = 0.0_dp
            else if (jj == kk) then
                t1 = 1.0_dp + eta1(ii) * eta1(ii)
                cdf1 = INV_PI * atan(eta1(ii)) + 0.5_dp
                dd(ii) = cdf1
                pdf1(ii) = INV_PI / t1
                dpdf1(ii) = -2.0_dp * eta1(ii) * INV_PI / (t1 * t1)
                pdf2(ii) = 0.0_dp
                dpdf2(ii) = 0.0_dp
            else
                t1 = 1.0_dp + eta1(ii) * eta1(ii)
                t2 = 1.0_dp + eta2(ii) * eta2(ii)
                cdf1 = INV_PI * atan(eta1(ii)) + 0.5_dp
                cdf2 = INV_PI * atan(eta2(ii)) + 0.5_dp
                dd(ii) = cdf1 - cdf2
                pdf1(ii) = INV_PI / t1
                pdf2(ii) = INV_PI / t2
                dpdf1(ii) = -2.0_dp * eta1(ii) * INV_PI / (t1 * t1)
                dpdf2(ii) = -2.0_dp * eta2(ii) * INV_PI / (t2 * t2)
            end if
        end do

    end select

end subroutine compute_probs_fast

subroutine solve_tridiag(nn, diag, offdiag, rhs, solution, status)
    integer(int32), intent(in) :: nn
    real(dp), intent(in) :: diag(nn), offdiag(nn-1), rhs(nn)
    real(dp), intent(out) :: solution(nn)
    integer(int32), intent(out) :: status
    real(dp) :: work_diag(nn), work_rhs(nn), m
    integer(int32) :: ii

    status = 0
    work_diag(1) = diag(1)
    work_rhs(1) = rhs(1)
    if (abs(work_diag(1)) < PIVOT_EPS) then
        status = 1
        return
    end if
    do ii = 2, nn
        m = offdiag(ii-1) / work_diag(ii-1)
        work_diag(ii) = diag(ii) - m * offdiag(ii-1)
        if (abs(work_diag(ii)) < PIVOT_EPS) then
            status = ii
            return
        end if
        work_rhs(ii) = rhs(ii) - m * work_rhs(ii-1)
    end do

    solution(nn) = work_rhs(nn) / work_diag(nn)
    do ii = nn-1, 1, -1
        solution(ii) = (work_rhs(ii) - offdiag(ii) * solution(ii+1)) / work_diag(ii)
    end do
end subroutine solve_tridiag

end subroutine ormll_path_enet
