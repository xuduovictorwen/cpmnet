! elastic net CPM path, one coefficient at a time
! build without -ffast-math or -march=native

subroutine ormll_path_enet(n, k, p, x, y, wt, link, &
                           nlambda, lambda_path, alpha_en, &
                           penalty_factor, &
                           alpha_init, beta_init, &
                           max_iter_outer, tol, &
                           alpha_mat, beta_mat, logL_vec, &
                           converged_vec, niter_vec, salloc, shuffle_seed)

use cpm_links
use, intrinsic :: ISO_FORTRAN_ENV, only: dp => real64, int32
implicit none

integer(int32), intent(in) :: n, k, p, nlambda, link, max_iter_outer, shuffle_seed
integer(int32), intent(in) :: y(n)
real(dp), intent(in) :: x(n, p), wt(n), lambda_path(nlambda), alpha_en, tol
real(dp), intent(in) :: penalty_factor(p)
real(dp), intent(in) :: alpha_init(k), beta_init(p)

real(dp), intent(out) :: alpha_mat(k, nlambda), beta_mat(p, nlambda), logL_vec(nlambda)
logical, intent(out) :: converged_vec(nlambda)
integer(int32), intent(out) :: niter_vec(nlambda), salloc

integer(int32) :: i_lambda, iter, i, j, l, tri_status, ls
integer(int32) :: m, jj, itmp, nseed
integer(int32), allocatable :: ord(:), seedv(:)
real(dp) :: u
logical :: accepted, dpos
real(dp) :: lambda, lambda1, lambda2, max_change, beta_old, s, w, z, z_j
real(dp) :: inv_n, inv_d, inv_d2, dev_old, dev_new, step_t

real(dp), allocatable :: alpha_curr(:), beta_curr(:), alpha_prev(:), beta_prev(:)
real(dp), allocatable :: lp(:), d(:), pdf1(:), pdf2(:), dpdf1(:), dpdf2(:)
real(dp), allocatable :: grad_alpha(:), hess_diag(:), hess_offdiag(:), alpha_step(:)

allocate(alpha_curr(k), beta_curr(p), alpha_prev(k), beta_prev(p), &
         lp(n), d(n), pdf1(n), pdf2(n), dpdf1(n), dpdf2(n), &
         grad_alpha(k), hess_diag(k), hess_offdiag(max(k-1,1)), alpha_step(k), &
         ord(p), stat=salloc)
if (salloc /= 0) return

! visiting order, redrawn every cycle
call random_seed(size = nseed)
allocate(seedv(nseed))
seedv = shuffle_seed + 37 * [(i, i = 0, nseed - 1)]
call random_seed(put = seedv)
ord = [(l, l = 1, p)]

alpha_curr = alpha_init
beta_curr = beta_init
inv_n = 1.0_dp / real(n, dp)

lambdas: do i_lambda = 1, nlambda

    lambda = lambda_path(i_lambda)
    lambda1 = lambda * alpha_en
    lambda2 = lambda * (1.0_dp - alpha_en)

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

        call compute_probs(n, k, alpha_curr, lp, y, link, d, pdf1, pdf2, dpdf1, dpdf2)

        ! intercept gradient and tridiagonal hessian
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
            ! singular hessian: keep the lambdas already fitted
            salloc = 1000 + tri_status
            niter_vec(i_lambda) = -1
            exit lambdas
        end if

        dev_old = 0.0_dp
        do i = 1, n
            if (d(i) <= 0.0_dp) then
                dev_old = huge(1.0_dp)
                exit
            end if
            dev_old = dev_old - wt(i) * log(d(i))
        end do

        ! damped newton step
        step_t = 1.0_dp
        accepted = .false.
        do ls = 1, MAX_HALVING
            alpha_curr = alpha_prev - step_t * alpha_step
            call compute_probs(n, k, alpha_curr, lp, y, link, d, pdf1, pdf2, dpdf1, dpdf2)
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
            alpha_curr = alpha_prev
            call compute_probs(n, k, alpha_curr, lp, y, link, d, pdf1, pdf2, dpdf1, dpdf2)
        end if

        ! coefficients, in a fisher-yates shuffled order
        do m = p, 2, -1
            call random_number(u)
            jj = 1 + int(u * real(m, dp))
            itmp = ord(m); ord(m) = ord(jj); ord(jj) = itmp
        end do
        do m = 1, p
            l = ord(m)
            beta_old = beta_curr(l)

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

            z = w * beta_old - z_j
            if (abs(z) <= lambda1 * penalty_factor(l)) then
                beta_curr(l) = 0.0_dp
            else
                beta_curr(l) = sign(1.0_dp, z) * (abs(z) - lambda1 * penalty_factor(l)) / (w + lambda2 * penalty_factor(l))
            end if

            s = beta_curr(l) - beta_old
            if (s /= 0.0_dp) then
                do i = 1, n
                    lp(i) = lp(i) + x(i, l) * s
                end do
                call compute_probs(n, k, alpha_curr, lp, y, link, d, pdf1, pdf2, dpdf1, dpdf2)
            end if
        end do

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

    call compute_probs(n, k, alpha_curr, lp, y, link, d, pdf1, pdf2, dpdf1, dpdf2)
    logL_vec(i_lambda) = 0.0_dp
    do i = 1, n
        logL_vec(i_lambda) = logL_vec(i_lambda) + wt(i) * log(d(i))
    end do
    logL_vec(i_lambda) = -2.0_dp * logL_vec(i_lambda)

    alpha_mat(:, i_lambda) = alpha_curr
    beta_mat(:, i_lambda) = beta_curr

end do lambdas

end subroutine ormll_path_enet
