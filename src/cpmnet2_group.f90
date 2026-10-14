! group elastic net CPM core, adaptively bounded gradient descent (Helwig 2025 JCGS)
! x arrives orthonormalized per group, beta is beta-tilde

! lambda_max at beta = 0
subroutine cpmnet2_lambda_max_group(n, k, p, ngroup, x, y, wt, link, &
                                    group_start, group_size, group_pf, &
                                    alpha_init, lmax, salloc)
use cpm_links
use, intrinsic :: ISO_FORTRAN_ENV, only: dp => real64, int32
implicit none
integer(int32), intent(in) :: n, k, p, ngroup, link, y(n)
integer(int32), intent(in) :: group_start(ngroup), group_size(ngroup)
real(dp), intent(in) :: x(n, p), wt(n), group_pf(ngroup), alpha_init(k)
real(dp), intent(out) :: lmax
integer(int32), intent(out) :: salloc
real(dp), allocatable :: lp(:), d(:), pdf1(:), pdf2(:), dpdf1(:), dpdf2(:), s(:), w(:), wres(:)
real(dp) :: gnorm, gg, acc, val, inv_n
integer(int32) :: g, gs, dg, jj, i
allocate(lp(n), d(n), pdf1(n), pdf2(n), dpdf1(n), dpdf2(n), s(n), w(n), wres(n), stat=salloc)
if (salloc /= 0) return
inv_n = 1.0_dp / real(n, dp); lp = 0.0_dp
call compute_probs(n, k, alpha_init, lp, y, link, d, pdf1, pdf2, dpdf1, dpdf2)
call score_weight(n, k, y, d, pdf1, pdf2, dpdf1, dpdf2, s, w)
do i = 1, n
    wres(i) = -wt(i) * s(i)
end do
lmax = 0.0_dp
do g = 1, ngroup
    gs = group_start(g); dg = group_size(g); gnorm = 0.0_dp
    do jj = 0, dg - 1
        acc = 0.0_dp
        do i = 1, n
            acc = acc + x(i, gs + jj) * wres(i)
        end do
        gg = acc * inv_n
        gnorm = gnorm + gg * gg
    end do
    val = sqrt(gnorm) / (group_pf(g) * sqrt(real(dg, dp)))
    if (val > lmax) lmax = val
end do
end subroutine cpmnet2_lambda_max_group


! full path fit
subroutine cpmnet2_path_group(n, k, p, ngroup, x, y, wt, link, &
                              group_start, group_size, group_pf, group_gamma, &
                              nlambda, lambda_path, alpha_en, &
                              alpha_init, beta_init, max_iter_outer, tol, &
                              alpha_mat, beta_mat, logL_vec, &
                              converged_vec, niter_vec, salloc, shuffle_seed)
use cpm_links
use, intrinsic :: ISO_FORTRAN_ENV, only: dp => real64, int32
implicit none
integer(int32), intent(in) :: n, k, p, ngroup, nlambda, link, max_iter_outer, y(n)
integer(int32), intent(in) :: group_start(ngroup), group_size(ngroup), shuffle_seed
real(dp), intent(in) :: x(n, p), wt(n), lambda_path(nlambda), alpha_en, tol
real(dp), intent(in) :: group_pf(ngroup), group_gamma(ngroup), alpha_init(k), beta_init(p)
real(dp), intent(out) :: alpha_mat(k, nlambda), beta_mat(p, nlambda), logL_vec(nlambda)
logical, intent(out) :: converged_vec(nlambda)
integer(int32), intent(out) :: niter_vec(nlambda), salloc

integer(int32) :: m, gtmp, nseed, gj
integer(int32), allocatable :: ord(:), seedv(:)
real(dp) :: u
integer(int32) :: il, iter, i, j, g, gs, dg, jj, ls, tri_status, col, maxdg
real(dp) :: lambda, l1g, l2g, delta, dgam, bnorm, scale, thr, bchg, acc, inv_n
real(dp) :: dev_old, dev_new, step_t, max_change
logical :: accepted, dpos
real(dp), allocatable :: alpha_c(:), beta_c(:), alpha_p(:), beta_p(:)
real(dp), allocatable :: lp(:), d(:), pdf1(:), pdf2(:), dpdf1(:), dpdf2(:), s(:), w(:), wres(:)
real(dp), allocatable :: grad_a(:), hdiag(:), hoff(:), astep(:), bblk(:), bnew(:)

maxdg = 1
do g = 1, ngroup
    if (group_size(g) > maxdg) maxdg = group_size(g)
end do
allocate(alpha_c(k), beta_c(p), alpha_p(k), beta_p(p), &
         lp(n), d(n), pdf1(n), pdf2(n), dpdf1(n), dpdf2(n), s(n), w(n), wres(n), &
         grad_a(k), hdiag(k), hoff(max(k-1,1)), astep(k), bblk(maxdg), bnew(maxdg), &
         ord(ngroup), stat=salloc)
if (salloc /= 0) return
! group order, redrawn every sweep
call random_seed(size = nseed)
allocate(seedv(nseed))
seedv = shuffle_seed + 37 * [(i, i = 0, nseed - 1)]
call random_seed(put = seedv)
ord = [(g, g = 1, ngroup)]
alpha_c = alpha_init; beta_c = beta_init; inv_n = 1.0_dp / real(n, dp)

lambdas: do il = 1, nlambda
    lambda = lambda_path(il)
    lp = 0.0_dp
    do col = 1, p
        if (beta_c(col) /= 0.0_dp) then
            do i = 1, n
                lp(i) = lp(i) + x(i, col) * beta_c(col)
            end do
        end if
    end do
    converged_vec(il) = .false.

    do iter = 1, max_iter_outer
        alpha_p = alpha_c; beta_p = beta_c
        call compute_probs(n, k, alpha_c, lp, y, link, d, pdf1, pdf2, dpdf1, dpdf2)

        ! intercepts, damped tridiagonal newton
        grad_a = 0.0_dp
        do i = 1, n
            j = y(i)
            if (j == 0) then
                grad_a(1) = grad_a(1) + wt(i) * pdf1(i) / d(i)
            else if (j == k) then
                grad_a(k) = grad_a(k) - wt(i) * pdf1(i) / d(i)
            else
                grad_a(j)   = grad_a(j)   - wt(i) * pdf1(i) / d(i)
                grad_a(j+1) = grad_a(j+1) + wt(i) * pdf2(i) / d(i)
            end if
        end do
        hdiag = 0.0_dp; hoff = 0.0_dp
        do i = 1, n
            j = y(i)
            if (j == 0) then
                hdiag(1) = hdiag(1) + wt(i) * (pdf1(i)*pdf1(i) + d(i)*dpdf1(i)) / (d(i)*d(i))
            else if (j == k) then
                hdiag(k) = hdiag(k) + wt(i) * (pdf1(i)*pdf1(i) - d(i)*dpdf1(i)) / (d(i)*d(i))
            else
                hdiag(j)   = hdiag(j)   + wt(i) * (pdf1(i)*pdf1(i) - d(i)*dpdf1(i)) / (d(i)*d(i))
                hdiag(j+1) = hdiag(j+1) + wt(i) * (pdf2(i)*pdf2(i) + d(i)*dpdf2(i)) / (d(i)*d(i))
                hoff(j)    = hoff(j)    - wt(i) * pdf1(i) * pdf2(i) / (d(i)*d(i))
            end if
        end do
        call solve_tridiag(k, hdiag, hoff, grad_a, astep, tri_status)
        if (tri_status /= 0) then
            ! singular hessian: keep the lambdas already fitted
            salloc = 1000 + tri_status
            niter_vec(il) = -1
            exit lambdas
        end if
        dev_old = 0.0_dp
        do i = 1, n
            if (d(i) <= 0.0_dp) then; dev_old = huge(1.0_dp); exit; end if
            dev_old = dev_old - wt(i) * log(d(i))
        end do
        step_t = 1.0_dp; accepted = .false.
        do ls = 1, MAX_HALVING
            alpha_c = alpha_p - step_t * astep
            call compute_probs(n, k, alpha_c, lp, y, link, d, pdf1, pdf2, dpdf1, dpdf2)
            dpos = .true.
            do i = 1, n
                if (d(i) <= 0.0_dp) then; dpos = .false.; exit; end if
            end do
            if (dpos) then
                dev_new = 0.0_dp
                do i = 1, n
                    dev_new = dev_new - wt(i) * log(d(i))
                end do
                if (dev_new <= dev_old) then; accepted = .true.; exit; end if
            end if
            step_t = 0.5_dp * step_t
        end do
        if (.not. accepted) then
            alpha_c = alpha_p
            call compute_probs(n, k, alpha_c, lp, y, link, d, pdf1, pdf2, dpdf1, dpdf2)
        end if

        ! coefficients, one sweep over the groups in a shuffled order
        call score_weight(n, k, y, d, pdf1, pdf2, dpdf1, dpdf2, s, w)
        delta = 0.0_dp
        do i = 1, n
            if (wt(i) * w(i) > delta) delta = wt(i) * w(i)
        end do
        if (delta > 0.0_dp) then
        do i = 1, n
            wres(i) = -wt(i) * s(i)
        end do
        do m = ngroup, 2, -1
            call random_number(u)
            gj = 1 + int(u * real(m, dp))
            gtmp = ord(m); ord(m) = ord(gj); ord(gj) = gtmp
        end do
        do m = 1, ngroup
            g = ord(m)
            gs = group_start(g); dg = group_size(g)
            dgam = delta * group_gamma(g)
            bnorm = 0.0_dp
            do jj = 0, dg - 1
                acc = 0.0_dp
                do i = 1, n
                    acc = acc + x(i, gs + jj) * wres(i)
                end do
                bblk(jj+1) = beta_c(gs + jj) + acc * inv_n / dgam
                bnorm = bnorm + bblk(jj+1) * bblk(jj+1)
            end do
            bnorm = sqrt(bnorm)
            l1g = lambda * alpha_en          * group_pf(g) * sqrt(real(dg, dp))
            l2g = lambda * (1.0_dp-alpha_en) * group_pf(g) * sqrt(real(dg, dp))
            if (dgam * bnorm <= l1g .or. bnorm <= 0.0_dp) then
                bnew(1:dg) = 0.0_dp
            else
                thr = 1.0_dp - l1g / (dgam * bnorm)
                scale = thr / (1.0_dp + l2g / dgam)
                bnew(1:dg) = scale * bblk(1:dg)
            end if
            do jj = 0, dg - 1
                col = gs + jj
                bchg = bnew(jj+1) - beta_c(col)
                if (bchg /= 0.0_dp) then
                    do i = 1, n
                        wres(i) = wres(i) - delta * x(i, col) * bchg
                        lp(i)   = lp(i)   + x(i, col) * bchg
                    end do
                    beta_c(col) = bnew(jj+1)
                end if
            end do
        end do
        end if

        max_change = 0.0_dp
        do j = 1, k
            if (abs(alpha_c(j) - alpha_p(j)) > max_change) max_change = abs(alpha_c(j) - alpha_p(j))
        end do
        do col = 1, p
            if (abs(beta_c(col) - beta_p(col)) > max_change) max_change = abs(beta_c(col) - beta_p(col))
        end do
        if (max_change < tol) then
            converged_vec(il) = .true.; niter_vec(il) = iter; exit
        end if
        if (iter == max_iter_outer) niter_vec(il) = max_iter_outer
    end do

    call compute_probs(n, k, alpha_c, lp, y, link, d, pdf1, pdf2, dpdf1, dpdf2)
    logL_vec(il) = 0.0_dp
    do i = 1, n
        logL_vec(il) = logL_vec(il) + wt(i) * log(d(i))
    end do
    logL_vec(il) = -2.0_dp * logL_vec(il)
    alpha_mat(:, il) = alpha_c
    beta_mat(:, il) = beta_c
end do lambdas

end subroutine cpmnet2_path_group
