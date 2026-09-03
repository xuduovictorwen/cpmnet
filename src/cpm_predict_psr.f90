! per-observation CDF for PSR
subroutine cpm_predict_psr(n, p, k, newx, alpha, beta, link_type, &
                           obs_idx, cdf_out, n_lambda)
  implicit none
  integer, intent(in) :: n, p, k, n_lambda, link_type
  double precision, intent(in) :: newx(n, p)
  double precision, intent(in) :: alpha(k, n_lambda)
  double precision, intent(in) :: beta(p, n_lambda)
  integer, intent(in) :: obs_idx(n)
  double precision, intent(out) :: cdf_out(n, n_lambda)
  integer :: i, l, idx
  double precision :: xb, eta
  double precision, parameter :: inv_sqrt2 = 0.7071067811865475d0
  double precision, parameter :: inv_pi = 0.31830988618379067154d0
  double precision :: xb_vec(n)
  cdf_out = 0.0d0
  do l = 1, n_lambda
    xb_vec = matmul(newx, beta(:, l))
    do i = 1, n
      idx = obs_idx(i)
      if (idx <= 0) then
        cdf_out(i, l) = 0.0d0
      else if (idx > k) then
        cdf_out(i, l) = 1.0d0
      else
        xb = xb_vec(i)
        eta = alpha(idx, l) + xb
        select case(link_type)
        case(1) ! logistic
          cdf_out(i, l) = 1.0d0 - 1.0d0 / (1.0d0 + exp(-eta))
        case(2) ! probit
          cdf_out(i, l) = 1.0d0 - 0.5d0 * (1.0d0 + erf(eta * inv_sqrt2))
        case(3) ! log-log
          cdf_out(i, l) = 1.0d0 - exp(-exp(-eta))   ! 1-F = P(Y<=c) [loglog]
        case(4) ! complementary log-log
          cdf_out(i, l) = exp(-exp(eta))   ! 1-F = P(Y<=c) [cloglog]
        case(5) ! cauchit
          cdf_out(i, l) = 0.5d0 - atan(eta) * inv_pi
        end select
      end if
    end do
  end do
end subroutine cpm_predict_psr
