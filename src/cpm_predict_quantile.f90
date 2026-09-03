! conditional quantiles at tau levels
subroutine cpm_predict_quantile(n, p, k, newx, alpha, beta, y_mapping, &
                                link_type, tau_levels, n_tau, predictions, &
                                n_lambda, tol)
  implicit none
  integer, intent(in) :: n, p, k, n_lambda, link_type, n_tau
  double precision, intent(in) :: newx(n, p)
  double precision, intent(in) :: alpha(k, n_lambda)
  double precision, intent(in) :: beta(p, n_lambda)
  double precision, intent(in) :: y_mapping(k+1)
  double precision, intent(in) :: tau_levels(n_tau)
  double precision, intent(in) :: tol
  double precision, intent(out) :: predictions(n, n_lambda, n_tau)
  integer :: i, j, l, t, q_idx
  double precision :: xb, tau
  double precision, parameter :: inv_sqrt2 = 0.7071067811865475d0
  double precision, parameter :: inv_pi = 0.31830988618379067154d0
  double precision :: xb_vec(n)
  double precision :: cdf_array(k+1)
  predictions = 0.0d0
  do l = 1, n_lambda
    xb_vec = matmul(newx, beta(:, l))
    do i = 1, n
      xb = xb_vec(i)
      do j = 1, k
        cdf_array(j) = link_cdf(alpha(j, l) + xb)
      end do
      cdf_array(k+1) = 1.0d0
      do t = 1, n_tau
        tau = tau_levels(t)
        q_idx = 1
        do j = 1, k
          if (cdf_array(j) >= tau) exit
          q_idx = j + 1
        end do
        if (q_idx <= k .and. abs(cdf_array(q_idx) - tau) < tol) then
          predictions(i, l, t) = (y_mapping(q_idx) + y_mapping(q_idx + 1)) * 0.5d0
        else
          predictions(i, l, t) = y_mapping(q_idx)
        end if
      end do
    end do
  end do

contains

  pure function link_cdf(eta) result(c)
    double precision, intent(in) :: eta
    double precision :: c
    select case(link_type)
    case(1) ! logistic
      c = 1.0d0 - 1.0d0 / (1.0d0 + exp(-eta))
    case(2) ! probit
      c = 1.0d0 - 0.5d0 * (1.0d0 + erf(eta * inv_sqrt2))
    case(3) ! log-log
      c = 1.0d0 - exp(-exp(-eta))   ! 1-F = P(Y<=c) [loglog]
    case(4) ! complementary log-log
      c = exp(-exp(eta))   ! 1-F = P(Y<=c) [cloglog]
    case(5) ! cauchit
      c = 0.5d0 - atan(eta) * inv_pi
    case default
      c = 0.0d0
    end select
  end function link_cdf

end subroutine cpm_predict_quantile
