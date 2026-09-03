! conditional median or mean
subroutine cpm_predict(n, p, k, newx, alpha, beta, y_mapping, &
                       link_type, pred_type, predictions, &
                       n_lambda, tol)
  implicit none
  integer, intent(in) :: n, p, k, n_lambda, link_type, pred_type
  double precision, intent(in) :: newx(n, p)
  double precision, intent(in) :: alpha(k, n_lambda)
  double precision, intent(in) :: beta(p, n_lambda)
  double precision, intent(in) :: y_mapping(k+1)
  double precision, intent(in) :: tol
  double precision, intent(out) :: predictions(n, n_lambda)
  integer :: i, j, l, median_idx
  double precision :: xb, cdf_curr, cdf_prev, pmf_j, mean_pred
  double precision, parameter :: inv_sqrt2 = 0.7071067811865475d0
  double precision, parameter :: inv_pi = 0.31830988618379067154d0
  double precision :: xb_vec(n)
  double precision :: cdf_array(k+1)  ! only needed for median
  predictions = 0.0d0
  do l = 1, n_lambda
    xb_vec = matmul(newx, beta(:, l))
    if (pred_type == 1) then

      do i = 1, n
        xb = xb_vec(i)
        do j = 1, k
          cdf_array(j) = link_cdf(alpha(j, l) + xb)
        end do
        cdf_array(k+1) = 1.0d0
        median_idx = 1
        do j = 1, k
          if (cdf_array(j) >= 0.5d0) exit
          median_idx = j + 1
        end do
        if (median_idx <= k .and. abs(cdf_array(median_idx) - 0.5d0) < tol) then
          predictions(i, l) = (y_mapping(median_idx) + y_mapping(median_idx + 1)) * 0.5d0
        else
          predictions(i, l) = y_mapping(median_idx)
        end if
      end do
    else
      ! mean
      do i = 1, n
        xb = xb_vec(i)
        mean_pred = 0.0d0
        cdf_prev = 0.0d0
        do j = 1, k
          cdf_curr = link_cdf(alpha(j, l) + xb)
          pmf_j = cdf_curr - cdf_prev
          mean_pred = mean_pred + y_mapping(j) * pmf_j
          cdf_prev = cdf_curr
        end do

        mean_pred = mean_pred + y_mapping(k+1) * (1.0d0 - cdf_prev)
        predictions(i, l) = mean_pred
      end do
    end if
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

end subroutine cpm_predict
