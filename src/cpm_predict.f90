! conditional median or mean
subroutine cpm_predict(n, p, k, newx, alpha, beta, y_mapping, &
                       link_type, pred_type, predictions, &
                       n_lambda, tol)
  use cpm_links
  use, intrinsic :: ISO_FORTRAN_ENV, only: dp => real64, int32
  implicit none
  integer(int32), intent(in) :: n, p, k, n_lambda, link_type, pred_type
  real(dp), intent(in) :: newx(n, p), alpha(k, n_lambda), beta(p, n_lambda)
  real(dp), intent(in) :: y_mapping(k+1), tol
  real(dp), intent(out) :: predictions(n, n_lambda)
  integer(int32) :: i, j, l, median_idx
  real(dp) :: xb, cdf_curr, cdf_prev, mean_pred
  real(dp) :: xb_vec(n), cdf_array(k+1)
  predictions = 0.0_dp
  do l = 1, n_lambda
    xb_vec = matmul(newx, beta(:, l))
    if (pred_type == 1) then
      do i = 1, n
        xb = xb_vec(i)
        do j = 1, k
          cdf_array(j) = link_cdf(alpha(j, l) + xb, link_type)
        end do
        cdf_array(k+1) = 1.0_dp
        median_idx = 1
        do j = 1, k
          if (cdf_array(j) >= 0.5_dp) exit
          median_idx = j + 1
        end do
        if (median_idx <= k .and. abs(cdf_array(median_idx) - 0.5_dp) < tol) then
          predictions(i, l) = (y_mapping(median_idx) + y_mapping(median_idx + 1)) * 0.5_dp
        else
          predictions(i, l) = y_mapping(median_idx)
        end if
      end do
    else
      do i = 1, n
        xb = xb_vec(i)
        mean_pred = 0.0_dp
        cdf_prev = 0.0_dp
        do j = 1, k
          cdf_curr = link_cdf(alpha(j, l) + xb, link_type)
          mean_pred = mean_pred + y_mapping(j) * (cdf_curr - cdf_prev)
          cdf_prev = cdf_curr
        end do
        predictions(i, l) = mean_pred + y_mapping(k+1) * (1.0_dp - cdf_prev)
      end do
    end if
  end do
end subroutine cpm_predict
