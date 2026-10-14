! conditional CDF at user thresholds
subroutine cpm_predict_cdf(n, p, k, newx, alpha, beta, y_mapping, &
                           link_type, thresholds, n_thresh, cdf_out, &
                           n_lambda)
  use cpm_links
  use, intrinsic :: ISO_FORTRAN_ENV, only: dp => real64, int32
  implicit none
  integer(int32), intent(in) :: n, p, k, n_lambda, link_type, n_thresh
  real(dp), intent(in) :: newx(n, p), alpha(k, n_lambda), beta(p, n_lambda)
  real(dp), intent(in) :: y_mapping(k+1), thresholds(n_thresh)
  real(dp), intent(out) :: cdf_out(n, n_lambda, n_thresh)
  integer(int32) :: i, j, l, t, idx(n_thresh)
  real(dp) :: xb_vec(n)
  ! the category below each threshold, 0 below the support, k+1 at or above its top
  do t = 1, n_thresh
    if (thresholds(t) < y_mapping(1)) then
      idx(t) = 0
    else if (thresholds(t) >= y_mapping(k+1)) then
      idx(t) = k + 1
    else
      idx(t) = 1
      do j = k, 1, -1
        if (y_mapping(j) <= thresholds(t)) then
          idx(t) = j
          exit
        end if
      end do
    end if
  end do
  do l = 1, n_lambda
    xb_vec = matmul(newx, beta(:, l))
    do t = 1, n_thresh
      if (idx(t) == 0) then
        cdf_out(:, l, t) = 0.0_dp
      else if (idx(t) == k + 1) then
        cdf_out(:, l, t) = 1.0_dp
      else
        do i = 1, n
          cdf_out(i, l, t) = link_cdf(alpha(idx(t), l) + xb_vec(i), link_type)
        end do
      end if
    end do
  end do
end subroutine cpm_predict_cdf
