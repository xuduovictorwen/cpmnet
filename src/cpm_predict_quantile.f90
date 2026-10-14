! conditional quantiles at tau levels
subroutine cpm_predict_quantile(n, p, k, newx, alpha, beta, y_mapping, &
                                link_type, tau_levels, n_tau, predictions, &
                                n_lambda, tol)
  use cpm_links
  use, intrinsic :: ISO_FORTRAN_ENV, only: dp => real64, int32
  implicit none
  integer(int32), intent(in) :: n, p, k, n_lambda, link_type, n_tau
  real(dp), intent(in) :: newx(n, p), alpha(k, n_lambda), beta(p, n_lambda)
  real(dp), intent(in) :: y_mapping(k+1), tau_levels(n_tau), tol
  real(dp), intent(out) :: predictions(n, n_lambda, n_tau)
  integer(int32) :: i, j, l, t, q_idx
  real(dp) :: xb, tau
  real(dp) :: xb_vec(n), cdf_array(k+1)
  predictions = 0.0_dp
  do l = 1, n_lambda
    xb_vec = matmul(newx, beta(:, l))
    do i = 1, n
      xb = xb_vec(i)
      do j = 1, k
        cdf_array(j) = link_cdf(alpha(j, l) + xb, link_type)
      end do
      cdf_array(k+1) = 1.0_dp
      do t = 1, n_tau
        tau = tau_levels(t)
        q_idx = 1
        do j = 1, k
          if (cdf_array(j) >= tau) exit
          q_idx = j + 1
        end do
        if (q_idx <= k .and. abs(cdf_array(q_idx) - tau) < tol) then
          predictions(i, l, t) = (y_mapping(q_idx) + y_mapping(q_idx + 1)) * 0.5_dp
        else
          predictions(i, l, t) = y_mapping(q_idx)
        end if
      end do
    end do
  end do
end subroutine cpm_predict_quantile
