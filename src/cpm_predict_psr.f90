! conditional CDF at each observation's own category
subroutine cpm_predict_psr(n, p, k, newx, alpha, beta, link_type, &
                           obs_idx, cdf_out, n_lambda)
  use cpm_links
  use, intrinsic :: ISO_FORTRAN_ENV, only: dp => real64, int32
  implicit none
  integer(int32), intent(in) :: n, p, k, n_lambda, link_type, obs_idx(n)
  real(dp), intent(in) :: newx(n, p), alpha(k, n_lambda), beta(p, n_lambda)
  real(dp), intent(out) :: cdf_out(n, n_lambda)
  integer(int32) :: i, l
  real(dp) :: xb_vec(n)
  do l = 1, n_lambda
    xb_vec = matmul(newx, beta(:, l))
    do i = 1, n
      if (obs_idx(i) <= 0) then
        cdf_out(i, l) = 0.0_dp
      else if (obs_idx(i) > k) then
        cdf_out(i, l) = 1.0_dp
      else
        cdf_out(i, l) = link_cdf(alpha(obs_idx(i), l) + xb_vec(i), link_type)
      end if
    end do
  end do
end subroutine cpm_predict_psr
