! gradient of the negative log-likelihood at beta = 0, for lambda_max

subroutine cpm_lambda_max_grad(n, k, p, x, y, wt, link, alpha, u, salloc)

use cpm_links
use, intrinsic :: ISO_FORTRAN_ENV, only: dp => real64, int32
implicit none

integer(int32), intent(in)  :: n, k, p, link, y(n)
real(dp),       intent(in)  :: x(n, p), wt(n), alpha(k)
real(dp),       intent(out) :: u(p)
integer(int32), intent(out) :: salloc

integer(int32) :: i, j, l
real(dp) :: inv_n, inv_d
real(dp), allocatable :: lp(:), d(:), pdf1(:), pdf2(:), dpdf1(:), dpdf2(:)

allocate(lp(n), d(n), pdf1(n), pdf2(n), dpdf1(n), dpdf2(n), stat=salloc)
if (salloc /= 0) return

inv_n = 1.0_dp / real(n, dp)
lp = 0.0_dp
call compute_probs(n, k, alpha, lp, y, link, d, pdf1, pdf2, dpdf1, dpdf2)

do l = 1, p
    u(l) = 0.0_dp
    do i = 1, n
        j = y(i)
        inv_d = wt(i) / d(i)
        if (j == 0) then
            u(l) = u(l) + inv_d * pdf1(i) * x(i, l)
        else if (j == k) then
            u(l) = u(l) - inv_d * pdf1(i) * x(i, l)
        else
            u(l) = u(l) - inv_d * x(i, l) * (pdf1(i) - pdf2(i))
        end if
    end do
    u(l) = u(l) * inv_n
end do

end subroutine cpm_lambda_max_grad
