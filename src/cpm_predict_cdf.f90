! gfortran -shared -fPIC -O2 -o cpm_predict_cdf.so cpm_predict_cdf.f90
subroutine cpm_predict_cdf(n, p, k, newx, alpha, beta, y_mapping, &
                           link_type, thresholds, n_thresh, cdf_out, &
                           n_lambda)
  implicit none
  integer, intent(in) :: n, p, k, n_lambda, link_type, n_thresh
  double precision, intent(in) :: newx(n, p)
  double precision, intent(in) :: alpha(k, n_lambda)
  double precision, intent(in) :: beta(p, n_lambda)
  double precision, intent(in) :: y_mapping(k+1)
  double precision, intent(in) :: thresholds(n_thresh)
  double precision, intent(out) :: cdf_out(n, n_lambda, n_thresh)
  integer :: i, j, l, t, idx
  double precision :: eta, xb, c_val
  double precision, parameter :: inv_sqrt2 = 0.7071067811865475d0
  double precision, parameter :: inv_pi = 0.31830988618379067154d0
  double precision :: xb_vec(n)
  double precision :: cdf_array(k)
  cdf_out = 0.0d0
  select case(link_type)
  case(1) ! logistic
    do l = 1, n_lambda
      xb_vec = matmul(newx, beta(:, l))
      do i = 1, n
        xb = xb_vec(i)
        do j = 1, k
          eta = alpha(j, l) + xb
          cdf_array(j) = 1.0d0 - 1.0d0 / (1.0d0 + exp(-eta))
        end do
        do t = 1, n_thresh
          c_val = thresholds(t)
          if (c_val < y_mapping(1)) then
            cdf_out(i, l, t) = 0.0d0
          else if (c_val >= y_mapping(k+1)) then
            cdf_out(i, l, t) = 1.0d0
          else
            idx = k
            do j = k, 1, -1
              if (y_mapping(j) <= c_val) then
                idx = j
                exit
              end if
              idx = j - 1
            end do
            if (idx < 1) then
              cdf_out(i, l, t) = 0.0d0
            else
              cdf_out(i, l, t) = cdf_array(idx)
            end if
          end if
        end do
      end do
    end do
  case(2) ! probit
    do l = 1, n_lambda
      xb_vec = matmul(newx, beta(:, l))
      do i = 1, n
        xb = xb_vec(i)
        do j = 1, k
          eta = alpha(j, l) + xb
          cdf_array(j) = 1.0d0 - 0.5d0 * (1.0d0 + erf(eta * inv_sqrt2))
        end do
        do t = 1, n_thresh
          c_val = thresholds(t)
          if (c_val < y_mapping(1)) then
            cdf_out(i, l, t) = 0.0d0
          else if (c_val >= y_mapping(k+1)) then
            cdf_out(i, l, t) = 1.0d0
          else
            idx = k
            do j = k, 1, -1
              if (y_mapping(j) <= c_val) then
                idx = j
                exit
              end if
              idx = j - 1
            end do
            if (idx < 1) then
              cdf_out(i, l, t) = 0.0d0
            else
              cdf_out(i, l, t) = cdf_array(idx)
            end if
          end if
        end do
      end do
    end do
  case(3) ! log-log
    do l = 1, n_lambda
      xb_vec = matmul(newx, beta(:, l))
      do i = 1, n
        xb = xb_vec(i)
        do j = 1, k
          eta = alpha(j, l) + xb
          cdf_array(j) = 1.0d0 - exp(-exp(-eta))   ! 1-F = P(Y<=c) [loglog]
        end do
        do t = 1, n_thresh
          c_val = thresholds(t)
          if (c_val < y_mapping(1)) then
            cdf_out(i, l, t) = 0.0d0
          else if (c_val >= y_mapping(k+1)) then
            cdf_out(i, l, t) = 1.0d0
          else
            idx = k
            do j = k, 1, -1
              if (y_mapping(j) <= c_val) then
                idx = j
                exit
              end if
              idx = j - 1
            end do
            if (idx < 1) then
              cdf_out(i, l, t) = 0.0d0
            else
              cdf_out(i, l, t) = cdf_array(idx)
            end if
          end if
        end do
      end do
    end do
  case(4) ! complementary log-log
    do l = 1, n_lambda
      xb_vec = matmul(newx, beta(:, l))
      do i = 1, n
        xb = xb_vec(i)
        do j = 1, k
          eta = alpha(j, l) + xb
          cdf_array(j) = exp(-exp(eta))   ! 1-F = P(Y<=c) [cloglog]
        end do
        do t = 1, n_thresh
          c_val = thresholds(t)
          if (c_val < y_mapping(1)) then
            cdf_out(i, l, t) = 0.0d0
          else if (c_val >= y_mapping(k+1)) then
            cdf_out(i, l, t) = 1.0d0
          else
            idx = k
            do j = k, 1, -1
              if (y_mapping(j) <= c_val) then
                idx = j
                exit
              end if
              idx = j - 1
            end do
            if (idx < 1) then
              cdf_out(i, l, t) = 0.0d0
            else
              cdf_out(i, l, t) = cdf_array(idx)
            end if
          end if
        end do
      end do
    end do
  case(5) ! cauchit
    do l = 1, n_lambda
      xb_vec = matmul(newx, beta(:, l))
      do i = 1, n
        xb = xb_vec(i)
        do j = 1, k
          eta = alpha(j, l) + xb
          cdf_array(j) = 0.5d0 - atan(eta) * inv_pi
        end do
        do t = 1, n_thresh
          c_val = thresholds(t)
          if (c_val < y_mapping(1)) then
            cdf_out(i, l, t) = 0.0d0
          else if (c_val >= y_mapping(k+1)) then
            cdf_out(i, l, t) = 1.0d0
          else
            idx = k
            do j = k, 1, -1
              if (y_mapping(j) <= c_val) then
                idx = j
                exit
              end if
              idx = j - 1
            end do
            if (idx < 1) then
              cdf_out(i, l, t) = 0.0d0
            else
              cdf_out(i, l, t) = cdf_array(idx)
            end if
          end if
        end do
      end do
    end do
  end select
end subroutine cpm_predict_cdf
