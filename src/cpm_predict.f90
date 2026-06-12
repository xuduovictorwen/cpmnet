! Compile: gfortran -shared -fPIC -O2 -o cpm_predict.so cpm_predict.f90

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
  double precision :: eta, xb, cdf_curr, cdf_prev, pmf_j, mean_pred
  double precision, parameter :: inv_sqrt2 = 0.7071067811865475d0
  double precision, parameter :: inv_pi = 0.31830988618379067154d0
  
  double precision :: xb_vec(n)
  double precision :: cdf_array(k+1)  ! Only needed for median
  
  predictions = 0.0d0
  
  ! Branch on link_type OUTSIDE all loops
  select case(link_type)
  
  
  case(1)  ! Logistic
  
    do l = 1, n_lambda
      ! Use matmul for matrix-vector product
      xb_vec = matmul(newx, beta(:, l))
      
      if (pred_type == 1) then
        ! MEDIAN
        do i = 1, n
          xb = xb_vec(i)
          do j = 1, k
            eta = alpha(j, l) + xb
            cdf_array(j) = 1.0d0 - 1.0d0 / (1.0d0 + exp(-eta))
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
        ! MEAN - compute incrementally without storing PMF
        do i = 1, n
          xb = xb_vec(i)
          mean_pred = 0.0d0
          cdf_prev = 0.0d0
          
          do j = 1, k
            eta = alpha(j, l) + xb
            cdf_curr = 1.0d0 - 1.0d0 / (1.0d0 + exp(-eta))
            pmf_j = cdf_curr - cdf_prev
            mean_pred = mean_pred + y_mapping(j) * pmf_j
            cdf_prev = cdf_curr
          end do
          ! Last category: pmf = 1 - cdf_prev
          mean_pred = mean_pred + y_mapping(k+1) * (1.0d0 - cdf_prev)
          predictions(i, l) = mean_pred
        end do
      end if
    end do
    
  
  case(2)  ! Probit
  
    do l = 1, n_lambda
      xb_vec = matmul(newx, beta(:, l))
      
      if (pred_type == 1) then
        do i = 1, n
          xb = xb_vec(i)
          do j = 1, k
            eta = alpha(j, l) + xb
            cdf_array(j) = 1.0d0 - 0.5d0 * (1.0d0 + erf(eta * inv_sqrt2))
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
        do i = 1, n
          xb = xb_vec(i)
          mean_pred = 0.0d0
          cdf_prev = 0.0d0
          
          do j = 1, k
            eta = alpha(j, l) + xb
            cdf_curr = 1.0d0 - 0.5d0 * (1.0d0 + erf(eta * inv_sqrt2))
            pmf_j = cdf_curr - cdf_prev
            mean_pred = mean_pred + y_mapping(j) * pmf_j
            cdf_prev = cdf_curr
          end do
          mean_pred = mean_pred + y_mapping(k+1) * (1.0d0 - cdf_prev)
          predictions(i, l) = mean_pred
        end do
      end if
    end do
    
  
  case(3)  ! Log-log
  
    do l = 1, n_lambda
      xb_vec = matmul(newx, beta(:, l))
      
      if (pred_type == 1) then
        do i = 1, n
          xb = xb_vec(i)
          do j = 1, k
            eta = alpha(j, l) + xb
            cdf_array(j) = 1.0d0 - exp(-exp(-eta))   ! 1-F = P(Y<=c) [loglog]
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
        do i = 1, n
          xb = xb_vec(i)
          mean_pred = 0.0d0
          cdf_prev = 0.0d0
          
          do j = 1, k
            eta = alpha(j, l) + xb
            cdf_curr = 1.0d0 - exp(-exp(-eta))   ! 1-F = P(Y<=c) [loglog]
            pmf_j = cdf_curr - cdf_prev
            mean_pred = mean_pred + y_mapping(j) * pmf_j
            cdf_prev = cdf_curr
          end do
          mean_pred = mean_pred + y_mapping(k+1) * (1.0d0 - cdf_prev)
          predictions(i, l) = mean_pred
        end do
      end if
    end do
    
  
  case(4)  ! Complementary log-log
  
    do l = 1, n_lambda
      xb_vec = matmul(newx, beta(:, l))
      
      if (pred_type == 1) then
        do i = 1, n
          xb = xb_vec(i)
          do j = 1, k
            eta = alpha(j, l) + xb
            cdf_array(j) = exp(-exp(eta))   ! 1-F = P(Y<=c) [cloglog]
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
        do i = 1, n
          xb = xb_vec(i)
          mean_pred = 0.0d0
          cdf_prev = 0.0d0
          
          do j = 1, k
            eta = alpha(j, l) + xb
            cdf_curr = exp(-exp(eta))   ! 1-F = P(Y<=c) [cloglog]
            pmf_j = cdf_curr - cdf_prev
            mean_pred = mean_pred + y_mapping(j) * pmf_j
            cdf_prev = cdf_curr
          end do
          mean_pred = mean_pred + y_mapping(k+1) * (1.0d0 - cdf_prev)
          predictions(i, l) = mean_pred
        end do
      end if
    end do
    
  
  case(5)  ! Cauchit
  
    do l = 1, n_lambda
      xb_vec = matmul(newx, beta(:, l))
      
      if (pred_type == 1) then
        do i = 1, n
          xb = xb_vec(i)
          do j = 1, k
            eta = alpha(j, l) + xb
            cdf_array(j) = 0.5d0 - atan(eta) * inv_pi
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
        do i = 1, n
          xb = xb_vec(i)
          mean_pred = 0.0d0
          cdf_prev = 0.0d0
          
          do j = 1, k
            eta = alpha(j, l) + xb
            cdf_curr = 0.5d0 - atan(eta) * inv_pi
            pmf_j = cdf_curr - cdf_prev
            mean_pred = mean_pred + y_mapping(j) * pmf_j
            cdf_prev = cdf_curr
          end do
          mean_pred = mean_pred + y_mapping(k+1) * (1.0d0 - cdf_prev)
          predictions(i, l) = mean_pred
        end do
      end if
    end do
    
  end select
  
end subroutine cpm_predict