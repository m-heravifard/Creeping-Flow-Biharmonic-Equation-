!Auhor: Mohammad E. Heravifard

Program sphere_Flow 
    implicit none 
     
    integer :: i, j, i_max, j_max, Itr 
    real(8) :: r_max, a_sphere, pi, U, d_r, d_theta, error_criteria, d2_dr2, &
               d2_dtheta2, d_dtheta, sai_exact, cot_theta, d_tau, Rel_Coef 
    real(8), allocatable :: x(:,:), y(:,:), sai_old(:,:), sai_new(:,:), E2(:,:), &
                            theta(:,:), r(:,:), Lower2(:), Lower1(:), Diag(:), &
                            Upper1(:), Upper2(:), RHS(:), Mat_Sol(:) 
    real(8) :: d4_s_dr4, d4_s_dt4, d4_s_dt2_ds2, d3_s_dt_dr2, d3_s_dt2_dr, &
               d3_s_dt3, d2_s_dt2, d2_s_dr_dt, d1_s_dt, cot_t, u_r 

    pi = acos(-1d0)
     
    i_max = 51 
    j_max = 51 
     
    allocate( x(i_max,j_max), y(i_max,j_max), sai_old(i_max,j_max), &
              sai_new(i_max,j_max), E2(i_max,j_max), theta(i_max,j_max), &
              r(i_max,j_max), Lower2(i_max), Lower1(i_max), Diag(i_max), &
              Upper1(i_max), Upper2(i_max), RHS(i_max), Mat_Sol(i_max) ) 
     
    a_sphere = 1d0 
    r_max    = 5d0 
    U        = 1d0 
     
    d_tau = 1d-4 
     
    !=== grid generation ===! 
    do i = 1, i_max 
        do j = 1, j_max 
            theta(i,j) = 0d0 + (i-1)*(pi-0d0)/(i_max-1) 
            r(i,j)     = a_sphere + (j-1)*(r_max-a_sphere)/(j_max-1) 
            x(i,j)     = r(i,j)*cos(theta(i,j)) 
            y(i,j)     = r(i,j)*sin(theta(i,j)) 
        end do 
    end do 
    d_r     = r(1,2)   - r(1,1) 
    d_theta = theta(2,1) - theta(1,1) 
    !=======================! 
     
    !=== Initial condition ===! 
    do i = 1, i_max 
        do j = 1, j_max 
            sai_old(i,j) = 0.25d0*U*a_sphere*a_sphere * &
                           sin(theta(i,j))*sin(theta(i,j)) * &
                           ( a_sphere/r(i,j) - 3d0*r(i,j)/a_sphere + &
                             2d0*r(i,j)*r(i,j)/(a_sphere*a_sphere) ) 
            !E2(i,j) = 6*sin(theta(i,j))*sin(theta(i,j))/r(i,j) 
        end do 
    end do 
    sai_old(3:i_max-2, 3:j_max-2) = 0d0 
    !E2(2:i_max-1,2:j_max-1) = 0d0 
    sai_new = sai_old 
    !=========================! 
     
    !=== Main loop ===! 
    error_criteria = 10d0 
    Itr            = 0 
    do while (error_criteria > 1d-8) 
        Itr = Itr + 1 
         
        !! Jacobi (commented out)
        !do i=3,i_max-2 
        !    do j=3,j_max-2 
        !        d4_s_dr4 = ( sai_old(i,j-2) - 4*sai_old(i,j-1) + 6*0d0 - &
        !                     4*sai_old(i,j+1) + sai_old(i,j+2) ) / (d_r**4)  !1 
        !        d4_s_dt4 = ( sai_old(i-2,j) - 4*sai_old(i-1,j) + 6*0d0 - &
        !                     4*sai_old(i+1,j) + sai_old(i+2,j) ) / (d_theta**4) !11 
        !        d4_s_dt2_ds2 = ( (sai_old(i-1,j+1)-2*sai_old(i-1,j)+sai_old(i-1,j-1)) - &
        !                         2*(sai_old(i,j+1)-2*0d0+sai_old(i,j-1)) + &
        !                         (sai_old(i+1,j+1)-2*sai_old(i+1,j)+sai_old(i+1,j-1)) ) / &
        !                         (d_r*d_r*d_theta*d_theta) !4 
        !        d3_s_dt_dr2 = ( (sai_old(i+1,j+1)-2*sai_old(i+1,j)+sai_old(i+1,j-1)) - &
        !                        (sai_old(i-1,j+1)-2*sai_old(i-1,j)+sai_old(i-1,j-1)) ) / &
        !                        (2*d_theta*d_r*d_r) !3 
        !        d3_s_dt2_dr = ( (sai_old(i+1,j+1)-2*sai_old(i,j+1)+sai_old(i-1,j+1)) - &
        !                        (sai_old(i+1,j-1)-2*sai_old(i,j-1)+sai_old(i-1,j-1)) ) / &
        !                        (2*d_theta*d_theta*d_r) !9 
        !        d3_s_dt3 = (-sai_old(i-2,j) + 2*sai_old(i-1,j) - &
        !                     2*sai_old(i+1,j) + sai_old(i+2,j)) / (2*d_theta**3) !8 
        !        d2_s_dt2 = (sai_old(i+1,j)-2*0d0+sai_old(i-1,j)) / (d_theta**2) !5,10 
        !        d2_s_dr_dt = ((sai_old(i+1,j+1)-sai_old(i+1,j-1)) - &
        !                      (sai_old(i-1,j+1)-sai_old(i-1,j-1))) / &
        !                      (4*d_r*d_theta) !6 
        !        d1_s_dt = (sai_old(i+1,j)-sai_old(i-1,j)) / (2*d_theta) !2,7 
        !         
        !        cot_t = cos(theta(i,j))/sin(theta(i,j)) 
        !         
        !        sai_new(i,j) = - ( d4_s_dr4*r(i,j)**4 - 3*cot_t**3*d1_s_dt - &
        !             2*cot_t*r(i,j)**2*d3_s_dt_dr2 + 2*r(i,j)**2*d4_s_dt2_ds2 + &
        !             3*cot_t**2*d2_s_dt2 + 4*cot_t*r(i,j)*d2_s_dr_dt - 9*cot_t*d1_s_dt - &
        !             2*cot_t*d3_s_dt3 - 4*r(i,j)*d3_s_dt2_dr + 8*d2_s_dt2 + d4_s_dt4 ) / &
        !             ( (6/(d_r**4))*r(i,j)**4 + 2*r(i,j)**2*(4/(d_r*d_r*d_theta*d_theta)) + &
        !               3*cot_t**2*(-2/(d_theta**2)) + 8*(-2/(d_theta**2)) + &
        !               (6/(d_theta**4)) ) 
        !    end do 
        !end do  
         
        !! Gauss-Seidel (commented out)
        !do i=3,i_max-2 
        !    do j=3,j_max-2 
        !        d4_s_dr4=( sai_new(i,j-2)-4*sai_new(i,j-1)+6*0d0 - &
        !                   4*sai_old(i,j+1)+sai_old(i,j+2) )/(d_r**4) !1 
        !        d4_s_dt4=( sai_new(i-2,j)-4*sai_new(i-1,j)+6*0d0 - &
        !                   4*sai_old(i+1,j)+sai_old(i+2,j) )/(d_theta**4) !11 
        !        d4_s_dt2_ds2=( (sai_old(i-1,j+1)-2*sai_new(i-1,j)+sai_new(i-1,j-1)) - &
        !                       2*(sai_old(i,j+1)-2*0d0+sai_new(i,j-1)) + &
        !                       (sai_old(i+1,j+1)-2*sai_old(i+1,j)+sai_old(i+1,j-1)) ) / &
        !                       (d_r*d_r*d_theta*d_theta)!4 
        !        d3_s_dt_dr2=( (sai_old(i+1,j+1)-2*sai_old(i+1,j)+sai_old(i+1,j-1)) - &
        !                      (sai_old(i-1,j+1)-2*sai_new(i-1,j)+sai_new(i-1,j-1)) ) / &
        !                      (2*d_theta*d_r*d_r)!3 
        !        d3_s_dt2_dr=( (sai_old(i+1,j+1)-2*sai_old(i,j+1)+sai_old(i-1,j+1)) - &
        !                      (sai_old(i+1,j-1)-2*sai_new(i,j-1)+sai_new(i-1,j-1)) ) / &
        !                      (2*d_theta*d_theta*d_r)!9 
        !        d3_s_dt3=(-sai_new(i-2,j)+2*sai_new(i-1,j) - &
        !                  2*sai_old(i+1,j)+sai_old(i+2,j))/(2*d_theta**3) !8 
        !        d2_s_dt2=(sai_old(i+1,j)-2*0d0+sai_new(i-1,j))/(d_theta**2)!5,10 
        !        d2_s_dr_dt=((sai_old(i+1,j+1)-sai_old(i+1,j-1)) - &
        !                    (sai_old(i-1,j+1)-sai_new(i-1,j-1)))/(4*d_r*d_theta)!6 
        !        d1_s_dt=(sai_old(i+1,j)-sai_new(i-1,j))/(2*d_theta)!2,7 
        !         
        !        cot_t = cos(theta(i,j))/sin(theta(i,j)) 
        !         
        !        sai_new(i,j) = - ( d4_s_dr4*r(i,j)**4 - 3*cot_t**3*d1_s_dt - &
        !             2*cot_t*r(i,j)**2*d3_s_dt_dr2 + 2*r(i,j)**2*d4_s_dt2_ds2 + &
        !             3*cot_t**2*d2_s_dt2 + 4*cot_t*r(i,j)*d2_s_dr_dt - 9*cot_t*d1_s_dt - &
        !             2*cot_t*d3_s_dt3 - 4*r(i,j)*d3_s_dt2_dr + 8*d2_s_dt2 + d4_s_dt4 ) / &
        !             ( (6/(d_r**4))*r(i,j)**4 + 2*r(i,j)**2*(4/(d_r*d_r*d_theta*d_theta)) + &
        !               3*cot_t**2*(-2/(d_theta**2)) + 8*(-2/(d_theta**2)) + &
        !               (6/(d_theta**4)) ) 
        !    end do 
        !end do  
         
        ! Line method (active) 
        do j = 3, j_max-2 
            do i = 3, i_max-2 
                d4_s_dr4 = ( sai_old(i,j-2) - 4*sai_old(i,j-1) - &
                             4*sai_old(i,j+1) + sai_old(i,j+2) ) / (d_r**4) !1 
                d4_s_dt4 = 0d0 / (d_theta**4) !11 
                d4_s_dt2_ds2 = ( (sai_old(i-1,j+1)+sai_old(i-1,j-1)) - &
                                 2*(sai_old(i,j+1)+sai_old(i,j-1)) + &
                                 (sai_old(i+1,j+1)-2*sai_old(i+1,j)+sai_old(i+1,j-1)) ) / &
                                 (d_r*d_r*d_theta*d_theta) !4 
                d3_s_dt_dr2 = ( (sai_old(i+1,j+1)+sai_old(i+1,j-1)) - &
                                (sai_old(i-1,j+1)+sai_old(i-1,j-1)) ) / &
                                (2*d_theta*d_r*d_r) !3 
                d3_s_dt2_dr = ( (sai_old(i+1,j+1)-2*sai_old(i,j+1)+sai_old(i-1,j+1)) - &
                                (sai_old(i+1,j-1)-2*sai_old(i,j-1)+sai_old(i-1,j-1)) ) / &
                                (2*d_theta*d_theta*d_r) !9 
                d3_s_dt3  = 0d0 / (2*d_theta**3) !8 
                d2_s_dt2  = 0d0 / (d_theta**2)   !5,10 
                d2_s_dr_dt = ((sai_old(i+1,j+1)-sai_old(i+1,j-1)) - &
                              (sai_old(i-1,j+1)-sai_old(i-1,j-1))) / &
                              (4*d_r*d_theta) !6 
                d1_s_dt   = 0d0 / (2*d_theta)    !2,7 
                 
                cot_t = cos(theta(i,j))/sin(theta(i,j)) 
                 
                RHS(i) = - ( d4_s_dr4*r(i,j)**4 - 3*cot_t**3*d1_s_dt - &
                     2*cot_t*r(i,j)**2*d3_s_dt_dr2 + 2*r(i,j)**2*d4_s_dt2_ds2 + &
                     3*cot_t**2*d2_s_dt2 + 4*cot_t*r(i,j)*d2_s_dr_dt - 9*cot_t*d1_s_dt - &
                     2*cot_t*d3_s_dt3 - 4*r(i,j)*d3_s_dt2_dr + 8*d2_s_dt2 + d4_s_dt4 ) 
                Lower2(i) = - ( 1/(d_theta**4) - 1*(-2*cot_t)/(2*d_theta**3) ) 
                Lower1(i) = - ( -4/(d_theta**4) - &
                     2*(2*r(i,j)**2)/(d_r*d_r*d_theta*d_theta) + &
                     2*(-2*cot_t*r(i,j)**2)/(2*d_theta*d_r*d_r) + &
                     2*(-2*cot_t)/(2*d_theta**3) + &
                     1*(3*cot_t**2)/(d_theta**2) - &
                     1*(-3*cot_t**3-9*cot_t)/(2*d_theta) ) 
                Diag(i) = ( (6/(d_r**4))*r(i,j)**4 + &
                     2*r(i,j)**2*(4/(d_r*d_r*d_theta*d_theta)) + &
                     3*cot_t**2*(-2/(d_theta**2)) + &
                     8*(-2/(d_theta**2)) + (6/(d_theta**4)) ) 
                Upper1(i) = - ( -4/(d_theta**4) - &
                     2*(2*r(i,j)**2)/(d_r*d_r*d_theta*d_theta) - &
                     2*(-2*cot_t*r(i,j)**2)/(2*d_theta*d_r*d_r) - &
                     2*(-2*cot_t)/(2*d_theta**3) + &
                     1*(3*cot_t**2)/(d_theta**2) + &
                     1*(-3*cot_t**3-9*cot_t)/(2*d_theta) ) 
                Upper2(i) = - ( 1/(d_theta**4) + &
                     1*(-2*cot_t)/(2*d_theta**3) ) 
            end do 
             
            i = 1 
            RHS(i)    = sai_new(i,j) 
            Lower2(i) = 0d0 
            Lower1(i) = 0d0 
            Diag(i)   = 1d0 
            Upper1(i) = 0d0 
            Upper2(i) = 0d0 
             
            i = 2 
            RHS(i)    = sai_new(i,j) 
            Lower2(i) = 0d0 
            Lower1(i) = 0d0 
            Diag(i)   = 1d0 
            Upper1(i) = 0d0 
            Upper2(i) = 0d0 
             
            i = i_max 
            RHS(i)    = sai_new(i,j) 
            Lower2(i) = 0d0 
            Lower1(i) = 0d0 
            Diag(i)   = 1d0 
            Upper1(i) = 0d0 
            Upper2(i) = 0d0 
             
            i = i_max-1 
            RHS(i)    = sai_new(i,j) 
            Lower2(i) = 0d0 
            Lower1(i) = 0d0 
            Diag(i)   = 1d0 
            Upper1(i) = 0d0 
            Upper2(i) = 0d0 

            call PENTA(i_max, Lower2, Lower1, Diag, Upper1, Upper2, RHS, Mat_Sol) 
            sai_new(:,j) = Mat_Sol(:) 
        end do  
 
        ! Relaxation 
        Rel_Coef = 0.1d0 
        sai_new  = Rel_Coef*sai_new + (1d0-Rel_Coef)*sai_old 
         
        ! Error 
        error_criteria = 0d0 
        do i = 3, i_max-2 
            do j = 3, j_max-2 
                !sai_exact = 0.25d0*U*a_sphere*a_sphere*sin(theta(i,j))*sin(theta(i,j))* &
                !            (a_sphere/r(i,j) - 3*r(i,j)/a_sphere + &
                !             2*r(i,j)*r(i,j)/(a_sphere*a_sphere)) 
                !error_criteria = error_criteria + (sai_new(i,j)-sai_exact)**2 
                error_criteria = error_criteria + (sai_new(i,j)-sai_old(i,j))**2 
            end do 
        end do 
        error_criteria = sqrt(error_criteria/(i_max*j_max)) 
         
        if (mod(Itr,1) == 0) then 
            write(*,*) Itr, error_criteria 
        end if 
         
        ! updating 
        sai_old = sai_new 
         
    end do 
    !=================! 
     
    !=== plot results ===! 
    open(101, file="Contour.plt") 
    write(101,*) 'variables = "x", "y", "sai"' 
    write(101,*) 'zone i = ', i_max, ', j = ', j_max 
    do j = 1, j_max 
        do i = 1, i_max 
            write(101,*) x(i,j), y(i,j), sai_new(i,j) 
        end do 
    end do     
    !====================! 
     
    pause 
end program sphere_Flow
