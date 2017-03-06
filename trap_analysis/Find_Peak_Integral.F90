    SUBROUTINE find_peak_integral(peak_start, peak_end, fft_magnitude, peak_integral)
    
    IMPLICIT NONE
    REAL, PARAMETER :: freq_res = 40.E6 / 2**21    ! Number of bits (2^16) in the ADC software
        
    INTEGER :: peak_start                     ! Pointer for peak start
    INTEGER :: peak_end                      ! Pointer for peak end 
    INTEGER :: i, j, k                            ! Loop variables

    REAL(4) :: start                    ! Dummy variable for beginning of peak   
    REAL(4) :: start_less_one           ! Dummy variable for beginning of peak minus one  
    REAL(4) :: last                     ! Dummy variable for end of peak    
    REAL(4) :: last_plus_one            ! Dummy variable for end of peak plus one 
    REAL(4) :: peak_integral            ! Area of peak
    
    REAL(4), DIMENSION (1:32768) :: fft_magnitude      ! DFT magnitude output
    
    
           
    
    start = fft_magnitude(peak_start)
    start_less_one = fft_magnitude(peak_start - 1)
    
    ! If first point falls on the right side of a local peak max, keep moving to the left until move past the local max
    IF  (start_less_one > start)   THEN                    
        DO WHILE    (start_less_one > start) 
                    peak_start = peak_start - 1
                    start = fft_magnitude(peak_start)
                    start_less_one = fft_magnitude(peak_start - 1)
        END DO
    END IF        

    ! Find start of peak           
    DO WHILE    (start >= start_less_one)
                peak_start = peak_start - 1
                start = fft_magnitude(peak_start)
                start_less_one = fft_magnitude(peak_start - 1)
    END DO
    
    peak_start = peak_start + 1 
                     
                     
    
    last = fft_magnitude(peak_end)
    last_plus_one = fft_magnitude(peak_end + 1)
    
    ! If first point falls on the left side of a local peak max, keep moving to the right until move past the local max
    IF  (last_plus_one > last)  THEN                    
        DO WHILE    (last_plus_one > last)
                    last = fft_magnitude(peak_end)
                    last_plus_one = fft_magnitude(peak_end + 1)
                    peak_end = peak_end + 1                
        END DO
    END IF                
    
    ! Find end of peak
    DO WHILE    (last >= last_plus_one)
                last = fft_magnitude(peak_end)
                last_plus_one = fft_magnitude(peak_end + 1)
                peak_end = peak_end + 1                
    END DO
             
    
    ! Integrate peak from peak_start to peak_end
    peak_integral = 0
    DO  i = peak_start, peak_end
        peak_integral = peak_integral + ( fft_magnitude(i) * freq_res )
    END DO 
               
                
    END SUBROUTINE find_peak_integral