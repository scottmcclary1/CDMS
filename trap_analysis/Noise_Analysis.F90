    SUBROUTINE noise_analysis(length_multiplier, load_time, dead_time, data_type)

    IMPLICIT NONE
    INTEGER :: i, j, k, m, n                        ! Counting variables for DO loops
    INTEGER :: unit                                 ! Dummy variable for opening files
    INTEGER :: stat                                 ! Error message when opening files
    INTEGER :: data_type                            ! Pointer for .dat or .txt files to analyze  
    INTEGER :: n_main                                ! Pointer for file section
    INTEGER :: channel                              ! Pointer for channel to analyze (A or B)
    INTEGER :: number_of_files                      ! Number of files to analyze
    INTEGER :: signal_end                           ! Length of test signal to add onto noise file 
    INTEGER :: dead_time                            ! Number of points of dead time due to A250 recovering from impulse noise of trap gate pulse
    INTEGER :: derivative_peak                      ! Inflection point of derivative
    INTEGER :: der_opt                              ! Derivative option for looking for positive or negative derivative
    
	INTEGER :: scan_increment           ! Length of window function when performing fft_scan_start
	INTEGER :: temp_signal_start        ! Starting point for window function when performing fft_scan_start to find end of signal
	INTEGER :: scan_incr_mult           ! Number of "scan_increments" in each cycle of length "2500. * REAL(length_multiplier - load_time)"
                    
    INTEGER :: length_multiplier                    ! Multiplier for length of each trapping event (multiple of 1 ms or 2500 datapoints)
    INTEGER :: load_time                            ! Time between end of one trappng event and beginning of next (multiple of 1 ms or 2500 datapoints)
    INTEGER :: expt_num                             ! Maximum number of experiments per file based on signal length (length_multiplier)
    
    INTEGER(2), DIMENSION (1:1000000) :: input_array 
    REAL(4), DIMENSION (1:1048576) :: y_input                ! Input datapoints from trap data 
    REAL(4), DIMENSION (1:1048576) :: adj_input              ! Input datapoints from trap data  
    REAL(4), DIMENSION (1:1048576) :: derivative             ! First derivative of y_input          
    REAL(4), DIMENSION (1:1048576) :: fft_filter_output      ! Trap data after passing through 10 kHz high-pass FFT filter 
        
    ! Parameters for FFT of length 1048576 
    ! ===================================
    REAL(4), DIMENSION (1:524288) :: fft_magnitude    
    REAL(4), DIMENSION (1:524288) :: fft_mag_sum
    REAL(4), DIMENSION (1:524288) :: fft_mag_avg
    REAL(4), DIMENSION (1:524288) :: fft_mag_diff_squared
    REAL(4), DIMENSION (1:524288) :: fft_mag_std_dev
    
    REAL(4), ALLOCATABLE :: y_input_array (:,:)
    
    CHARACTER(len=100) :: filename                      ! Dummy variable for opening files
    CHARACTER(len=8) :: fmt                             ! Format descriptor for file_number
    CHARACTER(len=10) :: file_number                     ! File number to read

    REAL(4), DIMENSION (1:1048576) :: input                ! Input datapoints from trap data
    INTEGER, DIMENSION (1:2,1:400) :: derivative_peaks    ! Array for location of derivative peaks from derivative of y_input
    INTEGER :: pos_peaks                ! Counter for positive peaks in derivative of y_input
    INTEGER :: neg_peaks                ! Counter for negative peaks in derivative of y_input
    INTEGER :: max_trap_time                        ! Length of period for trapping event (multiple of 1 ms or 2500 datapoints)
    


    ! User selected channel option or test signal
    PRINT *, "Press 1 for Channel A..."
    PRINT *, "Press 2 for Channel B..."
    READ *, channel
    WRITE (97,*) "Press 1 for Channel A..."
    WRITE (97,*) "Press 2 for Channel B..."
    WRITE (97,*) channel

    ! User selected option for number of files in each run
    PRINT *, "Enter number of noise files..."
    READ *, number_of_files 
    WRITE (97,*) "Enter number of noise files..."
    WRITE (97,*) number_of_files 
    
    CLOSE (UNIT=97, STATUS='KEEP')
          
    max_trap_time = length_multiplier - load_time 
!    expt_num = INT(400 / length_multiplier) - 1 
    
    
    ! This outer loop scans through window lengths
    DO  i = 1, INT(((2500. * REAL(length_multiplier - load_time)) - REAL(dead_time)) / 1000.)
    
        ! Initial length of signal on which to perform FFT to look for peaks 
        signal_end = ((i-1) * 1000) + 1000 + dead_time
        scan_increment = signal_end - dead_time    
        
        ! Find sum of all points to average
        fft_mag_sum = 0
        ! For each window length of the outer loop, this inner loop scans through all the noise files
        DO  j = 0, (number_of_files-1)
            
            WRITE (*,1000) "Now analyzing FFT of window size = ", scan_increment
            WRITE (*,2000) "File Number = ", (j+1), "/", number_of_files
            1000 FORMAT (A35, I6)
            2000 FORMAT (A14, I4, A1, I4)
            
            unit = 1        
            fmt = '(I5.5)'
            
            IF  (data_type == 1)    THEN
            
                IF  (channel == 1)  THEN
                    WRITE (file_number, fmt) j+10000
                    filename = 'chA'//TRIM(file_number)//'.dat'
                ELSE IF (channel == 2)  THEN
                    WRITE (file_number, fmt) j+10000
                    filename = 'chB'//TRIM(file_number)//'.dat'
                END IF
             
                PRINT *, "Now reading from ", filename    
                OPEN (UNIT=unit, FILE=filename, STATUS="old", ACTION="read", IOSTAT=stat, FORM="binary")             
            
                input_array = 0
                y_input = 0
                
                ! Read input file into array "y_input"
                READ (unit, END=8) input_array            ! Read filename values and store in array "y"
                8   CONTINUE   
            
                DO  k = 1, 1000000
                    y_input(k) = input_array(k)
                END DO
            
            ELSE IF (data_type == 2)   THEN
            
                IF  (channel == 1)  THEN
                    WRITE (file_number, fmt) j+10000
                    filename = 'chA'//TRIM(file_number)//'.txt'
                ELSE IF (channel == 2)  THEN
                    WRITE (file_number, fmt) j+10000
                    filename = 'chB'//TRIM(file_number)//'.txt'
                END IF         
             
                PRINT *, "Now reading from ", filename    
                OPEN (UNIT=unit, FILE=filename, STATUS="old", ACTION="read", IOSTAT=stat)           
            
                y_input = 0
                
                ! Read input file into array "y_input"
                READ (unit, *, END=9) y_input            ! Read filename values and store in array "y"
                9   CONTINUE              
             
            END IF
            
            CLOSE (unit)
            
            
            adj_input = 0
            
            DO  k = 1, 1000000
                IF  (MOD(k,2) == 1) THEN
                    adj_input(k) = y_input(k+1)
                ELSE IF (MOD(k,2) == 0) THEN
                    adj_input(k) = y_input(k-1)
                END IF
            END DO
            
            y_input = adj_input
          
            ! Calculate the derivative of y_input
            DO  k = 1, 1048576
                IF  (k == 1)    THEN
                    derivative (k) = (y_input(k+1) - y_input(k))/((k+1)-k)
                ELSE IF (k == 1048576)  THEN
                    derivative (k) = (y_input(k) - y_input(k-1))/(k-(k-1))
                ELSE           
                    derivative (k) = 0.5 * (((y_input(k+1) - y_input(k))/((k+1)-k)) + ((y_input(k) - y_input(k-1))/(k-(k-1)))) 
                END IF
            END DO        
                
            derivative_peaks = 0
            pos_peaks = 0
            neg_peaks = 0
            
            DO  k = 1, 1048576
                IF  (derivative(k) > 5000)  THEN
                    IF  (pos_peaks == 0)    THEN
                        pos_peaks = pos_peaks + 1
                        derivative_peaks(1,pos_peaks) = k
                    ELSE IF (k > (derivative_peaks(1,pos_peaks) + 50))  THEN
                        pos_peaks = pos_peaks + 1
                        derivative_peaks(1,pos_peaks) = k
                    END IF
                    
                    IF  (pos_peaks == 400)  THEN
                        GOTO 15
                    END IF
                    
                ELSE IF (derivative(k) < -5000)    THEN
                    IF  (neg_peaks == 0)    THEN                    
                        neg_peaks = neg_peaks + 1
                        derivative_peaks(2,neg_peaks) = k
                    ELSE IF (k > (derivative_peaks(2,neg_peaks) + 50))  THEN                 
                        neg_peaks = neg_peaks + 1
                        derivative_peaks(2,neg_peaks) = k
                    END IF
                    
                    IF  (neg_peaks == 400)  THEN
                        GOTO 15
                    END IF
                    
                END IF
            END DO
            
            
            IF  (derivative_peaks(1,1) <= (2500 * (400 - (INT(400 / length_multiplier) * length_multiplier)))) THEN        
                expt_num = INT(400 / length_multiplier)        
            ELSE IF (derivative_peaks(1,1) > (2500 * (400 - (INT(400 / length_multiplier) * length_multiplier)))) THEN  
                expt_num = INT(400 / length_multiplier) - 1
            END IF
                
                
            DO  k = 1, expt_num
                IF  (((derivative_peaks(2,(k+1)) - derivative_peaks(1,k)) > ((2500 * REAL(max_trap_time + (load_time / 2.))) - 50)) &
                    & .AND. ((derivative_peaks(2,(k+1)) - derivative_peaks(1,k)) < ((2500 * REAL(max_trap_time + (load_time / 2.))) + 50))) THEN
                    
                    derivative_peak = derivative_peaks(1,k) + (2500 * load_time / 2.)
                    input = 0
                    
                    DO  m = 0, ((2500 * max_trap_time) - 1)
                        input(m+1) = y_input(derivative_peak + m)
                    END DO
                      
                    CALL fft_filter(input, fft_filter_output)
                    input = fft_filter_output               
                    
                    n_main = 1           
        
                    CALL fftransform(input, n_main, signal_end, fft_magnitude, length_multiplier, dead_time) 
                                             
                    DO  m = 1, 524288 
                        fft_mag_sum(m) = fft_mag_sum(m) +  fft_magnitude(m)
                    END DO
                                 
                ELSE IF (((derivative_peaks(2,(k)) - derivative_peaks(1,k)) > ((2500 * REAL(max_trap_time + (load_time / 2.))) - 50)) &
                    & .AND. ((derivative_peaks(2,(k)) - derivative_peaks(1,k)) < ((2500 * REAL(max_trap_time + (load_time / 2.))) + 50))) THEN
                    
                    derivative_peak = derivative_peaks(1,k) + (2500 * load_time / 2.)
                    input = 0
                    
                    DO  m = 0, ((2500 * max_trap_time) - 1)
                        input(m+1) = y_input(derivative_peak + m)
                    END DO
                        
                    CALL fft_filter(input, fft_filter_output)
                    input = fft_filter_output               
                    
                    n_main = 1
                    
                    CALL fftransform(input, n_main, signal_end, fft_magnitude, length_multiplier, dead_time) 
                                             
                    DO  m = 1, 524288 
                        fft_mag_sum(m) = fft_mag_sum(m) +  fft_magnitude(m)
                    END DO
                                    
                END IF 
            END DO                  

        END DO
        
        
        ! Find average for the following variables over all noise files for this window length
        fft_mag_avg = 0
        DO  j = 1, 524288
            fft_mag_avg(j) = fft_mag_sum(j)/REAL(expt_num*number_of_files)
        END DO

        
        ! Find the sum of the square of the difference between average and individual value
        fft_mag_diff_squared = 0 
        DO  j = 0, (number_of_files-1)
        
            WRITE (*,1000) "Now analyzing FFT of window size = ", scan_increment
            WRITE (*,2000) "File Number = ", (j+1), "/", number_of_files
            
            unit = 1
        
            IF  (data_type == 1)    THEN
            
                IF  (channel == 1)  THEN
                    WRITE (file_number, fmt) j+10000
                    filename = 'chA'//TRIM(file_number)//'.dat'
                ELSE IF (channel == 2)  THEN
                    WRITE (file_number, fmt) j+10000
                    filename = 'chB'//TRIM(file_number)//'.dat'
                END IF
             
                PRINT *, "Now reading from ", filename    
                OPEN (UNIT=unit, FILE=filename, STATUS="old", ACTION="read", IOSTAT=stat, FORM="binary")             
            
                input_array = 0
                y_input = 0
                
                ! Read input file into array "y_input"
                READ (unit, END=10) input_array            ! Read filename values and store in array "y"
                10   CONTINUE   
            
                DO  k = 1, 1000000
                    y_input(k) = input_array(k)
                END DO
            
            ELSE IF (data_type == 2)   THEN
            
                IF  (channel == 1)  THEN
                    WRITE (file_number, fmt) j+10000
                    filename = 'chA'//TRIM(file_number)//'.txt'
                ELSE IF (channel == 2)  THEN
                    WRITE (file_number, fmt) j+10000
                    filename = 'chB'//TRIM(file_number)//'.txt'
                END IF         
             
                PRINT *, "Now reading from ", filename    
                OPEN (UNIT=unit, FILE=filename, STATUS="old", ACTION="read", IOSTAT=stat)           
            
                y_input = 0
                
                ! Read input file into array "y_input"
                READ (unit, *, END=11) y_input            ! Read filename values and store in array "y"
                11   CONTINUE              
             
            END IF
            
            CLOSE (unit)
            
            
            adj_input = 0
            
            DO  k = 1, 1000000
                IF  (MOD(k,2) == 1) THEN
                    adj_input(k) = y_input(k+1)
                ELSE IF (MOD(k,2) == 0) THEN
                    adj_input(k) = y_input(k-1)
                END IF
            END DO
            
            y_input = adj_input
          
            ! Calculate the derivative of y_input
            DO  k = 1, 1048576
                IF  (k == 1)    THEN
                    derivative (k) = (y_input(k+1) - y_input(k))/((k+1)-k)
                ELSE IF (k == 1048576)  THEN
                    derivative (k) = (y_input(k) - y_input(k-1))/(k-(k-1))
                ELSE           
                    derivative (k) = 0.5 * (((y_input(k+1) - y_input(k))/((k+1)-k)) + ((y_input(k) - y_input(k-1))/(k-(k-1)))) 
                END IF
            END DO        
            
            derivative_peaks = 0
            pos_peaks = 0
            neg_peaks = 0
            
            DO  k = 1, 1048576
                IF  (derivative(k) > 5000)  THEN
                    IF  (pos_peaks == 0)    THEN
                        pos_peaks = pos_peaks + 1
                        derivative_peaks(1,pos_peaks) = k
                    ELSE IF (k > (derivative_peaks(1,pos_peaks) + 50))  THEN
                        pos_peaks = pos_peaks + 1
                        derivative_peaks(1,pos_peaks) = k
                    END IF
                    
                    IF  (pos_peaks == 400)  THEN
                        GOTO 15
                    END IF
                    
                ELSE IF (derivative(k) < -5000)    THEN
                    IF  (neg_peaks == 0)    THEN                    
                        neg_peaks = neg_peaks + 1
                        derivative_peaks(2,neg_peaks) = k
                    ELSE IF (k > (derivative_peaks(2,neg_peaks) + 50))  THEN                 
                        neg_peaks = neg_peaks + 1
                        derivative_peaks(2,neg_peaks) = k
                    END IF
                    
                    IF  (neg_peaks == 400)  THEN
                        GOTO 15
                    END IF
                    
                END IF
            END DO
            
            
            IF  (derivative_peaks(1,1) <= (2500 * (400 - (INT(400 / length_multiplier) * length_multiplier)))) THEN        
                expt_num = INT(400 / length_multiplier)        
            ELSE IF (derivative_peaks(1,1) > (2500 * (400 - (INT(400 / length_multiplier) * length_multiplier)))) THEN  
                expt_num = INT(400 / length_multiplier) - 1
            END IF
                
                
            DO  k = 1, expt_num
                IF  (((derivative_peaks(2,(k+1)) - derivative_peaks(1,k)) > ((2500 * REAL(max_trap_time + (load_time / 2.))) - 50)) &
                    & .AND. ((derivative_peaks(2,(k+1)) - derivative_peaks(1,k)) < ((2500 * REAL(max_trap_time + (load_time / 2.))) + 50))) THEN
                    
                    derivative_peak = derivative_peaks(1,k) + (2500 * load_time / 2.)
                    input = 0
                    
                    DO  m = 0, ((2500 * max_trap_time) - 1)
                        input(m+1) = y_input(derivative_peak + m)
                    END DO
                           
                    CALL fft_filter(input, fft_filter_output)
                    input = fft_filter_output               
                    
                    n_main = 1           
        
                    CALL fftransform(input, n_main, signal_end, fft_magnitude, length_multiplier, dead_time) 
                                                
                    DO  m = 1, 524288                  
                        fft_mag_diff_squared(m) = fft_mag_diff_squared(m) + ((fft_magnitude(m) - fft_mag_avg(m))**2)
                    END DO
                             
                ELSE IF (((derivative_peaks(2,(k)) - derivative_peaks(1,k)) > ((2500 * REAL(max_trap_time + (load_time / 2.))) - 50)) &
                    & .AND. ((derivative_peaks(2,(k)) - derivative_peaks(1,k)) < ((2500 * REAL(max_trap_time + (load_time / 2.))) + 50))) THEN
                    
                    derivative_peak = derivative_peaks(1,k) + (2500 * load_time / 2.)
                    input = 0
                    
                    DO  m = 0, ((2500 * max_trap_time) - 1)
                        input(m+1) = y_input(derivative_peak + m)
                    END DO
                        
                    CALL fft_filter(input, fft_filter_output)
                    input = fft_filter_output               
                    
                    n_main = 1
                    
                    CALL fftransform(input, n_main, signal_end, fft_magnitude, length_multiplier, dead_time) 
                                           
                    DO  m = 1, 524288                  
                        fft_mag_diff_squared(m) = fft_mag_diff_squared(m) + ((fft_magnitude(m) - fft_mag_avg(m))**2)
                    END DO
                             
                END IF 
            END DO      
            
        END DO
        
        ! Find standard deviation for the following variables
        fft_mag_std_dev = 0
        DO  j = 1, 524288
            fft_mag_std_dev(j) = SQRT(fft_mag_diff_squared(j)/REAL((expt_num*number_of_files)-1))  
        END DO   
        
        IF  (scan_increment < 10000) THEN
            fmt = '(I4.4)'
        ELSE IF (scan_increment >= 10000 .AND. scan_increment < 100000)    THEN
            fmt = '(I5.5)'
        ELSE IF (scan_increment >= 100000)   THEN
            fmt = '(I6.6)'
        END IF
        
        unit = (i+100) 
        WRITE (file_number, fmt) scan_increment
        filename = 'FFT_avg_'//TRIM(file_number)//'_pt_win.txt' 
        CALL open_new_file (unit, filename)        
        
        WRITE (unit,777) "Frequency", "FFT_Avg_Mag", "FFT_Mag_Std_Dev"
        777 FORMAT (TR4, A9, TR19, A11, TR18, A15)
             
!        DO  j = 2098, 71304     ! Write from 5 kHz to 170 kHz 
        DO  j = 840, 71304     ! Write from 2 kHz to 170 kHz
            WRITE (unit,888) REAL((j-1)*40.E6/2**24), fft_mag_avg(j), fft_mag_std_dev(j)
            888 FORMAT (F12.5, TR10, F20.5 TR10, F20.5)
        END DO        
        
        CLOSE (unit, STATUS = "KEEP")
        
    END DO
        
    15  CONTINUE
    
    END SUBROUTINE noise_analysis