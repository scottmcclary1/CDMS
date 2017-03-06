    SUBROUTINE fft_of_selected_section(channel, length_multiplier, load_time, data_type)
    
    IMPLICIT NONE    
    REAL(4) :: sigma_main                           ! Sigma value for Gaussian window function 
    REAL(4), PARAMETER :: pi = 3.14159265                      
       
    INTEGER(2), DIMENSION (1:1000000) :: input_array 
    REAL(4), DIMENSION (1:1048576) :: y_input            ! I/O y-axis datapoints 
    REAL(4), DIMENSION (1:1048576) :: adj_input         ! Input datapoints from trap data
    REAL(4), DIMENSION (1:1048576) :: derivative        ! First derivative of y_input
    REAL(4), DIMENSION (1:1048576) :: raw_input          ! Data of selected file & section moved to beginning of file  
    REAL(4), DIMENSION (1:1048576) :: raw_input_window   ! Windowed data of selected file & section moved to beginning of file  
    REAL(4), DIMENSION (1:1048576) :: filtered_input     ! Filtered data of selected file & section moved to beginning of file  
    REAL(4), DIMENSION (1:1048576) :: y_filter           ! Output from 10 kHz high-pass filter                         
    REAL(4), DIMENSION (1:1048576) :: fft_filter_output  ! Output from 1 kHz FFT filter
    REAL(4), DIMENSION (1:524288) :: mag_full            ! Magnitude of DFT for complete dataset
          
    INTEGER :: i, j, k
    INTEGER :: unit                     ! Dummy variable for opening files
    INTEGER :: stat                     ! Error message when opening files
    INTEGER :: channel                  ! Pointer for channel to analyze (A, B, or test signal) 
    INTEGER :: data_type                ! Pointer for .dat or .txt files to analyze      
    INTEGER :: n_main                   ! Pointer for file section
    INTEGER :: section
    INTEGER :: window_function          ! Indicator for application of apodization function 
    INTEGER :: start_of_signal          ! Length to cut off at beginning of signal before FFT
    INTEGER :: end_of_signal            ! Length of signal on which to perform FFT
    INTEGER :: derivative_peak          ! Inflection point of derivative
    INTEGER :: der_opt                              ! Derivative option for looking for positive or negative derivative
    
    REAL(4) :: charge                   ! Average charge value (windowed FFT)
    
    INTEGER :: length_multiplier        ! Multiplier for length of each trapping event (multiple of 1 ms or 2500 datapoints)
    INTEGER :: load_time                ! Time between end of one trappng event and beginning of next (multiple of 1 ms or 2500 datapoints)
      
    INTEGER :: single_event             ! Counter for number of trapping events with a single ion trapped    
    INTEGER :: multiple_event           ! Counter for number of trapping events with multiple ions trapped
    
    CHARACTER(len=50) :: filename       ! Dummy variable for opening files
    
    REAL(4) :: low_mass_cutoff          ! Low mass cutoff for frequency information
    REAL(4) :: high_mass_cutoff         ! High mass cutoff for frequency information
    INTEGER :: sig_cutoff               ! Signal length cutoff (in datapoints) for frequency information
              
    REAL(4), DIMENSION (1:1048576) :: input                ! Input datapoints from trap data
    INTEGER, DIMENSION (1:2,1:400) :: derivative_peaks    ! Array for location of derivative peaks from derivative of y_input
    INTEGER :: pos_peaks                ! Counter for positive peaks in derivative of y_input
    INTEGER :: neg_peaks                ! Counter for negative peaks in derivative of y_input
    INTEGER :: expt_num                             ! Maximum number of experiments per file based on signal length (length_multiplier)
    INTEGER :: max_trap_time                        ! Length of period for trapping event (multiple of 1 ms or 2500 datapoints)
    
    
      
          
          
!    ! Open file to store peak data for all signals found
!    unit = 31
!    filename = "test_output.txt"
!    CALL open_new_file (unit, filename)
                
    ! Open output files
    unit = 101
    filename = "mag_full_windowed.txt"
    CALL open_new_file (unit, filename)
    
    unit = 102
    filename = "raw_input.txt"
    CALL open_new_file (unit, filename)
    
    unit = 103
    filename = "filtered_input.txt"
    CALL open_new_file (unit, filename)
      
    unit = 104
    filename = "raw_window_input.txt"
    CALL open_new_file (unit, filename)
    
    unit = 105
    filename = "filter_window_input.txt"
    CALL open_new_file (unit, filename)
        
    unit = 107
    filename = "mag_full_no_window.txt"
    CALL open_new_file (unit, filename)
    
    unit = 108
    filename = "harmonic_analysis.txt"
    CALL open_new_file (unit, filename)
        
    ! Open input file
    unit = 1
    PRINT *, "Enter file name (chA#####.txt or chA#####.dat)..."
    READ *, filename        
    WRITE (97,*) "Enter file name (chA#####.txt or chA#####.dat)..."
    WRITE (97,*) filename
    
    ! Select file section        
    PRINT *, "Enter section # on which to perform DFT..."
    READ *, n_main 
    WRITE (97,*) "Enter section # on which to perform DFT..."
    WRITE (97,*) n_main
    
    section = n_main
    
!    PRINT *, "Enter beginning of signal..."
!    PRINT *, "Minimum should be 5000..."
!    READ *, start_of_signal 
!    WRITE (97,*) "Enter beginning of signal..."
!    WRITE (97,*) "Minimum should be 5000..."
!    WRITE (97,*) start_of_signal 
    start_of_signal = 5000                
    
!    PRINT *, "Enter end of signal..."
!    PRINT *, ""
!    PRINT *, "Total # Datapoints"
!    PRINT *, "=================="
!    PRINT *, " 2500 =  1 ms"
!    PRINT *, " 5000 =  2 ms"
!    PRINT *, "20000 =  8 ms"
!    PRINT *, "40000 = 16 ms"
!    PRINT *, "72500 = 29 ms"
!    READ *, end_of_signal
!    WRITE (97,*) "Enter end of signal..."
!    WRITE (97,*) ""
!    WRITE (97,*) "Total # Datapoints"
!    WRITE (97,*) "=================="
!    WRITE (97,*) " 2500 =  1 ms"
!    WRITE (97,*) " 5000 =  2 ms"
!    WRITE (97,*) "20000 =  8 ms"
!    WRITE (97,*) "40000 = 16 ms"
!    WRITE (97,*) "72500 = 29 ms"
!    WRITE (97,*) end_of_signal
    end_of_signal = length_multiplier*2500
    
    ! Write name of file currently being analzyed to output files
    WRITE (26,300) filename
    300 FORMAT (A12)         
    
    max_trap_time = length_multiplier - load_time

    IF  (data_type == 1)    THEN
    
        PRINT *, "Now reading from ", filename    
        OPEN (UNIT=unit, FILE=filename, STATUS="old", ACTION="read", IOSTAT=stat, FORM="binary")             
    
        input_array = 0
        y_input = 0
        
        ! Read input file into array "y_input"
        READ (unit, END=8) input_array            ! Read filename values and store in array "y"
        8   CONTINUE   
    
        DO  j = 1, 1000000
            y_input(j) = input_array(j)
        END DO
    
    ELSE IF (data_type == 2)   THEN
                 
        PRINT *, "Now reading from ", filename    
        OPEN (UNIT=unit, FILE=filename, STATUS="old", ACTION="read", IOSTAT=stat)           
    
        y_input = 0
        
        ! Read input file into array "y_input"
        READ (unit, *, END=9) y_input            ! Read filename values and store in array "y"
        9   CONTINUE              
     
    END IF
    
    CLOSE(1)
                           
    adj_input = 0
    
    DO  j = 1, 1000000
        IF  (MOD(j,2) == 1) THEN
            adj_input(j) = y_input(j+1)
        ELSE IF (MOD(j,2) == 0) THEN
            adj_input(j) = y_input(j-1)
        END IF
    END DO
    
    y_input = adj_input
        
    ! Calculate the derivative of y_input
    DO  j = 1, 1048576
        IF  (j == 1)    THEN
            derivative (j) = (y_input(j+1) - y_input(j))/((j+1)-j)
        ELSE IF (j == 1048576)  THEN
            derivative (j) = (y_input(j) - y_input(j-1))/(j-(j-1))
        ELSE           
            derivative (j) = 0.5 * (((y_input(j+1) - y_input(j))/((j+1)-j)) + ((y_input(j) - y_input(j-1))/(j-(j-1)))) 
        END IF
    END DO        


   derivative_peaks = 0 
   pos_peaks = 0
   neg_peaks = 0
    
    DO  j = 1, 1048576
        IF  (derivative(j) > 5000)  THEN
            IF  (pos_peaks == 0)    THEN
                pos_peaks = pos_peaks + 1
                derivative_peaks(1,pos_peaks) = j
            ELSE IF (j > (derivative_peaks(1,pos_peaks) + 50))  THEN
                pos_peaks = pos_peaks + 1
                derivative_peaks(1,pos_peaks) = j
            END IF
            
            IF  (pos_peaks == 400)  THEN
                GOTO 15
            END IF
                    
        ELSE IF (derivative(j) < -5000)    THEN
            IF  (neg_peaks == 0)    THEN                    
                neg_peaks = neg_peaks + 1
                derivative_peaks(2,neg_peaks) = j
            ELSE IF (j > (derivative_peaks(2,neg_peaks) + 50))  THEN                 
                neg_peaks = neg_peaks + 1
                derivative_peaks(2,neg_peaks) = j
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

    IF  (((derivative_peaks(2,(n_main+1)) - derivative_peaks(1,n_main)) > ((2500 * REAL(max_trap_time + (load_time / 2.))) - 50)) &
        & .AND. ((derivative_peaks(2,(n_main+1)) - derivative_peaks(1,n_main)) < ((2500 * REAL(max_trap_time + (load_time / 2.))) + 50))) THEN
        
        derivative_peak = derivative_peaks(1,n_main) + (2500 * load_time / 2.)
        input = 0
        
        DO  k = 0, ((2500 * max_trap_time) - 1)
            input(k+1) = y_input(derivative_peak + k)
        END DO
        
    ELSE IF (((derivative_peaks(2,(n_main)) - derivative_peaks(1,n_main)) > ((2500 * REAL(max_trap_time + (load_time / 2.))) - 50)) &
        & .AND. ((derivative_peaks(2,(n_main)) - derivative_peaks(1,n_main)) < ((2500 * REAL(max_trap_time + (load_time / 2.))) + 50))) THEN
        
        derivative_peak = derivative_peaks(1,n_main) + (2500 * load_time / 2.)
        input = 0
        
        DO  k = 0, ((2500 * max_trap_time) - 1)
            input(k+1) = y_input(derivative_peak + k)
        END DO
                           
    END IF 
    
    y_input = input 
     
    n_main = 1            
                    
    ! Move file section to beginning of new array                
    DO  i = (1 + ((n_main-1)*2500*length_multiplier) + start_of_signal), (((n_main-1)*2500*length_multiplier) + end_of_signal)
        raw_input(i - ((n_main-1)*2500*length_multiplier) - start_of_signal) = y_input(i)
    END DO
 
    ! Make sure everything after signal of interest is zero
    DO  i = (end_of_signal - start_of_signal + 1), 1048576
        raw_input(i) = 0
    END DO 

    window_function = 1    
    ! Apply a window function if selected
    IF  (window_function == 1) THEN            
        !   Gaussian Window
        sigma_main = 0.45
        DO  i = 0, (end_of_signal - (start_of_signal + 1))
            raw_input_window(i+1) = ( EXP( -0.5 * ((  (i - ((REAL(end_of_signal-(start_of_signal + 1)))/2)) / (sigma_main * ((REAL(end_of_signal-(start_of_signal + 1)))/2))  )**2)   )  ) * raw_input(i+1)
        END DO                
    ELSE IF (window_function == 2) THEN                
        !   Hann Window
        DO  i = 0, (end_of_signal - (start_of_signal + 1))
            IF  (i == 0)   THEN
                raw_input_window(i+1) = ( 0.5 + (0.5 * cos ((2*pi*i)/REAL(end_of_signal-(start_of_signal + 1)))) ) * raw_input(i+1)  
            ELSE
                raw_input_window(i+1) = ( 0.5 - (0.5 * cos ((2*pi*i)/REAL(end_of_signal-(start_of_signal + 1)))) ) * raw_input(i+1)
            END IF
        END DO                        
    ELSE IF (window_function == 3) THEN                
        !   Hamming Window
        DO  i = 0, (end_of_signal - (start_of_signal + 1))
            IF  (i == 0)   THEN
                raw_input_window(i+1) = ( 0.54 + (0.46 * cos ((2*pi*i)/REAL(end_of_signal-(start_of_signal + 1)))) ) * raw_input(i+1)  
            ELSE
                raw_input_window(i+1) = ( 0.54 - (0.46 * cos ((2*pi*i)/REAL(end_of_signal-(start_of_signal + 1)))) ) * raw_input(i+1)
            END IF
        END DO  
    ELSE IF  (window_function == 4) THEN 
        raw_input_window = raw_input                
    END IF
    

    ! Apply a filter in accordance with user-defined input
    CALL fft_filter(y_input, fft_filter_output)
    y_input = fft_filter_output
    
    ! Move filtered file section to new array                
    DO  i = (1 + ((n_main-1)*2500*length_multiplier) + start_of_signal), (((n_main-1)*2500*length_multiplier) + end_of_signal)
        filtered_input(i - ((n_main-1)*2500*length_multiplier) - start_of_signal) = y_input(i)
    END DO
 
    ! Make sure everything after signal of interest is zero
    DO  i = (end_of_signal - start_of_signal + 1), 1048576
        filtered_input(i) = 0
    END DO 
    
    ! Analyze data via peakfinder_VI subroutine
     CALL peakfinder_VI(filename, y_input, n_main, channel, charge, length_multiplier, (length_multiplier-load_time), start_of_signal, single_event, multiple_event, low_mass_cutoff, high_mass_cutoff, sig_cutoff, section)
    
    ! Perform FFT on selected windowed section 
!    start_of_signal = 150*2500     ! To do FT of particular part of particular section, enable and customize these two lines
!    end_of_signal = 245*2500       ! Also, see commented out stuff near beginning of section_fft subroutine
    CALL section_fft(filtered_input, n_main, start_of_signal, end_of_signal, mag_full, window_function)  
    WRITE (101, 3000)  mag_full 
    
    ! Perform FFT on selected section w/o a window
    window_function = 4
    CALL section_fft(filtered_input, n_main, start_of_signal, end_of_signal, mag_full, window_function)
    
    ! Write output data to files        
    WRITE (102, 3000) raw_input
    WRITE (103, 3000) filtered_input
    WRITE (104, 3000) raw_input_window 
    WRITE (107, 3000) mag_full
    3000 FORMAT (F25.5)              
    
    15 CONTINUE           
    
    END SUBROUTINE fft_of_selected_section
    
    
    
    
    
    ! FFT subroutine for looking at specific file section
    SUBROUTINE section_fft(fft_input, n_fft, signal_start, signal_end, fft_mag, window_fn)
    
    USE mkl_dfti
    
    IMPLICIT NONE
    
    TYPE (DFTI_DESCRIPTOR), POINTER :: My_Desc1_Handle  
    INTEGER :: Status                                       ! Sets DFT Parameters
    INTEGER :: istat                                        ! Allocate status
    INTEGER, INTENT(IN) :: window_fn                        ! User supplied input for applying a window function or not
    INTEGER :: i,k
    INTEGER, INTENT(IN) :: n_fft                            ! Pointer for file section
    INTEGER, INTENT(IN) :: signal_start                     ! Pointer for start of section                     
    INTEGER, INTENT(IN) :: signal_end                       ! Length of signal on which to perform DFT
    REAL(4) :: filter_sum                                   ! Rolling average sum to subtract from datapoints
    REAL(4) :: sigma_fft                                    ! Sigma value for Gaussian window function 
!    REAL(4) :: alpha_fft                                    ! Alpha value for Hann-Poisson window function
    REAL(4), PARAMETER :: pi = 3.14159265   
    REAL(4), DIMENSION (1:1048576), INTENT(IN) :: fft_input ! Input values after applying rolling average filter  
    CHARACTER :: filename
    INTEGER :: unit
        
    ! Parameters for FFT of length 1048576 
    ! ===================================
    INTEGER, PARAMETER :: npts = 1048576                    ! Number of datapoints    
    REAL(4), DIMENSION (1:524288), INTENT(OUT) :: fft_mag  ! DFT magnitude values   
    REAL(4), DIMENSION (1:1048576) :: fft_input_window      ! Input values after applying rolling average filter & windowing 
    REAL(4), DIMENSION (1:1048576) :: fft_output            ! DFT output values   
!    REAL(4), DIMENSION (1:524288) :: phase                  ! DFT phase values 

    
    IF  (window_fn == 1) THEN        
        !   Gaussian Window
        sigma_fft = 0.45
        DO  i = 0, (signal_end - (signal_start + 1))
            fft_input_window(i+1) = ( EXP( -0.5 * ((  (i - ((REAL(signal_end-(signal_start + 1)))/2)) / (sigma_fft * ((REAL(signal_end-(signal_start + 1)))/2))  )**2)   )  ) * fft_input(i+1)
        END DO 
!        k = 0   ! To do FFT of particular part of the section, enable this line and following loop, and disable above loop
!                ! Also enable and customize two commented out lines starting with "start_of_signal" and "end_of_signal" just before section_fft called
!        DO  i = signal_start, (signal_end - 1)     ! Get rid of this loop when done testing. The above one is correct.
!            k = k + 1
!            fft_input_window(k) = ( EXP( -0.5 * ((  (k - ((REAL(signal_end-(signal_start + 1)))/2)) / (sigma_fft * ((REAL(signal_end-(signal_start + 1)))/2))  )**2)   )  ) * fft_input(i+1)
!        END DO  
              
!        DO  i = 0, (signal_end - (signal_start + 1))
!            raw_input_window(i+1) =  ( EXP( -0.5 * ((  (i - ((REAL(signal_end-(signal_start + 1)))/2)) / (sigma_fft * ((REAL(signal_end-(signal_start + 1)))/2))  )**2)   )  ) * raw_input(i+1)
!        END DO                           
    ELSE IF (window_fn == 2) THEN            
        !   Hann Window
        DO  i = 0, (signal_end - (signal_start + 1))
            IF  (i == 0)   THEN
                fft_input_window(i+1) = ( 0.5 + (0.5 * cos ((2*pi*i)/REAL(signal_end-(signal_start + 1)))) ) * fft_input(i+1)  
            ELSE
                fft_input_window(i+1) = ( 0.5 - (0.5 * cos ((2*pi*i)/REAL(signal_end-(signal_start + 1)))) ) * fft_input(i+1)
            END IF
        END DO        
!        !   Hann Window
!        DO  i = 0, (signal_end - (signal_start + 1))
!            IF  (i == 0)   THEN
!                raw_input_window(i+1) = ( 0.5 + (0.5 * cos ((2*pi*i)/REAL(signal_end-(signal_start + 1)))) ) * raw_input(i+1)  
!            ELSE
!                raw_input_window(i+1) = ( 0.5 - (0.5 * cos ((2*pi*i)/REAL(signal_end-(signal_start + 1)))) ) * raw_input(i+1)
!            END IF
!        END DO        
    ELSE IF (window_fn == 3) THEN            
        !   Hamming Window
        DO  i = 0, (signal_end - (signal_start + 1))
            IF  (i == 0)   THEN
                fft_input_window(i+1) = ( 0.54 + (0.46 * cos ((2*pi*i)/REAL(signal_end-(signal_start + 1)))) ) * fft_input(i+1)  
            ELSE
                fft_input_window(i+1) = ( 0.54 - (0.46 * cos ((2*pi*i)/REAL(signal_end-(signal_start + 1)))) ) * fft_input(i+1)
            END IF
        END DO             
!        !   Hamming Window
!        DO  i = 0, (signal_end - (signal_start + 1))
!            IF  (i == 0)   THEN
!                raw_input_window(i+1) = ( 0.54 + (0.46 * cos ((2*pi*i)/REAL(signal_end-(signal_start + 1)))) ) * raw_input(i+1)  
!            ELSE
!                raw_input_window(i+1) = ( 0.54 - (0.46 * cos ((2*pi*i)/REAL(signal_end-(signal_start + 1)))) ) * raw_input(i+1)
!            END IF
!        END DO           
    ELSE IF (window_fn == 4) THEN
        DO  i = 1, 1048576
            fft_input_window(i) = fft_input(i)       
        END DO                   
    END IF
    
    
    DO  i = (signal_end - signal_start + 1), 1048576
        fft_input_window(i) = 0
    END DO        
        
    IF  (window_fn /= 4)    THEN
        WRITE (105, 6000) fft_input_window
        6000 FORMAT (F25.5)
    END IF
         
!    !   Hann-Poisson Window
!    alpha_fft = 3
!    DO  i = 0, (signal_end - (signal_start + 1))
!        IF  (i == 0)   THEN
!            fft_input(i+1) = (( 0.5 + (0.5 * cos ((2*pi*i)/(REAL(signal_end - (signal_start + 1))/2))) )) * ( exp(-alpha_fft * (i/(REAL(signal_end - (signal_start + 1))/2))) ) * fft_input(i+1)  
!            ELSE
!            fft_input(i+1) = ( 0.5 - (0.5 * cos ((2*pi*i)/(REAL(signal_end - (signal_start + 1))/2))) ) * ( exp(-alpha_fft * (i/(REAL(signal_end - (signal_start + 1))/2))) ) * fft_input(i+1)
!        END IF
!    END DO

    
    Status = DftiCreateDescriptor(My_Desc1_Handle, DFTI_SINGLE, DFTI_REAL, 1, npts)   
    Status = DftiSetValue(My_Desc1_Handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
    Status = DftiSetValue(My_Desc1_Handle, DFTI_NUMBER_OF_USER_THREADS, 1)
    Status = DftiSetValue(My_Desc1_Handle, DFTI_FORWARD_SCALE, 1.0)
    Status = DftiCommitDescriptor(My_Desc1_Handle)
    Status = DftiComputeForward(My_Desc1_Handle, fft_input_window, fft_output)
    Status = DftiFreeDescriptor(My_Desc1_Handle)
        
    ! Calculate magnitude of DFT from real and imaginary parts
    ! Real parts of DFT are in odd rows [(2*i)-1], Imaginary parts are in even rows (2*i) {for i = 1, npts}
    ! Real parts of DFT are in even rows [(2*i)-1], Imaginary parts are in odd rows (2*i) {for i = 0, (npts-1)}
    DO  i = 1, 524288
        fft_mag(i) = sqrt( ((fft_output ((2*i)-1))**2) + ((fft_output (2*i))**2) )
!        phase (i) = atan( (fft_output (2*i)) / (fft_output ((2*i)-1)) )    ! arctan (imag/real)
    END DO

     
    END SUBROUTINE section_fft