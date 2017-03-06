    SUBROUTINE trap_charge_calibration(length_multiplier, load_time, expt_num, data_type)
    
    IMPLICIT NONE
    CHARACTER (len = 100) :: filename                ! Dummy variable for opening files
    CHARACTER(len=8) :: fmt                          ! Format descriptor for file_number
    CHARACTER(len=8) :: file_number                  ! File number to read
    CHARACTER(len=8) :: file_freq                    ! File frequency to read
    CHARACTER(len=8) :: file_attenuation             ! File attenuation to read
    
    INTEGER :: nvals = 0                            ! Number of datapoints
    INTEGER :: calibration_channel                  ! Integer value to determine calibration option
    INTEGER :: noise_channel                        ! Integer value to determine channel of noise files
    INTEGER :: filter_option                        ! Indicator for application of a filter 
    INTEGER :: first_file_number                    ! Number of the first file to be analyzed
    INTEGER :: last_file_number                     ! Number of the last file to be analyzed
    INTEGER :: number_of_files                      ! Number of files to analyze (max = 100)              
    INTEGER :: i, j, k, l, m, n, y                  ! Loop variables
    INTEGER :: unit                                 ! Unit number for r/w files
    INTEGER :: stat                                 ! Error message when opening files
    INTEGER :: data_type                            ! Pointer for .dat or .txt files to analyze
    INTEGER :: signal_end                           ! Length of signal on which to apply FFT
    INTEGER :: window_function                      ! Indicator for application of apodization function
    INTEGER :: mult                                 ! Number of cycles in all of datafile (1048576/square_total_real)
    INTEGER :: start                                ! Start of cycle as determined by random number generator
    INTEGER :: n_main                               ! Pointer for file section
    INTEGER :: attenuation                          ! Attenuation of input signal in decibels
    INTEGER :: derivative_peak                      ! Inflection point of derivative
    INTEGER :: der_opt                              ! Derivative option for looking for positive or negative derivative
!    INTEGER :: cycle_start                          ! Cycle number to begin analyzing
!    INTEGER :: cycle_end                            ! Cycle number on which to end analysis
    
    INTEGER :: length_multiplier                    ! Multiplier for length of each trapping event (multiple of 1 ms or 2500 datapoints)
    INTEGER :: load_time                            ! Time between end of one trappng event and beginning of next (multiple of 1 ms or 2500 datapoints)
    INTEGER :: expt_num                             ! Maximum number of experiments per file based on signal length (length_multiplier)
     
    REAL(4) :: input_freq                           ! Frequency of input calibration signal
    REAL(4) :: freq_shift                           ! Frequency shift of test signal (in percent)
    REAl(4) :: avg_sum                              ! Sum of all magnitude values (at a given time length) for averaging purposes
    REAL(4) :: avg                                  ! Average magnitude value at a given time length
    REAL(4) :: diff_squared                         ! Sum of square of difference between avg and peak_magnitude(j)
    REAL(4) :: std_dev                              ! Standard deviation of avg 
    REAL(4) :: rel_std_dev                          ! Relative standard deviation of avg   
    REAl(4) :: freq_avg_sum                         ! Sum of all frequency values (at a given time length) for averaging purposes
    REAL(4) :: freq_avg                             ! Average frequency value at a given time length
    REAL(4) :: freq_diff_squared                    ! Sum of square of difference between freq_avg and freq_value(j)
    REAL(4) :: freq_std_dev                         ! Standard deviation of freq_avg  
    REAL(4) :: freq_rel_std_dev                     ! Relative standard deviation of freq_avg   
    REAl(4) :: integral_avg_sum                     ! Sum of all frequency values (at a given time length) for averaging purposes
    REAL(4) :: integral_avg                         ! Average frequency value at a given time length
    REAL(4) :: integral_diff_squared                ! Sum of square of difference between freq_avg and freq_value(j)
    REAL(4) :: integral_std_dev                     ! Standard deviation of freq_avg  
    REAL(4) :: integral_rel_std_dev                 ! Relative standard deviation of freq_avg     
!    REAl(4) :: weighted_freq_avg_sum                ! Sum of all weighted frequency values (at a given time length) for averaging purposes
!    REAL(4) :: weighted_freq_avg                    ! Average weighted frequency value at a given time length
!    REAL(4) :: weighted_freq_diff_squared           ! Sum of square of difference between weighted_freq_avg and weighted_freq_value(j)
!    REAL(4) :: weighted_freq_std_dev                ! Standard deviation of weighted_freq_avg  
!    REAL(4) :: weighted_freq_rel_std_dev            ! Relative standard deviation of weighted_freq_avg    
    REAL(4) :: time                                 ! Length of sampled signal based on number of datapoints selected
    REAL(4) :: cycles                               ! Number of cycles based on length of sampled signal and frequency
    REAL(4) :: millivolts                           ! Voltage of mock input signal (in mV)
    REAL(4) :: electrons                            ! Number of charges on mock input signal (in electrons)
    REAL(4) :: peak_magnitude                       ! FFT magnitude value returned from subroutine "calibration_fft"
    REAL(4) :: harm2_rel_mag                        ! Ratio of 2nd harmonic magnitude to fundamental
    REAL(4) :: freq_value                           ! Frequency value returned from subroutine "calibration_fft"
    REAL(4) :: integral_value                       ! Integral value returned from subroutine "calibration_fft"
    
    REAL(8) :: A1                                   ! Variable to define test signal (obtained from Origin)
    REAL(8) :: A2                                   ! Variable to define test signal (obtained from Origin)
    REAL(8) :: x0                                   ! Variable to define test signal (obtained from Origin)
    REAL(8) :: dx                                   ! Variable to define test signal (obtained from Origin)
    REAL(8) :: x                                    ! Variable to define test signal (obtained from Origin)
    
    REAL(8) :: square_total_real                    ! Length of test signal
    REAL(8) :: old_square_total_real                ! Length of previous test signal before frequency shift
    REAL(8) :: remainder                            ! Offset for next cycle in test signal
    
    INTEGER :: array_pointer                        ! Counter for creating test signal with frequency shift
    REAL(4) :: rand                                 ! Random number
    
    INTEGER(2), DIMENSION (1:1000000) :: input_array 
    REAL(4), DIMENSION (1:1048576) :: sq_wave_array     ! Square wave array
    REAL(4), DIMENSION (1:1048576) :: sq_wave           ! Square wave array after random start point selected
    REAL(4), DIMENSION (1:1048576) :: y_input           ! Input array
    REAL(4), DIMENSION (1:1048576) :: adj_input         ! Input datapoints from trap data
    REAL(4), DIMENSION (1:1048576) :: derivative        ! First derivative of y_input
    REAL(4), DIMENSION (1:1048576) :: y_filter          ! Output of 10 kHz high-pass filter                     
    REAL(4), DIMENSION (1:1048576) :: fft_filter_output ! Output of 10 kHz FFT filter
    REAL(4), DIMENSION (1:100) :: peaks                 ! Array for square wave peak magnitudes 
    REAL(4), DIMENSION (1:100) :: freq                  ! Array for frequency values of square wave max peak magnitudes  
    REAL(4), DIMENSION (1:100) :: peak_int              ! Array for integral values of square wave peak
    REAL(4), DIMENSION (1:100) :: rnd                   ! Random number array
    
    
    REAL(4), ALLOCATABLE :: peak_mag(:)             ! Array for peak magnitudes 
    REAL(4), ALLOCATABLE :: frequency(:)            ! Array for frequency values of max peak magnitudes
    REAL(4), ALLOCATABLE :: integral(:)             ! Array for integral values of peak
    
    
    
    
    INTEGER :: file_channel
    INTEGER :: number_of_cycles                      ! Number of cycles per file based on length datafile and frequency
    INTEGER :: cycle_length                          ! Number of datapoints per cycle based on frequency
    
    REAl(4) :: sum_max                              ! Sum of all local maxima for averaging purposes
    REAl(4) :: sum_min                              ! Sum of all local minima for averaging purposes
    REAL(4) :: max_avg                              ! Average of maxima values
    REAL(4) :: min_avg                              ! Average of minima values
    REAL(4) :: max_diff_squared                     ! Sum of square of difference between max_avg and max_val_array(i)
    REAL(4) :: min_diff_squared                     ! Sum of square of difference between min_avg and min_val_array(i)
    REAL(4) :: max_std_dev                          ! Standard deviation of max_avg 
    REAL(4) :: min_std_dev                          ! Standard deviation of min_avg
    REAL(4) :: amplitude                            ! Difference between max_avg and min_avg
    REAL(4) :: amp_error                            ! Amplitude error 
    
    
    REAL(4), ALLOCATABLE :: max_val_array(:)         ! Array for local maxima values 
    REAL(4), ALLOCATABLE :: min_val_array(:)         ! Array for local minima values
    REAL(4), ALLOCATABLE :: cycle_array(:)           ! Array for amplitude values in one cycle
    REAL(4), ALLOCATABLE :: y_input_array (:,:)
       
    
    
    
    ! User selected calibration options
    PRINT *, "Enter 1 for channel A calibration..."
    PRINT *, "Enter 2 for channel B calibration..."
    PRINT *, "Enter 3 for noiseless test signal (to see FFT magnitude uncertainty)..."
    PRINT *, "Enter 4 for test signal added onto noise files..."
    PRINT *, "Enter 5 for calibration of input charge to ADC bits..."
    READ *, calibration_channel
    WRITE (97,*) "Enter 1 for channel A calibration..."
    WRITE (97,*) "Enter 2 for channel B calibration..."
    WRITE (97,*) "Enter 3 for noiseless test signal (to see FFT magnitude uncertainty)..."
    WRITE (97,*) "Enter 4 for test signal added onto noise files..."
    WRITE (97,*) "Enter 5 for calibration of input charge to ADC bits..."
    WRITE (97,*) calibration_channel
    
    
    expt_num = INT(400 / length_multiplier)   
    
        
    ! Actions to take for calibrating channel "A" or "B"
    IF (calibration_channel == 1 .OR. calibration_channel == 2) THEN
        
        PRINT *, "Enter attenuation in decibels..."
        READ *, attenuation
        WRITE (97,*) "Enter attenuation in decibels..."
        WRITE (97,*) attenuation
        
        ! Select frequency of calbration signal        
        PRINT *, "Enter frequency of input calibration signal in kHz..."
        READ *, input_freq        
        WRITE (97,*) "Enter frequency of input calibration signal in kHz..."
        WRITE (97,*) input_freq 
        
        CLOSE (UNIT=97, STATUS='KEEP')
        
        fmt = '(I2.2)'
        WRITE (file_attenuation, fmt) attenuation
        
        fmt = '(I3.3)'
        WRITE (file_freq, fmt) INT(input_freq)
        filename = TRIM(file_attenuation)//'dB_'//TRIM(file_freq)//'kHz.txt'                        
        
        input_freq = input_freq * 1000
                        
!        PRINT *, "Enter output file name (100 characters max)..."
!        READ *, filename
!        WRITE (97,*) "Enter output file name (100 characters max)..."
!        WRITE (97,*) filename
        unit = 101    
        CALL open_new_file (unit, filename)
        
        unit = 102
        filename = "peak_magnitudes.txt"  
        CALL open_new_file (unit, filename)
        
        ! Headers for output files 
!        WRITE (101,1002) "Freq", "Data_Pts", "Seconds", "Cycles", "Magnitude", "Std_Dev", "Rel_Std_Dev", "Freq_Avg", "Freq_Std_Dev", "Freq_Rel_Std_Dev", &
!                         & "Integral_Avg", "Integral_Std_Dev", "Integral_Rel_Std_Dev"
!        1002 FORMAT  (T5, A5, T13, A9, T27, A8, T40, A7, T56, A10, T80, A8, T93, A11, T115, A8, T128, A12, T145, A16, TR10, A12, TR10, A16, TR5, A20)
        WRITE (101,1002) "Freq", "Data_Pts", "Seconds", "Cycles", "Harm1+2_Magnitude", "Std_Dev", "Rel_Std_Dev", "Freq_Avg", "Freq_Std_Dev", "Freq_Rel_Std_Dev", "Harm2_Rel_Mag"
        1002 FORMAT  (T5, A5, T13, A9, T27, A8, T40, A7, T52, A18, T80, A8, T93, A11, T115, A8, T128, A12, T145, A16, T166, A13)
                
        
        ! User selected option for number of files in each run
        ! Maximum number of files is only limited by how the code is currently written
        ! Can be increased by adding extra code to associated "DO" loop below
!        PRINT *, "Enter number of files..."
!        READ *, number_of_files  
!        WRITE (97,*) "Enter number of files..."
!        WRITE (97,*) number_of_files 
        PRINT *, "Enter number of first file..."
        READ *, first_file_number
        PRINT *, "Enter number of last file..."
        READ *, last_file_number
        WRITE (97,*) "Enter number of first file..."
        WRITE (97,*) first_file_number
        WRITE (97,*) "Enter number of last file..."
        WRITE (97,*) last_file_number
        
        
        number_of_files = last_file_number - first_file_number + 1
        ALLOCATE (y_input_array(1:1048576, 1:number_of_files))   

        ! This DO loop assigns name of datafile to open
        ! This streamlines data acquisition and analysis
        DO  j = 0, (number_of_files - 1)
         
            IF (last_file_number < 10000) fmt = '(I4.4)'
            IF (last_file_number >= 10000) fmt = '(I5.5)'
!            fmt = '(I5.5)'
!            unit = j + 10002 
!            ! No need to delete first two files with triggered boards, so...
!            unit = j + 10000
            unit = j + first_file_number
            
            IF  (data_type == 1)    THEN
            
                IF  (calibration_channel == 1)  THEN
!                    WRITE (file_number, fmt) (j+10002)
                    WRITE (file_number, fmt) (j+first_file_number)
                    filename = 'chA'//TRIM(file_number)//'.dat'
                ELSE IF (calibration_channel == 2)  THEN
!                    WRITE (file_number, fmt) (j+10002)
                    WRITE (file_number, fmt) (j+first_file_number)
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
            
                IF  (calibration_channel == 1)  THEN
!                    WRITE (file_number, fmt) (j+10002)
                    WRITE (file_number, fmt) (j+first_file_number)
                    filename = 'chA'//TRIM(file_number)//'.txt'
                ELSE IF (calibration_channel == 2)  THEN
!                    WRITE (file_number, fmt) (j+10002)
                    WRITE (file_number, fmt) (j+first_file_number)
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
            
            y_input_array(:,(j+1)) = adj_input
            
        END DO
        
                    
        
        ! Performs calibration for both filtered and non-filtered
!        DO  y = 1, 2
        DO  y = 1, 1
        
            IF  (y == 1)    THEN
                filter_option = 1
            ELSE
                filter_option = 2
            END IF
                
                    
            ! Performs calibration for both Gaussian windowed data and non-windowed data
!            DO  n = 1, 2
            DO  n = 1, 1
                
                ! Writes type of window in output file
                IF  (y == 1 .AND. n == 1)    THEN
                    window_function = 1
                    WRITE (101,201) "Gaussian_Window_Data_Filtered"
                    WRITE (102,201) "Gaussian_Window_Data_Filtered"
                    201 FORMAT (A30)
                ELSE IF (y == 1 .AND. n == 2)   THEN
                    window_function = 4
                    WRITE (101,201) "Non-Windowed_Data_Filtered"
                    WRITE (102,201) "Non-Windowed_Data_Filtered"
                ELSE IF (y == 2 .AND. n == 1)   THEN
                    window_function = 1
                    WRITE (101,201) "Gaussian_Window_Data_No_Filter"
                    WRITE (102,201) "Gaussian_Window_Data_No_Filter"
                ELSE IF (y == 2 .AND. n == 2)   THEN
                    window_function = 4
                    WRITE (101,201) "Non-Windowed_Data_No_Filter"
                    WRITE (102,201) "Non-Windowed_Data_No_Filter"
                END IF
                
            
                ! Allocate size of arrays based on number of calibration files (13 sections per file)
                ALLOCATE(peak_mag(1:(number_of_files*expt_num)), frequency(1:(number_of_files*expt_num)), integral(1:(number_of_files*expt_num)))
                
!                number_of_files = 10
                
                ! Looks at different lengths of signals
                ! Takes FFT at given length of many files and finds avg FFT magnitude
                DO  i = 1, INT(2500. * REAL(length_multiplier) / 2000.)
!                    signal_end = i  
                    signal_end = ((i-1) * 2000) + 2500
                    PRINT *, signal_end  
                    
                    ! This DO loop assigns name of datafile to open
                    ! This streamlines data acquisition and analysis
                    DO  j = 1, number_of_files
                     
                        y_input = y_input_array(:,j)
                        
                        ! Apply a filter in accordance with user-defined input            
                        IF  (filter_option == 1)    THEN
                            CALL fft_filter(y_input, fft_filter_output)
                            y_input = fft_filter_output
                        END IF
                        
                        ! Finds FFT magnitude for each section at given length                
                        DO  k = 1, expt_num
                            n_main = k
                            
                            CALL calibration_fft(y_input, signal_end, n_main, window_function, input_freq, j, peak_magnitude, freq_value, integral_value, number_of_files, length_multiplier, harm2_rel_mag)
                            
                            ! Assign mag, freq, and integral values to new array
                            peak_mag(k + ((j-1)*expt_num)) = peak_magnitude
                            frequency(k + ((j-1)*expt_num)) = freq_value   
                            integral(k + ((j-1)*expt_num)) = integral_value   
                        END DO
                           
                    END DO    
                    
                    
                    ! Sum all peak magnitude values at a given frequency
                    ! Sum all frequency values at a given signal length
                    ! Sum all integral values at a given signal length
                    avg_sum = 0
                    freq_avg_sum = 0
                    integral_avg_sum = 0
                        DO  j = 1, (number_of_files*expt_num)
                            avg_sum = avg_sum + peak_mag(j)
                            freq_avg_sum = freq_avg_sum + frequency(j)
                            integral_avg_sum = integral_avg_sum + integral(j)
                        END DO
                    
                    ! Find the average magnitude of a given frequency
                    ! Find the average frequency at a given signal length
                    ! Find the average integral at a given signal length
                    avg = avg_sum/(number_of_files*expt_num)
                    freq_avg = freq_avg_sum/(number_of_files*expt_num)
                    integral_avg = integral_avg_sum/(number_of_files*expt_num)
                    
                    ! Sum of the square of the difference betweeen each value and the average
                    diff_squared = 0
                    freq_diff_squared = 0
                    integral_diff_squared = 0
                        DO  j = 1, (number_of_files*expt_num)
                            diff_squared = diff_squared + ((peak_mag(j) - avg)**2)
                            freq_diff_squared = freq_diff_squared + ((frequency(j) - freq_avg)**2)
                            integral_diff_squared = integral_diff_squared + ((integral(j) - integral_avg)**2)
                        END DO
                    
                    ! Standard deviation of magnitude average
                    ! Standard deviation of frequency average
                    ! Standard deviation of integral average
                    std_dev = SQRT(diff_squared/((number_of_files*expt_num)-1))
                    freq_std_dev = SQRT(freq_diff_squared/((number_of_files*expt_num)-1))
                    integral_std_dev = SQRT(integral_diff_squared/((number_of_files*expt_num)-1))
                    
                    ! Relative standard deviation of magnitude average
                    ! Relative standard deviation of frequency average
                    ! Relative standard deviation of integral average
                    rel_std_dev = std_dev / avg * 100
                    freq_rel_std_dev = freq_std_dev / freq_avg * 100
                    integral_rel_std_dev = integral_std_dev / integral_avg * 100
                               
                    time = REAL(signal_end - 1) * 16. / (40.E6)   
                    cycles = time * input_freq          
                    
                    ! Write data to output files
!                    WRITE (101,200) input_freq, signal_end, time, cycles, avg, std_dev, rel_std_dev, freq_avg, freq_std_dev, freq_rel_std_dev, &
!                                    & integral_avg, integral_std_dev, integral_rel_std_dev
!                    200 FORMAT  (T2, F8.1, T14, I6, T25, F11.9, T37, F10.3, T48, F20.5, T69, F20.5, T95, F9.6, T113, F11.4, T128, F11.4, T151, F9.6, &
!                                    & TR2, F20.0, TR2, F20.0, TR15, F9.6)
                    WRITE (101,200) input_freq, signal_end, time, cycles, avg, std_dev, rel_std_dev, freq_avg, freq_std_dev, freq_rel_std_dev, harm2_rel_mag
                    200 FORMAT  (T2, F8.1, T14, I6, T25, F11.9, T37, F10.3, T48, F20.5, T69, F20.5, T95, F9.6, T113, F11.4, T128, F11.4, T151, F9.6, T170, F9.6)
                    
                    WRITE (102,400) signal_end
                    400 FORMAT  (I6)
                    
                    WRITE (102,600) (peak_mag(k), k=1,(number_of_files*expt_num))
                    600 FORMAT (T10, F20.5)
               
                END DO     
                    
                DEALLOCATE(peak_mag, frequency, integral)  
                
            END DO            
        END DO     
        
        DEALLOCATE (y_input_array)
         
        ! Close all output files and keep them (i.e., don't delete them)
        CLOSE (101, STATUS = "KEEP")  
        CLOSE (102, STATUS = "KEEP")  
            
     
       
    ! Actions to take for observing FFT magnitude uncertainty with noiseless test signal    
    ELSE IF (calibration_channel == 3) THEN                    
        
        ! User defined frequency of test signal            
        PRINT *, "Enter desired frequency of test signal..."
        READ *, input_freq  
        WRITE (97,*) "Enter desired frequency of test signal..."
        WRITE (97,*) input_freq
                  
        ! User chosen frequency shift of test signal
        PRINT *, "By what percent would you like the frequency to shift?"
        PRINT *, "(0 = no frequency shift, 5 = 5%, etc.)"
        READ *, freq_shift
        WRITE (97,*) "By what percent would you like the frequency to shift?"
        WRITE (97,*) "(0 = no frequency shift, 5 = 5%, etc.)"
        WRITE (97,*) freq_shift
        
!        ! Option for looking at a selected number of cycles
!        PRINT *, "Enter starting number of cycles value..."
!        READ *, cycle_start
!        WRITE (97,*) "Enter starting number of cycles value..."
!        WRITE (97,*) cycle_start
!                    
!        PRINT *, "Enter ending number of cycles value (total must be ~5 cycles)..."
!        READ *, cycle_end
!        WRITE (97,*) "Enter ending number of cycles value (total must be ~5 cycles)..."
!        WRITE (97,*) cycle_end

        ! User defined ion charge
        PRINT *, "Enter desired input charge (in electrons) for test signal..."
        READ *, electrons        
        WRITE (97,*) "Enter desired input charge (in electrons) for test signal..."
        WRITE (97,*) electrons
        
        ! User defined signal length  
        PRINT *, "Enter length of test signal..."
        PRINT *, " 2500 =  1 ms"
        PRINT *, " 5000 =  2 ms"
        PRINT *, "20000 =  8 ms"
        PRINT *, "40000 = 16 ms"
        PRINT *, "72500 = 29 ms"
        READ *, signal_end
        WRITE (97,*) "Enter length of test signal..."
        WRITE (97,*) " 2500 =  1 ms"
        WRITE (97,*) " 5000 =  2 ms"
        WRITE (97,*) "20000 =  8 ms"
        WRITE (97,*) "40000 = 16 ms"
        WRITE (97,*) "72500 = 29 ms"
        WRITE (97,*) signal_end
        
        CLOSE (UNIT=97, STATUS='KEEP')
        
        ! Open output files
        unit = 301
        filename = "test_signal_array.txt"  
        CALL open_new_file (unit, filename)
        
        unit = 302
        filename = "test_signal_peak_data.txt"  
        CALL open_new_file (unit, filename)
        
        unit = 303
        filename = "test_signal_peak_mag.txt"  
        CALL open_new_file (unit, filename)
        
        ! Headers for output files 
        WRITE (302,3101) "Freq", "Data_Pts", "Seconds", "Cycles", "Magnitude", "Std_Dev", "Rel_Std_Dev", "Freq_Avg", "Freq_Std_Dev", "Freq_Rel_Std_Dev", &
                         & "Integral_Avg", "Integral_Std_Dev", "Integral_Rel_Std_Dev"
        3101 FORMAT  (T5, A5, T13, A9, T27, A8, T40, A7, T56, A10, T80, A8, T93, A11, T115, A8, T128, A12, T145, A16, TR10, A12, TR10, A16, TR5, A20)
                
        
        ! Variable names and values obtained from Origin for creating test signal
        ! Test signal approximates a real signal and follows a Boltzmann distribution
        ! Details of test signal and fitting function obtained from Simion and Origin, respectively
        A1 = -3.35668061258415E-4
        A2 = 0.994194949645035
        x0 = -0.505420267974969
        dx = 0.040649051968342
        
        
        square_total_real = 40.E6 / (input_freq * 16.)
        old_square_total_real = square_total_real
        mult = INT((2500*length_multiplier) / square_total_real)
        remainder = 0      
        array_pointer = 0
        
        ! Create test signal based on an 18.4% duty cycle (based off of Simion simulations)
        ! Also must take into account the offset of each cycle based on limited time resolution
        ! (i.e., not every frequency will have an integer number of datapoints per cycle due to a finite time
        ! resolution in data collection)
        ! This offset changes the center of the Boltzmann distribution slightly so that not every cycle looks exactly the same       
        DO  i = 1, mult-1              
            DO  j = 1, INT(square_total_real)
                array_pointer = array_pointer + 1
                x = -1.715265866 + ((j-1) * 3.430531732 / square_total_real) - (remainder / square_total_real * 100 * (3.430531732 / square_total_real))
                IF (x <= 0)    THEN
                    sq_wave_array (array_pointer) = (A2 + (A1-A2)/(1 + EXP((x-x0)/dx))) * NINT(2.31 * electrons)
                ELSE IF (x > 0) THEN
                    sq_wave_array (array_pointer) = (A2 + (A1-A2)/(1 + EXP((-x-x0)/dx))) * NINT(2.31 * electrons)
                END IF                                 
            END DO
            
            ! Non-integer portion of each cycle length is added onto frequency shift to create a steadily changing frequency
            remainder = square_total_real - INT(square_total_real)
            square_total_real = (old_square_total_real * (1 + ((freq_shift/100)*i/((mult-1)*signal_end/((2500*length_multiplier)))))) + remainder
            
        END DO       
        
        ! The rest of the array is filled with zeroes
        DO  i = (mult * INT(square_total_real)), 1048576
            sq_wave_array (i) = 0
            sq_wave (i) = 0
        END DO        
                
                 
        unit = 1
        n_main = 1   
                    
        CALL RANDOM_SEED()
        CALL RANDOM_NUMBER(rnd)
            
            ! 100 files are averaged with a random starting position within the cycle
            DO  j = 1, 100    
                start = NINT(rnd(j) * old_square_total_real)                
                DO  k = 1, (2500*length_multiplier)
                    sq_wave (k) = sq_wave_array (start + k)
                END DO
                
                CALL calibration_fft(sq_wave, signal_end, n_main, window_function, input_freq, unit, peak_magnitude, freq_value, integral_value, number_of_files, length_multiplier, harm2_rel_mag)
            
                ! Assign mag, freq, and integral values to new array
                peaks(j) = peak_magnitude
                freq(j) = freq_value
                peak_int(j) = integral_value                
            END DO
       

        ! Sum all peak magnitude values at a given frequency
        ! Sum all frequency values at a given signal length
        ! Sum all integral values at a given signal length
        avg_sum = 0
        freq_avg_sum = 0
        integral_avg_sum = 0
            DO  j = 1, 100
                avg_sum = avg_sum + peaks(j)
                freq_avg_sum = freq_avg_sum + freq(j)
                integral_avg_sum = integral_avg_sum + peak_int(j)
            END DO
        
        ! Find the average magnitude of a given frequency
        ! Find the average frequency at a given signal length
        ! Find the average integral at a given signal length
        avg = avg_sum/100
        freq_avg = freq_avg_sum/100
        integral_avg = integral_avg_sum/100
        
        ! Sum of the square of the difference betweeen each value and the average
        diff_squared = 0
        freq_diff_squared = 0
        integral_diff_squared = 0
            DO  j = 1, 100
                diff_squared = diff_squared + ((peaks(j) - avg)**2)
                freq_diff_squared = freq_diff_squared + ((freq(j) - freq_avg)**2)
                integral_diff_squared = integral_diff_squared + ((peak_int(j) - integral_avg)**2)
            END DO
        
        ! Standard deviation of magnitude average
        ! Standard deviation of frequency average
        ! Standard deviation of integral average
        std_dev = SQRT(diff_squared/99)
        freq_std_dev = SQRT(freq_diff_squared/99) 
        integral_std_dev = SQRT(integral_diff_squared/99) 
                
        ! Relative standard deviation of magnitude average
        ! Relative standard deviation of frequency average
        ! Relative standard deviation of integral average
        rel_std_dev = std_dev / avg * 100
        freq_rel_std_dev = freq_std_dev / freq_avg * 100
        integral_rel_std_dev = integral_std_dev / integral_avg * 100
                    
    
        time = REAL(signal_end - 1) * 16. / (40.E6)   
        cycles = time * input_freq  
        
        ! Write data to output files
        WRITE (302,3200) input_freq, signal_end, time, cycles, avg, std_dev, rel_std_dev, freq_avg, freq_std_dev, freq_rel_std_dev, &
                        & integral_avg, integral_std_dev, integral_rel_std_dev
        3200 FORMAT  (T2, F8.1, T14, I6, T25, F11.9, T37, F10.3, T48, F20.5, T69, F20.5, T95, F9.6, T113, F11.4, T128, F11.4, T151, F9.6, &
                        & TR2, F20.0, TR2, F20.0, TR15, F9.6)
        
        WRITE (303,3600) signal_end, cycles, peak_magnitude
        3600 FORMAT (I6, T6, F8.4, T20, F20.5)
                           
        WRITE (301,3000) sq_wave
        3000 FORMAT (F15.5)        
         
        ! Close all output files and keep them (i.e., don't delete them)
        CLOSE (301, STATUS = "KEEP")  
        CLOSE (302, STATUS = "KEEP")    
        CLOSE (303, STATUS = "KEEP")
              
     
    ! Actions to take for calibrating to a test signal overlaid with noise           
    ELSE IF (calibration_channel == 4) THEN                      
        
        ! User defined input frequency       
        PRINT *, "Enter desired frequency of test signal..."
        READ *, input_freq
        WRITE (97,*) "Enter desired frequency of test signal..."
        WRITE (97,*) input_freq
            
!        ! Option for looking at a selected number of cycles      
!        PRINT *, "Enter starting number of cycles value..."
!        READ *, cycle_start
!        WRITE (97,*) "Enter starting number of cycles value..."
!        WRITE (97,*) cycle_start
!                    
!        PRINT *, "Enter ending number of cycles value (total must be ~5 cycles)..."
!        READ *, cycle_end
!        WRITE (97,*) "Enter ending number of cycles value (total must be ~5 cycles)..."
!        WRITE (97,*) cycle_end
        
        ! User defined ion charge
        PRINT *, "Enter desired input charge (in electrons) for test signal..."
        READ *, electrons      
        WRITE (97,*) "Enter desired input charge (in electrons) for test signal..."
        WRITE (97,*) electrons 
        
        PRINT *, "Are noise files from channel A or channel B?"
        PRINT *, "Enter 1 for A..."
        PRINT *, "Enter 2 for B..."
        READ *, noise_channel
        WRITE (97,*) "Are noise files from channel A or channel B?"
        WRITE (97,*) "Enter 1 for A..."
        WRITE (97,*) "Enter 2 for B..."
        WRITE (97,*) noise_channel
          
        ! User selected option for number of files in each run
        ! Maximum number of files is only limited by how the code is currently written
        ! Can be increased by adding extra code to associated "DO" loop below
        PRINT *, "Enter number of files..."
        READ *, number_of_files  
        WRITE (97,*) "Enter number of files..."
        WRITE (97,*) number_of_files
        
        CLOSE (UNIT=97, STATUS='KEEP')


        ! Open output files
        unit = 301
        filename = "mag_freq_integral_values.txt"  
        CALL open_new_file (unit, filename)
        
        unit = 302
        filename = "test_signal_peak_data.txt"  
        CALL open_new_file (unit, filename)
        
        unit = 303
        filename = "test_signal_peak_mag.txt"  
        CALL open_new_file (unit, filename)
        
        unit = 304
        filename = "freq_data_all_lengths.txt"
        CALL open_new_file (unit, filename)
        
        ! Headers for output files 
        WRITE (302,5101) "Freq", "Data_Pts", "Seconds", "Cycles", "Magnitude", "Std_Dev", "Rel_Std_Dev", "Freq_Avg", "Freq_Std_Dev", "Freq_Rel_Std_Dev", &
                         & "Integral_Avg", "Integral_Std_Dev", "Integral_Rel_Std_Dev"
        5101 FORMAT  (T5, A5, T13, A9, T27, A8, T40, A7, T56, A10, T80, A8, T93, A11, T115, A8, T128, A12, T145, A16, TR10, A12, TR10, A16, TR5, A20)
                
        
        ! Variable names and values obtained from Origin for creating test signal
        ! Test signal approximates a real signal and follows a Boltzmann distribution
        ! Details of test signal and fitting function obtained from Simion and Origin, respectively
        A1 = -3.35668061258415E-4
        A2 = 0.994194949645035
        x0 = -0.505420267974969
        dx = 0.040649051968342
        
        
        square_total_real = 40.E6 / (input_freq * 16.)
        old_square_total_real = square_total_real
        mult = INT((2500*length_multiplier) / square_total_real)
        remainder = 0      
        array_pointer = 0
        
        ! Create test signal based on an 18.4% duty cycle (based off of Simion simulations)
        ! Also must take into account the offset of each cycle based on limited time resolution
        ! (i.e., not every frequency will have an integer number of datapoints per cycle due to a finite time
        ! resolution in data collection)
        ! This offset changes the center of the Boltzmann distribution slightly so that not every cycle looks exactly the same       
        DO  i = 1, mult-1              
            DO  j = 1, INT(square_total_real)
                array_pointer = array_pointer + 1
                x = -1.715265866 + ((j-1) * 3.430531732 / square_total_real) - (remainder / square_total_real * 100 * (3.430531732 / square_total_real))
                IF (x <= 0)    THEN
                    sq_wave_array (array_pointer) = (A2 + (A1-A2)/(1 + EXP((x-x0)/dx))) * NINT(2.31 * electrons)
                ELSE IF (x > 0) THEN
                    sq_wave_array (array_pointer) = (A2 + (A1-A2)/(1 + EXP((-x-x0)/dx))) * NINT(2.31 * electrons)
                END IF                  
            END DO
            
            ! Non-integer portion of each cycle length is added onto frequency shift to create a steadily changing frequency
            remainder = square_total_real - INT(square_total_real)
            square_total_real = (old_square_total_real * (1 + ((freq_shift/100)*i/((mult-1)*signal_end/((2500*length_multiplier)))))) + remainder
            
        END DO       
        
        ! The rest of the array is filled with zeroes
        DO  i = (mult * INT(square_total_real)), 1048576
            sq_wave_array (i) = 0
        END DO   
          
        ! Performs calibration for both filtered and non-filtered
        DO  y = 1, 2
        
            IF  (y == 1)    THEN
                filter_option = 1
            ELSE
                filter_option = 2
            END IF
                        
           
            ! Performs calibration for both Gaussian windowed data and non-windowed data
            DO  n = 1, 2
                
                ! Writes type of window in output file
                IF  (y == 1 .AND. n == 1)    THEN
                    window_function = 1
                    WRITE (101,201) "Gaussian_Window_Data_Filtered"
                    WRITE (102,201) "Gaussian_Window_Data_Filtered"
                ELSE IF (y == 1 .AND. n == 2)   THEN
                    window_function = 4
                    WRITE (101,201) "Non-Windowed_Data_Filtered"
                    WRITE (102,201) "Non-Windowed_Data_Filtered"
                ELSE IF (y == 2 .AND. n == 1)   THEN
                    window_function = 1
                    WRITE (101,201) "Gaussian_Window_Data_No_Filter"
                    WRITE (102,201) "Gaussian_Window_Data_No_Filter"
                ELSE IF (y == 2 .AND. n == 2)   THEN
                    window_function = 4
                    WRITE (101,201) "Non-Windowed_Data_No_Filter"
                    WRITE (102,201) "Non-Windowed_Data_No_Filter"
                END IF
            
            
                ! Allocate size of arrays based on number of calibration files (13 sections per file)
                ALLOCATE(peak_mag(1:(number_of_files*(expt_num-1))), frequency(1:(number_of_files*(expt_num-1))), integral(1:(number_of_files*(expt_num-1))))
                
!                number_of_files = 10        
                
                ! Looks at different lengths of signals
                ! Takes FFT at given length of many files and finds avg FFT magnitude 
                DO  i = 1, INT(2500. * REAL(length_multiplier) / 2000.)
                    signal_end = ((i-1) * 2000) + 2500
                    
                    PRINT *, signal_end  
                
                    ! This DO loop assigns name of datafile to open
                    ! This streamlines data acquisition and analysis
                    DO  j = 0, (number_of_files - 1)
                    
                        fmt = '(I5.5)'  
                        unit = j + 10000 
                                         
                        IF  (data_type == 1)    THEN
                        
                            IF  (noise_channel == 1)  THEN
                                WRITE (file_number, fmt) j+10000
                                filename = 'chA'//TRIM(file_number)//'.dat'
                            ELSE IF (noise_channel == 2)  THEN
                                WRITE (file_number, fmt) j+10000
                                filename = 'chB'//TRIM(file_number)//'.dat'
                            END IF
                         
                            PRINT *, "Now reading from ", filename    
                            OPEN (UNIT=unit, FILE=filename, STATUS="old", ACTION="read", IOSTAT=stat, FORM="binary")             
                                                    
                            input_array = 0
                            y_input = 0
                            
                            ! Read input file into array "y_input"
                            READ (unit, END=81) input_array            ! Read filename values and store in array "y"
                            81   CONTINUE   
                        
                            DO  k = 1, 1000000
                                y_input(k) = input_array(k)
                            END DO
                        
                        ELSE IF (data_type == 2)   THEN
                        
                            IF  (noise_channel == 1)  THEN
                                WRITE (file_number, fmt) j+10000
                                filename = 'chA'//TRIM(file_number)//'.txt'
                            ELSE IF (noise_channel == 2)  THEN
                                WRITE (file_number, fmt) j+10000
                                filename = 'chB'//TRIM(file_number)//'.txt'
                            END IF                         
                         
                            PRINT *, "Now reading from ", filename    
                            OPEN (UNIT=unit, FILE=filename, STATUS="old", ACTION="read", IOSTAT=stat)           
                                                    
                            y_input = 0
                            
                            ! Read input file into array "y_input"
                            READ (unit, *, END=91) y_input            ! Read filename values and store in array "y"
                            91   CONTINUE              
                         
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
                            ELSE IF (j == 1048576)  THEN
                                derivative (k) = (y_input(k) - y_input(k-1))/(k-(k-1))
                            ELSE           
                                derivative (k) = 0.5 * (((y_input(k+1) - y_input(k))/((k+1)-k)) + ((y_input(k) - y_input(k-1))/(k-(k-1)))) 
                            END IF
                        END DO        
                        
                                    
                        der_opt = 0
                        22  CONTINUE
                        k = 1
                        
                        IF  (der_opt == 0)  THEN
                            DO WHILE    (derivative(k) > -5000) ! Looks for negative going derivative from end pulse @ 29 ms
                                k = k + 1
                                derivative_peak = k + (2500 * load_time)  ! Adding 2500 datapoints brings us to the beginning of the start pulse
                                
                                IF  (k == (2500 * (length_multiplier + 1)))  THEN
                                    der_opt = 1
                                    GOTO 22
                                END IF                
                            END DO 
                       
                        ELSE IF (der_opt == 1)  THEN 
                            
                            DO WHILE    (derivative(k) < 5000) ! Looks for positive going derivative from start pulse @ 29
                                k = k + 1
                                derivative_peak = k 
                                
                                IF  (k == (2500 * (length_multiplier + 1)))  THEN
                                    PRINT *, "Problem taking derivative..."
                                    PRINT *, "Is the BNC Pulser on?"
                                    GOTO 15
                                END IF                
                            END DO             
                        END IF 
                                               
                        
                        DO  k = 0, (1048576 - derivative_peak - 1)
                            y_input(k+1) = y_input(derivative_peak + k)
                        END DO
                        
                        DO  k = (1048576 - derivative_peak), 1048576
                            y_input(k) = 0
                        END DO            
                         
                         
                        ! Apply a filter in accordance with user-defined input 
                        IF  (filter_option == 1)    THEN
                            CALL fft_filter(y_input, fft_filter_output)
                            y_input = fft_filter_output
                        END IF
                        
                        WRITE (301,6998) i
                        WRITE (301,6999) j
                        6998 FORMAT (I6)
                        6999 FORMAT (I2)
                       
                            ! Finds FFT magnitude for each section at given length 
                            ! Each test signal has a random starting position within the cycle
                            DO  m = 1, (expt_num - 1)
                                n_main = m
                                CALL RANDOM_SEED()
                                CALL RANDOM_NUMBER(rand)
                                start = NINT(rand * old_square_total_real)
                                        
                                    DO  k = 1, (2500 * (length_multiplier - load_time))
                                        y_input(k + ((m-1)*(2500*length_multiplier))) = y_input(k + ((m-1)*(2500*length_multiplier))) + sq_wave_array (start + k)
                                    END DO
                                            
                                CALL calibration_fft(y_input, signal_end, n_main, window_function, input_freq, j, peak_magnitude, freq_value, integral_value, number_of_files, length_multiplier, harm2_rel_mag)
                                
                                ! Assign mag, freq & integral values to new array
                                peak_mag(m + ((j-1)*(expt_num - 1))) = peak_magnitude
                                frequency(m + ((j-1)*(expt_num - 1))) = freq_value   
                                integral(m + ((j-1)*(expt_num - 1))) = integral_value  
                                    
                                ! Write data to output file
                                WRITE (301,7000) peak_mag(m + ((j-1)*(expt_num - 1))), frequency(m + ((j-1)*(expt_num - 1))), integral(m + ((j-1)*(expt_num - 1)))
                                7000 FORMAT (T5, F15.2, T30, F9.2, T45, F9.2, TR5, F20.0)
                                
                            END DO
                        
                        15 CONTINUE
                                                       
                    END DO
                        
                                
                    
                    ! Sum all peak magnitude values at a given frequency
                    ! Sum all frequency values at a given signal length
                    avg_sum = 0
                    freq_avg_sum = 0
                    integral_avg_sum = 0
                        DO  j = 1, (number_of_files*(expt_num - 1))
                            avg_sum = avg_sum + peak_mag(j)
                            freq_avg_sum = freq_avg_sum + frequency(j)
                            integral_avg_sum = integral_avg_sum + integral(j)
                            
                            WRITE (304,7001) frequency(j)
                            7001 FORMAT (F9.2)
                        END DO
                    
                    ! Find the average magnitude of a given frequency
                    ! Find the average frequency at a given signal length
                    avg = avg_sum/(number_of_files*(expt_num - 1))
                    freq_avg = freq_avg_sum/(number_of_files*(expt_num - 1))
                    integral_avg = integral_avg_sum/(number_of_files*(expt_num - 1))
                                
                    ! Sum of the square of the difference betweeen each value and the average
                    diff_squared = 0
                    freq_diff_squared = 0
                    integral_diff_squared = 0
                        DO  j = 1, (number_of_files*(expt_num - 1))
                            diff_squared = diff_squared + ((peak_mag(j) - avg)**2)
                            freq_diff_squared = freq_diff_squared + ((frequency(j) - freq_avg)**2)
                            integral_diff_squared = integral_diff_squared + ((integral(j) - integral_avg)**2)
                        END DO
                    
                    ! Standard deviation of magnitude average
                    ! Standard deviation of frequency average
                    std_dev = SQRT(diff_squared/((number_of_files*(expt_num - 1))-1))
                    freq_std_dev = SQRT(freq_diff_squared/((number_of_files*(expt_num - 1))-1)) 
                    integral_std_dev = SQRT(integral_diff_squared/((number_of_files*(expt_num - 1))-1)) 
                    
                    ! Relative standard deviation of magnitude average
                    ! Relative standard deviation of frequency average
                    rel_std_dev = std_dev / avg * 100
                    freq_rel_std_dev = freq_std_dev / freq_avg * 100
                    integral_rel_std_dev = integral_std_dev / integral_avg * 100
                               
                    time = REAL(signal_end - 1) * 16. / (40.E6)   
                    cycles = time * input_freq          
                    
                    ! Write data to output files
                    WRITE (302,6200) input_freq, signal_end, time, cycles, avg, std_dev, rel_std_dev, freq_avg, freq_std_dev, freq_rel_std_dev, &
                                    & integral_avg, integral_std_dev, integral_rel_std_dev
                    6200 FORMAT  (T2, F8.1, T15, I5, T25, F11.9, T37, F10.3, T48, F20.5, T69, F20.5, T95, F9.6, T113, F11.4, T128, F11.4, T151, F9.6, &
                                    & TR2, F20.0, TR2, F20.0, TR15, F9.6)
                    
                    WRITE (303,6600) signal_end, cycles
                    6600 FORMAT (I6, T6, F8.4)
                    
                    WRITE (303,6800) (peak_mag(j), j=1,(number_of_files*(expt_num - 1)))
                    6800 FORMAT (T20, F20.5)
                    
                END DO           
                   
                DEALLOCATE(peak_mag, frequency, integral) 
            
            END DO 
        END DO     
             
        ! Close all output files and keep them (i.e., don't delete them)
        CLOSE (301, STATUS = "KEEP")  
        CLOSE (302, STATUS = "KEEP")  
        CLOSE (303, STATUS = "KEEP")  
        CLOSE (304, STATUS = "KEEP")
        
    
    ! Calibration of input charge to ADC bits    
    ELSE IF (calibration_channel == 5) THEN
       
        ! Select frequency of calbration signal        
        PRINT *, "Enter frequency of input calibration signal (in kHz)..."
        READ *, input_freq 
        WRITE (97,*) "Enter frequency of input calibration signal (in kHz)..."
        WRITE (97,*) input_freq 
        input_freq = input_freq * 1000
        
        ! Select channel for files
        PRINT *, ""
        PRINT *, "Are files from channel A or channel B?"
        PRINT *, "Enter 1 for A..."
        PRINT *, "Enter 2 for B..."
        READ *, file_channel
        WRITE (97,*) ""
        WRITE (97,*) "Are files from channel A or channel B?"
        WRITE (97,*) "Enter 1 for A..."
        WRITE (97,*) "Enter 2 for B..."
        WRITE (97,*) file_channel
          
        ! User selected option for number of files in each run
        ! Maximum number of files is only limited by how the code is currently written
        ! Can be increased by adding extra code to associated "DO" loop below
        PRINT *, "Enter number of files..."
        READ *, number_of_files  
        WRITE (97,*) "Enter number of files..."
        WRITE (97,*) number_of_files
        
        CLOSE (UNIT=97, STATUS='KEEP')
                
        
        number_of_cycles = INT(1.E6 * input_freq * 16. / 40.E6)     ! = input_freq Hz * 0.4 s
        
        ALLOCATE(max_val_array(1:((number_of_cycles-1200)*number_of_files)), min_val_array(1:((number_of_cycles-1200)*number_of_files)))
        
        ! This DO loop assigns name of datafile to open
        ! This streamlines data acquisition and analysis
        DO  j = 0, (number_of_files - 1)
    
            fmt = '(I5.5)'
            unit = j + 10000
            
            IF  (data_type == 1)    THEN
            
                IF  (file_channel == 1)  THEN
                    WRITE (file_number, fmt) j+10002
                    filename = 'chA'//TRIM(file_number)//'.dat'
                ELSE IF (file_channel == 2)  THEN
                    WRITE (file_number, fmt) j+10002
                    filename = 'chB'//TRIM(file_number)//'.dat'
                END IF
             
                PRINT *, "Now reading from ", filename    
                OPEN (UNIT=unit, FILE=filename, STATUS="old", ACTION="read", IOSTAT=stat, FORM="binary")             
            
                input_array = 0
                y_input = 0
                
                ! Read input file into array "y_input"
                READ (unit, END=84) input_array            ! Read filename values and store in array "y"
                84   CONTINUE   
            
                DO  k = 1, 1000000
                    y_input(k) = input_array(k)
                END DO
            
            ELSE IF (data_type == 2)   THEN
            
                IF  (file_channel == 1)  THEN
                    WRITE (file_number, fmt) j+10002
                    filename = 'chA'//TRIM(file_number)//'.txt'
                ELSE IF (file_channel == 2)  THEN
                    WRITE (file_number, fmt) j+10002
                    filename = 'chB'//TRIM(file_number)//'.txt'
                END IF             
             
                PRINT *, "Now reading from ", filename    
                OPEN (UNIT=unit, FILE=filename, STATUS="old", ACTION="read", IOSTAT=stat)           
            
                y_input = 0
                
                ! Read input file into array "y_input"
                READ (unit, *, END=94) y_input            ! Read filename values and store in array "y"
                94   CONTINUE              
             
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
             
            CALL band_pass_fft_filter(y_input, input_freq, fft_filter_output)
            
            
            i = 0    
            square_total_real = (40.E6 / (16. * input_freq))        
            old_square_total_real = square_total_real
            
            DO WHILE (i < (number_of_cycles-1200))
                
                i = i + 1
                cycle_length = INT(square_total_real)  
                                
                ALLOCATE(cycle_array(1:cycle_length))
                
                DO  k = 1, cycle_length
                    cycle_array(k) = fft_filter_output((cycle_length*(i-1)) + k)
                END DO             
                
                max_val_array(i + (j*(number_of_cycles-1200))) = MAXVAL(cycle_array, DIM=1)
                min_val_array(i + (j*(number_of_cycles-1200))) = MINVAL(cycle_array, DIM=1)
                
                
                remainder = square_total_real - INT(square_total_real)
                square_total_real = old_square_total_real + remainder
                
!                IF  (remainder == 0) THEN
!                    square_total_real = old_square_total_real
!                END IF                    
                
                DEALLOCATE(cycle_array)
                
            END DO            
            
        END DO
        
        sum_max = 0
        sum_min = 0
        DO  k = 1, ((number_of_cycles-1200)*number_of_files)
            sum_max = sum_max + max_val_array(k)
            sum_min = sum_min + min_val_array(k)
        END DO
        
        max_avg = sum_max / ((number_of_cycles-1200)*number_of_files)
        min_avg = sum_min / ((number_of_cycles-1200)*number_of_files)
        
        ! Sum of the square of the difference betweeen each value and the average
        max_diff_squared = 0
        min_diff_squared = 0
        DO  k = 1, ((number_of_cycles-1200)*number_of_files)
            max_diff_squared = max_diff_squared + ((max_val_array(k) - max_avg)**2)
            min_diff_squared = min_diff_squared + ((min_val_array(k) - min_avg)**2)
        END DO
        
        ! Standard deviation of max average
        ! Standard deviation of max average
        max_std_dev = SQRT(max_diff_squared/(((number_of_cycles-1200)*number_of_files)-1)) 
        min_std_dev = SQRT(min_diff_squared/(((number_of_cycles-1200)*number_of_files)-1))
        
        amplitude = max_avg - min_avg
        amp_error = SQRT((max_std_dev)**2 + (min_std_dev)**2)
        
        PRINT *, "Frequency: ", input_freq
        PRINT *, "Amplitude: ", amplitude
        PRINT *, "Error:     ", amp_error                
            
        
        DEALLOCATE(max_val_array, min_val_array)
        
        
    END IF 
          
    END SUBROUTINE trap_charge_calibration
       
    
  
  
    
    ! Subroutine for FFT of calibration data
    SUBROUTINE calibration_fft (data_in, signal_end, n_fft, window_fn, freq_input, file_number, peak_magnitudes, freq_values, integral_values, number_of_files, len_mult, harm2_rel_mag)
    
    USE mkl_dfti
    
    IMPLICIT NONE
        
    REAL, PARAMETER :: freq_res = 40.E6 / 2**24   
             
    TYPE (DFTI_DESCRIPTOR), POINTER :: My_Desc1_Handle  
    INTEGER :: Status                                   ! Sets DFT Parameters
    INTEGER :: istat                                    ! Allocate status
    INTEGER :: i, j
    INTEGER :: maxfind_start = 0                        ! Pointer for local maximum               
    INTEGER :: temp_signal_end                          ! Temporary pointer for end of signal
    INTEGER :: win_fn                                   ! Indicator for application of apodization function (1 = Gaussian, 2 = none)
    INTEGER :: file_number                              ! Pointer for input file number
    INTEGER :: number_of_files                          ! Number of files to analyze (max = 100) 
    INTEGER :: peak_start = 0                           ! Pointer for peak start
    INTEGER :: peak_end = 0                             ! Pointer for peak end   
    INTEGER :: len_mult                                 ! Multiplier for length of each trapping event (multiple of 1 ms or 2,500 datapoints)
    
    INTEGER, INTENT(IN) :: n_fft                            ! Pointer for file section
    INTEGER, INTENT(IN) :: signal_end                       ! Length of signal on which to perform DFT
    INTEGER, INTENT(IN) :: window_fn                        ! User supplied input for applying a window function or not
    REAL(4) :: freq_input                                   ! Frequency of input calibration signal
    REAL(4) :: sigma_fft                                    ! Sigma value for Gaussian window function 
!    REAL(4) :: alpha_fft                                    ! Alpha value for Hann-Poisson window function
    REAL(4) :: start                                        ! Dummy variable for beginning of peak   
    REAL(4) :: start_less_one                               ! Dummy variable for beginning of peak minus one  
    REAL(4) :: last                                         ! Dummy variable for end of peak    
    REAL(4) :: last_plus_one                                ! Dummy variable for end of peak plus one 
    REAL(4) :: peak_magnitudes                              ! Peak magnitudes 
    REAL(4) :: freq_values                                  ! Frequency values of max peak magnitudes 
    REAL(4) :: integral_values                               ! Integral value returned from subroutine "calibration_fft"
              
    REAL(4), PARAMETER :: pi = 3.14159265                      
    REAL(4), DIMENSION (1:1048576), INTENT(IN) :: data_in   ! Raw input values 
    
    REAL(4), DIMENSION (1:1048576) :: sig_in            ! Input array
!    REAL(4), DIMENSION (1:180) :: peak_magnitudes       ! Array for peak magnitudes 
!    REAL(4), DIMENSION (1:180) :: freq_values           ! Array for frequency values of max peak magnitudes 
    REAL(4), DIMENSION (1:501) :: maxfind_array         ! Array to find max fft magnitude
    	 
    	 
!    ! Parameters for FFT of length 1048576 
!    ! ===================================
    INTEGER, PARAMETER :: npts = 1048576                    ! Number of datapoints  
    REAL(4), DIMENSION (1:524288) :: freq               ! Frequency x-value after FFT  
    REAL(4), DIMENSION (1:524288) :: fft_mag                ! DFT magnitude values                  
    REAL(4), DIMENSION (1:1048576) :: fft_input             ! Input values after applying rolling average filter 
    REAL(4), DIMENSION (1:1048576) :: fft_output            ! DFT output values   
!    REAL(4), DIMENSION (1:524288) :: phase                 ! DFT phase values


    ! Extra variables needed for "Harmonic Analysis"
    INTEGER :: channel                  ! Pointer for channel to analyze (A, B, or test signal)
    INTEGER :: freq_mult                ! Multiplier for length of FFT (base is 2^16)
    INTEGER :: t                        ! Flag used in user-selected options for channels 3 and 5
    REAL(4) :: fund_mag_win_1           ! No use here; needs declared so "Harmonic_Analysis" will work
    REAL(4), DIMENSION(1:10) :: harm_freq
                                        ! Frequency of first 10 harmonics
    REAL(4), DIMENSION(1:10) :: harm_mag
                                        ! Magnitude of first 10 harmonics
    REAL(4) :: harm2_rel_mag
    
    
    
                 
              
    DO  i = (1 + ((n_fft-1)*(2500*len_mult))), (((n_fft-1)*(2500*len_mult)) + signal_end)
        fft_input(i - ((n_fft-1)*(2500*len_mult))) = data_in(i)
    END DO
    
    !   Make sure everything after signal of interest is zero
    DO  i = (signal_end + 1), 1048576
        fft_input(i) = 0
    END DO 
          
         
    IF  (window_fn == 1) THEN        
        !   Gaussian Window
        sigma_fft = 0.45
        DO  i = 0, (signal_end - 1)
            fft_input(i+1) = ( EXP( -0.5 * ((  (i - ((REAL(signal_end-1))/2)) / (sigma_fft * ((REAL(signal_end-1))/2))  )**2)   )  ) * fft_input(i+1) 
        END DO    
    ELSE IF (window_fn == 2) THEN            
        !   Hann Window
        DO  i = 0, (signal_end - 1)
            IF  (i == 0)   THEN
                fft_input(i+1) = ( 0.5 + (0.5 * cos ((2*pi*i)/REAL(signal_end-1))) ) * fft_input(i+1)  
            ELSE
                fft_input(i+1) = ( 0.5 - (0.5 * cos ((2*pi*i)/REAL(signal_end-1))) ) * fft_input(i+1)
            END IF
        END DO        
    ELSE IF (window_fn == 3) THEN 
        !   Hamming Window
        DO  i = 0, (signal_end - 1)
            IF  (i == 0)   THEN
                fft_input(i+1) = ( 0.54 + (0.46 * cos ((2*pi*i)/REAL(signal_end-1))) ) * fft_input(i+1)  
            ELSE
                fft_input(i+1) = ( 0.54 - (0.46 * cos ((2*pi*i)/REAL(signal_end-1))) ) * fft_input(i+1)
            END IF
        END DO 
    END IF
    
!    !   Hann-Poisson Window
!    alpha_fft = 3
!    DO  i = 0, (signal_end - 1)
!        IF  (i == 0)   THEN
!            fft_input(i+1) = (( 0.5 + (0.5 * cos ((2*pi*i)/(REAL(signal_end-1)/2))) )) * ( exp(-alpha_fft * (i/(REAL(signal_end-1)/2))) ) * fft_input(i+1)  
!        ELSE
!            fft_input(i+1) = ( 0.5 - (0.5 * cos ((2*pi*i)/(REAL(signal_end-1)/2))) ) * ( exp(-alpha_fft * (i/(REAL(signal_end-1)/2))) ) * fft_input(i+1)
!        END IF
!    END DO

    
        
    Status = DftiCreateDescriptor(My_Desc1_Handle, DFTI_SINGLE, DFTI_REAL, 1, npts)   
    Status = DftiSetValue(My_Desc1_Handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
    Status = DftiSetValue(My_Desc1_Handle, DFTI_NUMBER_OF_USER_THREADS, 1)
    Status = DftiSetValue(My_Desc1_Handle, DFTI_FORWARD_SCALE, 1.0)
    Status = DftiCommitDescriptor(My_Desc1_Handle)
    Status = DftiComputeForward(My_Desc1_Handle, fft_input, fft_output)
    Status = DftiFreeDescriptor(My_Desc1_Handle)
        
    ! Calculate magnitude of DFT from real and imaginary parts
    ! Real parts of DFT are in odd rows [(2*i)-1], Imaginary parts are in even rows (2*i) {for i = 1, npts}
    ! Real parts of DFT are in even rows [(2*i)-1], Imaginary parts are in odd rows (2*i) {for i = 0, (npts-1)}
    DO  i = 1, 524288
        fft_mag(i) = sqrt( ((fft_output ((2*i)-1))**2) + ((fft_output (2*i))**2) )
!        phase (i) = atan( (fft_output (2*i)) / (fft_output ((2*i)-1)) )    ! arctan (imag/real)
    END DO
    
    ! Use "Harmonic_Analysis" to determine fundamental frequency and to sum intensities of first two harmonics
    fund_mag_win_1 = 0.; channel = 9; t = 0; freq_mult = 16
    CALL Harmonic_Analysis(channel, freq_mult, fft_mag, freq_input, fund_mag_win_1, signal_end, harm_freq, harm_mag, t)
    
    peak_magnitudes = SUM(harm_mag(1:2))
    freq_values = harm_freq(1)
    harm2_rel_mag = harm_mag(2)/harm_mag(1)
        
    
!    DO  i = 1, 524288
!        freq (i) = (i-1) * 40.E6 / 2**24
!    END DO
!    
!    ! Starting point to look for frequency is 251 points before expected peak value
!    maxfind_start = NINT(freq_input/(40.E6 / 2**24)) - 251
!      
!!    PRINT *, maxfind_start 
!                                 
!        ! Make array of all fft_magnitude values from peak_start to peak_end
!        ! MAXLOC(maxfind_array, DIM=1) +  maxfind_start equals the datapoint of peak max
!        ! Therefore freq(MAXLOC(maxfind_array, DIM=1) +  maxfind_start) equals freq of max magnitude
!        ! And fft_magnitude(MAXLOC(maxfind_array, DIM=1) +  maxfind_start) equals max fft magnitude
!        DO  i = 1, 501
!            maxfind_array(i) = fft_mag(maxfind_start + i)
!        END DO        
!    
!!    PRINT *, MAXLOC(maxfind_array, DIM=1)
!    
!    ! Write peak magnitude and frequency values to array
!    peak_magnitudes = fft_mag(MAXLOC(maxfind_array, DIM=1) + maxfind_start)
!    freq_values = freq(MAXLOC(maxfind_array, DIM=1) + maxfind_start)  
!    
!!    peak_magnitudes( ((file_number-1)*18) + n_fft) = fft_mag(MAXLOC(maxfind_array, DIM=1) + maxfind_start)
!!    freq_values( ((file_number-1)*18) + n_fft) = freq(MAXLOC(maxfind_array, DIM=1) + maxfind_start)  
!    
!
!    ! Find beginning and end of range to integrate for windowed data   
!    peak_start = MAXLOC(maxfind_array, DIM=1) + maxfind_start - 1
!    start = fft_mag(peak_start)
!    start_less_one = fft_mag(peak_start - 1)
!    
!    ! Find start of peak           
!    DO WHILE    (start >= start_less_one)
!                peak_start = peak_start - 1
!                start = fft_mag(peak_start)
!                start_less_one = fft_mag(peak_start - 1)
!    END DO
!    
!    peak_start = peak_start + 1 
!                                
!    
!    peak_end = MAXLOC(maxfind_array, DIM=1) + maxfind_start + 1            
!    last = fft_mag(peak_end)
!    last_plus_one = fft_mag(peak_end + 1)
!    
!    ! Find end of peak
!    DO WHILE    (last >= last_plus_one)           
!                last = fft_mag(peak_end)
!                last_plus_one = fft_mag(peak_end + 1)
!                peak_end = peak_end + 1     
!    END DO
!             
!    
!    ! Integrate peak from peak_start to peak_end
!    integral_values = 0
!    DO  j = peak_start, peak_end
!        integral_values = integral_values + (fft_mag(j) * 40.E6 / 2**24)
!    END DO 
    
    
    END SUBROUTINE calibration_fft