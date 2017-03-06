    SUBROUTINE find_frequency_shift(channel, length_multiplier, load_time, dead_time)
     
    IMPLICIT NONE    
    REAL(4) :: sigma_main                           ! Sigma value for Gaussian window function 
    REAL(4), PARAMETER :: pi = 3.14159265                      
       
    INTEGER(2), DIMENSION (1:1000000) :: input_array 
    REAL(4), DIMENSION (1:1048576) :: y_input            ! I/O y-axis datapoints 
    REAL(4), DIMENSION (1:1048576) :: adj_input          ! Input datapoints from trap data
    REAL(4), DIMENSION (1:1048576) :: derivative         ! First derivative of y_input  
    REAL(4), DIMENSION (1:1048576) :: filtered_input     ! Filtered data of selected file & section moved to beginning of file  
    REAL(4), DIMENSION (1:1048576) :: y_filter           ! Output from 10 kHz high-pass filter                         
    REAL(4), DIMENSION (1:1048576) :: fft_filter_output  ! Output from 1 kHz FFT filter
          
    INTEGER :: i, j, k
    INTEGER :: unit                     ! Dummy variable for opening files
    INTEGER :: stat                     ! Error message when opening files
    INTEGER :: channel                  ! Pointer for channel to analyze (A, B, or test signal)
    INTEGER :: number_of_files          ! Number of files to analyze   
    INTEGER :: n_main                   ! Pointer for file section
    INTEGER :: section
    INTEGER :: window_function          ! Indicator for application of apodization function 
    INTEGER :: start_of_signal          ! Length to cut off at beginning of signal before FFT
    INTEGER :: end_of_signal            ! Length of signal on which to perform FFT
    INTEGER :: files_analyzed           ! Counter for number of files analyzed
    INTEGER :: derivative_peak          ! Inflection point of derivative
    INTEGER :: der_opt                  ! Derivative option for looking for positive or negative derivative
    INTEGER :: expt_num                 ! Maximum number of experiments per file based on signal length (length_multiplier)
    INTEGER :: dead_time                ! Number of points of dead time due to A250 recovering from impulse noise of trap gate pulse
    
    REAL(4) :: charge                   ! Average charge value (windowed FFT)
    
    INTEGER :: length_multiplier        ! Multiplier for length of each trapping event (multiple of 1 ms or 2500 datapoints)
    INTEGER :: load_time                ! Time between end of one trappng event and beginning of next (multiple of 1 ms or 2500 datapoints)
      
    INTEGER :: single_event             ! Counter for number of trapping events with a single ion trapped    
    INTEGER :: multiple_event           ! Counter for number of trapping events with multiple ions trapped
    
    LOGICAL :: od                       ! Logical variable for inquiring whether a unit is opened or not
    
    CHARACTER(len=8) :: fmt             ! Format descriptor for file_number
    CHARACTER(len=10) :: file_number    ! File number to read
    CHARACTER(len=50) :: filename       ! Dummy variable for opening files
    
    REAL(4) :: low_mass_cutoff = 0         ! Low mass cutoff for frequency information
    REAL(4) :: high_mass_cutoff = 0         ! High mass cutoff for frequency information
    INTEGER :: sig_cutoff = 0               ! Signal length cutoff (in datapoints) for frequency information
    
    INTEGER :: trap_time_cutoff         ! Time length cutoff (in ms) for frequency information
     
    REAL(4), DIMENSION (1:1048576) :: input                ! Input datapoints from trap data
    INTEGER, DIMENSION (1:2,1:400) :: derivative_peaks    ! Array for location of derivative peaks from derivative of y_input
    INTEGER :: pos_peaks                ! Counter for positive peaks in derivative of y_input
    INTEGER :: neg_peaks                ! Counter for negative peaks in derivative of y_input
    INTEGER :: max_trap_time                        ! Length of period for trapping event (multiple of 1 ms or 2500 datapoints)
    
    
          
          
    PRINT *, "Enter last file number..."
    READ *, number_of_files
    WRITE (97,*) "Enter last file number..."
    WRITE (97,*) number_of_files
    number_of_files = number_of_files + 1 - 10000            
    expt_num = INT(400 / length_multiplier) - 1 
    
    CLOSE (UNIT=97, STATUS='KEEP')       

    fmt = '(I5.5)' 
    files_analyzed = 0        
    
    max_trap_time = length_multiplier - load_time

        
    !$OMP PARALLEL NUM_THREADS(8) DEFAULT(none) COPYIN(files_analyzed, number_of_files, fmt, channel, dead_time) &
    !$OMP SHARED(files_analyzed, number_of_files, fmt, channel, length_multiplier, load_time, expt_num, dead_time) &
    !$OMP SHARED(low_mass_cutoff, high_mass_cutoff, sig_cutoff, max_trap_time)
        
    !$OMP DO PRIVATE(filename, file_number, unit, stat, y_input, input_array, adj_input, derivative, n_main) &
    !$OMP PRIVATE(derivative_peak, charge, der_opt, fft_filter_output, i, j, single_event, multiple_event) &
    !$OMP PRIVATE(derivative_peaks, od, pos_peaks, neg_peaks, section, input) &
    
    !$OMP SCHEDULE(guided)
    
    DO  i = 0, (number_of_files - 1)
                        
        ! Open input file
        unit = (i + 10000)                
              
        input_array = 0
        y_input = 0
          
        WRITE (file_number, fmt) unit
        filename = 'chA'//TRIM(file_number)//'.dat'
         
        OPEN (UNIT=unit, FILE=filename, STATUS="old", ACTION="read", IOSTAT=stat, FORM="binary")             
        
        INQUIRE (UNIT=unit, OPENED=od)
        
        IF  (.NOT. od)  THEN
            GOTO 15
        ELSE
            ! Read input file into array "input_array"
            READ (unit, END=8) input_array            ! Read filename values and store in array "y"
            8   CONTINUE
        END IF         
        
        DO  j = 1, 1000000
            y_input(j) = input_array(j)
        END DO
        
        CLOSE(unit)
    
        files_analyzed = files_analyzed + 1
        PRINT *, "Now reading from ", filename
        PRINT *, files_analyzed, "/", number_of_files 
            
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
            
            
        DO  j = 1, expt_num
            IF  (((derivative_peaks(2,(j+1)) - derivative_peaks(1,j)) > ((2500 * REAL(max_trap_time + (load_time / 2.))) - 50)) &
                & .AND. ((derivative_peaks(2,(j+1)) - derivative_peaks(1,j)) < ((2500 * REAL(max_trap_time + (load_time / 2.))) + 50))) THEN
                
                derivative_peak = derivative_peaks(1,j) + (2500 * load_time / 2.)
                input = 0
                
                DO  k = 0, ((2500 * max_trap_time) - 1)
                    input(k+1) = y_input(derivative_peak + k)
                END DO
                
                DO  k = ((2500 * max_trap_time) + 1), 1048576
                    input(k) = 0
                END DO                        
                           
                CALL fft_filter(input, fft_filter_output)
                input = fft_filter_output               
                
                n_main = 1
                section = j            
    
                CALL peakfinder_VI(filename, input, n_main, channel, charge, length_multiplier, max_trap_time, dead_time, single_event, multiple_event, low_mass_cutoff, high_mass_cutoff, sig_cutoff, section)
                               
            ELSE IF (((derivative_peaks(2,(j)) - derivative_peaks(1,j)) > ((2500 * REAL(max_trap_time + (load_time / 2.))) - 50)) &
                & .AND. ((derivative_peaks(2,(j)) - derivative_peaks(1,j)) < ((2500 * REAL(max_trap_time + (load_time / 2.))) + 50))) THEN
                
                derivative_peak = derivative_peaks(1,j) + (2500 * load_time / 2.)
                input = 0
                
                DO  k = 0, ((2500 * max_trap_time) - 1)
                    input(k+1) = y_input(derivative_peak + k)
                END DO
                
                DO  k = ((2500 * max_trap_time) + 1), 1048576
                    input(k) = 0
                END DO                        
                           
                CALL fft_filter(input, fft_filter_output)
                input = fft_filter_output               
                
                n_main = 1
                section = j            
    
                CALL peakfinder_VI(filename, input, n_main, channel, charge, length_multiplier, max_trap_time, dead_time, single_event, multiple_event, low_mass_cutoff, high_mass_cutoff, sig_cutoff, section)
                                     
            END IF 
        END DO  
            
            
        15 CONTINUE 
    
    END DO
    !$OMP END DO
    !$OMP END PARALLEL 
    
    END SUBROUTINE find_frequency_shift