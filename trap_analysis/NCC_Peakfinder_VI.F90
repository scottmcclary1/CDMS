    SUBROUTINE  peakfinder_VI(filename, sig_in, n_ptr, channel, charge, length_multiplier, max_trap_time, dead_time, single_event, multiple_event, low_mass_cutoff, high_mass_cutoff, sig_cutoff, section)
        ! sig_in corresponds to input from the main code and n_ptr to n_main. The rest of the variables have the same name.
           
    IMPLICIT NONE
    CHARACTER(*) :: filename            ! Dummy variable for opening files
    CHARACTER(len=100) :: file          ! Dummy variable for opening files
    INTEGER :: unit                     ! Dummy variable for opening files
    
    REAL(4) :: rnd                      ! Random number
    REAL(4) :: max                      ! Dummy variable for peak maximum   
    REAL(4) :: max_less_one             ! Dummy variable for peak maximum minus one  
    REAL(4) :: start                    ! Dummy variable for beginning of peak   
    REAL(4) :: start_less_one           ! Dummy variable for beginning of peak minus one  
    REAL(4) :: last                     ! Dummy variable for end of peak    
    REAL(4) :: last_plus_one            ! Dummy variable for end of peak plus one 
    
    REAL(4) :: freq_max                 ! Frequency value of max FFT magnitude
    REAL(4) :: mag_max                  ! Max FFT magntitude value
    
    REAL(4) :: end_time                 ! Trapping time
    REAL(4) :: full_end_time            ! Trapping time for taking FFT of whole signal (using full_signal_end)
    REAL(4) :: midpoint                 ! Midpoint between find_length_mag MAXVAL and MINVAL
    
    REAL(4) :: temp_freq                ! Temporary frequency value        
    REAL(4) :: freq_input               ! Frequency of peak determined from windowed FFT 
    REAL(4) :: cycles                   ! Number of cycles ion was trapped (windowed FFT)
    REAL(4) :: mass_to_charge           ! Average m/z value (windowed FFT)
    REAL(4) :: final_freq               ! Final frequency value (windowed FFT)
    REAL(4) :: final_mag                ! Final magnitude value (windowed FFT)
    REAL(4) :: charge                   ! Average charge value (windowed FFT)
    REAL(4) :: mass                     ! Average mass value (windowed FFT)  
    
    REAL(4) :: avg
    REAL(4) :: stdev  
    
    REAL(4) :: mass_to_charge_full      ! m/z value over full trapping time
    REAL(4) :: charge_full              ! Charge value over full trapping time
    REAL(4) :: mass_full                ! Mass value over full trapping time    
    
    REAL(4) :: freq_sum                 ! Summation of frequency for averaging purposes
    REAL(4) :: freq_avg                 ! Frequency average
    REAL(4) :: freq_diff_squared        ! The square of (find_length_freq(i) - freq_avg)
    REAL(4) :: freq_std_dev             ! Standard deviation of freq_avg
    REAL(4) :: freq_rel_std_dev         ! Relative standard deviation of freq_avg
    REAL(4) :: freq_sum_1                 ! Summation of frequency for averaging purposes
    REAL(4) :: freq_avg_1                 ! Frequency average
    REAL(4) :: freq_diff_squared_1        ! The square of (find_length_freq(i) - freq_avg)
    REAL(4) :: freq_std_dev_1             ! Standard deviation of freq_avg
    REAL(4) :: freq_sum_2                 ! Summation of frequency for averaging purposes
    REAL(4) :: freq_avg_2                 ! Frequency average
    REAL(4) :: freq_diff_squared_2        ! The square of (find_length_freq(i) - freq_avg)
    REAL(4) :: freq_std_dev_2             ! Standard deviation of freq_avg
    REAL(4) :: freq_shift                 ! Frequency shift from freq_avg_1 & freq_avg_2
        
    REAL(4) :: mass_sum                 ! Summation of mass for averaging purposes
    REAL(4) :: mass_diff_squared        ! The square of (mass_array(i) - mass)
    REAL(4) :: mass_std_dev             ! Standard deviation of mass
    REAL(4) :: mass_sum_1                 ! Summation of frequency for averaging purposes
    REAL(4) :: mass_avg_1                 ! Frequency average
    REAL(4) :: mass_diff_squared_1        ! The square of (find_length_freq(i) - freq_avg)
    REAL(4) :: mass_stdev_1             ! Standard deviation of freq_avg
    REAL(4) :: mass_sum_2                 ! Summation of frequency for averaging purposes
    REAL(4) :: mass_avg_2                 ! Frequency average
    REAL(4) :: mass_diff_squared_2        ! The square of (find_length_freq(i) - freq_avg)
    REAL(4) :: mass_stdev_2             ! Standard deviation of freq_avg
    REAL(4) :: mass_shift                 ! Frequency shift from freq_avg_1 & freq_avg_2
    
    REAL(4) :: mass_to_charge_sum       ! Summation of m/z for averaging purposes
    REAL(4) :: mass_to_charge_diff_squared      ! The square of (mass_to_charge_array(i) - mass_to_charge)
    REAL(4) :: mass_to_charge_std_dev   ! Standard deviation of m/z
    REAL(4) :: m_z_sum_1                 ! Summation of frequency for averaging purposes
    REAL(4) :: m_z_avg_1                 ! Frequency average
    REAL(4) :: m_z_diff_squared_1        ! The square of (find_length_freq(i) - freq_avg)
    REAL(4) :: m_z_stdev_1             ! Standard deviation of freq_avg
    REAL(4) :: m_z_sum_2                 ! Summation of frequency for averaging purposes
    REAL(4) :: m_z_avg_2                 ! Frequency average
    REAL(4) :: m_z_diff_squared_2        ! The square of (find_length_freq(i) - freq_avg)
    REAL(4) :: m_z_stdev_2             ! Standard deviation of freq_avg
    REAL(4) :: m_z_shift                 ! Frequency shift from freq_avg_1 & freq_avg_2
    
    REAL(4) :: charge_sum               ! Summation of charge for averaging purposes
    REAL(4) :: charge_diff_squared      ! The square of (charge_array(i) - charge)
    REAL(4) :: charge_std_dev           ! Standard deviation of charge
    REAL(4) :: charge_sum_1                 ! Summation of frequency for averaging purposes
    REAL(4) :: charge_avg_1                 ! Frequency average
    REAL(4) :: charge_diff_squared_1        ! The square of (find_length_freq(i) - freq_avg)
    REAL(4) :: charge_stdev_1             ! Standard deviation of freq_avg
    REAL(4) :: charge_sum_2                 ! Summation of frequency for averaging purposes
    REAL(4) :: charge_avg_2                 ! Frequency average
    REAL(4) :: charge_diff_squared_2        ! The square of (find_length_freq(i) - freq_avg)
    REAL(4) :: charge_stdev_2             ! Standard deviation of freq_avg
    REAL(4) :: charge_shift                 ! Frequency shift from freq_avg_1 & freq_avg_2
    
    INTEGER :: peak_start               ! Pointer for peak start
    INTEGER :: peak_end                 ! Pointer for peak end  
    INTEGER :: loop_end                 ! Limits the "DO WHILE" loop to not overlap the subsequent section
    INTEGER :: exit_loop                ! Exits the "DO WHILE" loop when the frequency change is greater than 5%
    INTEGER :: cycle_length             ! Length of cycle based on frequency (windowed data)
    INTEGER :: i, j, k                  ! Loop variables
    INTEGER :: j_low, j_high            ! When "peak" found, look in narrow range to see if it's actually a maximum
    INTEGER :: max_in_range             ! Maximum in fft_magnitude between j_low and j_high
    INTEGER :: n_ptr                    ! Pointer for file section
    INTEGER :: channel                  ! Pointer for channel to analyze (A, B, or test signal)       
    INTEGER :: write_fft_window         ! Pointer for deciding whether to write detailed fft info for channel 5
    INTEGER :: number_of_peaks          ! Counter for number of peaks found in each section of datafile            
    INTEGER :: local_max                ! Pointer for local maximum  
    INTEGER :: DFT_max                  ! Pointer for peak maximum for parabolic interpolation
    INTEGER :: temp_signal_end          ! Temporary pointer for end of signal
    INTEGER :: final_signal_end         ! Pointer for end of signal
    INTEGER :: full_signal_end          ! Pointer for end of signal for taking FFT of whole signal
    INTEGER :: scan_increment           ! Length of window function when performing fft_scan_start
    INTEGER :: temp_signal_start        ! Starting point for window function when performing fft_scan_start to find end of signal
    INTEGER :: points_to_average        ! Number of points averaged to find mass_to_charge, charge, and mass
    INTEGER :: inflection_point         ! Inflection point of magnitude from fft_scan_start
    INTEGER :: section_loop             ! Re-does "Scan Start" loop if channel = 3 or 4
    INTEGER :: dead_time                ! Number of points of dead time due to A250 recovering from impulse noise of trap gate pulse
    INTEGER :: freq_mult                ! Multiplier for length of FFT (base is 2^16)
    INTEGER :: inflection_pt            ! Inflection point for frequency shift analysis

    INTEGER :: length_multiplier        ! Multiplier for length of each trapping event (multiple of 1 ms or 2500 datapoints)
    INTEGER :: max_trap_time            ! Length of period for trapping event (multiple of 1 ms or 2500 datapoints)
    INTEGER :: max_datapoints           ! Length window to scan for signal (multiple of 1 ms or 2500 datapoints)
    INTEGER :: single_event             ! Counter for number of trapping events with a single ion trapped    
    INTEGER :: multiple_event           ! Counter for number of trapping events with multiple ions trapped    
    INTEGER :: array_length

    REAL(4), DIMENSION (1:1048576) :: sig_in            ! Input array
    REAL(4), DIMENSION (1:4, 1:100) :: peak_magnitudes  ! Array for peak magnitudes 
    REAL(4) :: fund_freq_win_1          ! Fundamental frequency as determined by the first window
    REAL(4) :: fund_mag_win_1           ! Magnitude of fundamental frequency as determined by the first window
    REAL(4), DIMENSION (1:4, 1:100) :: freq_values      ! Array for frequency values of max peak magnitudes
    !REAL(4), DIMENSION (1:4, 1:100) :: peak_raw         ! Array for datapoint of peak values
    INTEGER, DIMENSION (1:4, 1:100) :: peak_raw         ! Array for datapoint of peak values 
    REAL(4), DIMENSION (1:4, 1:100) :: initial_peak_magnitudes   
    REAL(4), DIMENSION (1:4, 1:100) :: initial_freq_values
    INTEGER, DIMENSION (1:4, 1:100) :: initial_peak_raw ! Analogues of above arrays for initial FT of whole trapping event
        
    CHARACTER(len=4) :: window_number                   ! Number of window being scanned across ion's trapping event
    CHARACTER(len=2) :: freq_mult_char                  ! Converting freq_mult into character format    
    REAL(4), DIMENSION (1:524288) :: freq_array         ! Array of frequencies corresponding to fft_magnitude
    INTEGER :: t                                        ! Flag used in user-selected options for channels 3 and 5
    REAL(4), DIMENSION (1:524288) :: fft_magnitude      ! DFT magnitude output
    REAL(4), DIMENSION (1:75000) :: max_fftmag          ! Array for values of max fft magnitudes (~52,000 m/z cutoff)
   ! Is the ~52,000 m/z cutoff problematic? Also, how does this determine what the cutoff is? Also, P22 w/scaffold has a higher m/z than 52,000
    REAL(4) :: DFT_mag, DFT_mag_less_one, DFT_mag_plus_one ! Values of fft_magnitude corresponding to max of peak and adjacent peaks
    REAL(4) :: interp_max                               ! Position of peak max determined by parabolic interpolation
    REAL(4) :: interp_mag                               ! Magnitude of peak max determined by parabolic interpolation
    REAL(4), DIMENSION(1:10) :: harm_freq
                                        ! Frequency of first 10 harmonics
    REAL(4), DIMENSION(1:10) :: harm_mag
                                        ! Magnitude of first 10 harmonics
            
    ! Arrays for variable max trapping time
    ! =====================================
    REAL(4), ALLOCATABLE :: find_length_mag (:)        ! Array for peak magnitudes (variable trapping time)
    REAL(4), ALLOCATABLE :: find_length_freq (:)       ! Array for frequency values of max peak magnitudes (variable trapping time)     
    REAL(4), ALLOCATABLE :: charge_array (:)           ! Array for charge values (variable trapping time)      
    REAL(4), ALLOCATABLE :: mass_to_charge_array (:)   ! Array for m/z values (variable trapping time)
    REAL(4), ALLOCATABLE :: mass_array (:)             ! Array for mass values (variable trapping time) 
    REAL(4), ALLOCATABLE :: freq_avg_array (:)         ! Array for freq_sum values
    REAL(4), ALLOCATABLE :: freq_stdev_array (:)         ! Array for freq_sum values
    INTEGER, ALLOCATABLE :: final_signal_end_array (:) ! Pointer for end of signal (variable trapping time)   
    REAL(4), ALLOCATABLE :: harm2_rel_mag (:)
    REAL(4) :: final_signal_end_array_avg
    REAL(4) :: harm2_rel_mag_avg
    REAL(4) :: harm2_rel_mag_std_dev
    REAL(4) :: harm2_rel_mag_init_avg
    REAL(4) :: harm2_rel_mag_end_avg
    
    REAL(4) :: frequency_slope                      ! Slope of best fit line through frequency vs time plot (units of Hz/data point)
    REAL(4) :: slope_top, slope_bottom              ! Numerator and denominator used to compute frequency_slope
    REAL(4) :: frequency_relative_slope             ! frequency_slope divided by average frequency of ion
    REAL(4) :: frequency_intercept                  ! Intercept of frequency's best fit line at time 0
    REAL(4) :: frequency_sum_sq                     ! Square-rooted sum of the squares of difference between measured frequency and best fit
    REAL(4) :: fit                                  ! Used to compute frequency_sum_sq
    
    REAL(8) :: square_total_real                    ! Length of test signal
    REAL(8) :: old_square_total_real                ! Length of previous test signal before frequency shift
    REAL(8) :: remainder                            ! Offset for next cycle in test signal
            
    INTEGER :: window_length 
    INTEGER :: start_of_signal          ! Length to cut off at beginning of signal before FFT
    INTEGER :: max_fftmag_start
    INTEGER :: max_fftmag_end
    REAL(4) :: mag_cutoff               ! FFT magnitude analysis cutoff value
    REAL(4), PARAMETER :: pi = 4*ATAN(1.)
        
    REAL(4) :: y_sum                    ! Sum of y_values for linear regression analysis
    REAL(4) :: x_sum                    ! Sum of x_values for linear regression analysis
    REAL(4) :: y_bar                    ! Average of y_values
    REAL(4) :: x_bar                    ! Average of x_values
    REAL(4) :: cov                      ! Covariance of x_values and y_values
    REAL(4) :: var                      ! Variance of x_values and y_values
    REAL(4) :: beta                     ! Slope of linear regression line
    REAL(4) :: alpha                    ! Y-intercept of linear regression line
    REAL(4) :: x_int                    ! X-intercept of linear regression line   
    REAL(4) :: mass_to_charge_lin_reg   ! m/z calculated from linear regression
    REAL(4) :: mass_lin_reg             ! Mass calculated from linear regression m/z & averaged charge   
    REAL(4) :: ss_err                   ! Residual sum of squares
    REAL(4) :: ss_tot                   ! Total sum of squares
    REAL(4) :: r_squared                ! R_squared value
    REAL(4) :: fit_rms                  ! rms value for linear fit
    
    REAL(4) :: low_mass_cutoff          ! Low mass cutoff for frequency information
    REAL(4) :: high_mass_cutoff         ! High mass cutoff for frequency information
    INTEGER :: sig_cutoff               ! Signal length cutoff (in datapoints) for frequency information
      
    CHARACTER(len=8) :: fmt                             ! Format descriptor for file_section
    CHARACTER(len=10) :: file_section                   ! File section to read
    
    INTEGER :: section
    
    
!    ! Create array of FFT frequency values based on frequency resolution (sampling rate/number of datapoints)
!    ! Sampling rate is 2.5 MHz (40.E6/16). Number of datapoints is 1E6, but zero-padded to 2^20 (1048576).
!    ! First frequency is 0 Hz which results in the (i - 1) on the right side.
!    DO  i = 1, 524288
!        freq(i) = (i-1) * 40.E6 / 2**24
!    END DO
    
    final_signal_end = 0
    end_time = 0
    cycles = 0
    points_to_average = 0
    freq_avg = 0
    freq_std_dev = 0
    mass_to_charge = 0
    mass_to_charge_std_dev = 0
    charge = 0
    charge_std_dev = 0
    mass = 0
    mass_std_dev = 0
    loop_end = 0
    single_event = 0
    multiple_event = 0
    t = 0

    ! Take initial FT of entire trapping event to resolve out multiple ions better. Count as multiple ion event and skip analysis     
IF (channel == 1 .OR. channel == 2 .OR. channel == 5 .OR. channel == 8 .OR. channel == 12) THEN
    initial_peak_magnitudes = 0
    initial_freq_values = 0
    initial_peak_raw = 0       
    peak_start = 0
    peak_end = 0
    local_max = 0
    
    temp_signal_end = length_multiplier*2500
    
    ! Since this initial transform gets zero-filled to 1048576 points, use freq_mult = 16. When windows begin further down, this changes
    freq_mult = 16
!    IF  (temp_signal_end <= 131072)   THEN          
!        freq_mult = 2 
!    ELSE IF (temp_signal_end > 131072 .AND. temp_signal_end <= 262144)   THEN                       
!        freq_mult = 4  
!    ELSE IF (temp_signal_end > 262144 .AND. temp_signal_end <= 524288)   THEN                     
!        freq_mult = 8 
!    ELSE IF (temp_signal_end > 524288 .AND. temp_signal_end <= 1048576)   THEN                  
!        freq_mult = 16 
!    END IF    
    ! 131072=2^17; 262144=2^18; 524288=2^19; 1048576=2^20. Base length of FFT is 2^16. When the signal exceeds 2^16, the length of the FFT
    ! must be multiplied by freq_mult to reach the required power of 2.
    
    CALL fftransform(sig_in, n_ptr, temp_signal_end, fft_magnitude, length_multiplier, dead_time)         
!    CALL fft_scan_start (sig_in, n_ptr, dead_time, (temp_signal_end-dead_time), fft_magnitude, length_multiplier, freq_mult)

!    ! The higher-resolution, initial FT of the whole event requires smoothing or else false "peaks" are found which are just jitters in a large peak
!    DO j = 5, 524284
!        fft_magnitude(j) = SUM(fft_magnitude(j-4:j+4))/9
!    END DO
    
!    j = NINT(freq_mult * 108.)   ! Set starting value of magnitude at ~4 kHz (~150,000 m/z) 
        j = NINT(freq_mult*53.)     ! Set starting value of magnitude at ~2 kHz (~627,000 m/z)
        number_of_peaks = 0
        
        ! Scan fft_magnitude file for freqencies between 2 kHz and 50 kHz (to avoid noise from the quadrupole)
        DO WHILE (j < NINT(freq_mult * 1310.8))               
!        DO WHILE (j < NINT(freq_mult * 2621.5))
            239 CONTINUE
            j = j + 1 
            
            temp_freq = (j-1) * 40.E6 / (2**20 * freq_mult)
            ! Search "frequency resolution" for explanation of frequency calculation
            
            ! Since noise increases as time^0.5, extrapolate what cutoff should be based on cutoff values for window of 1000
            ! This is approximate, but it works pretty well (adding 2e5 instead of 2e4 because why not?)
            IF (temp_freq < 2500.) THEN
!                mag_cutoff = ((temp_signal_end - dead_time)/1000)**0.5*(-476400.00485 + 805.44827*temp_freq - 0.1634*temp_freq**2) + 100000.
                mag_cutoff = ((temp_signal_end - dead_time)/1000)**0.5*(-625795.09913 + 899.86954*temp_freq - 0.17706*temp_freq**2) + 100000.
            END IF
            IF (temp_freq >= 2500. .AND. temp_freq < 5300.) THEN
!                mag_cutoff = ((temp_signal_end - dead_time)/1000)**0.5*(295659.50889 + 225114.57697*SIN(pi*(temp_freq - 911.16232)/2933.17311)) + 100000.
                mag_cutoff = ((temp_signal_end - dead_time)/1000)**0.5*(308401.28991 + 214789.20359*SIN(pi*(temp_freq - 973.19009)/2900.72165)) + 100000.
            END IF
            IF (temp_freq >= 5300.) THEN
!                mag_cutoff = ((temp_signal_end - dead_time)/1000)**0.5*(13061.96812*EXP(39151.37476/(temp_freq + 18700.51715))) + 100000.
                mag_cutoff = ((temp_signal_end - dead_time)/1000)**0.5*(13295.0054*EXP(54306.12804/(temp_freq + 23835.00743))) + 100000.
            END IF
!            mag_cutoff = ((temp_signal_end - dead_time)/1000)**0.5*16376.49992*EXP(20130.15946/(temp_freq + 7560.51439)) + 200000
!            CALL mag_cutoff_69k ((temp_signal_end-dead_time), temp_freq, mag_cutoff) 
!            CALL mag_cutoff_317k (window_length, temp_freq, mag_cutoff)
            
            ! If fft_magnitude is greater than "mag_cutoff" then find peak max, its frequency, and FFT magnitude
            IF  (fft_magnitude(j) >= mag_cutoff) THEN 
            
                IF  (number_of_peaks == 100)    THEN
                    IF (channel /= 5) GOTO 11
                END IF
                
                ! If peak found, see if it's the maximum within a 0.5% range. If not, skip. This bypasses jitter on peaks
                j_low = CEILING(0.995*j)
                j_high = FLOOR(1.005*j)
                mag_max = 0
                DO i = j_low, j_high
                    IF (fft_magnitude(i) > mag_max) THEN
                        mag_max = fft_magnitude(i)
                        max_in_range = i
                    END IF
                END DO
                IF (j /= max_in_range) GOTO 239
                
                
                number_of_peaks = number_of_peaks + 1
                local_max = j
                max = fft_magnitude(local_max)
                max_less_one = fft_magnitude(local_max - 1)
                
                ! Find local peak maximum            
                DO WHILE    (max >= max_less_one)
                            local_max = local_max + 1
                            max = fft_magnitude(local_max)
                            max_less_one = fft_magnitude(local_max - 1)
                END DO
            
                ! Peak max info stored in first column of arrays
                local_max = local_max - 1
                DFT_mag = fft_magnitude(local_max)
                DFT_mag_less_one = fft_magnitude(local_max - 1)
                DFT_mag_plus_one = fft_magnitude(local_max + 1)
                
                ! Do parabolic interpolation with maximum and adjacent points to obtain more accurate peak max
                CALL parabolic_interpolation (local_max,DFT_mag,DFT_mag_less_one,DFT_mag_plus_one,interp_mag,interp_max)
                
                initial_peak_raw(1, number_of_peaks) = interp_max
                initial_peak_magnitudes(1, number_of_peaks) = interp_mag
                initial_freq_values(1, number_of_peaks) = (interp_max-1) * 40.E6 / (2**20 * freq_mult)
!                initial_peak_raw(1, number_of_peaks) = local_max
!                initial_peak_magnitudes(1, number_of_peaks) = fft_magnitude(local_max)
!                initial_freq_values(1, number_of_peaks) = (local_max-1) * 40.E6 / (2**20 * freq_mult)
                ! Search "frequency resolution" for explanation of frequency calculation
                        
                peak_end = local_max + 1
                last = fft_magnitude(peak_end)
                last_plus_one = fft_magnitude(peak_end + 1)
                
                ! Find end of peak
                DO WHILE    (last >= last_plus_one)
                            last = fft_magnitude(peak_end)
                            last_plus_one = fft_magnitude(peak_end + 1)
                            peak_end = peak_end + 1                
                END DO          
                
                IF  (channel == 1 .OR. channel == 2 .OR. channel == 5 .OR. channel == 8 .OR. channel == 12)    THEN
                    ! Write frequency, FFT magnitude, and section number to file "Complete_Peak_Data"                     
                    WRITE (21,1000) filename, section, (temp_signal_end-dead_time), initial_freq_values(1, number_of_peaks), initial_peak_magnitudes(1, number_of_peaks)
                    1000 FORMAT (A12, TR4, I5, TR12, I6, T44, F15.2, T64, F20.5)
                END IF
                
                ! Insure that "j" continues to move forward and peaks are not counted twice
                IF  (j < peak_end) THEN
                    j = peak_end
                ELSE
                    j = j + 1
                END IF

            END IF
        END DO
                             
        ! If the number of peaks > 1, make sure the second peak is a harmonic of the first
        ! If not, discard this peak and move on to next file section
        IF  (number_of_peaks == 1)  THEN
            single_event = 1  ! Counter for number of trapping events with a single ion trapped                    
            GOTO 33
        ELSE IF (number_of_peaks > 1)   THEN
            IF (NINT((initial_freq_values(1,2)/initial_freq_values(1,1))) > 1)   THEN
        
                DO  j = 2, number_of_peaks
                    IF  (ABS(NINT((initial_freq_values(1,j)/initial_freq_values(1,1))) - (initial_freq_values(1,j)/initial_freq_values(1,1))) < 0.05) THEN
                        CONTINUE            ! Check that peaks beyond the first are harmoncis
                    ELSE
                        multiple_event = 1  ! Counter for number of trapping events with multiple ions trapped 
                        IF (channel /= 5) GOTO 11
                    END IF                   
                
                    IF (j == number_of_peaks)   THEN
                        single_event = 1  ! Counter for number of trapping events with a single ion trapped
                        GOTO 33
                    END IF
                END DO
                
            ELSE
                multiple_event = 1  !Counter for number of trapping events with multiple ions trapped 
                IF (channel /= 5) GOTO 11
            END IF  
        END IF    
        
    ! If no peaks found if full event FT, this may be because it was lost quickly and didn't rise above all that noise. Analyze anyway.    
!    IF  (number_of_peaks == 0)  THEN
!        IF (channel /= 5) GOTO 11
!    END IF    
    
    
    33 CONTINUE  
END IF     
    
    
    
    
    

    
!    ! Use this for "adaptive window" DO LOOP below.
!    ! This should only be used for looking at charges <15 or else it is just time-consuming.
!    IF  (max_trap_time > 129)   THEN
!        max_datapoints = (2500 * 129) - dead_time
!    ELSE
!        max_datapoints = (2500 * max_trap_time) - dead_time
!    END IF
    
    
    ! Start with small window of time and search for peaks in the FT. If none found, increase window length.
    ! Once peak found, move to next loop which scans that window across the trapping event, doing an FT each time.
!    DO  i = 1, NINT((72500. - 1500. - dead_time) / 1000.)        !Original
!    DO  i = 3, INT((REAL(max_datapoints) / 1000.) - 2)           !Adaptive window
!    DO  i = 3, INT(((72500. - REAL(dead_time)) / 1000.) - 7)     !For x16 Program
    DO  i = 3, INT(((72500. - REAL(dead_time)) / 1000.) - 0)     !For x8 Program
    ! How were the end points of these loops determined?
        
        
        ! Reset all values of peak-finding arrrays to zero
        peak_magnitudes = 0
        freq_values = 0
        peak_raw = 0        
        peak_start = 0
        peak_end = 0
        local_max = 0
    
        ! Initial length of signal on which to perform FFT to look for peaks 
        window_length = ((i-1) * 1000) + 1000
        temp_signal_end = ((i-1) * 1000) + 1000 + dead_time
        
        
        IF  (temp_signal_end <= 131072)   THEN          
            freq_mult = 2 
        ELSE IF (temp_signal_end > 131072 .AND. temp_signal_end <= 262144)   THEN                       
            freq_mult = 4  
        ELSE IF (temp_signal_end > 262144 .AND. temp_signal_end <= 524288)   THEN                     
            freq_mult = 8 
        ELSE IF (temp_signal_end > 524288 .AND. temp_signal_end <= 1048576)   THEN                  
            freq_mult = 16 
        END IF    
        ! 131072=2^17; 262144=2^18; 524288=2^19; 1048576=2^20. Base length of FFT is 2^16. When the signal exceeds 2^16, the length of the FFT
        ! must be multiplied by freq_mult to reach the required power of 2.
                 
        CALL fft_scan_start (sig_in, n_ptr, dead_time, (temp_signal_end-dead_time), fft_magnitude, length_multiplier, freq_mult)  
        
!       j = NINT(freq_mult * 108.)   ! Set starting value of magnitude at ~4 kHz (~150,000 m/z) 
        j = NINT(freq_mult*53.)     ! Set starting value of magnitude at ~2 kHz (~627,000 m/z)
        number_of_peaks = 0
        
        ! Scan fft_magnitude file for freqencies between 4 kHz and 50 kHz (to avoid quadrupole noise)
        DO WHILE (j < NINT(freq_mult * 1310.8))
!        DO WHILE (j < NINT(freq_mult * 2621.5))               
            j = j + 1 
            
            temp_freq = (j-1) * 40.E6 / (2**20 * freq_mult)
            ! Search "frequency resolution" for explanation of frequency calculation
            
            
            IF (temp_freq < 2500.) THEN
!                mag_cutoff = (window_length/1000)**0.5*(-476400.00485 + 805.44827*temp_freq - 0.1634*temp_freq**2) + 100000.
                mag_cutoff = (window_length/1000)**0.5*(-625795.09913 + 899.86954*temp_freq - 0.17706*temp_freq**2) + 100000.
            END IF
            IF (temp_freq >= 2500. .AND. temp_freq < 5300.) THEN
!                mag_cutoff = (window_length/1000)**0.5*(295659.50889 + 225114.57697*SIN(pi*(temp_freq - 911.16232)/2933.17311)) + 100000.
                mag_cutoff = (window_length/1000)**0.5*(308401.28991 + 214789.20359*SIN(pi*(temp_freq - 973.19009)/2900.72165)) + 100000.
            END IF
            IF (temp_freq >= 5300.) THEN
!                mag_cutoff = (window_length/1000)**0.5*(13061.96812*EXP(39151.37476/(temp_freq + 18700.51715))) + 100000.
                mag_cutoff = (window_length/1000)**0.5*(13295.0054*EXP(54306.12804/(temp_freq + 23835.00743))) + 100000.
            END IF
!            CALL mag_cutoff_69k (window_length, temp_freq, mag_cutoff)
!            CALL mag_cutoff_317k (window_length, temp_freq, mag_cutoff)
            
            ! If fft_magnitude is greater than "mag_cutoff" then find peak max, its frequency, and FFT magnitude
            IF  (fft_magnitude(j) >= mag_cutoff) THEN           
                
                IF  (number_of_peaks == 100)    THEN
                    IF (channel /= 5) GOTO 11
                END IF
                
                number_of_peaks = number_of_peaks + 1
                local_max = j
                max = fft_magnitude(local_max)
                max_less_one = fft_magnitude(local_max - 1)
                
                ! Find local peak maximum            
                DO WHILE    (max >= max_less_one)
                            local_max = local_max + 1
                            max = fft_magnitude(local_max)
                            max_less_one = fft_magnitude(local_max - 1)
                END DO
            
                ! Peak max info stored in first column of arrays
                local_max = local_max - 1
                DFT_mag = fft_magnitude(local_max)
                DFT_mag_less_one = fft_magnitude(local_max - 1)
                DFT_mag_plus_one = fft_magnitude(local_max + 1)
                
                ! Do parabolic interpolation with maximum and adjacent points to obtain more accurate peak max
                CALL parabolic_interpolation (local_max,DFT_mag,DFT_mag_less_one,DFT_mag_plus_one,interp_mag,interp_max)
                
                peak_raw(1, number_of_peaks) = interp_max
                peak_magnitudes(1, number_of_peaks) = interp_mag
                freq_values(1, number_of_peaks) = (interp_max-1) * 40.E6 / (2**20 * freq_mult)
!                peak_raw(1, number_of_peaks) = local_max
!                peak_magnitudes(1, number_of_peaks) = fft_magnitude(local_max)
!                freq_values(1, number_of_peaks) = (local_max-1) * 40.E6 / (2**20 * freq_mult)
                ! Search "frequency resolution" for explanation of frequency calculation
                       
                ! If "2nd harmonic" < 40% of "fundamental," then what program thinks is "fundamental" is probably the 2nd harmonic 
                ! rising above the threshold first. In this case, go to longer window until fundamental correctly determined.
                IF (number_of_peaks == 1) THEN
                    fund_freq_win_1 = freq_values(1, number_of_peaks)
                    fund_mag_win_1 = peak_magnitudes(1, number_of_peaks)
                    CALL Harmonic_Analysis(channel, freq_mult, fft_magnitude, fund_freq_win_1, fund_mag_win_1, temp_signal_end, harm_freq, harm_mag, t)
                    IF (harm_mag(2)/harm_mag(1) < 0.4) GOTO 382
                END IF
                
                peak_start = local_max
                start = fft_magnitude(peak_start)
                start_less_one = fft_magnitude(peak_start - 1)
                
                ! Find start of peak           
                DO WHILE    (start >= start_less_one)
                            peak_start = peak_start - 1
                            start = fft_magnitude(peak_start)
                            start_less_one = fft_magnitude(peak_start - 1)
                END DO
                
                ! Peak start info stored in second column of arrays
                peak_start = peak_start + 1 
                    ! This puts peak_start one point above the actual minimum (so on the peak)
                peak_raw(2, number_of_peaks) = peak_start
                peak_magnitudes(2, number_of_peaks) = fft_magnitude(peak_start)
                freq_values(2, number_of_peaks) = (peak_start-1) * 40.E6 / (2**20 * freq_mult)  
                        
                peak_end = local_max + 1
                last = fft_magnitude(peak_end)
                last_plus_one = fft_magnitude(peak_end + 1)
                
                ! Find end of peak
                DO WHILE    (last >= last_plus_one)
                            last = fft_magnitude(peak_end)
                            last_plus_one = fft_magnitude(peak_end + 1)
                            peak_end = peak_end + 1                
                END DO
                
                ! Peak end info stored in third column of arrays
                peak_raw(3, number_of_peaks) = peak_end
                    ! This puts peak_end one point above the actual minimum (so not on the peak)
                peak_magnitudes(3, number_of_peaks) = fft_magnitude(peak_end)
                freq_values(3, number_of_peaks) = (peak_end-1) * 40.E6 / (2**20 * freq_mult)          
                
!                IF  (channel == 1 .OR. channel == 2 .OR. channel == 5 .OR. channel == 8 .OR. channel == 12)    THEN
!                    ! Write frequency, FFT magnitude, and section number to file "Complete_Peak_Data"                     
!                    WRITE (21,1000) filename, section, window_length, freq_values(1, number_of_peaks), peak_magnitudes(1, number_of_peaks)
!                    1000 FORMAT (A12, TR4, I5, TR12, I6, T44, F15.2, T64, F20.5)
!                END IF
                
                ! Insure that "j" continues to move forward and peaks are not counted twice
                IF  (j < peak_end) THEN
                    j = peak_end
                ELSE
                    j = j + 1
                END IF

            END IF
        END DO
                             
!        ! If the number of peaks > 1, make sure the second peak is a harmonic of the first
!        ! If not, discard this peak and move on to next file section
!        IF  (number_of_peaks == 1)  THEN
!            single_event = 1  !Counter for number of trapping events with a single ion trapped                    
!            GOTO 42
!        ELSE IF (number_of_peaks > 1)   THEN
!            IF (NINT((freq_values(1,2)/freq_values(1,1))) > 1)   THEN
!        
!                DO  j = 2, number_of_peaks
!                    IF  (ABS(NINT((freq_values(1,j)/freq_values(1,1))) - (freq_values(1,j)/freq_values(1,1))) < 0.05) THEN
!                        CONTINUE
!                    ELSE
!                        multiple_event = 1  !Counter for number of trapping events with multiple ions trapped 
!                        single_event = 0
!                        GOTO 11
!                    END IF                   
!                
!                    IF (j == number_of_peaks)   THEN
!                        single_event = 1  !Counter for number of trapping events with a single ion trapped
!                        GOTO 42
!                    END IF
!                END DO
!                
!            ELSE
!                multiple_event = 1  !Counter for number of trapping events with multiple ions trapped 
!                single_event = 0
!                GOTO 11
!            END IF  
!        END IF
        ! Rely on initial, full-event FT to determine multiple ion events
        IF (number_of_peaks >= 1) GOTO 42
        382 CONTINUE
    END DO     
        
        
    IF  (number_of_peaks == 0)  THEN
        GOTO 11
    END IF    
    
    
    42 CONTINUE    
            
            
    ! Find maximum number of cycles an ion oscillationg at 178 kHz could have in maximum trapping time
    array_length = NINT(2500. * REAL(length_multiplier) * 178000. * (16. / 40.E6))
    
    ! Allocate size of arrays based on maximum number of cycles
!    ALLOCATE (find_length_freq(1:array_length), find_length_mag(1:array_length), charge_array(1:array_length), &
!             & mass_to_charge_array(1:array_length), mass_array(1:array_length), final_signal_end_array(1:array_length))
    ALLOCATE (find_length_freq(1:array_length), find_length_mag(1:array_length), charge_array(1:array_length), &
             & mass_to_charge_array(1:array_length), mass_array(1:array_length), final_signal_end_array(1:array_length), harm2_rel_mag(1:array_length))
    
    section_loop = 0
    
    108 CONTINUE

    ! Initialize all arrays to zero
    max_fftmag = 0
    find_length_freq = 0
    find_length_mag = 0
    charge_array = 0
    mass_to_charge_array = 0
    mass_array = 0
    final_signal_end_array = 0
    harm2_rel_mag = 0
    
    
    IF  (window_length < 67000) THEN
        scan_increment = window_length + 1500
    ELSE
!        scan_increment = window_length
        scan_increment = 67000
    END IF
          
        IF  (scan_increment <= 131072)   THEN 
            freq_mult = 2 
        ELSE IF (scan_increment > 131072 .AND. scan_increment <= 262144)   THEN                
            freq_mult = 4  
        ELSE IF (scan_increment > 262144 .AND. scan_increment <= 524288)   THEN                   
            freq_mult = 8 
        ELSE IF (scan_increment > 524288 .AND. scan_increment <= 1048576)   THEN                       
            freq_mult = 16 
        END IF    
        
    temp_signal_start = dead_time
    cycle_length = NINT(40.E6 / (16. * freq_values(1,1)))        
    end_time = REAL(scan_increment - 1) * REAL(16. / 40.E6)
    
        
    IF  (freq_values(1,1) > 16500)    THEN
        max_fftmag_start = INT((((freq_values(1,1)/2)-1000) * 2**20 * freq_mult / 40.E6) + 1)
    ELSE        
        max_fftmag_start = INT(((freq_values(1,1)-1000) * 2**20 * freq_mult / 40.E6) + 1)
    END IF    
    
    IF  (max_fftmag_start <= NINT(freq_mult * 108.))   THEN
        max_fftmag_start = NINT(freq_mult * 108.) + 1
    END IF
    
    ! Adjust the starting position of the FFT window until the starting point reaches max_trap_time
    ! This insures that the scan does not bleed over into the next file section
    exit_loop = 0
    i = 0
     
    IF (t == 1) t = 2
    ! Watch out for parentheses with .OR. and .AND.; the .AND. is taken into account first
    IF ((channel == 3 .OR. channel == 5) .AND. (t /= 2)) THEN
        
        t = 1       ! Prevent this construct from being accessed again. Search for "108 CONTINUE" and "GOTO 108" to see why this would be accessed twice
        PRINT *, 'Would you like file of FFT for each window of scan?'
        PRINT *, 'Press 1 for no; press 2 for yes...'
        READ *, write_fft_window
        WRITE (97,*) 'Would you like file of FFT for each window of scan?'
        WRITE (97,*) 'Press 1 for no; press 2 for yes...'
        WRITE (97,*) write_fft_window
        
        CLOSE (UNIT=97, STATUS='KEEP')
        
        WRITE (108, 345)'End_of_Signal','Harm_1_Freq','Harm_1_Mag','Harm_2_Freq','Harm2_Mag','Harm_3_Freq','Harm3_Mag','Harm_4_Freq','Harm4_Mag', &
        'Harm_5_Freq','Harm5_Mag','Harm_6_Freq','Harm6_Mag','Harm_7_Freq','Harm7_Mag','Harm_8_Freq','Harm8_Mag','Harm_9_Freq','Harm9_Mag', &
        'Harm_10_Freq','Harm10_Mag'
        345 FORMAT (TR5, A13, TR3, A11, TR4, A10, TR5, A11, TR4, A10, TR5, A11, TR4, A10, TR5, A11, TR4, A10, TR5, A11, TR4, A10, TR5, A11, TR4, A10, TR5, &
            A11, TR4, A10, TR5, A11, TR4, A10, TR5, A11, TR4, A10, TR5, A12, TR3, A11)
    END IF
    
    ! Store initial, fundamental frequency for reference in harmonic analysis subroutine.
    fund_freq_win_1 = freq_values(1,1)
    fund_mag_win_1 = peak_magnitudes(1,1)
    
    
    ! Find m/z, charge, & mass by scanning an FFT window of length "scan_increment" across the data in increments of "cycle_length"
    DO WHILE (temp_signal_start <= ((2500*max_trap_time) - scan_increment))
        
        i = i + 1
        loop_end = i 
           
        temp_signal_end = scan_increment + temp_signal_start
                                  
        CALL fft_scan_start (sig_in, n_ptr, temp_signal_start, scan_increment, fft_magnitude, length_multiplier, freq_mult)  
                      
        IF ((channel == 3 .OR. channel == 5) .AND. (t == 1)) THEN
            
!            CALL Harmonic_Analysis(channel, freq_mult, fft_magnitude, fund_freq_win_1, fund_mag_win_1, temp_signal_end, harm_freq, harm_mag, t)
          
            IF (write_fft_window == 2) THEN
                WRITE (window_number, "(I4)") i
                WRITE (freq_mult_char, '(I2)') freq_mult
                unit = 201
                unit = unit + 1
                filename = 'mag_windowed_window'//TRIM(window_number)//'.txt'
    !            filename = 'mag_window_#'//TRIM(window_number)//'freq_mult'//TRIM(freq_mult_char)//'.txt'
                CALL open_new_file (unit, filename)
            
                ! Write frequency and magnitude data for each window stepped across trapping event, if user-selected.
                freq_array(1) = 0
                WRITE (unit, '(2F25.10)') freq_array(1), fft_magnitude(1)
                DO  j = 2, (freq_mult * 32768)
                    freq_array(j) = (j - 1)*40.E6/(2**20*freq_mult)
                    WRITE (unit, '(2F25.10)') freq_array(j), fft_magnitude(j)
                END DO
                
                CLOSE (unit, STATUS = "KEEP")
            END IF
                        
        END IF     
                     
                       
!        max_fftmag = 0      
!              
!        IF  (i == 1)    THEN   
!            DO  j = (max_fftmag_start + 1), NINT(freq_mult * 2621.5)  
!                max_fftmag(j) = fft_magnitude(j)  
!            END DO
!            
!            temp_freq = (MAXLOC(max_fftmag, DIM=1)-1) * 40.E6 / (2**20 * freq_mult)
!            
!            max_fftmag_start = INT(((temp_freq * 0.9) * 2**20 * freq_mult / 40.E6) + 1)
!            max_fftmag_end = INT(((temp_freq * 1.1) * 2**20  * freq_mult / 40.E6) + 1)
!            
!            IF  (max_fftmag_start <= NINT(freq_mult * 108.))   THEN
!                max_fftmag_start = NINT(freq_mult * 108.) + 1
!            END IF
!            
!            IF  (max_fftmag_end > NINT(freq_mult * 4281.5))   THEN
!                max_fftmag_end = NINT(freq_mult * 4281.5)
!            END IF
!           
!        ELSE        
!            DO  j = (max_fftmag_start + 1), max_fftmag_end  
!                max_fftmag(j) = fft_magnitude(j)  
!            END DO
!        END IF
!         
!                                
!        ! Store freq & magnitude info into array for each change in starting point 
!        DFT_max = MAXLOC(max_fftmag, DIM=1)
!        DFT_mag = MAXVAL(max_fftmag, DIM=1)
!        DFT_mag_less_one = max_fftmag(DFT_max - 1)
!        DFT_mag_plus_one = max_fftmag(DFT_max + 1)
!        
!        ! Do parabolic interpolation with maximum and adjacent points to obtain more accurate peak max
!        CALL parabolic_interpolation (DFT_max,DFT_mag,DFT_mag_less_one,DFT_mag_plus_one,interp_mag,interp_max)

        ! Use harmonic analysis to get peak magnitudes of first two harmonics.
        CALL Harmonic_Analysis(channel, freq_mult, fft_magnitude, fund_freq_win_1, fund_mag_win_1, temp_signal_end, harm_freq, harm_mag, t)
        
        
!        find_length_freq(i) = (interp_max-1) * 40.E6 / (2**20 * freq_mult)
!        find_length_freq(i) = (MAXLOC(max_fftmag, DIM=1)-1) * 40.E6 / (2**20 * freq_mult)
!        find_length_mag(i) = interp_mag
!        find_length_mag(i) = MAXVAL(max_fftmag, DIM=1)
        find_length_freq(i) = harm_freq(1)
!        find_length_mag(i) = harm_mag(1)
        find_length_mag(i) = SUM(harm_mag(1:2))
        freq_input = find_length_freq(i)            
        cycle_length = NINT(40.E6 / (16. * freq_input))
        final_mag = find_length_mag(i)
        harm2_rel_mag(i) = harm_mag(2)/harm_mag(1)
        
        ! 2nd Generation trap (New 40 MHz Boards, Channel A, new 2SK152 FET, -171 C)
!        charge_array(i) = final_mag / (end_time * 797420.358600896)      ! Calibration factor based on time length of signal
!        charge_array(i) = final_mag / (end_time * 1264190)      ! Calibration factor based on time length of signal and addition of 2nd harmonic
!        charge_array(i) = 0.9832*final_mag / (end_time * 1264190) 
        charge_array(i) = (4.13843/4.57135)*0.9832*final_mag / (end_time * 1264190)
            ! Calibration based on 1st two harmonics, with added correction (0.9832) so average z the same as with normal calibration
            ! The fraction in front is a correction for the board installed 2014-09-12 which has a different gain. Used masses of 
            ! 2014-04-14 sample #2 of WHV measured on 2014-06-16 with old board and 2014-09-15 with new board (same FET).
                                               
        ! 2nd Generation trap
        ! Calculate m/z based on calibration factor obtained from Simion simulations
!        mass_to_charge_array(i) = (2540452974954.05) * (1/(find_length_freq(i))**2)    ! NCC TRAP 44 - Conetrap, detector tube at 0.0 V
        mass_to_charge_array(i) = (2510080000000) * (1/(find_length_freq(i))**2)    ! NCC TRAP 44 - Conetrap, detector tube at 0.6 V
        mass_array(i) = mass_to_charge_array(i) * charge_array(i)
        final_signal_end_array(i) = scan_increment + temp_signal_start
               
        ! If analyzing real data, then exit loop when frequency varies by more than 2% (sign that end of signal has been reached)
        IF  (i >= 2)    THEN
            IF  (channel == 1 .OR. channel == 2 .OR. channel == 4 .OR. channel == 8 .OR. channel == 12) THEN
                IF  ((find_length_freq(i) > (find_length_freq(1) * 1.02)) .OR. (find_length_freq(i) < (find_length_freq(1) * 0.98))) THEN
                    exit_loop = 1
                END IF
            
            ! For analysis of test signal or specific section...   
            ! If this is the first pass of this loop, it will analyze entire signal 
            ! If this the second pass of this loop, then it will analyze to find m/z, charge, & mass
            ELSE IF (channel == 3 .OR. channel == 5)    THEN
                IF  (section_loop == 1)  THEN
                    IF  ((find_length_freq(i) > (find_length_freq(1) * 1.02)) .OR. (find_length_freq(i) < (find_length_freq(1) * 0.98))) THEN
                        exit_loop = 1
                    END IF
                END IF
            END IF
        END IF
        
        IF  (exit_loop == 1)   EXIT
        
        temp_signal_start = temp_signal_start + cycle_length   
        
    END DO   
    
    
    
    IF  (channel == 7)  THEN

        points_to_average = loop_end - 10
        
        IF  (points_to_average < 100)   THEN
            GOTO 11
        END IF

        ALLOCATE(freq_avg_array(1:INT(points_to_average/50)), freq_stdev_array(1:INT(points_to_average/50)))
                
        DO  i = 1, INT(points_to_average/50)
                            
            freq_sum = 0
            DO  j = (((i-1)*50) + 1), (((i-1)*50) + 50)
                freq_sum = freq_sum + find_length_freq(j)
            END DO     
            
            freq_avg = freq_sum / 50.  
            
            freq_diff_squared = 0
            DO  j = (((i-1)*50) + 1), (((i-1)*50) + 50)
                freq_diff_squared = freq_diff_squared + ((find_length_freq(j) - freq_avg)**2)
            END DO
            
            freq_std_dev = SQRT(freq_diff_squared/49.)
            freq_avg_array(i) = freq_avg
            freq_stdev_array(i) = freq_std_dev
            
        END DO
        
        IF  (freq_stdev_array(1) > (freq_avg_array(1)/200.)) THEN        
            exit_loop = 1            
        ELSE IF (freq_stdev_array(1) <= (freq_avg_array(1)/200.) .AND. freq_stdev_array(points_to_average/50) >= (freq_avg_array(points_to_average/50)/200.))   THEN
            exit_loop = 1                
        ELSE IF (freq_stdev_array(1) < (freq_avg_array(1)/2000.) .OR. freq_stdev_array(points_to_average/50) < (freq_avg_array(points_to_average/50)/2000.))    THEN
            exit_loop = 1            
        ELSE IF (freq_stdev_array(1) <= (freq_avg_array(1)/200.) .AND. freq_stdev_array(points_to_average/50) <= (freq_avg_array(points_to_average/50)/200.))   THEN
                       
            DO  i = 1, (INT(points_to_average/50) - 1)            
                IF  (freq_avg_array(i+1) < (freq_avg_array(1) - 20)) THEN                                    
                    inflection_pt = i
                    GOTO 88                    
                END IF                 
            END DO  
        
            exit_loop = 1
            
        END IF 
        
        IF  (exit_loop == 1)    THEN
            GOTO 11
        END IF
        
        88  CONTINUE
                
                
        inflection_pt = inflection_pt * 50         

        CALL avg_stdev (find_length_freq, 1, inflection_pt, avg, stdev)
            freq_avg_1 = avg
            freq_std_dev_1 = stdev
        CALL avg_stdev (mass_array, 1, inflection_pt, avg, stdev)
            mass_avg_1 = avg
            mass_stdev_1 = stdev
        CALL avg_stdev (mass_to_charge_array, 1, inflection_pt, avg, stdev)
            m_z_avg_1 = avg
            m_z_stdev_1 = stdev
        CALL avg_stdev (charge_array, 1, inflection_pt, avg, stdev)
            charge_avg_1 = avg
            charge_stdev_1 = stdev
            
        CALL avg_stdev (find_length_freq, 1, (inflection_pt + 51), avg, stdev)
            freq_avg_2 = avg
            freq_std_dev_2 = stdev
        CALL avg_stdev (mass_array, 1, (inflection_pt + 51), avg, stdev)
            mass_avg_2 = avg
            mass_stdev_2 = stdev
        CALL avg_stdev (mass_to_charge_array, 1, (inflection_pt + 51), avg, stdev)
            m_z_avg_2 = avg
            m_z_stdev_2 = stdev
        CALL avg_stdev (charge_array, 1, (inflection_pt + 51), avg, stdev)
            charge_avg_2 = avg
            charge_stdev_2 = stdev
                        
        freq_shift = freq_avg_1 - freq_avg_2
        mass_shift = mass_avg_1 - mass_avg_2
        m_z_shift = m_z_avg_2 - m_z_avg_1
        charge_shift = charge_avg_1 - charge_avg_2
        
        
        IF  (freq_std_dev_1 > 50. .OR. freq_std_dev_2 > 50.) THEN
            GOTO 11
        END IF
        
        IF  (freq_shift < 10.)   THEN
            GOTO 11
        END IF    
                  
        WRITE (27,242)  filename, section, freq_avg_1, freq_std_dev_1, freq_avg_2, freq_std_dev_2, inflection_pt, freq_shift, &
                        & mass_avg_1, mass_stdev_1, mass_avg_2, mass_stdev_2, mass_shift, &
                        & m_z_avg_1, m_z_stdev_1, m_z_avg_2, m_z_stdev_2, m_z_shift, &
                        & charge_avg_1, charge_stdev_1, charge_avg_2, charge_stdev_2, charge_shift                       
        242 FORMAT (A12, T19, I2, TR7, F11.2, TR8, F10.3, TR7, F11.2, TR5, F10.3, TR15, I5, TR9, F9.3, TR5, &
                    & F15.2, TR5, F15.2, TR5, F15.2, TR5, F15.2, TR6, F12.3, TR5, &
                    & F11.3, TR4, F11.3, TR8, F11.3, TR4, F11.3, TR9, F9.3, TR8, &
                    & F7.2, TR13, F7.2, TR13, F7.2, TR13, F7.2, TR10, F9.3)    
        
        CALL RANDOM_SEED()
        CALL RANDOM_NUMBER(rnd)
        unit = 100 + (rnd * 1000)
        fmt = '(I2.2)'
        
        WRITE (file_section, fmt) section
        file = 'Freq_Data-'//TRIM(filename)//'-#'//TRIM(file_section)//'.txt'
             
        CALL open_new_file (unit, file)            
        
        WRITE (unit,241) "Section_#", "End_of_Signal", "Max_Freq_Win", "Mag_Win", "m/z_Win", "Charge_Win", "Mass_Win"
        WRITE (unit,300) filename
        
        DO  j = 1, loop_end
            WRITE (unit,202) section, final_signal_end_array(j), find_length_freq(j), find_length_mag(j), mass_to_charge_array(j), & 
                   & charge_array(j), mass_array(j)
        END DO        
        
        DEALLOCATE(freq_avg_array, freq_stdev_array)
           
        CLOSE (unit, STATUS = "KEEP")     
        GOTO 11  
              
    END IF    
    
    
    
    ! If program is set to analyze test signal or specific section, then it will re-run previous "DO WHILE" loop
    ! This allows a full analysis of the file, as well as an accurate output of the m/z, charge, & mass
    
    section_loop = section_loop + 1
    IF  (channel == 3 .OR. channel == 5 .OR. channel == 6)    THEN        
        
        IF  (channel == 6)  THEN 
        
            ! Open file to store info from FFT scan across data section
            CALL RANDOM_SEED()
            CALL RANDOM_NUMBER(rnd)
            unit = 100 + (rnd * 1000)
        
            fmt = '(I2.2)'
            
            WRITE (file_section, fmt) section
            file = 'Freq_Data-'//TRIM(filename)//'-#'//TRIM(file_section)//'.txt'
                 
            CALL open_new_file (unit, file)            
            
            WRITE (unit,241) "Section_#", "End_of_Signal", "Max_Freq_Win", "Mag_Win", "m/z_Win", "Charge_Win", "Mass_Win"
            241 FORMAT (T15, A10, TR5, A16, TR5, A12, TR10, A7, TR10, A7, TR10, A10, TR10, A8)
            
            ! Write name of file currently being analzyed to output files
            WRITE (unit,300) filename
            300 FORMAT (A12)         
        
            DO  i = 1, loop_end
                WRITE (unit,202) section, final_signal_end_array(i), find_length_freq(i), find_length_mag(i), mass_to_charge_array(i), & 
                       & charge_array(i), mass_array(i)
                202 FORMAT (T19, I2, TR12, I6, TR3, F20.5, F20.5, TR5, F11.3, TR10, F7.2, TR2, F20.5) 
            END DO        
                  
            CLOSE (unit, STATUS = "KEEP")   
            GOTO 11
        
        ELSE IF (channel /= 6 .AND. section_loop == 1)   THEN    
        
            DO  i = 1, loop_end
                WRITE (26,202) section, final_signal_end_array(i), find_length_freq(i), find_length_mag(i), mass_to_charge_array(i), & 
                       & charge_array(i), mass_array(i)
            END DO        
            
!            IF (channel /= 3 .OR. channel /= 5) GOTO 108
            GOTO 108
                
        END IF     
        
    END IF
    
            
    points_to_average = loop_end - 1
    
    IF  (points_to_average < 1) THEN 
        GOTO 11
    END IF
    
    final_signal_end = final_signal_end_array(points_to_average)
    end_time = (final_signal_end - (dead_time + 1)) * REAL(16. / 40.E6) 
    
    ! For ions not trapped the whole time, subtract off last 5 ms of time it was "trapped." This is because the FT will still see the
    ! ion even after it crashes since it is looking at a width of time, causing the FT magnitude to drop precipitously before the 
    ! program "loses" the ion from its 2% frequency shift (which usually happens when noise exceeds the ion's shrinking peak)
    IF (end_time < 0.005) GOTO 11
    IF (end_time < (REAL(length_multiplier) - 4)/1000) THEN
        DO i = 3, points_to_average     ! Start at 3 so if first point triggers IF block, points_to_average > 1 (can't do std dev of 1)
            IF (final_signal_end_array(i) > final_signal_end - 12500) THEN
                final_signal_end = final_signal_end_array(i - 1)
                end_time = (final_signal_end - (dead_time + 1)) * REAL(16. / 40.E6)
                points_to_average = i - 1
                EXIT
            END IF
        END DO
    END IF
    IF (points_to_average < 1) GOTO 11
      
        
    ! Calculate values to output to normal data files
    CALL avg_stdev (find_length_freq, 1, points_to_average, avg, stdev)
        freq_avg = avg
        freq_std_dev = stdev
    CALL avg_stdev (mass_array, 1, points_to_average, avg, stdev)
        mass = avg
        mass_std_dev = stdev
    CALL avg_stdev (mass_to_charge_array, 1, points_to_average, avg, stdev)
        mass_to_charge = avg
        mass_to_charge_std_dev = stdev
    CALL avg_stdev (charge_array, 1, points_to_average, avg, stdev)
        charge = avg
        charge_std_dev = stdev
    CALL avg_stdev (harm2_rel_mag, 1, points_to_average, avg, stdev)
        harm2_rel_mag_avg = avg
        harm2_rel_mag_std_dev = stdev
    
!    IF (loop_end >= 5) THEN    
!        harm2_rel_mag_init_avg = SUM(harm2_rel_mag(1:5))/5
!        harm2_rel_mag_end_avg = SUM(harm2_rel_mag((loop_end - 4):loop_end))/5
!    ELSE
!        harm2_rel_mag_init_avg = SUM(harm2_rel_mag(1:loop_end))/loop_end
!        harm2_rel_mag_end_avg = harm2_rel_mag_init_avg 
!    END IF
    
    final_signal_end_array_avg = 0.
    DO i = 1, points_to_average
        final_signal_end_array_avg = final_signal_end_array_avg + REAL(final_signal_end_array(i))
    END DO
    final_signal_end_array_avg = final_signal_end_array_avg/REAL(points_to_average)
    
    slope_top = 0.
    slope_bottom = 0.
    DO i = 1, points_to_average
        slope_top = slope_top + (REAL(final_signal_end_array(i)) - final_signal_end_array_avg)*(find_length_freq(i) - freq_avg)
        slope_bottom = slope_bottom + (REAL(final_signal_end_array(i)) - final_signal_end_array_avg)**2
    END DO
    frequency_slope = slope_top/slope_bottom
        ! This is equivalent to standard slope of best fit line equation: slope = [sum(x*y) - sum(x)*sum(y)/n] / [sum(x^2)-(sum(x))^2/n]
    frequency_relative_slope = frequency_slope/freq_avg
    frequency_intercept = freq_avg - frequency_slope*final_signal_end_array_avg
    
    frequency_sum_sq = 0.
    DO i = 1, points_to_average
        fit = frequency_intercept + frequency_slope*REAL(final_signal_end_array(i))
        frequency_sum_sq = frequency_sum_sq + (find_length_freq(i) - fit)**2
    END DO
    
    frequency_sum_sq = SQRT(frequency_sum_sq/REAL(points_to_average))
        
        
!    ! The following is for linear regression analysis
!    ! =============================================== 
!
!    ! Find sum of all points to average
!    y_sum = 0
!    x_sum = 0
!    DO  i = 1, points_to_average
!        y_sum = y_sum + find_length_freq(i)      
!        x_sum = x_sum + ((final_signal_end_array(i) - (dead_time + 1)) * REAL(16. / 40.E6))
!!        x_sum = x_sum + final_signal_end_array(i)  
!    END DO
!    
!    
!    ! Find average for the following variables
!    y_bar = y_sum / REAL(points_to_average)
!    x_bar = x_sum / REAL(points_to_average)        
!    
!    
!    ! Find the sum of the square of the difference between average and individual value
!    cov = 0
!    var = 0  
!    DO  i = 1, points_to_average
!        cov = cov + ( (((final_signal_end_array(i) - (dead_time + 1)) * REAL(16. / 40.E6)) - x_bar) * (find_length_freq(i) - y_bar) )
!        var = var + ( (((final_signal_end_array(i) - (dead_time + 1)) * REAL(16. / 40.E6)) - x_bar)**2 )  
!!        cov = cov + ( (final_signal_end_array(i) - x_bar) * (find_length_freq(i) - y_bar) )
!!        var = var + ( (final_signal_end_array(i) - x_bar)**2 )
!    END DO
!    
!    
!    !   Beta is slope
!    !   Alpha is y-intercept
!    !   x_int is where y = 0        
!    beta = cov/var
!    alpha = y_bar - (beta * x_bar)
!    x_int = -alpha/beta
!    
!    ss_err = 0
!    ss_tot = 0 
!    DO  i = 1, points_to_average
!!        ss_err = ss_err + ((find_length_freq(i) - ((beta * final_signal_end_array(i)) + alpha))**2)
!        ss_err = ss_err + ((find_length_freq(i) - ((beta * ((final_signal_end_array(i) - (dead_time + 1)) * REAL(16. / 40.E6))) + alpha))**2)
!        ss_tot = ss_tot + ((find_length_freq(i) - y_bar)**2)
!    END DO
!    
!    r_squared = 1 - (ss_err / ss_tot) 
!    fit_rms = SQRT(ss_err / points_to_average) 
!    
!!    mass_to_charge_lin_reg = (2540452974954.05) * (1/alpha**2)   ! NCC TRAP 44 - Conetrap
!    mass_to_charge_lin_reg = (2510080000000) * (1/alpha**2)   ! NCC TRAP 44 - Conetrap
!!    mass_to_charge_lin_reg = (7505193258391.12) * (1/alpha**2)   ! NCC TRAP 82 - Conetrap 
!    mass_lin_reg = mass_to_charge_lin_reg * charge 
!
!    ! ==============================
!    ! End linear regression analysis
     
     
            
    !  The following IF statement outputs freq data files for ions trapped for the max trap time
    !  Must undo multithreading in Trap_Analysis main program code or errors will occur 
    !  ==========================================================================================   
    IF  (channel == 8)  THEN                     
        IF  (mass > low_mass_cutoff .AND. mass < high_mass_cutoff .AND. final_signal_end > sig_cutoff)  THEN   ! Variable
            
            ! Open file to store info from FFT scan across data section
            CALL RANDOM_SEED()
            CALL RANDOM_NUMBER(rnd)
            unit = 100 + (rnd * 100000)
        
            fmt = '(I2.2)'
            
            WRITE (file_section, fmt) section
            file = 'Freq_Data-'//TRIM(filename)//'-#'//TRIM(file_section)//'.txt'
                 
            CALL open_new_file (unit, file)            
            
            WRITE (unit,241) "Section_#", "End_of_Signal", "Max_Freq_Win", "Mag_Win", "m/z_Win", "Charge_Win", "Mass_Win"
            
            ! Write name of file currently being analyzed to output files
            WRITE (unit,300) filename
        
!            DO  i = 1, loop_end
            DO i = 1, points_to_average
                WRITE (unit,202) section, final_signal_end_array(i), find_length_freq(i), find_length_mag(i), mass_to_charge_array(i), & 
                       & charge_array(i), mass_array(i)
            END DO        
                  
            CLOSE (unit, STATUS = "KEEP")   
        END IF
    END IF
    !  ==========================================================================================
      
                  
    cycles = freq_avg * end_time
            
    DEALLOCATE (find_length_freq, find_length_mag, charge_array, mass_to_charge_array, mass_array, final_signal_end_array, harm2_rel_mag)
           
         
    ! Write out info to output files    
    IF  (channel == 1 .OR. channel == 2 .OR. channel == 4 .OR. channel == 8 .OR. channel == 12)    THEN    
!        WRITE (22,2500)  filename, section, final_signal_end, end_time, cycles, points_to_average, freq_avg, freq_std_dev, &
!                         & mass_to_charge, mass_to_charge_std_dev, charge, charge_std_dev, mass, mass_std_dev
!        2500 FORMAT (A12, T19, I2, TR11, I6, TR10, F11.9, TR4, F8.2, TR5, I6, TR9, F11.2, TR5, F10.3, TR5, &
!                    & F11.3, TR4, F11.3, TR8, F7.2, TR5, F7.2, TR5, F15.2, TR2, F15.2, TR5)
!        WRITE (22,2500)  filename, section, end_time, points_to_average, freq_avg, &
!                         & mass_to_charge, mass_to_charge_std_dev, charge, charge_std_dev, mass, &
!                         & harm2_rel_mag_avg, harm2_rel_mag_std_dev
!        2500 FORMAT (A12, T19, I2, TR8, F11.9, TR7, I6, TR9, F11.2, TR3, &
!                    & F11.3, TR4, F11.3, TR6, F7.2, TR5, F7.2, TR5, F15.2, TR2, &
!                    & F22.6, TR2, F22.6, TR2)
!        WRITE (22,2500)  filename, section, end_time, points_to_average, freq_avg, frequency_slope, frequency_intercept, frequency_sum_sq, &
!                         & mass_to_charge, mass_to_charge_std_dev, charge, charge_std_dev, mass, &
!                         & harm2_rel_mag_avg, harm2_rel_mag_std_dev
!        2500 FORMAT (A12, T19, I2, TR8, F11.9, TR7, I6, TR9, F11.2, TR7, ES12.5, TR2, F11.2, TR3, F11.5, TR6 &
!                    & F11.3, TR4, F11.3, TR6, F7.2, TR5, F7.2, TR5, F15.2, TR2, &
!                    & F22.6, TR2, F22.6, TR2)
        WRITE (22,2500)  filename, section, end_time, points_to_average, freq_avg, frequency_slope, frequency_relative_slope, frequency_sum_sq, &
                         & mass_to_charge, mass_to_charge_std_dev, charge, charge_std_dev, mass, &
                         & harm2_rel_mag_avg, harm2_rel_mag_std_dev
        2500 FORMAT (A12, T19, I2, TR8, F11.9, TR7, I6, TR9, F11.2, TR7, ES12.5, TR7, ES12.5, TR3, F11.5, TR6 &
                    & F11.3, TR4, F11.3, TR6, F7.2, TR5, F7.2, TR5, F15.2, TR2, &
                    & F22.6, TR2, F22.6, TR2)
                    
        ! Only write to "filtered" data file if ion trapped the whole time and its 2nd harmonic ratio is palatable            
        IF (end_time > (REAL(length_multiplier) - 4)/1000 .AND. harm2_rel_mag_avg > 0.5 .AND. harm2_rel_mag_avg < 0.7) THEN
            WRITE (23,2501)  filename, section, end_time, points_to_average, freq_avg, frequency_slope, frequency_relative_slope, frequency_sum_sq, &
                         & mass_to_charge, mass_to_charge_std_dev, charge, charge_std_dev, mass, &
                         & harm2_rel_mag_avg, harm2_rel_mag_std_dev
            2501 FORMAT (A12, T19, I2, TR8, F11.9, TR7, I6, TR9, F11.2, TR7, ES12.5, TR7, ES12.5, TR3, F11.5, TR6 &
                    & F11.3, TR4, F11.3, TR6, F7.2, TR5, F7.2, TR5, F15.2, TR2, &
                    & F22.6, TR2, F22.6, TR2)
        END IF
    END IF
              
    ! Print select values to screen
    WRITE (*,2700) TRIM(filename), "#", section, "m/z:", mass_to_charge, "z:", charge, "m:", mass, "time:", end_time
    2700 FORMAT (A12, TR1, A1, I2, TR2, A4, TR1, F9.2, TR2, A2, TR1, F7.2, TR2, A2, F15.2, TR2, A5, TR1, F8.6)
     
    11  CONTINUE         
        
    END SUBROUTINE  peakfinder_VI
