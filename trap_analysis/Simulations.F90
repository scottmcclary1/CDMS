    SUBROUTINE simulations(channel, data_type, fmt, length_multiplier, max_trap_time, load_time, dead_time, expt_num)
    
    IMPLICIT NONE
            
    INTEGER :: i, j, k, m, n                        ! Counting variables for DO loops
    INTEGER :: unit                                 ! Dummy variable for opening files
    INTEGER :: stat                                 ! Error message when opening files
    
    INTEGER :: data_type                            ! Pointer for .dat or .txt files to analyze
    INTEGER :: noise_channel                        ! Pointer for noise channel
    INTEGER :: number_of_files                      ! Number of files to analyze
    INTEGER :: test_signal_option                   ! Indicator for SIMION test signal or generic test signal 
    INTEGER :: number_of_ions                       ! Number of ions to analyze for each charge state in simulations
    INTEGER :: files_analyzed                       ! Counter for number of files analyzed
    INTEGER :: channel                              ! Pointer for channel to analyze (A, B, or test signal)  
    INTEGER :: n_main                               ! Pointer for file section
    INTEGER :: section
    INTEGER :: simulation_option
    INTEGER :: number_of_charge_states
    INTEGER :: expt_num                             ! Maximum number of experiments per file based on signal length (length_multiplier)
    INTEGER :: array_pointer                        ! Counter for creating test signal with frequency shift
    INTEGER :: signal_end                           ! Length of test signal to add onto noise file 
    INTEGER :: single_event                 ! Counter for number of trapping events with a single ion trapped    
    INTEGER :: multiple_event               ! Counter for number of trapping events with multiple ions trapped    
    INTEGER :: dead_time                            ! Number of points of dead time due to A250 recovering from impulse noise of trap gate pulse
    
    
    REAL(8) :: y_value                                      ! y-values for test signal
    REAL(8) :: z_value                                      ! z-values for test signal
    REAL(8) :: tof                                          ! TOF for test signal
    REAL(8) :: potential                                    ! Potential values for test signal
    REAL(8) :: R1                                           ! Bi-linear interpolation variable
    REAL(8) :: R2                                           ! Bi-linear interpolation variable
    
    INTEGER :: index_z1y1                                   ! Index for interpolation of potential values for test signal
    INTEGER :: index_z2y1                                   ! Index for interpolation of potential values for test signal
    INTEGER :: index_z1y2                                   ! Index for interpolation of potential values for test signal
    INTEGER :: index_z2y2                                   ! Index for interpolation of potential values for test signal
    
    REAL(8), DIMENSION(1:3, 1:15311) :: potential_table     ! Array for y-values, z-values, and potential values
    
    REAL(4) :: mass_to_charge           ! Average m/z value (windowed FFT)
    REAL(4) :: mass                     ! Average mass value (windowed FFT)    
    REAL(4) :: charge                   ! Average charge value (windowed FFT)
    REAL(4) :: charge_sum
    REAL(4) :: charge_avg
    REAL(4) :: charge_diff_squared
    REAL(4) :: charge_exp_diff_squared
    REAL(4) :: charge_std_dev           ! Standard deviation of charge (windowed FFT)
    REAL(4) :: charge_SD_exp_mean
    
    INTEGER :: derivative_peak                      ! Inflection point of derivative
    INTEGER :: der_opt                              ! Derivative option for looking for positive or negative derivative 
    INTEGER :: length_multiplier                    ! Multiplier for length of each trapping event (multiple of 1 ms or 2500 datapoints)
    INTEGER :: max_trap_time                        ! Length of period for trapping event (multiple of 1 ms or 2500 datapoints)
    INTEGER :: load_time                            ! Time between end of one trappng event and beginning of next (multiple of 1 ms or 2500 datapoints)
       
    REAL(4) :: max_cycles                           ! Maximum number of cycles based on frequency
    REAL(4) :: input_freq                           ! Frequency of test signal
    REAL(4) :: freq_shift                           ! Frequency shift of test signal (in percent)
    REAL(4) :: electrons                            ! Number of charges on test signal (in electrons)
    REAL(4) :: rand_1                               ! Random number for random starting point of test signal
    REAL(4) :: rand_2                               ! Random number for random starting point of test signal
    REAL(4) :: rand_3                               ! Random number for random starting point of test signal
    REAL(4) :: rand_4                               ! Random number for random starting point of test signal
    REAL(4) :: rand_5                               ! Random number for random starting point of test signal
    REAL(4) :: rand_6                               ! Random number for random starting point of test signal
    REAL(4) :: rand_7                               ! Random number for random starting point of test signal
    REAL(4) :: rand_8                               ! Random number for random starting point of test signal
    REAL(8) :: A1                                   ! Variable to define test signal (obtained from Origin)
    REAL(8) :: A2                                   ! Variable to define test signal (obtained from Origin)
    REAL(8) :: x0                                   ! Variable to define test signal (obtained from Origin)
    REAL(8) :: dx                                   ! Variable to define test signal (obtained from Origin)
    REAL(8) :: x                                    ! Variable to define test signal (obtained from Origin)
    
    REAL(8) :: square_total_real                    ! Length of test signal
    REAL(8) :: old_square_total_real                ! Length of previous test signal before frequency shift
    REAL(8) :: remainder                            ! Offset for next cycle in test signal
    INTEGER :: mult                                 ! Number of cycles a section of datafile (56250/square_total_real)
    INTEGER :: start                                ! Random start of cycle for test signal as determined by random number generator
     
    INTEGER(2), DIMENSION (1:1000000) :: input_array  
    REAL(4), DIMENSION (1:1048576) :: y_input                ! Input datapoints from trap data
    REAL(4), DIMENSION (1:1048576) :: adj_input              ! Input datapoints from trap data
    REAL(4), DIMENSION (1:1048576) :: derivative             ! First derivative of y_input
    REAL(4), DIMENSION (1:1048576) :: y_filter               ! Trap data after passing through 10 kHz high-pass filter               
    REAL(4), DIMENSION (1:1048576) :: fft_filter_output      ! Trap data after passing through 10 kHz high-pass FFT filter     
    REAL(4), DIMENSION (1:1048576) :: sq_wave_array     ! Test signal array (not actually a square wave)
    REAL(4), DIMENSION (1:8) :: random_array            ! Array of 5 random numbers
      
    REAL(4),ALLOCATABLE :: charge_avg_array (:)         ! Array for charge values to determine average
    REAL(4), ALLOCATABLE :: sq_wave_array_allocatable(:,:)  ! Allocatable test signal array (not actually a square wave)
    REAL(4), ALLOCATABLE :: output(:)                       ! Output for optimization purposes   
    
    REAL(4) :: low_mass_cutoff          ! Low mass cutoff for frequency information
    REAL(4) :: high_mass_cutoff         ! High mass cutoff for frequency information
    INTEGER :: sig_cutoff               ! Signal length cutoff (in datapoints) for frequency information
        
    CHARACTER(len=100) :: filename                      ! Dummy variable for opening files
    CHARACTER(len=100) :: scratch_char                  ! Dummy variable to dicard header lines
    CHARACTER(len=8) :: fmt                             ! Format descriptor for file_number
    CHARACTER(len=10) :: file_number                     ! File number to read
    
    REAL(4), DIMENSION (1:1048576) :: input                ! Input datapoints from trap data
    INTEGER, DIMENSION (1:2,1:400) :: derivative_peaks    ! Array for location of derivative peaks from derivative of y_input
    INTEGER :: pos_peaks                ! Counter for positive peaks in derivative of y_input
    INTEGER :: neg_peaks                ! Counter for negative peaks in derivative of y_input
    
    
    ! Variables for "streaking simulation"
    REAL(4) :: low_mz, high_mz, low_m, high_m           ! User-input ranges of m/z and m for streaking simulation
    REAL(4), ALLOCATABLE :: rnd_m(:), rnd_mz(:)         ! Random arrays for determining m and m/z of each ion
    REAL(4), ALLOCATABLE :: m_array(:), mz_array(:)     ! Arrays holding input m and m/z of each ion in streaking simulation
    REAL(4), ALLOCATABLE :: charge_in(:), freq_in(:)    ! Arrays holding input z and freq of each ion, determined from m_array and mz_array
    
        ! Select channel for noise files
        PRINT *, ""
        PRINT *, "Are noise files from channel A or channel B?"
        PRINT *, "Enter 1 for A..."
        PRINT *, "Enter 2 for B..."
        PRINT *, "Enter 3 for no noise..."
        READ *, noise_channel
        WRITE (97,*) ""
        WRITE (97,*) "Are noise files from channel A or channel B?"
        WRITE (97,*) "Enter 1 for A..."
        WRITE (97,*) "Enter 2 for B..."
        WRITE (97,*) "Enter 3 for no noise..."
        WRITE (97,*) noise_channel
        
        IF  (noise_channel < 3) THEN  
            ! User selected option for number of files in each run
            PRINT *, "Enter number of noise files..."
            READ *, number_of_files 
            WRITE (97,*) "Enter number of noise files..."
            WRITE (97,*) number_of_files
        END IF 
    
        ! User selected option for SIMION test signal or generica test signal
        PRINT *, "Do you want to use a SIMION file or create a generic test signal?"
        PRINT *, "Enter 1 for SIMION file (entitled test_signal.txt) (NOT SUPPORTED)..."
        PRINT *, "Enter 2 for generic test signal..."
        READ *, test_signal_option
        WRITE (97,*) "Do you want to use a SIMION file or create a generic test signal?"
        WRITE (97,*) "Enter 1 for SIMION file (entitled test_signal.txt) (NOT SUPPORTED)..."
        WRITE (97,*) "Enter 2 for generic test signal..."
        WRITE (97,*) test_signal_option
        
        ! Select type of simulation to perfom
        PRINT *, "What do you want to simulate?"
        PRINT *, "Enter 1 for charge states 1-100 (random length & frequency)..."
        PRINT *, "Enter 2 for BSA monomer..."
        PRINT *, "Enter 3 for BSA dimer..."
        PRINT *, "Enter 4 for Pyruvate Kinase..."
        PRINT *, "Enter 5 for Cytochrome C..."
        PRINT *, "Enter 6 for Myoglobin..."
        PRINT *, "Enter 7 for Alcohol Dehydrogenase..."
        PRINT *, "Enter 8 for Ubiquitin..."
        PRINT *, "Enter 9 for Streaking Simulation..."
        READ *, simulation_option   
        WRITE (97,*) "What do you want to simulate?"
        WRITE (97,*) "Enter 1 for charge states 1-100 (random length & frequency)..."
        WRITE (97,*) "Enter 2 for BSA monomer..."
        WRITE (97,*) "Enter 3 for BSA dimer..."
        WRITE (97,*) "Enter 4 for Pyruvate Kinase..."
        WRITE (97,*) "Enter 5 for Cytochrome C..."
        WRITE (97,*) "Enter 6 for Myoglobin..."
        WRITE (97,*) "Enter 7 for Alcohol Dehydrogenase..."
        WRITE (97,*) "Enter 8 for Ubiquitin..."
        WRITE (97,*) "Enter 9 for Streaking Simulation..."
        WRITE (97,*) simulation_option      
        
                
!        ! Opens output file for test signal
!        unit = 46
!        filename = "test_signal only.txt"
!        CALL open_new_file (unit, filename)
        
!        unit = 101
!        filename = "mag_full_windowed.txt"
!        CALL open_new_file (unit, filename)
        
        
        IF  (simulation_option == 1)    THEN
            number_of_charge_states = 100
            
            unit = 227
            filename = "charge_states_1-100_(42)_accuracy.txt"
            CALL open_new_file (unit, filename)
                               
            ! Headers for file #227
            WRITE (227,243) "Charge_Input", "Charge_Avg", "Charge_Std_Dev", "Charge_Std_Dev_Exp_Mean", "_n_"
            243 FORMAT (TR5, A12, TR6, A10, TR10, A14, TR10, A23, TR5, A3)     
    
            filename = "charge_states_1-100_(42)_verbose.txt"
            
        ELSE IF (simulation_option == 2)    THEN
            number_of_charge_states = 25
            filename = "bsa_monomer_simulation.txt"
        ELSE IF (simulation_option == 3)    THEN
            number_of_charge_states = 49
            filename = "bsa_dimer_simulation.txt"
        ELSE IF (simulation_option == 4)    THEN
            number_of_charge_states = 8
            filename = "pyruvate_kinase_simulation.txt"
        ELSE IF (simulation_option == 5)    THEN
            number_of_charge_states = 12
            filename = "cytochrome_C_simulation.txt"
        ELSE IF (simulation_option == 6)    THEN
            number_of_charge_states = 16 
            filename = "myoglobin_simulation.txt"
        ELSE IF (simulation_option == 7)    THEN
            number_of_charge_states = 17
            filename = "alcohol_dehydrogenase_simulation.txt"  
        ELSE IF (simulation_option == 8)    THEN
            number_of_charge_states = 7
            filename = "ubiquitin_simulation.txt"  
        ELSE IF (simulation_option == 9)    THEN
            number_of_charge_states = 1
            filename = "streaking_simulation.txt"  
        END IF
        
        unit = 22
        CALL open_new_file (unit, filename)
        
        ! Headers for file #22
!        WRITE (22,230) "File", "Section_#", " Signal_End ", "End_Time_(s)", "Cycles", "#_Points_Avgd", "Freq_Avg", &
!                       & "Freq_Std_Dev", "m/z", "m/z_Std_Dev", "Charge", "Charge_Std_Dev", "Mass", "Mass_Std_Dev"
!!                       & "Freq_(y-int)", "Slope", "R_Squared", "Fit_RMS", "m/z_Lin_Reg", "Mass_Lin_Reg", &
!!                       & "m/z_full", "Charge_full", "Mass_full"
!        230 FORMAT (T3, A4, T15, A10, TR4, A13, TR6, A12, TR5, A6, TR3, A13, TR7, A8, TR6, &
!                   & A12, TR8, A3, TR9, A11, TR6, A6, TR3, A14, TR8, A4, TR8, A13, TR5)
!!                   & A12, TR7, A5, TR8, A9, TR7, A7, TR5, A11, TR7, A13, &
!!                   & A8, TR8, A11, TR7, A9) 
        WRITE (22,230) "File", "Section_#", "End_Time_(s)", "#_Points_Avgd", "Freq_Avg", &
           & "m/z", "m/z_Std_Dev", "Charge", "Charge_Std_Dev", "Mass", &
           & "Harm2_Rel_Mag_Avg", "Harm2_Rel_Mag_Std_Dev"
        230 FORMAT (T3, A4, T15, A10, TR4, A12, TR5, A13, TR7, A8, TR6, &
               & A3, TR9, A11, TR6, A6, TR3, A14, TR8, A4, TR8, A22, TR2, A22, TR2)       
        
        ! Options for creating mock signal to evaluate accuracy & precision
        DO  k = 1, number_of_charge_states
        
            freq_shift = 0
            
            IF  (simulation_option == 1)    THEN
                electrons = k
            ELSE IF (simulation_option == 2)    THEN
                electrons = 53 - k
!                input_freq = SQRT(2540452974954.05/(66766./electrons))
                input_freq = SQRT(2510080000000/(66766./electrons))
            ELSE IF (simulation_option == 3)    THEN
                electrons = 105 - k
!                input_freq = SQRT(2540452974954.05/(66766. *2./electrons))
                input_freq = SQRT(2510080000000/(66766. *2./electrons))
            ELSE IF (simulation_option == 4)    THEN
                electrons = 36 - k
!                input_freq = SQRT(2540452974954.05/(232000./electrons))
                input_freq = SQRT(2510080000000/(232000./electrons))
            ELSE IF (simulation_option == 5)    THEN
                electrons = 20 - k
!                input_freq = SQRT(2540452974954.05/(12384./electrons))
                input_freq = SQRT(2510080000000/(12384./electrons))
            ELSE IF (simulation_option == 6)    THEN
                electrons = 28 - k
!                input_freq = SQRT(2540452974954.05/(17083./electrons))
                input_freq = SQRT(2510080000000/(17083./electrons))
            ELSE IF (simulation_option == 7)    THEN
                electrons = 47 - k
!                input_freq = SQRT(2540452974954.05/(35250./electrons))
                input_freq = SQRT(2510080000000/(35250./electrons))
            ELSE IF (simulation_option == 8)    THEN
                electrons = 13 - k
!                input_freq = SQRT(2540452974954.05/(8565./electrons))
                input_freq = SQRT(2510080000000/(8565./electrons))
            END IF
                            
            IF (simulation_option /= 9) THEN 
                WRITE (22, 245) "Charge_Input:", electrons 
                245 FORMAT (A13, TR5, F5.0)
            END IF
            
            
            IF  (simulation_option == 1)    THEN
                number_of_ions = 42
                ALLOCATE (charge_avg_array(1:number_of_ions))
                charge_avg_array = 0
            ELSE IF (simulation_option == 2)    THEN
            
                ! BSA monomer ion count
                IF  (k == 1)    THEN
                    number_of_ions = 19
                ELSE IF (k == 2)    THEN
                    number_of_ions = 22
                ELSE IF (k == 3)    THEN
                    number_of_ions = 25
                ELSE IF (k == 4)    THEN
                    number_of_ions = 33
                ELSE IF (k == 5)    THEN
                    number_of_ions = 41
                ELSE IF (k == 6)    THEN
                    number_of_ions = 51
                ELSE IF (k == 7)    THEN
                    number_of_ions = 58
                ELSE IF (k == 8)    THEN
                    number_of_ions = 66
                ELSE IF (k == 9)    THEN
                    number_of_ions = 77
                ELSE IF (k == 10)    THEN
                    number_of_ions = 88
                ELSE IF (k == 11)    THEN
                    number_of_ions = 94
                ELSE IF (k == 12)    THEN
                    number_of_ions = 98
                ELSE IF (k == 13)    THEN
                    number_of_ions = 100
                ELSE IF (k == 14)    THEN
                    number_of_ions = 100
                ELSE IF (k == 15)    THEN
                    number_of_ions = 93
                ELSE IF (k == 16)    THEN
                    number_of_ions = 89
                ELSE IF (k == 17)    THEN
                    number_of_ions = 82
                ELSE IF (k == 18)    THEN
                    number_of_ions = 76
                ELSE IF (k == 19)    THEN
                    number_of_ions = 73
                ELSE IF (k == 20)    THEN
                    number_of_ions = 67
                ELSE IF (k == 21)    THEN
                    number_of_ions = 59
                ELSE IF (k == 22)    THEN
                    number_of_ions = 53
                ELSE IF (k == 23)    THEN
                    number_of_ions = 44
                ELSE IF (k == 24)    THEN
                    number_of_ions = 38
                ELSE IF (k == 25)    THEN
                    number_of_ions = 32
                END IF
               
                number_of_ions = number_of_ions * 10
                
!                number_of_ions = 200
                            
            ELSE IF (simulation_option == 3)    THEN
               
                ! BSA dimer ion count 
                IF  (k == 1)    THEN
                    number_of_ions = 19
                ELSE IF (k == 2)    THEN
                    number_of_ions = 20
                ELSE IF (k == 3)    THEN
                    number_of_ions = 22
                ELSE IF (k == 4)    THEN
                    number_of_ions = 23
                ELSE IF (k == 5)    THEN
                    number_of_ions = 25
                ELSE IF (k == 6)    THEN
                    number_of_ions = 28
                ELSE IF (k == 7)    THEN
                    number_of_ions = 33
                ELSE IF (k == 8)    THEN
                    number_of_ions = 37
                ELSE IF (k == 9)    THEN
                    number_of_ions = 41
                ELSE IF (k == 10)    THEN
                    number_of_ions = 46
                ELSE IF (k == 11)    THEN
                    number_of_ions = 51
                ELSE IF (k == 12)    THEN
                    number_of_ions = 54
                ELSE IF (k == 13)    THEN
                    number_of_ions = 58
                ELSE IF (k == 14)    THEN
                    number_of_ions = 62
                ELSE IF (k == 15)    THEN
                    number_of_ions = 66
                ELSE IF (k == 16)    THEN
                    number_of_ions = 71
                ELSE IF (k == 17)    THEN
                    number_of_ions = 77
                ELSE IF (k == 18)    THEN
                    number_of_ions = 83
                ELSE IF (k == 19)    THEN
                    number_of_ions = 88
                ELSE IF (k == 20)    THEN
                    number_of_ions = 91
                ELSE IF (k == 21)    THEN
                    number_of_ions = 94
                ELSE IF (k == 22)    THEN
                    number_of_ions = 96
                ELSE IF (k == 23)    THEN
                    number_of_ions = 98
                ELSE IF (k == 24)    THEN
                    number_of_ions = 99
                ELSE IF (k == 25)    THEN
                    number_of_ions = 100
                ELSE IF (k == 26)    THEN
                    number_of_ions = 100
                ELSE IF (k == 27)    THEN
                    number_of_ions = 100
                ELSE IF (k == 28)    THEN
                    number_of_ions = 96
                ELSE IF (k == 29)    THEN
                    number_of_ions = 93
                ELSE IF (k == 30)    THEN
                    number_of_ions = 91
                ELSE IF (k == 31)    THEN
                    number_of_ions = 89
                ELSE IF (k == 32)    THEN
                    number_of_ions = 86
                ELSE IF (k == 33)    THEN
                    number_of_ions = 82
                ELSE IF (k == 34)    THEN
                    number_of_ions = 79
                ELSE IF (k == 35)    THEN
                    number_of_ions = 76
                ELSE IF (k == 36)    THEN
                    number_of_ions = 75
                ELSE IF (k == 37)    THEN
                    number_of_ions = 73
                ELSE IF (k == 38)    THEN
                    number_of_ions = 70
                ELSE IF (k == 39)    THEN
                    number_of_ions = 67
                ELSE IF (k == 40)    THEN
                    number_of_ions = 63
                ELSE IF (k == 41)    THEN
                    number_of_ions = 59
                ELSE IF (k == 42)    THEN
                    number_of_ions = 56
                ELSE IF (k == 43)    THEN
                    number_of_ions = 53
                ELSE IF (k == 44)    THEN
                    number_of_ions = 49
                ELSE IF (k == 45)    THEN
                    number_of_ions = 44
                ELSE IF (k == 46)    THEN
                    number_of_ions = 41
                ELSE IF (k == 47)    THEN
                    number_of_ions = 38
                ELSE IF (k == 48)    THEN
                    number_of_ions = 35
                ELSE IF (k == 49)    THEN
                    number_of_ions = 32
                END IF  
                
                number_of_ions = number_of_ions * 10
                
!                number_of_ions = 200
                            
            ELSE IF (simulation_option == 4)    THEN  
            
                ! Pyruvate Kinase tetramer ion count
                IF  (k == 1)    THEN
                    number_of_ions = 12
                ELSE IF (k == 2)    THEN
                    number_of_ions = 53
                ELSE IF (k == 3)    THEN
                    number_of_ions = 100
                ELSE IF (k == 4)    THEN
                    number_of_ions = 88
                ELSE IF (k == 5)    THEN
                    number_of_ions = 47
                ELSE IF (k == 6)    THEN
                    number_of_ions = 35
                ELSE IF (k == 7)    THEN
                    number_of_ions = 18
                ELSE IF (k == 8)    THEN
                    number_of_ions = 9
                END IF
            
                number_of_ions = number_of_ions * 10
                
!                number_of_ions = 200
                            
            ELSE IF (simulation_option == 5)    THEN  
            
                ! Cytochrome C ion count 
                IF  (k == 1)    THEN
                    number_of_ions = 8
                ELSE IF (k == 2)    THEN
                    number_of_ions = 10
                ELSE IF (k == 3)    THEN
                    number_of_ions = 25
                ELSE IF (k == 4)    THEN
                    number_of_ions = 40
                ELSE IF (k == 5)    THEN
                    number_of_ions = 63
                ELSE IF (k == 6)    THEN
                    number_of_ions = 100
                ELSE IF (k == 7)    THEN
                    number_of_ions = 100
                ELSE IF (k == 8)    THEN
                    number_of_ions = 60
                ELSE IF (k == 9)    THEN
                    number_of_ions = 40
                ELSE IF (k == 10)    THEN
                    number_of_ions = 32
                ELSE IF (k == 11)    THEN
                    number_of_ions = 23
                ELSE IF (k == 12)    THEN
                    number_of_ions = 8
                END IF
            
                number_of_ions = number_of_ions * 10
                
!                number_of_ions = 200
                            
            ELSE IF (simulation_option == 6)    THEN
            
                ! Myoglobin ion count
                IF  (k == 1)    THEN
                    number_of_ions = 28
                ELSE IF (k == 2)    THEN
                    number_of_ions = 36
                ELSE IF (k == 3)    THEN
                    number_of_ions = 47
                ELSE IF (k == 4)    THEN
                    number_of_ions = 60
                ELSE IF (k == 5)    THEN
                    number_of_ions = 68
                ELSE IF (k == 6)    THEN
                    number_of_ions = 87
                ELSE IF (k == 7)    THEN
                    number_of_ions = 90
                ELSE IF (k == 8)    THEN
                    number_of_ions = 100
                ELSE IF (k == 9)    THEN
                    number_of_ions = 99
                ELSE IF (k == 10)    THEN
                    number_of_ions = 92
                ELSE IF (k == 11)    THEN
                    number_of_ions = 70
                ELSE IF (k == 12)    THEN
                    number_of_ions = 44
                ELSE IF (k == 13)    THEN
                    number_of_ions = 33
                ELSE IF (k == 14)    THEN
                    number_of_ions = 21
                ELSE IF (k == 15)    THEN
                    number_of_ions = 11
                ELSE IF (k == 16)    THEN
                    number_of_ions = 8
                END IF
                
                number_of_ions = number_of_ions * 10
                
!                number_of_ions = 200
                            
            ELSE IF (simulation_option == 7)    THEN
            
                ! Alcohol Dehydrogenase ion count
                IF  (k == 1)    THEN
                    number_of_ions = 29
                ELSE IF (k == 2)    THEN
                    number_of_ions = 40
                ELSE IF (k == 3)    THEN
                    number_of_ions = 56
                ELSE IF (k == 4)    THEN
                    number_of_ions = 70
                ELSE IF (k == 5)    THEN
                    number_of_ions = 87
                ELSE IF (k == 6)    THEN
                    number_of_ions = 100
                ELSE IF (k == 7)    THEN
                    number_of_ions = 99
                ELSE IF (k == 8)    THEN
                    number_of_ions = 97
                ELSE IF (k == 9)    THEN
                    number_of_ions = 95
                ELSE IF (k == 10)    THEN
                    number_of_ions = 93
                ELSE IF (k == 11)    THEN
                    number_of_ions = 90
                ELSE IF (k == 12)    THEN
                    number_of_ions = 81
                ELSE IF (k == 13)    THEN
                    number_of_ions = 73
                ELSE IF (k == 14)    THEN
                    number_of_ions = 63
                ELSE IF (k == 15)    THEN
                    number_of_ions = 55
                ELSE IF (k == 16)    THEN
                    number_of_ions = 48
                ELSE IF (k == 17)    THEN
                    number_of_ions = 45
                END IF
                
                number_of_ions = number_of_ions * 10
                
!                number_of_ions = 200
                        
            ELSE IF (simulation_option == 8)    THEN
            
                ! Alcohol Dehydrogenase ion count
                IF  (k == 1)    THEN
                    number_of_ions = 50
                ELSE IF (k == 2)    THEN
                    number_of_ions = 86
                ELSE IF (k == 3)    THEN
                    number_of_ions = 100
                ELSE IF (k == 4)    THEN
                    number_of_ions = 96
                ELSE IF (k == 5)    THEN
                    number_of_ions = 92
                ELSE IF (k == 6)    THEN
                    number_of_ions = 77
                ELSE IF (k == 7)    THEN
                    number_of_ions = 23
                END IF
                
                number_of_ions = number_of_ions * 10
                
!                number_of_ions = 200
                       
            ELSE IF (simulation_option == 9)    THEN
                PRINT *, "Enter number of ions to simulate..."
                READ *, number_of_ions        
                WRITE (97,*) "Enter number of ions to simulate..."
                WRITE (97,*) number_of_ions 
!                number_of_ions = 10000

                ! Determine m/z and m of each ion based on ranges set by user
                PRINT *, "Enter lower m/z of interest..."
                READ *, low_mz
                PRINT *, "Enter upper m/z of interest..."
                READ *, high_mz
                PRINT *, "Enter lower mass of interest..."
                READ *, low_m
                PRINT *, "Enter upper mass of interest..."
                READ *, high_m
                WRITE (97,*) "Enter lower m/z of interest..."
                WRITE (97,*) low_mz
                WRITE (97,*) "Enter upper m/z of interest..."
                WRITE (97,*) high_mz
                WRITE (97,*) "Enter lower mass of interest..."
                WRITE (97,*) low_m
                WRITE (97,*) "Enter upper mass of interest..."
                WRITE (97,*) high_m
                
                CLOSE (UNIT=97, STATUS='KEEP')
                
                ALLOCATE (rnd_m(1:number_of_ions), rnd_mz(1:number_of_ions))
                CALL RANDOM_SEED()
                CALL RANDOM_NUMBER(rnd_m)
                CALL RANDOM_NUMBER(rnd_mz) 
                
                ALLOCATE (m_array(1:number_of_ions), mz_array(1:number_of_ions), charge_in(1:number_of_ions), freq_in(1:number_of_ions))
                
                unit = 434
                filename = "input_ions_Sep15.txt"
                CALL open_new_file (unit, filename) 
                
                WRITE (434, 435) "Input_Freq", "Input_Charge", "Input_m/z", "Input_Mass"
                435 FORMAT (TR5, A10, TR4, A12, TR2, A9, TR4, A10)
                
                DO i = 1, number_of_ions
                    m_array(i) = low_m + rnd_m(i)*(high_m - low_m)
                    mz_array(i) = low_mz + rnd_mz(i)*(high_mz - low_mz)
                    
                    charge_in(i) = m_array(i) / mz_array(i)
!                    freq_in(i) = SQRT(2540452974954.05 / mz_array(i)) 
                    freq_in(i) = SQRT(2510080000000 / mz_array(i))
                    
                    WRITE (434, 436) freq_in(i), charge_in(i), mz_array(i), m_array(i)
                    436 FORMAT (4F14.2) 
                END DO    
            END IF
            
            
            files_analyzed = 0     
            
!!!            !$OMP PARALLEL NUM_THREADS(8) DEFAULT(private) COPYIN(simulation_option, number_of_files, data_type) &
!!!            !$OMP COPYIN(test_signal_option, noise_channel, channel, fmt, input_freq, electrons, files_analyzed) &
!!!            !$OMP SHARED(charge_avg_array, number_of_ions, electrons, number_of_files, simulation_option, files_analyzed, data_type) &
!!!            !$OMP SHARED(input_freq, test_signal_option, noise_channel, channel, fmt, length_multiplier, max_trap_time)
!!!            
!!!            !$OMP DO SCHEDULE(guided) PRIVATE(sq_wave_array, array_pointer, square_total_real, old_square_total_real) &
!!!            !$OMP PRIVATE(mult, remainder, random_array, rand_1, rand_2, rand_3, rand_4, rand_5, charge, m, unit, j)
                  
            DO  m = 1, number_of_ions
                
                IF  (simulation_option == 1)    THEN
                    input_freq = 29000 + (rand_1 * 36000)
                ELSE IF (simulation_option == 9)    THEN 
!                    mass = (rand_6 * 9.E5) + 1.E5
!                    mass_to_charge = (rand_7 * 7500.) + 7500.
!                    electrons = mass / mass_to_charge
!!                    input_freq = SQRT(2540452974954.05 / mass_to_charge)
!                    input_freq = SQRT(2510080000000 / mass_to_charge)
                    electrons = charge_in(m)
                    input_freq = freq_in(m)
                END IF
                
                CALL RANDOM_SEED()
                CALL RANDOM_NUMBER(random_array)
                rand_1 = random_array(1)
                rand_2 = random_array(2)
                rand_3 = random_array(3)
                rand_4 = random_array(4)
                rand_5 = random_array(5)
                rand_6 = random_array(6) 
                rand_7 = random_array(7) 
                rand_8 = random_array(8)         
        
                files_analyzed = files_analyzed + 1
                PRINT *, "Analyzing Charge State: ", electrons, files_analyzed, "/", number_of_ions  
                
                
                max_cycles = (2500 * max_trap_time * input_freq * 16) / 40E6
!                signal_end = INT(((200 + (rand_2 * (max_cycles - 200))) * 40.E6) / (input_freq * 16.))
                signal_end = 2500 * max_trap_time
                
                    
                ! Actions to perform if test signal is created from SIMION file
                IF (test_signal_option == 1)    THEN                         
                
                    !   Open input file
                    unit = 801
                    filename = "potential_table(y_z_E).txt"
                    
                    PRINT *, "Now reading from ", filename    
                    OPEN (UNIT=unit, FILE=filename, STATUS="old", ACTION="read", IOSTAT=stat)   
            
                    DO  i = 1, 15311
                        READ (801, *) y_value, z_value, potential
                        potential_table(1,i) = y_value
                        potential_table(2,i) = z_value
                        potential_table(3,i) = potential/10000
                    END DO        
                    
                    CLOSE(801)                    
                                              
                    !   Open input file
                    unit = 802
                    filename = "test_signal.txt"
                    
                    PRINT *, "Now reading from ", filename    
                    OPEN (UNIT=unit, FILE=filename, STATUS="old", ACTION="read", IOSTAT=stat)   
                            
                    !   Read input file
                    READ (802, *) scratch_char  ! Read header
                    
                    ! Find number of entries from file "test_signal.txt"
                    n = 0
                    DO
                        READ (802, *, IOSTAT=stat) tof, y_value, z_value
                        IF (stat /= 0) EXIT                   ! Stop reading values when end of file is reached
                        n = n + 1
                    END DO        
                    
                    CLOSE(802)
                                                
                    PRINT *, "Now reading from ", filename    
                    OPEN (UNIT=unit, FILE=filename, STATUS="old", ACTION="read", IOSTAT=stat)           
                               
                    ALLOCATE(sq_wave_array_allocatable(1:3, 1:n), output(1:n))
                    
                    !   Read input file
                    READ (802, *) scratch_char  ! Read header
                    
                    DO  i = 1, n
                        READ (802, *, IOSTAT=stat) tof, y_value, z_value
                        sq_wave_array_allocatable(1,i) = ABS(y_value/0.127)
                        sq_wave_array_allocatable(2,i) = ABS(z_value/0.127)
!                        sq_wave_array_allocatable(1,i) = ABS(y_value/0.0635)
!                        sq_wave_array_allocatable(2,i) = ABS(z_value/0.0635)
                        sq_wave_array_allocatable(3,i) = tof
                    END DO  
                    
                    CLOSE(802)
                                
                    ! index = (z_value * 61) + (y_value + 1)
                    DO  i = 1, n   
                        index_z1y1 = (INT(sq_wave_array_allocatable(2,i)) * 61) + (INT(sq_wave_array_allocatable(1,i)) + 1)
                        index_z2y1 = ((INT(sq_wave_array_allocatable(2,i)) + 1) * 61) + (INT(sq_wave_array_allocatable(1,i)) + 1)
                        index_z1y2 = (INT(sq_wave_array_allocatable(2,i)) * 61) + (INT(sq_wave_array_allocatable(1,i) + 1 + 1))
                        index_z2y2 = ((INT(sq_wave_array_allocatable(2,i)) + 1) * 61) + (INT(sq_wave_array_allocatable(1,i) + 1 + 1))
                        
                        R1 = (((INT(sq_wave_array_allocatable(2,i)) + 1) - sq_wave_array_allocatable(2,i)) / ((INT(sq_wave_array_allocatable(2,i)) + 1) - INT(sq_wave_array_allocatable(2,i))) * potential_table(3,index_z1y1) ) + &
                             & ( (sq_wave_array_allocatable(2,i) - INT(sq_wave_array_allocatable(2,i))) / ((INT(sq_wave_array_allocatable(2,i)) + 1) - INT(sq_wave_array_allocatable(2,i))) * potential_table(3,index_z2y1) )
                        
                        R2 = (((INT(sq_wave_array_allocatable(2,i)) + 1) - sq_wave_array_allocatable(2,i)) / ((INT(sq_wave_array_allocatable(2,i)) + 1) - INT(sq_wave_array_allocatable(2,i))) * potential_table(3,index_z1y2) ) + &
                             & ( (sq_wave_array_allocatable(2,i) - INT(sq_wave_array_allocatable(2,i))) / ((INT(sq_wave_array_allocatable(2,i)) + 1) - INT(sq_wave_array_allocatable(2,i))) * potential_table(3,index_z2y2) )
                        
                        output(i) = (((INT(sq_wave_array_allocatable(1,i)) + 1) - sq_wave_array_allocatable(1,i)) / ((INT(sq_wave_array_allocatable(1,i)) + 1) - INT(sq_wave_array_allocatable(1,i))) * R1 ) + &
                             & ( (sq_wave_array_allocatable(1,i) - INT(sq_wave_array_allocatable(1,i))) / ((INT(sq_wave_array_allocatable(1,i)) + 1) - INT(sq_wave_array_allocatable(1,i))) * R2 )
                    END DO
                        
                
                    DO i = 1, signal_end
!                        sq_wave_array(i) = output(((i-1)*100)+1) * NINT(2 * electrons)
                        sq_wave_array(i) = output(((i-1)*100)+1) * NINT(2.31 * electrons)
                    END DO
                    
                    ! The rest of the array is filled with zeroes
                    DO  i = (signal_end + 1), 1048576
                        sq_wave_array (i) = 0
                    END DO        
                    
                    DEALLOCATE(sq_wave_array_allocatable, output)
                                
                                
                ! Actions to perform if test signal is based on generic equation
                ELSE IF (test_signal_option == 2)   THEN
                
                    ! Variable names and values obtained from Origin for creating test signal
                    ! Test signal approximates a real signal and follows a Boltzmann distribution
                    ! Details of test signal and fitting function obtained from Simion and Origin, respectively
                    A1 = -3.35668061258415E-4   ! NCC TRAP 44 - CONETRAP
                    A2 = 0.994194949645035      ! NCC TRAP 44 - CONETRAP
                    x0 = -0.505420267974969     ! NCC TRAP 44 - CONETRAP
                    dx = 0.040649051968342      ! NCC TRAP 44 - CONETRAP
                    
                    square_total_real = 40.e6 / (input_freq * 16.)
                    old_square_total_real = square_total_real
                    mult = INT((2500*max_trap_time) / square_total_real)
                    remainder = 0      
                    array_pointer = 0
                    sq_wave_array = 0
                    
                    ! Create test signal based on an 29.15% duty cycle (based on Simion simulations for NCC TRAP 44 - CONETRAP)
                    ! Also must take into account the offset of each cycle based on limited time resolution
                    ! (i.e., not every frequency will have an integer number of datapoints per cycle due to a finite time
                    ! resolution in data collection)
                    ! This offset changes the center of the Boltzmann distribution slightly so that not every cycle looks exactly the same       
                    DO  i = 1, mult-1         
                        DO  j = 1, INT(square_total_real)
                            array_pointer = array_pointer + 1
                            
                            ! NCC TRAP 44 - CONETRAP
                            x = -1.715265866 + ((j-1) * 3.430531732 / square_total_real) - (remainder / square_total_real * 100 * (3.430531732 / square_total_real))
                            
                            IF (x <= 0)    THEN
                                sq_wave_array (array_pointer) = (A2 + (A1-A2)/(1 + EXP((x-x0)/dx))) * NINT(2.31 * electrons)
                            ELSE IF (x > 0) THEN
                                sq_wave_array (array_pointer) = (A2 + (A1-A2)/(1 + EXP((-x-x0)/dx))) * NINT(2.31 * electrons)
                            END IF              
                        END DO
                        
                        ! Non-integer portion of each cycle length is added onto frequency shift to create a steadily changing frequency
                        remainder = square_total_real - INT(square_total_real)
                        square_total_real = (old_square_total_real * (1 - ((freq_shift/100)*i/((mult-1)*signal_end/(2500*length_multiplier))))) + remainder
                        
                    END DO                         
                    
                END IF
                
            
                IF  (noise_channel < 3) THEN
                    
                    n_main = CEILING(rand_4 * expt_num)
                    IF  (n_main == 0)   THEN
                        n_main = 1
                    ELSE IF (n_main == (expt_num + 1))  THEN
                        n_main = expt_num
                    END IF
                          
                    ! Randomly select a noise file
                    j = FLOOR(rand_3 * number_of_files) 
                    IF  (j ==  number_of_files) THEN
                        j = number_of_files
                    END IF
                            
                
                    ! Open input noise file
                    unit = (j + 10000)
                            
                    IF  (data_type == 1)    THEN
                    
                        IF  (noise_channel == 1)    THEN                
                            WRITE (file_number, fmt) j+10000
                            filename = 'chA'//TRIM(file_number)//'.dat' 
                        ELSE IF (noise_channel == 2)    THEN                
                            WRITE (file_number, fmt) j+10000
                            filename = 'chB'//TRIM(file_number)//'.dat'            
                        END IF          
                    
                        PRINT *, "Now reading from ", filename    
                        OPEN (UNIT=unit, FILE=filename, STATUS="old", ACTION="read", IOSTAT=stat, FORM="binary")             
                    
                        input_array = 0
                        y_input = 0
                        
                        ! Read input file into array "y_input"
                        READ (unit, END=11) input_array            ! Read filename values and store in array "y"
                        11   CONTINUE 
                    
                        DO  n = 1, 1000000
                            y_input(n) = input_array(n)
                        END DO
                    
                    ELSE IF (data_type == 2)   THEN
                           
                        IF  (noise_channel == 1)    THEN                
                            WRITE (file_number, fmt) j+10000
                            filename = 'chA'//TRIM(file_number)//'.txt' 
                        ELSE IF (noise_channel == 2)    THEN                
                            WRITE (file_number, fmt) j+10000
                            filename = 'chB'//TRIM(file_number)//'.txt'            
                        END IF    
                     
                        PRINT *, "Now reading from ", filename    
                        OPEN (UNIT=unit, FILE=filename, STATUS="old", ACTION="read", IOSTAT=stat)           
                              
                        y_input = 0
                        
                        ! Read input file into array "y_input"
                        READ (unit, *, END=82) y_input            ! Read filename values and store in array "y"
                        82   CONTINUE 
                     
                    END IF
                    
                    CLOSE (unit) 
                 
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
                                GOTO 17
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
                                GOTO 17
                            END IF
                            
                        END IF
                    END DO                        
                      
                            
                    IF  (((derivative_peaks(2,(n_main+1)) - derivative_peaks(1,n_main)) > ((2500 * REAL(max_trap_time + (load_time / 2.))) - 50)) &
                        & .AND. ((derivative_peaks(2,(n_main+1)) - derivative_peaks(1,n_main)) < ((2500 * REAL(max_trap_time + (load_time / 2.))) + 50))) THEN
                        
                        derivative_peak = derivative_peaks(1,n_main) + (2500 * load_time / 2.)
                        input = 0
                        
                        DO  i = 0, ((2500 * max_trap_time) - 1)
                            input(i+1) = y_input(derivative_peak + i)
                        END DO
                                                
                        CALL fft_filter(input, fft_filter_output)
                        input = fft_filter_output
                        
                    ELSE IF (((derivative_peaks(2,(n_main)) - derivative_peaks(1,n_main)) > ((2500 * REAL(max_trap_time + (load_time / 2.))) - 50)) &
                        & .AND. ((derivative_peaks(2,(n_main)) - derivative_peaks(1,n_main)) < ((2500 * REAL(max_trap_time + (load_time / 2.))) + 50))) THEN
                        
                        derivative_peak = derivative_peaks(1,n_main) + (2500 * load_time / 2.)
                        input = 0
                        
                        DO  i = 0, ((2500 * max_trap_time) - 1)
                            input(i+1) = y_input(derivative_peak + i)
                        END DO
                           
                        CALL fft_filter(input, fft_filter_output)
                        input = fft_filter_output
                                         
                    END IF 
                                                   
                ELSE
                    n_main = 1
                    y_input = 0
                    input = 0
                END IF
                
                
!                WRITE (22, 244) filename, "True_Signal_End:", signal_end
!                244 FORMAT (A12, TR4, A16, TR5, I6)                  
                
                
                section = n_main
                n_main = 1 
                           
                ! Add test signal onto noise file for troubleshooting/optimization purposes    
                DO  i = 1, signal_end
                    input(i) = input(i) + sq_wave_array (NINT(rand_5 * old_square_total_real) + i)
                END DO
                
                ! Analyze test signal (added to noise file) to determine accuracy and limits of program 
                CALL peakfinder_VI(filename, input, n_main, channel, charge, length_multiplier, max_trap_time, dead_time, single_event, multiple_event, low_mass_cutoff, high_mass_cutoff, sig_cutoff, section)
                
                IF  (simulation_option == 1)    THEN
                    charge_avg_array(m) = charge
                END IF
            
                17 CONTINUE 
                                   
            END DO
            
!!!            !$OMP END DO
!!!            !$OMP BARRIER
!!!            !$OMP END PARALLEL   
                       
            
            
            IF  (simulation_option == 1)    THEN
            
                n = 0
                charge_sum = 0
                DO  i = 1, number_of_ions            
                    IF  (charge_avg_array(i) /= 0)  THEN
                        n = n + 1
                        charge_sum = charge_sum + charge_avg_array(i)
                    END IF                
                END DO
                
                charge_avg = charge_sum/n
                
                ! Find the sum of the square of the difference between average and individual value
                charge_diff_squared = 0  
                charge_exp_diff_squared = 0        
                DO  i = 1, number_of_ions
                    IF  (charge_avg_array(i) /= 0)  THEN
                        charge_diff_squared = charge_diff_squared + ((charge_avg_array(i) - charge_avg)**2) 
                        charge_exp_diff_squared = charge_exp_diff_squared + ((charge_avg_array(i) - REAL(electrons))**2)
                    END IF  
                END DO
                                
                ! Find standard deviation for the following variables
                charge_std_dev = SQRT(charge_diff_squared/REAL(n-1))
                charge_SD_exp_mean = SQRT(charge_exp_diff_squared/REAL(n-1))    
                
                WRITE (227, 3100) electrons, charge_avg, charge_std_dev, charge_SD_exp_mean, n
                3100 FORMAT (TR5, F5.1, TR13, F8.3, TR10, F11.6, TR19, F11.6, TR11, I3) 
                
                DEALLOCATE (charge_avg_array)
                
            END IF         
            
        END DO
        
!        ! Close all output files and keep them (i.e., don't delete them)
        CLOSE (22, STATUS = "KEEP")  
        
        IF  (simulation_option == 1)    THEN   
            CLOSE (227, STATUS = "KEEP") 
        END IF
                 
    END SUBROUTINE simulations    