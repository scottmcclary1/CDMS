    PROGRAM FFT
    USE mkl_dfti
    USE trap_mod
    USE omp_lib
    IMPLICIT NONE
     
   !REAL(4), PARAMETER :: PI = 3.14159265                            
   !REAL(4), PARAMETER :: NPTS = 1048576         !Number of points for FFT
   !REAL(4), PARAMETER :: NPTS10 = 1000000       !Max number of data points from input file
       
    REAL(8) :: start_run_time          !Time when program starts
    REAL(8) :: end_run_time            !Time when program ends
    REAL(8) :: total_time              !Total time program took to run
                  
    INTEGER :: i, j, k, m, n           !Counting variables for DO loops
    INTEGER :: status                  !I/O variable for status of file
    INTEGER :: unit                    !Dummy variable for opening files
    INTEGER :: stat                    !Error message when opening files
    INTEGER :: n_main                  !Pointer for file section
    INTEGER :: number_of_sections      !Number of trapping events per file assuming triggered trapping
    INTEGER :: channel                 !Pointer for channel to analyze (A, B, or test signal)
    INTEGER :: data_type               !Pointer for .dat or .txt files to analyze
    INTEGER :: noise_channel           !Pointer for noise channel
    INTEGER :: number_of_files         !Number of files to analyze
    INTEGER :: number_of_events        !Number of trapping events in files being analyzed
    INTEGER :: signal_end              !Length of test signal to add onto noise file 
    INTEGER :: test_signal_option      !Indicator for SIMION test signal or generic test signal  
    INTEGER :: number_of_ions          !Number of ions to analyze for each charge state in mock BSA signal
    INTEGER :: files_analyzed          !Counter for number of files analyzed
    INTEGER :: derivative_peak         !Inflection point of derivative
    INTEGER :: der_opt                 !Derivative option for looking for positive or negative derivative 
    INTEGER :: length_multiplier       !Multiplier for length of each trapping event (multiple of 1 ms or
                                       !   2500 datapoints) This is period of trapping event in ms,
                                       !   including the loading time for ions.
    INTEGER :: max_trap_time           !Length of each trapping event (multiple of 1 ms or 2500 datapoints)
                                       !   1 ms is exactly equal to 2500 data points because the boards
                                       !   collect at 2.5 MHz.
    INTEGER :: run                     !Run number used for naming data and peaks files
    CHARACTER(len=2) :: run_number     !Converting run to character format
    INTEGER :: load_time               !Time between end of one trapping event and beginning of next
                                       !   (multiple of 1 ms or 2500 datapoints)
    INTEGER :: expt_num                !Maximum number of experiments per file based on signal length
                                       !   (length_multiplier), i.e. # of trapping events.
    INTEGER :: dead_time               !Number of points of dead time due to A250 recovering from impulse
                                       !   noise of trap gate pulse
        
    INTEGER :: final_signal_end        !Pointer for end of signal
    REAL(4) :: end_time                !Trapping time
    REAL(4) :: cycles                  !Number of cycles ion was trapped (windowed FFT)
    INTEGER :: points_to_average       !Number of points averaged to find mass_to_charge, charge, and mass
    REAL(4) :: freq_avg                !Frequency average
    REAL(4) :: freq_std_dev            !Standard deviation of freq_avg
    REAL(4) :: mass_to_charge          !Average m/z value (windowed FFT)
    REAL(4) :: mass_to_charge_std_dev  !Standard deviation of m/z (windowed FFT)
    REAL(4) :: charge                  !Average charge value (windowed FFT)
    REAL(4) :: charge_std_dev          !Standard deviation of charge (windowed FFT)
    REAL(4) :: mass                    !Average mass value (windowed FFT)
    REAL(4) :: mass_std_dev            !Standard deviation of mass (windowed FFT)
    
    
    CHARACTER(len=10), ALLOCATABLE :: filename_array (:)     !Array for filename
    INTEGER, ALLOCATABLE :: n_main_array (:)                 !Array for file section
    INTEGER, ALLOCATABLE :: final_signal_end_array (:)       !Pointer for end of signal
    REAL(4), ALLOCATABLE :: end_time_array (:)               !Trapping time
    REAL(4), ALLOCATABLE :: cycles_array (:)                 !Number of cycles ion was trapped (windowed FFT)
    INTEGER, ALLOCATABLE :: points_to_average_array (:)      !Number of points averaged to find mass_to_charge, charge, and mass
    REAL(4), ALLOCATABLE :: freq_avg_array (:)               !Frequency average
    REAL(4), ALLOCATABLE :: freq_std_dev_array  (:)          !Standard deviation of freq_avg
    REAL(4), ALLOCATABLE :: mass_to_charge_array (:)         !Average m/z value (windowed FFT)
    REAL(4), ALLOCATABLE :: mass_to_charge_std_dev_array (:) !Standard deviation of m/z (windowed FFT)
    REAL(4), ALLOCATABLE :: charge_array (:)                 !Average charge value (windowed FFT)
    REAL(4), ALLOCATABLE :: charge_std_dev_array (:)         !Standard deviation of charge (windowed FFT)
    REAL(4), ALLOCATABLE :: mass_array (:)                   !Average mass value (windowed FFT)
    REAL(4), ALLOCATABLE :: mass_std_dev_array (:)           !Standard deviation of mass (windowed FFT)
                ! Why are there arrays which seem to duplicate the previous block of definitions?    
    
    INTEGER(2) :: input_array(NPTS10)       !Data read into here from raw binary files.
                                            !   1000000 values because files are 400 ms, and
                                            !   collection rate is 2.5 MHz. 2.5 MHz * 0.4 s = 1000000.
    REAL(4)    :: y_input(NPTS)             !Input datapoints from trap data
                                            !   This goes to 2^20. It must be in this format for the FFT, right?    
    REAL(4)    :: derivative(NPTS)          !First derivative of y_input
    REAL(4)    :: y_filter(NPTS)            !Trap data after passing through 10 kHz high-pass filter               
                                            !   So this filters low frequencies out of the raw data?    
    REAL(4)    :: fft_filter_output(NPTS)   !Trap data after passing through 10 kHz high-pass FFT filter 
                                            !   And this does an FFT and then chops off low frequencies?    
    REAL(4)    :: fft_magnitude(NPTS/2)     !FFT magnitude output
    
    REAL(4) :: ytmp              !Used to correct order of points in y_input (initially, every pair of points switched).  
    REAL(8) :: y_value           !y-values for test signal
    REAL(8) :: z_value           !z-values for test signal
    REAL(8) :: tof               !TOF for test signal
    REAL(8) :: potential         !Potential values for test signal
    REAL(8) :: R1                !Bi-linear interpolation variable
    REAL(8) :: R2                !Bi-linear interpolation variable
               ! For which subroutine are these variables? And what are these variables?    
    
    INTEGER :: index_z1y1        !Index for interpolation of potential values for test signal
    INTEGER :: index_z2y1        !Index for interpolation of potential values for test signal
    INTEGER :: index_z1y2        !Index for interpolation of potential values for test signal
    INTEGER :: index_z2y2        !Index for interpolation of potential values for test signal
    
    REAL(8), DIMENSION(1:3, 1:15311) :: potential_table     !Array for y-values, z-values, and potential values
    
    REAL(4), ALLOCATABLE :: sq_wave_array_allocatable(:,:)  !Allocatable test signal array (not actually a square wave)
    REAL(4), ALLOCATABLE :: output(:)                       !Output for optimization purposes   
    
    
    REAL(4) :: max_cycles                !Maximum number of cycles based on frequency
    REAL(4) :: input_freq                !Frequency of test signal
    INTEGER :: input_freq_int            !Converting input_freq to integer format
    CHARACTER(len=7) :: input_freq_char  !Converting input_freq_int to character format
    REAL(4) :: input_freq_2              !Second frequency of test signal    
    REAL(4) :: freq_shift                !Frequency shift of test signal (in percent)
    REAL(4) :: electrons                 !Number of charges on test signal (in electrons)
    INTEGER :: electrons_int             !Converting electrons to integer format
    CHARACTER(len=7) :: electrons_char   !Converting electrons_int to character format
    REAL(4) :: rand_1                    !Random number for random starting point of test signal
    REAL(4) :: rand_2                    !Random number for random starting point of test signal
    REAL(4) :: rand_3                    !Random number for random starting point of test signal
    REAL(4) :: rand_4                    !Random number for random starting point of test signal
    REAL(4) :: rand_5                    !Random number for random starting point of test signal
    REAL(8) :: A1                        !Variable to define test signal (obtained from Origin)
    REAL(8) :: A2                        !Variable to define test signal (obtained from Origin)
    REAL(8) :: x0                        !Variable to define test signal (obtained from Origin)
    REAL(8) :: dx                        !Variable to define test signal (obtained from Origin)
    REAL(8) :: x                         !Variable to define test signal (obtained from Origin)
    
    REAL(8) :: square_total_real         !Length of test signal
    REAL(8) :: old_square_total_real     !Length of previous test signal before frequency shift
    REAL(8) :: remainder                 !Offset for next cycle in test signal
    
    INTEGER :: array_pointer             !Counter for creating test signal with frequency shift
    INTEGER :: mult                      !Number of cycles a section of datafile (56250/square_total_real)
    INTEGER :: start                     !Random start of cycle for test signal as determined by random number generator
     
    REAL(4) :: sq_wave_array(NPTS)       !Test signal array (not actually a square wave)
    REAL(4) :: random_array(5)           !Array of 5 random numbers
      
    CHARACTER(len=100) :: filename       !Dummy variable for opening files
    CHARACTER(len=100) :: scratch_char   !Dummy variable to dicard header lines
    CHARACTER(len=8)   :: fmt            !Format descriptor for file_number
    CHARACTER(len=10)  :: file_number    !File number to read
    
    
    INTEGER :: simulation_option
    INTEGER :: number_of_charge_states
    REAL(4) :: charge_sum
    REAL(4) :: charge_avg
    REAL(4) :: charge_diff_squared
    REAL(4) :: charge_exp_diff_squared
    REAL(4) :: charge_SD_exp_mean
    
    REAL(4),ALLOCATABLE :: charge_avg_array(:)     !Array for charge values to determine average
    
    INTEGER :: single_event                 !Counter for number of trapping events with a single ion trapped    
    INTEGER :: multiple_event               !Counter for number of trapping events with multiple ions trapped
    CHARACTER(len=12),ALLOCATABLE :: single_ion_filename(:)   !filenames of types of trapping events
    CHARACTER(len=12),ALLOCATABLE :: multiple_ion_filename(:) !filenames of types of trapping events
    CHARACTER(len=12),ALLOCATABLE :: no_ion_filename(:)       !filenames of types of trapping events

    INTEGER,ALLOCATABLE :: single_ion_section(:)   !sections of types of trapping events
    INTEGER,ALLOCATABLE :: multiple_ion_section(:) !sections of types of trapping events
    INTEGER,ALLOCATABLE :: no_ion_section(:)       !sections of types of trapping events

    CHARACTER(len=2),ALLOCATABLE :: single_ion_section_char(:)   !Character versions of above arrays
    CHARACTER(len=2),ALLOCATABLE :: multiple_ion_section_char(:) !Character versions of above arrays
    CHARACTER(len=2),ALLOCATABLE :: no_ion_section_char(:)       !Character versions of above arrays

    INTEGER :: total_events = 0            !Counter for total number of potential trapping events   
    INTEGER :: total_single_events = 0     !Counter for number of trapping events with a single ion trapped    
    INTEGER :: total_multiple_events  = 0  !Counter for number of trapping events with multiple ions trapped
    INTEGER :: total_empty_events = 0      !Counter for number of trapping events with no ions
    
    
    LOGICAL :: od                         !Logical variable for inquiring whether a unit is opened or not
        
    REAL(4) :: low_mass_cutoff  = 0       !Low mass cutoff for frequency information
    REAL(4) :: high_mass_cutoff = 0       !High mass cutoff for frequency information
    INTEGER :: sig_cutoff = 0             !Signal length cutoff (in datapoints) for frequency information
    
    REAL(4) :: input(NPTS)                !Input datapoints from trap data; contains data from just one
                                          !   trapping event from a file for individual analysis.
    INTEGER :: derivative_peaks(2,400)    ! Array for location of derivative peaks from derivative of y_input
    INTEGER :: pos_peaks                  !Counter for positive peaks in derivative of y_input
    INTEGER :: neg_peaks                  !Counter for negative peaks in derivative of y_input
    
    INTEGER :: section

    INTEGER :: mythrd = 0            !DKB-debug
    CHARACTER(12), PARAMETER :: procname="TrapAn03"  !DKB-debug


!=========================================================================================

    PRINT *, 'Start'
    
!Open file to store user input
    unit = 97
    filename = 'user_input_Sep15_2014.txt'
    OPEN (UNIT=unit, FILE=filename, STATUS="replace", ACTION="write", IOSTAT=status)
    
!User selected channel option or test signal
    PRINT *, "Press 1 for Channel A..."
    PRINT *, "Press 2 for Channel B (NOT SUPPORTED)..."
    PRINT *, "Press 3 to add a test signal onto a noise file and analyze..."
    PRINT *, "Press 4 for simulated data options..."    
    PRINT *, "Press 5 to perform FFT on selected file and section..."
    PRINT *, "Press 6 to output freq data for all files..."
    PRINT *, "Press 7 to scan all files for frequency shift > 20 Hz..."
    PRINT *, "Press 8 to output freq data for selected mass range and time trapped..."
    PRINT *, "Press 9 for trap calibration options..."
    PRINT *, "Press 10 for noise analysis..."
    PRINT *, "Press 11 for voltage calibration..."
    PRINT *, "Press 12 to add a test signal onto multiple noise files and analyze..."
    READ *, channel
    WRITE (97,*) "Press 1 for Channel A..."
    WRITE (97,*) "Press 2 for Channel B (NOT SUPPORTED)..."
    WRITE (97,*) "Press 3 to add a test signal onto a noise file and analyze..."
    WRITE (97,*) "Press 4 for simulated data options..."    
    WRITE (97,*) "Press 5 to perform FFT on selected file and section..."
    WRITE (97,*) "Press 6 to output freq data for all files..."
    WRITE (97,*) "Press 7 to scan all files for frequency shift > 20 Hz..."
    WRITE (97,*) "Press 8 to output freq data for selected mass range and time trapped..."
    WRITE (97,*) "Press 9 for trap calibration options..."
    WRITE (97,*) "Press 10 for noise analysis..."
    WRITE (97,*) "Press 11 for voltage calibration..."
    WRITE (97,*) "Press 12 to add a test signal onto multiple noise files and analyze..."
    WRITE (97,*) channel
    
!Delete user input file for channel 11 since no user input is required
    IF (channel == 11) CLOSE (UNIT=97, STATUS='DELETE')
          
      
!Only ask for the following information for the following options.
    IF  (channel == 3 .OR. channel == 4 .OR. channel == 5 .OR. channel == 9 .OR. channel == 10 .OR. channel == 12)   THEN 
        
        PRINT *, "Are files .dat or .txt?"
        PRINT *, "Press 1 for .dat (binary)..."
        PRINT *, "Press 2 for .txt (text)..."
        READ *, data_type
        WRITE (97,*) "Are files .dat or .txt?"
        WRITE (97,*) "Press 1 for .dat (binary)..."
        WRITE (97,*) "Press 2 for .txt (text)..."
        WRITE (97,*) data_type
        
    END IF         
        
      
!Only ask for the following information if a non-voltage calibration option is selected above
    IF  (channel < 11 .OR. channel == 12)   THEN          
        
        ! User selected option for maximum potential trapping time (determined by experiment) 
        PRINT *, "How long is one period (in ms)?"
        READ *, length_multiplier
        WRITE (97,*) "How long is one period (in ms)?"
        WRITE (97,*), length_multiplier
        
        PRINT *, "How long is maximum trapping time (in ms)?"
        READ *, max_trap_time
        WRITE (97,*) "How long is maximum trapping time (in ms)?"
        WRITE (97,*) max_trap_time   
        
        load_time = length_multiplier - max_trap_time
                
    END IF       
    
    dead_time = 5000 
        ! dead_time allows for A250 recovery after impulse.       
        
!File names are different for channel 4 (simulated data options, see below) so as not to
!overwrite data analysis files.
    IF  (channel == 1 .OR. channel == 2 .OR. channel == 8 .OR. channel == 12)    THEN
                          
        IF (channel /= 12) THEN  
            PRINT *, "Enter run number..."
            READ *, run
            WRITE (97,*) "Enter run number..."
            WRITE (97,*) run
            WRITE (run_number, "(I2.2)") run
            
            ! Open file to output identity of single, multiple, and no ion trapping events
            IF (channel /= 8) THEN      ! For channel 8, open this in "Frequency_Data"
                unit = 20
                filename = 'event_types_Sep15_2014_run'//TRIM(run_number)//'.txt'
                CALL open_new_file (unit, filename)
            END IF
               
            ! Open file to store peak data for all signals found
            unit = 21
            filename = 'peaks_Sep15_2014_run'//TRIM(run_number)//'.txt'
            CALL open_new_file (unit, filename)
            
            ! Open file to store verbose data info for windowed data
            unit = 22
            filename = 'data_Sep15_2014_run'//TRIM(run_number)//'.txt'
            CALL open_new_file (unit, filename)
            
            ! Open file to store good ion data info for windowed data
            unit = 23
            filename = 'data_Sep15_2014_filtered_run'//TRIM(run_number)//'.txt'
            CALL open_new_file (unit, filename)
            
            ! Headers for verbose output files
            IF (channel /= 8) THEN
                WRITE (20,210) "Single_Ion_File","Section_#","Multiple_Ion_File","Section_#","No_Ion_File","Section_#"
                210 FORMAT (A15, TR3, A9, TR3, A17, TR5, A9, TR3, A17, TR5, A9, TR3)
            END IF
            
            WRITE (21,220) "File", "Section_#", "Signal_Length", "Max_Freq", "Magnitude"
            220 FORMAT (A4, TR10, A10, TR5, A13, T49, A9, T71, A10)
            
!            WRITE (22,230) "File", "Section_#", " Signal_End ", "End_Time_(s)", "Cycles", "#_Points_Avgd", "Freq_Avg", &
!                           & "Freq_Std_Dev", "m/z", "m/z_Std_Dev", "Charge", "Charge_Std_Dev", "Mass", "Mass_Std_Dev"
!            230 FORMAT (T3, A4, T15, A10, TR4, A13, TR6, A12, TR5, A6, TR3, A13, TR7, A8, TR6, &
!                       & A12, TR8, A3, TR9, A11, TR6, A6, TR3, A14, TR8, A4, TR8, A13, TR5)
!            WRITE (22,230) "File", "Section_#", "End_Time_(s)", "#_Points_Avgd", "Freq_Avg", &
!               & "m/z", "m/z_Std_Dev", "Charge", "Charge_Std_Dev", "Mass", &
!               & "Harm2_Rel_Mag_Avg", "Harm2_Rel_Mag_Std_Dev"
!            230 FORMAT (T3, A4, T15, A10, TR4, A12, TR5, A13, TR7, A8, TR6, &
!                   & A3, TR9, A11, TR6, A6, TR3, A14, TR8, A4, TR8, A22, TR2, A22, TR2)
!            WRITE (22,230) "File", "Section_#", "End_Time_(s)", "#_Points_Avgd", "Freq_Avg", "Freq_Slope", "Freq_Int", "Freq_Sum_Sq", &
!               & "m/z", "m/z_Std_Dev", "Charge", "Charge_Std_Dev", "Mass", &
!               & "Harm2_Rel_Mag_Avg", "Harm2_Rel_Mag_Std_Dev"
!            230 FORMAT (T3, A4, T15, A10, TR4, A12, TR5, A13, TR6, A8, TR7, A10, TR7, A8, TR4, A11, TR8, &
!                   & A3, TR9, A11, TR6, A6, TR3, A14, TR8, A4, TR8, A22, TR2, A22, TR2)
            WRITE (22,230) "File", "Section_#", "End_Time_(s)", "#_Points_Avgd", "Freq_Avg", "Freq_Slope", "Freq_Rel_Slope", "Freq_Sum_Sq", &
               & "m/z", "m/z_Std_Dev", "Charge", "Charge_Std_Dev", "Mass", &
               & "Harm2_Rel_Mag_Avg", "Harm2_Rel_Mag_Std_Dev"
            230 FORMAT (T3, A4, T15, A10, TR4, A12, TR5, A13, TR6, A8, TR7, A10, TR7, A14, TR4, A11, TR8, &
                   & A3, TR9, A11, TR6, A6, TR3, A14, TR8, A4, TR8, A22, TR2, A22, TR2)
                   
                   
            WRITE (23,280) "File", "Section_#", "End_Time_(s)", "#_Points_Avgd", "Freq_Avg", "Freq_Slope", "Freq_Rel_Slope", "Freq_Sum_Sq", &
               & "m/z", "m/z_Std_Dev", "Charge", "Charge_Std_Dev", "Mass", &
               & "Harm2_Rel_Mag_Avg", "Harm2_Rel_Mag_Std_Dev"
            280 FORMAT (T3, A4, T15, A10, TR4, A12, TR5, A13, TR6, A8, TR7, A10, TR7, A14, TR4, A11, TR8, &
                   & A3, TR9, A11, TR6, A6, TR3, A14, TR8, A4, TR8, A22, TR2, A22, TR2)
        END IF
               
!        ! Open file to store peak data for all signals found
!        unit = 31
!        filename = "test_output.txt"
!        CALL open_new_file (unit, filename)        
    END IF
    
    fmt = '(I5.5)'      !This is format descriptor for the file number. Format is cIw.m,
                        !   where c is # of repeats, I is integer, w is total width, and
                        !   m is minimum # of digits.
    
    
!=========================================================================================
!=========================================================================================
    
!Actions to take for analyzing channel "A" or channel "B"
    IF  (channel == 1 .OR. channel == 2 .OR. channel == 12) THEN
        
        ! User selected option for number of files in each run (maximum = 100)  ! That statement is obsolete.
        ! Maximum number of files is only limited by how the code is currently written
        ! Can be increased by adding extra code to associated "DO" loop below
        PRINT *, "Enter last file number..."
        READ *, number_of_files  
        WRITE (97,*) "Enter last file number..."
        WRITE (97,*) number_of_files
        number_of_files = number_of_files + 1 - 10000       ! + 1 - 10000 because files start at 10000.
        
        ! Allocate each of following arrays with maximum number of trapping events in run
        number_of_events = number_of_files*INT(400/length_multiplier)
        ALLOCATE(single_ion_filename(1:number_of_events), multiple_ion_filename(1:number_of_events),no_ion_filename(1:number_of_events))
        ALLOCATE(single_ion_section(1:number_of_events), multiple_ion_section(1:number_of_events),no_ion_section(1:number_of_events))
        ALLOCATE(single_ion_section_char(1:number_of_events), multiple_ion_section_char(1:number_of_events),no_ion_section_char(1:number_of_events))
        single_ion_filename = '0'; multiple_ion_filename = '0'; no_ion_filename = '0'
        single_ion_section = 0; multiple_ion_section = 0; no_ion_section = 0
        single_ion_section_char = '0'; multiple_ion_section_char = '0'; no_ion_section_char = '0'
        
        IF (channel == 1 .OR. channel == 2) CLOSE (UNIT=97, STATUS='KEEP')
       !write(6,"(a,' :  375 : thread=',i4)") procname,mythrd        !DKB-debug
        
        IF (channel == 12) THEN     ! See channel == 3 code for notes comparison and explanation of variables.
            ! User selected option for SIMION test signal or generica test signal
            PRINT *, "Do you want to use a SIMION file or create a generic test signal?"
           !PRINT *, "Enter 1 for SIMION file (entitled test_signal.txt) (NOT SUPPORTED!)..."
            PRINT *, "Enter 1 for SIMION file (entitled signal_data.txt)..."
            PRINT *, "Enter 2 for generic test signal..."
            READ *, test_signal_option    
            WRITE (97,*) "Do you want to use a SIMION file or create a generic test signal?"
           !WRITE (97,*) "Enter 1 for SIMION file (entitled test_signal.txt) (NOT SUPPORTED!)..."
            WRITE (97,*) "Enter 1 for SIMION file (entitled signal_data.txt)..."
            WRITE (97,*) "Enter 2 for generic test signal..."
            WRITE (97,*) test_signal_option
            
            IF (test_signal_option == 2) THEN
                PRINT *, "Enter desired frequency of test signal in Hz..."
                WRITE (97,*) "Enter desired frequency of test signal in Hz..."
            ELSE IF (test_signal_option == 1) THEN
                PRINT *, "Enter expected frequency of test signal in Hz..."
                WRITE (97,*) "Enter expected frequency of test signal in Hz..."
            END IF
            READ *, input_freq
            WRITE (97,*) input_freq
            input_freq_int = input_freq
            WRITE (input_freq_char, "(I6)") input_freq_int
            
            IF (test_signal_option == 2) THEN
                PRINT *, "By what percent would you like the frequency to shift?"
                PRINT *, "(0 = no frequency shift, 5 = 5%, etc.)"
                READ *, freq_shift
                WRITE (97,*) "By what percent would you like the frequency to shift?"
                WRITE (97,*) "(0 = no frequency shift, 5 = 5%, etc.)"
                WRITE (97,*) freq_shift
            END IF
            
            ! Apply signal to entire noise file so each section can be analyzed.
            signal_end = NPTS10
           !PRINT *, "Enter length of test signal in ms..."
           !READ *, signal_end
            
            PRINT *, "Enter desired input charge (in electrons) for test signal..."
            READ *, electrons
            WRITE (97,*) "Enter desired input charge (in electrons) for test signal..."
            WRITE (97,*) electrons
            electrons_int = electrons
            WRITE (electrons_char, "(I4)") electrons_int
            
            CLOSE (UNIT=97, STATUS='KEEP')
            
            ! Open file to output identity of single, multiple, and no ion trapping events
            unit = 20
            filename = 'test_signal_event_types_Sep15_2014'//TRIM(input_freq_char)//'Hz'//TRIM(electrons_char)//'e.txt'
            CALL open_new_file (unit, filename)
            
            ! Open file to store peak data for all signals found
            unit = 21
           !filename = 'test_signal_peaks_Sep15_2014'//TRIM(input_freq_char)//'Hz'//TRIM(electrons_char)//'e.txt'
            filename = 'test_signal_peaks_Feb03_2014_run'//TRIM(run_number)//'.txt'
            CALL open_new_file (unit, filename)
            
            ! Open file to store verbose data info for windowed data
            unit = 22
           !filename = 'test_signal_data_Sep15_2014'//TRIM(input_freq_char)//'Hz'//TRIM(electrons_char)//'e.txt'
            filename = 'test_signal_data_Feb03_2014_run'//TRIM(run_number)//'.txt'
            CALL open_new_file (unit, filename)
            
            ! Open file to store good data info for windowed data
            unit = 23
            filename = 'test_signal_data_Sep15_2014_filtered'//TRIM(input_freq_char)//'Hz'//TRIM(electrons_char)//'e.txt'
            CALL open_new_file (unit, filename)
            
            ! Headers for verbose output files
            WRITE (20,260) "Single_Ion_File","Section_#","Multiple_Ion_File","Section_#","No_Ion_File","Section_#"
            260 FORMAT (A15, TR3, A9, TR3, A17, TR5, A9, TR3, A17, TR5, A9, TR3)
            
            WRITE (21,240) "File", "Section_#", "Signal_Length", "Max_Freq", "Magnitude"
            240 FORMAT (A4, TR10, A10, TR5, A13, T49, A9, T71, A10)
            
!            WRITE (22,250) "File", "Section_#", " Signal_End ", "End_Time_(s)", "Cycles", "#_Points_Avgd", "Freq_Avg", &
!                           & "Freq_Std_Dev", "m/z", "m/z_Std_Dev", "Charge", "Charge_Std_Dev", "Mass", "Mass_Std_Dev"
!            250 FORMAT (T3, A4, T15, A10, TR4, A13, TR6, A12, TR5, A6, TR3, A13, TR7, A8, TR6, &
!                       & A12, TR8, A3, TR9, A11, TR6, A6, TR3, A14, TR8, A4, TR8, A13, TR5)
!            WRITE (22,250) "File", "Section_#", "End_Time_(s)", "#_Points_Avgd", "Freq_Avg", &
!               & "m/z", "m/z_Std_Dev", "Charge", "Charge_Std_Dev", "Mass", &
!               & "Harm2_Rel_Mag_Avg", "Harm2_Rel_Mag_Std_Dev"
!            250 FORMAT (T3, A4, T15, A10, TR4, A12, TR5, A13, TR7, A8, TR6, &
!                   & A3, TR9, A11, TR6, A6, TR3, A14, TR8, A4, TR8, A22, TR2, A22, TR2)
!            WRITE (22,250) "File", "Section_#", "End_Time_(s)", "#_Points_Avgd", "Freq_Avg", "Freq_Slope", "Freq_Int", "Freq_Sum_Sq", &
!               & "m/z", "m/z_Std_Dev", "Charge", "Charge_Std_Dev", "Mass", &
!               & "Harm2_Rel_Mag_Avg", "Harm2_Rel_Mag_Std_Dev"
!            250 FORMAT (T3, A4, T15, A10, TR4, A12, TR5, A13, TR6, A8, TR7, A10, TR7, A8, TR4, A11, TR8, &
!                   & A3, TR9, A11, TR6, A6, TR3, A14, TR8, A4, TR8, A22, TR2, A22, TR2)
           WRITE (22,250) "File", "Section_#", "End_Time_(s)", "#_Points_Avgd", "Freq_Avg", "Freq_Slope", "Freq_Rel_Slope", "Freq_Sum_Sq", &
               & "m/z", "m/z_Std_Dev", "Charge", "Charge_Std_Dev", "Mass", &
               & "Harm2_Rel_Mag_Avg", "Harm2_Rel_Mag_Std_Dev"
            250 FORMAT (T3, A4, T15, A10, TR4, A12, TR5, A13, TR6, A8, TR7, A10, TR7, A14, TR4, A11, TR8, &
                   & A3, TR9, A11, TR6, A6, TR3, A14, TR8, A4, TR8, A22, TR2, A22, TR2)
                   
           
           WRITE (23,290) "File", "Section_#", "End_Time_(s)", "#_Points_Avgd", "Freq_Avg", "Freq_Slope", "Freq_Rel_Slope", "Freq_Sum_Sq", &
               & "m/z", "m/z_Std_Dev", "Charge", "Charge_Std_Dev", "Mass", &
               & "Harm2_Rel_Mag_Avg", "Harm2_Rel_Mag_Std_Dev"
            290 FORMAT (T3, A4, T15, A10, TR4, A12, TR5, A13, TR6, A8, TR7, A10, TR7, A14, TR4, A11, TR8, &
                   & A3, TR9, A11, TR6, A6, TR3, A14, TR8, A4, TR8, A22, TR2, A22, TR2)
        END IF
        
        ! Call initial clock time (to measure program running time)        
        CALL CLOCKX(start_run_time)
        
        IF (channel == 12) THEN
            IF (test_signal_option ==2 ) THEN
                CALL RANDOM_SEED()
                CALL RANDOM_NUMBER(random_array)
                rand_1 = random_array(1)
                rand_2 = random_array(2)
                rand_3 = random_array(3)
                rand_4 = random_array(4)
                rand_5 = random_array(5)
                
                ! Generic test signal used
                A1 = -3.35668061258415E-4   ! NCC TRAP 44 - CONETRAP
                A2 = 0.994194949645035      ! NCC TRAP 44 - CONETRAP
                x0 = -0.505420267974969     ! NCC TRAP 44 - CONETRAP
                dx = 0.040649051968342      ! NCC TRAP 44 - CONETRAP
                    ! The values here work with the Boltzmann function in Origin (search for "Boltzmann" in Origin) to 
                    ! produce a sigmoid. A1 and A2 are in close agreement with equation 2.5 of Nathan Contino's thesis.
                    ! x0 and dx differ by the same factor from values in equation 2.5, producing a narrower sigmoid.
                
                square_total_real = 40.e6/(input_freq*16.)
                old_square_total_real = square_total_real
                mult = INT(signal_end/square_total_real)
               !mult = INT((2500*length_multiplier)/square_total_real)
                remainder = 0      
                array_pointer = 0
                
                ! Create test signal based on an 29.15% duty cycle (based off of Simion simulations for NCC TRAP 44 - CONETRAP)
                ! Also must take into account the offset of each cycle based on limited time resolution
                ! (i.e., not every frequency will have an integer number of datapoints per cycle due to a finite time
                ! resolution in data collection)
                ! This offset changes the center of the Boltzmann distribution slightly so that not every cycle looks exactly the same       
                expt_num = INT(400/length_multiplier)     ! This definition of expt_num is only valid for triggered data files
                DO k = 1, expt_num  ! Nested DO loops allow for frequency shifts to be written fresh for each trapping event
                    i = 0           ! i is counting the number of test signal cycles written for current trapping event
                    DO WHILE (array_pointer >= (k-1)*2500*length_multiplier .AND. array_pointer < k*2500*length_multiplier)
                        i = i + 1
                   !DO  i = 1, mult-1       ! Each iteration of this outer loop writes data for one cycle of test signal
                       !DO  i = 1, INT((mult-1)/2.)       
                        DO  j = 1, INT(square_total_real)
                            array_pointer = array_pointer + 1
                            ! NCC TRAP 44 - CONETRAP
                            x = -1.715265866 + ((j-1)*3.430531732/square_total_real) - (remainder/square_total_real * 100*(3.430531732/square_total_real))
                            IF (x <= 0)    THEN
                                sq_wave_array(array_pointer) = (A2 + (A1-A2)/(1+EXP((x-x0)/dx))) * NINT(2.31*electrons)
                            ELSE IF (x > 0) THEN
                                sq_wave_array(array_pointer) = (A2 + (A1-A2)/(1+EXP((-x-x0)/dx))) * NINT(2.31*electrons)
                            END IF  
                        END DO
                        
                        !Non-integer portion of each cycle length is added onto frequency
                        !shift to create a steadily changing frequency
                        remainder = square_total_real - INT(square_total_real)
                        square_total_real = (old_square_total_real * (1 - (expt_num*freq_shift/100)*i/(mult-1))) + remainder
                            !freq_shift is shift across NPTS10 data points (nominally 1000000). Multiply by
                            !   expt_num to cause that shift over each trapping event
                       !square_total_real = (old_square_total_real * (1 - (freq_shift/100)*i/(mult-1))) + remainder
                       !square_total_real = (old_square_total_real * (1 - ((freq_shift/100)*i/((mult-1)*signal_end/(2500*length_multiplier))))) + remainder 
                    END DO   
                    square_total_real = old_square_total_real
                END DO
            END IF
            
            IF (test_signal_option == 1) THEN
                ! Open input file
                unit = 801
                filename = "signal_data.txt"            
                PRINT *, "Now reading from ", TRIM(filename)    
                OPEN (UNIT=unit, FILE=filename, STATUS="old", ACTION="read", IOSTAT=stat)
                DO i = 1, signal_end
                    READ (801,*) sq_wave_array(i)
                END DO 
                sq_wave_array = sq_wave_array
                sq_wave_array = sq_wave_array*NINT(2.31*electrons)
            END IF
            
           !DO  i = mult*INT(square_total_real), NPTS
           !    sq_wave_array (i) = 0
           !END DO
            
        END IF

        files_analyzed = 0
       !write(6,"(a,' :  564 : thread=',i4)") procname,mythrd        !DKB-debug

        ! OpenMP (multithreading) lines begin with !$. Any other lines beginning with ! are comments.     
        !$OMP PARALLEL DEFAULT(none) &
        !$OMP SHARED(files_analyzed, number_of_files, fmt, channel, length_multiplier, max_trap_time, load_time, expt_num, dead_time) &
        !$OMP SHARED(total_events, total_single_events, total_multiple_events, total_empty_events, low_mass_cutoff, high_mass_cutoff, sig_cutoff) &
        !$OMP SHARED(old_square_total_real, signal_end, sq_wave_array) &
        !$OMP SHARED(single_ion_filename, single_ion_section, multiple_ion_filename, multiple_ion_section, no_ion_filename, no_ion_section) &
        !$OMP PRIVATE(mythrd)

        !DKB-debug: The following shared variables are allocatable arrays:
        !DKB-debug:    single_ion_filename,    single_ion_section,
        !DKB-debug:    multiple_ion_filename,  multiple_ion_section,
        !DKB-debug:    no_ion_filename,        no_ion_section

        !xxxSHARED(filename_array, n_main_array, final_signal_end_array, end_time_array, cycles_array, points_to_average_array) & 
        !xxxSHARED(freq_avg_array, freq_std_dev_array, mass_to_charge_array, mass_to_charge_std_dev_array) &
        !xxxSHARED(charge_array, charge_std_dev_array, mass_array, mass_std_dev_array, i)

!#ifdef _OPENMP                         !DKB-deug
!        mythrd = omp_get_thread_num()  !DKB-deug
!#endif                                 !DKB-deug
        !$OMP DO PRIVATE(filename, file_number, unit, stat, ytmp, y_input, input_array, derivative, n_main) &
        !$OMP PRIVATE (derivative_peak, charge, der_opt, fft_filter_output, i, j, single_event, multiple_event, od) &
        !$OMP PRIVATE (derivative_peaks, pos_peaks, neg_peaks, section, input) &
        !$OMP PRIVATE (rand_5) &
        
        !xxxPRIVATE(final_signal_end, end_time, cycles, points_to_average, freq_avg, freq_std_dev) &
        !xxxPRIVATE(mass_to_charge, mass_to_charge_std_dev, charge, charge_std_dev, mass, mass_std_dev) &
        
        !$OMP SCHEDULE(guided)
        !xxxx SCHEDULE(static)
    
        ! This DO loop assigns name of datafile to open
        ! This streamlines data acquisition and analysis
        DO  i = 0, number_of_files-1
                        
            ! Open input file
            unit = i + 10000               
        
            ! If files are .dat then use the following to open
            ! ================================================
            input_array = 0
              
            !$OMP CRITICAL (READ_FILE)
              
           !IF  (channel == 1)  THEN
           !    WRITE (file_number, fmt) i+10000
           !    filename = 'chA'//TRIM(file_number)//'.dat'
           !ELSE IF (channel == 2)  THEN
           !    WRITE (file_number, fmt) i+10000
           !    filename = 'chB'//TRIM(file_number)//'.dat'
           !END IF
         
            WRITE (file_number, fmt) unit
            filename = 'chA'//TRIM(file_number)//'.dat'
            !write(6,"(a,' :  619 : thread=',i4,'  file=',a)") procname,mythrd,TRIM(filename)   !DKB-debug
             
            OPEN (UNIT=unit, FILE=filename, STATUS="old", ACTION="read", IOSTAT=stat, FORM="binary") 
                ! FORM="unformatted" allows Fortran to interpret binary files properly.
            
           !INQUIRE (UNIT=unit, OPENED=od)
           !IF  (.NOT. od)  THEN
           !   GOTO 15
           !ELSE
               READ (unit, END=8) input_array
    8          CONTINUE
           !END IF  
            y_input(1:NPTS10) = input_array
            y_input(NPTS10+1:NPTS) = 0
            
            files_analyzed = files_analyzed + 1
            PRINT *, "Now reading from ", TRIM(filename)
            PRINT *, files_analyzed, "/", number_of_files            

            !$OMP END CRITICAL (READ_FILE)
            ! ================================================   
                        
            !Switch each pair of points since each pair gets switched in the data file
            !(e.g. the 4th point in the file is actually the 3rd point chronologically,
            !and vice versa)
            DO  j = 2, NPTS10, 2
               ytmp = y_input(j-1)
               y_input(j-1) = y_input(j)
               y_input(j) = ytmp
            END DO
            
            !At this point, y_input is a zero-padded array with raw data points in correct chronological order.
                
            IF (channel == 12) THEN
               ! Choose a random starting point in the cycle to realistically approximate unknown aspect of real signal
               !start = FLOOR(rand_5*(old_square_total_real + 1))
               !IF (start == (old_square_total_real + 1)) THEN
               !    start = old_square_total_real
               !END IF
                
               !Choose a random starting point in the cycle to realistically approximate unknown aspect of real signal
               !CALL RANDOM_NUMBER(rand_5)
                DO j = 1, signal_end
                   !y_input(j) = y_input(j) + sq_wave_array(NINT(rand_5*old_square_total_real)+j)
                    y_input(j) = y_input(j) + sq_wave_array(j)
                END DO
            END IF
            
            
           !Calculate the derivative of y_input
           !It would be faster only to do this for NPTS10 points (nominally 1000000) since points higher than that are 0.
            DO  j = 1, NPTS
                IF (j == 1)    THEN
                    derivative (j) = (y_input(j+1) - y_input(j))/((j+1)-j)
                    ! Denominator is 1, but writing it this way highlights that derivative = rise / run
                ELSE IF (j == NPTS)  THEN
                    derivative(j) = (y_input(j) - y_input(j-1))/(j-(j-1))
                ELSE           
                    derivative(j) = 0.5 * (((y_input(j+1) - y_input(j))/((j+1)-j)) + ((y_input(j) - y_input(j-1))/(j-(j-1)))) 
                    ! This is the average slope around the neighboring two points.
                END IF
            END DO        

            derivative_peaks = 0 
            pos_peaks = 0
            neg_peaks = 0

            DO  j = 1, NPTS
                IF  (derivative(j) > 5000)  THEN
                    IF  (pos_peaks == 0)    THEN
                        pos_peaks = pos_peaks + 1
                        derivative_peaks(1,pos_peaks) = j
                        ! First row of derivative_peaks contains position of positive-going peaks.
                    ELSE IF (j > (derivative_peaks(1,pos_peaks) + 50))  THEN
                        pos_peaks = pos_peaks + 1
                        derivative_peaks(1,pos_peaks) = j
                        ! The ELSE IF prevents every point of an impulse from being labeled a peak.
                    END IF
                    
                    IF  (pos_peaks == 400)  THEN
                        GOTO 15
                        ! 15 brings control to the end of the overall DO loop for this particular file.
                    END IF
                    
                ELSE IF (derivative(j) < -5000)    THEN
                    IF  (neg_peaks == 0)    THEN                    
                        neg_peaks = neg_peaks + 1
                        derivative_peaks(2,neg_peaks) = j
                        ! Second row of derivative_peaks contains position of negative-going peaks.
                    ELSE IF (j > (derivative_peaks(2,neg_peaks) + 50))  THEN                 
                        neg_peaks = neg_peaks + 1
                        derivative_peaks(2,neg_peaks) = j
                    END IF
                    
                    IF  (neg_peaks == 400)  THEN
                        GOTO 15
                    END IF
                    
                END IF
            END DO
            
            
            !Calculate the number of full trapping events within the file. Multiply by 2500 to
            !bring from units of ms to units of data points (since 2500 data points per 1 ms).
            !Note that this method is compatible with or without triggering between the BNC
            !pulsers and the FPGA boards.
            IF (derivative_peaks(1,1) <= 2500*(400-INT(400/length_multiplier)*length_multiplier)) THEN        
                expt_num = INT(400/length_multiplier)        
            ELSE IF (derivative_peaks(1,1) > 2500*(400-INT(400/length_multiplier)*length_multiplier)) THEN  
                expt_num = INT(400/length_multiplier)-1
            END IF
            
        
            DO  j = 1, expt_num
                !A trapping event is bracketed by a positive beginning pulse and a negative
                !ending pulse. The first pulse in an untriggered file may be either positive
                !or negative. The first of the following IF statements analyzes trapping
                !events for files in which the first pulse is positive. The second is for
                !files beginning with a negative pulse.
                IF ((derivative_peaks(2,j+1)-derivative_peaks(1,j) > 2500*REAL(max_trap_time+load_time/2.)-50) &
                    & .AND. (derivative_peaks(2,j+1)-derivative_peaks(1,j) < 2500*REAL(max_trap_time+load_time/2.)+50)) THEN
                   !How do this and the following IF statement work? Mathematically, I
                   !don't understand why the load_time/2 is there; it seems to make it
                   !impossible to fall within the right range.
                    
                    derivative_peak = derivative_peaks(1,j) + 2500*load_time/2.
                    !Why is this shifted?                    
                    
                    !Write an array which is composed of data just from the trapping event of interest.
                    input = 0
                    DO  k = 0, 2500*max_trap_time-1
                        input(k+1) = y_input(derivative_peak+k)
                    END DO
                    !With the above shift of derivative_peak, doesn't this take the FFT of
                    !points after the negative end pulse?
                     
                   !CALL fft_filter(input, NPTS, fft_filter_output)
                    CALL fft_filter(input, fft_filter_output)
                    input = fft_filter_output     
                    !At this point, input is still in the time domain but has had a 1 kHz
                    !high-pass filter passed over it.          
                    
                    n_main = 1
                    section = j            
        
                    CALL peakfinder_VI(filename, input, n_main, channel, charge, length_multiplier,  &
                           max_trap_time, dead_time, single_event, multiple_event, low_mass_cutoff,  &
                           high_mass_cutoff, sig_cutoff, section)
                    total_events = total_events + 1
                    total_single_events = total_single_events + single_event  
                    total_multiple_events = total_multiple_events + multiple_event
                    IF (single_event == 0 .AND. multiple_event == 0) total_empty_events = total_empty_events + 1
                    IF (single_event == 1) THEN
                        single_ion_filename(total_single_events) = filename
                        single_ion_section(total_single_events) = section
                    ELSE IF (multiple_event == 1) THEN
                        multiple_ion_filename(total_multiple_events) = filename
                        multiple_ion_section(total_multiple_events) = section
                    ELSE IF (single_event /= 1 .AND. multiple_event /= 1) THEN
                        no_ion_filename(total_events - total_single_events - total_multiple_events) = filename
                        no_ion_section(total_events - total_single_events - total_multiple_events) = section
                    END IF
                        
                    
                ELSE IF ((derivative_peaks(2,j)-derivative_peaks(1,j) > 2500*REAL(max_trap_time+load_time/2.)-50) &
                    & .AND. (derivative_peaks(2,j)-derivative_peaks(1,j) < 2500*REAL(max_trap_time+load_time/2.)+50)) THEN

                    derivative_peak = derivative_peaks(1,j) + 2500*load_time/2.
                    input = 0
 
                    DO  k = 0, 2500*max_trap_time-1
                        input(k+1) = y_input(derivative_peak+k)
                    END DO

                   !CALL fft_filter(input, NPTS, fft_filter_output)
                    CALL fft_filter(input, fft_filter_output)
                    !write(6,"(a,' :  795 : thread=',i4)") procname,mythrd        !DKB-debug
                    input = fft_filter_output               

                    n_main = 1
                    section = j            

                    CALL peakfinder_VI(filename, input, n_main, channel, charge, length_multiplier,  &
                           max_trap_time, dead_time, single_event, multiple_event, low_mass_cutoff,  &
                           high_mass_cutoff, sig_cutoff, section)
                    total_events = total_events + 1
                    total_single_events = total_single_events + single_event  
                    total_multiple_events = total_multiple_events + multiple_event 
                    IF (single_event == 0 .AND. multiple_event == 0) total_empty_events = total_empty_events + 1
                    IF (single_event == 1) THEN
                        single_ion_filename(total_single_events) = filename
                        single_ion_section(total_single_events) = section
                    ELSE IF (multiple_event == 1) THEN
                        multiple_ion_filename(total_multiple_events) = filename
                        multiple_ion_section(total_multiple_events) = section
                    ELSE IF (single_event /= 1 .AND. multiple_event /= 1) THEN
                        no_ion_filename(total_events - total_single_events - total_multiple_events) = filename
                        no_ion_section(total_events - total_single_events - total_multiple_events) = section
                    END IF                      
                END IF
            END DO

   15       CONTINUE
            
            !$OMP CRITICAL (CLOSE_FILE)
            CLOSE (unit)
            !$OMP END CRITICAL (CLOSE_FILE)
            
        END DO
        !$OMP END DO
        !$OMP END PARALLEL    
        
        
        !Call end clock time to determine program running time
        !total_time is in seconds
        CALL CLOCKX(end_run_time)
        total_time = (end_run_time - start_run_time) / 1.E6

        !Write program running time to "complete peak data" file
        WRITE(21,222) "End: Program running time =", total_time, "seconds"
        WRITE(21,223) "Single ion trapping events =", total_single_events, "/", total_events
        WRITE(21,224) "Multiple ion trapping events =", total_multiple_events, "/", total_events
        WRITE(21,225) "No ion trapping events =", (total_events - total_multiple_events - total_single_events), "/", total_events
  222   FORMAT(/, A27, TR5, F11.2, TR5, A7)
  223   FORMAT(/, A28, TR5, I7, A1, I7)
  224   FORMAT(/, A30, TR5, I7, A1, I7)
  225   FORMAT(/, A24, TR5, I7, A1, I7)
    
        !Close all output files and keep them (i.e., don't delete them)
        CLOSE (21, STATUS = "KEEP")
        CLOSE (22, STATUS = "KEEP")  
        
            
!===============================================================================
!Actions to take for creating a test signal, adding onto a noise file, and analyzing

    ELSE IF (channel == 3)  THEN    
    
        !Open file to store info from FFT scan across data section
        unit = 26
        filename = "full_test_signal.txt"
        CALL open_new_file (unit, filename)
        
        !Open output file for FFT magnitude of windowed data
        unit = 31
        filename = "mag_full_windowed.txt"
        CALL open_new_file (unit, filename)
        
        !Open output file for FFT magnitude of non-windowed data
        unit = 32
        filename = "mag_full_non-windowed.txt"
        CALL open_new_file (unit, filename)        
        
        !Open output file for test signal
        unit = 46
        filename = "test_signal only.txt"
        CALL open_new_file (unit, filename)
        
        !Open output file for test signal
        unit = 47
        filename = "test_signal_&_noise.txt"
        CALL open_new_file (unit, filename)
        
        unit = 108
        filename = "harmonic_analysis.txt"
        CALL open_new_file (unit, filename)
        
        !Headers for file #26
        WRITE (26,242) "Section_#", "#_Datapoints_Win", "Max_Freq_Win", "Mag_Win", "m/z_Win", "Charge_Win", "Mass_Win"
        242 FORMAT (T15, A10, TR5, A16, TR5, A12, TR10, A7, TR10, A7, TR10, A10, TR10, A8)
            
        !Select channel for noise files
        PRINT *, ""
        PRINT *, "Are noise files from channel A or channel B?"
        PRINT *, "Enter 1 for A..."
        PRINT *, "Enter 2 for B..."
        PRINT *, "Enter 3 for no noise..."
        READ *, noise_channel
        WRITE (97,*) "Are noise files from channel A or channel B?"
        WRITE (97,*) "Enter 1 for A..."
        WRITE (97,*) "Enter 2 for B..."
        WRITE (97,*) "Enter 3 for no noise..."
        WRITE (97,*) noise_channel
        
       !IF  (noise_channel < 3) THEN  
       !    ! User selected option for number of files in each run
       !    PRINT *, "Enter number of noise files..."
       !    READ *, number_of_files
       !    WRITE (97,*) "Enter number of noise files..."
       !    WRITE (97,*) number_of_files 
       !END IF 
    
        ! User selected option for SIMION test signal or generica test signal
        PRINT *, "Do you want to use a SIMION file or create a generic test signal?"
       !PRINT *, "Enter 1 for SIMION file (entitled test_signal.txt) (NOT SUPPORTED!)..."
        PRINT *, "Enter 1 for SIMION file (entitled signal_data.txt)..."
        PRINT *, "Enter 2 for generic test signal..."
        READ *, test_signal_option    
        WRITE (97,*) "Do you want to use a SIMION file or create a generic test signal?"
       !WRITE (97,*) "Enter 1 for SIMION file (entitled test_signal.txt) (NOT SUPPORTED!)..."
        WRITE (97,*) "Enter 1 for SIMION file (entitled signal_data.txt)..."
        WRITE (97,*) "Enter 2 for generic test signal..."
        WRITE (97,*) test_signal_option
            
        ! If generic test signal, must choose frequency, frequency shift, and length of signal
        IF (test_signal_option == 2)    THEN    
            ! Create test signal for troubleshooting/optimization purposes      
            PRINT *, "Enter desired frequency of test signal in Hz..."
            READ *, input_freq
            WRITE (97,*) "Enter desired frequency of test signal in Hz..."
            WRITE (97,*) input_freq
                    
           !PRINT *, "Enter second desired frequency of test signal..."
           !READ *, input_freq_2                
                
            ! User chosen frequency shift of test signal
            PRINT *, "By what percent would you like the frequency to shift?"
            PRINT *, "(0 = no frequency shift, 5 = 5%, etc.)"
            READ *, freq_shift    
            WRITE (97,*) "By what percent would you like the frequency to shift?"
            WRITE (97,*) "(0 = no frequency shift, 5 = 5%, etc.)"
            WRITE (97,*) freq_shift   
        END IF
                
        ! User selects length of test signal
       !PRINT *, "Enter length of test signal in ms..."
      !!PRINT *, "length < 65536"
      !!PRINT *, " 5000 =  2 ms"
      !!PRINT *, "20000 =  8 ms"
      !!PRINT *, "40000 = 16 ms"
      !!PRINT *, "65536 = 26.2 ms"
      !!READ *, signal_end
      !!WRITE (97,*) "Enter length of test signal in ms..."
      !!WRITE (97,*) signal_end
      !!signal_end = signal_end*2500  
        signal_end = NPTS10            !nominally 400*2500=1000000
        
        ! User selects charge of test signal                
        PRINT *, "Enter desired input charge (in electrons) for test signal..."
        READ *, electrons               
        WRITE (97,*) "Enter desired input charge (in electrons) for test signal..."
        WRITE (97,*) electrons
        
       !CLOSE (UNIT=97, STATUS='KEEP')
        
        ! Call initial clock time (to measure program running time)
        CALL CLOCKX(start_run_time)
        
          
        CALL RANDOM_SEED()
        CALL RANDOM_NUMBER(random_array)
        rand_1 = random_array(1)
        rand_2 = random_array(2)
        rand_3 = random_array(3)
        rand_4 = random_array(4)
        rand_5 = random_array(5)        
            
        ! Actions to perform if test signal is created from SIMION file
        IF (test_signal_option == 1)    THEN 
        
            ! Open input file
            unit = 801
            filename = "signal_data.txt"            
            
            PRINT *, "Now reading from ", TRIM(filename)
            OPEN (UNIT=unit, FILE=filename, STATUS="old", ACTION="read", IOSTAT=stat)
            
            DO i = 1, signal_end
                READ (801, *) sq_wave_array(i)
            END DO 
            
            sq_wave_array = sq_wave_array
            sq_wave_array = sq_wave_array*NINT(2.31*electrons)         
            
            ! The rest of the array is filled with zeroes
            DO  i = (signal_end + 1), NPTS
                sq_wave_array (i) = 0
            END DO        
            
            WRITE (46,777) sq_wave_array
  777       FORMAT (F25.10)
            
           !DEALLOCATE(sq_wave_array_allocatable, output)
                        
        ! Actions to perform if test signal is based on generic equation
        ELSE IF (test_signal_option == 2)   THEN
            ! Variable names and values obtained from Origin for creating test signal
            ! Test signal approximates a real signal and follows a Boltzmann distribution
            ! Details of test signal and fitting function obtained from Simion and Origin, respectively
            A1 = -3.35668061258415E-4   ! NCC TRAP 44 - CONETRAP
            A2 = 0.994194949645035      ! NCC TRAP 44 - CONETRAP
            x0 = -0.505420267974969     ! NCC TRAP 44 - CONETRAP
            dx = 0.040649051968342      ! NCC TRAP 44 - CONETRAP
                ! The values here work with the Boltzmann function in Origin (search for "Boltzmann" in Origin) to 
                ! produce a sigmoid. A1 and A2 are in close agreement with equation 2.5 of Nathan Contino's thesis.
                ! x0 and dx differ by the same factor from values in equation 2.5, producing a narrower sigmoid.
            
            square_total_real = 40.e6 / (input_freq * 16.)      ! So this is number of points per cycle of test signal
            old_square_total_real = square_total_real
            mult = INT((signal_end) / square_total_real)        ! So this is number of cycles of test signal
            !mult = INT((2500*length_multiplier) / square_total_real)
            remainder = 0      
            array_pointer = 0
            
            !Create test signal based on an 29.15% duty cycle (based off of Simion simulations
            !for NCC TRAP 44 - CONETRAP).  Also must take into account the offset of each cycle
            !based on limited time resolution (i.e., not every frequency will have an integer
            !number of datapoints per cycle due to a finite time resolution in data collection).
            !This offset changes the center of the Boltzmann distribution slightly so that not
            !every cycle looks exactly the same.
            expt_num = INT(400/length_multiplier)   !This definition of expt_num is only valid for triggered data files
            DO k = 1, expt_num
                i = 0           ! i is counting the number of test signal cycles written for current trapping event
                DO WHILE (array_pointer >= (k-1)*2500*length_multiplier .AND. array_pointer < k*2500*length_multiplier)
                    i = i + 1
               !DO  i = 1, mult-1       ! Each iteration of this outer loop writes data for one cycle of test signal
                   !DO  i = 1, INT((mult-1)/2.)       
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

                    !Non-integer portion of each cycle length is added onto frequency shift
                    !to create a steadily changing frequency.
                    remainder = square_total_real - INT(square_total_real)
                    square_total_real = (old_square_total_real * (1 - ((expt_num*freq_shift/100)*i/((mult-1))))) + remainder
                        !freq_shift is shift across NPTS10 data points (nominally 1000000). Multiply by
                        !expt_num to cause that shift over each trapping event
                   !square_total_real = (old_square_total_real * (1 - ((freq_shift/100)*i/((mult-1))))) + remainder
                   !square_total_real = (old_square_total_real * (1 - ((freq_shift/100)*i/((mult-1)*signal_end/(2500*length_multiplier))))) + remainder 
                END DO   
                square_total_real = old_square_total_real
            END DO
            
            
!            ! Second half of signal for simulated signal where the frequency shifts suddenly
!            square_total_real = 40.e6 / (input_freq_2 * 16.)
!            old_square_total_real = square_total_real
!            remainder = 0             
!             
!            DO  i = (INT((mult-1)/2.) + 1), mult-1         
!                DO  j = 1, INT(square_total_real)
!                    array_pointer = array_pointer + 1
!!                    ! NCC TRAP 44 - CONETRAP
!!                    x = -1.715265866 + ((j-1) * 3.430531732 / square_total_real) - (remainder / square_total_real * 100 * (3.430531732 / square_total_real))
!                    
!                    ! NCC TRAP 82 - CONETRAP
!                    x = -0.7745617799 + ((j-1) * 1.54912356 / square_total_real) - (remainder / square_total_real * 100 * (1.54912356 / square_total_real))
!                    IF (x <= 0)    THEN
!                        sq_wave_array (array_pointer) = (A2 + (A1-A2)/(1 + EXP((x-x0)/dx))) * NINT(2.31 * electrons)
!                    ELSE IF (x > 0) THEN
!                        sq_wave_array (array_pointer) = (A2 + (A1-A2)/(1 + EXP((-x-x0)/dx))) * NINT(2.31 * electrons)
!                    END IF  
!                    
!                    WRITE (46,777) sq_wave_array(array_pointer), INT(square_total_real), square_total_real
!                    777 FORMAT (F25.10, TR5, I5, TR5, F25.10)
!                    
!                END DO
!                
!                ! Non-integer portion of each cycle length is added onto frequency shift to create a steadily changing frequency
!                remainder = square_total_real - INT(square_total_real)
!                square_total_real = (old_square_total_real * (1 - ((freq_shift/100)*i/((mult-1)*signal_end/(2500*length_multiplier))))) + remainder
!                
!            END DO     
            
            
            
!            ! Sine Wave Test Signal
!            square_total_real = 40.e6 / (input_freq * 16.)
!            old_square_total_real = square_total_real
!            mult = INT((2500*length_multiplier) / square_total_real)
!            remainder = 0      
!            array_pointer = 0
!            
!            DO  i = 1, (mult-1) 
!!            DO  i = 1, INT((mult-1)/2.)          
!                DO  j = 1, INT(square_total_real)
!                    array_pointer = array_pointer + 1
!                    sq_wave_array (array_pointer) = NINT((2.31/2) * electrons) * SIN(2 * PI * input_freq * j * 16./40.e6) + NINT((2.31/2) * electrons)           
!                END DO
!                
!                ! Non-integer portion of each cycle length is added onto frequency shift to create a steadily changing frequency
!                remainder = square_total_real - INT(square_total_real)
!                square_total_real = (old_square_total_real * (1 - ((freq_shift/100)*i/((mult-1)*signal_end/(2500*length_multiplier))))) + remainder
!                
!            END DO       
!            
!            
!            square_total_real = 40.e6 / (input_freq_2 * 16.)
!            old_square_total_real = square_total_real
!            remainder = 0     
!            array_pointer = 0             
!             
!            DO  i = 1, (mult-1) 
!!            DO  i = (INT((mult-1)/2.) + 1), mult-1         
!                DO  j = 1, INT(square_total_real)
!                    array_pointer = array_pointer + 1
!                    sq_wave_array (array_pointer) = NINT((2.31/2) * electrons) * SIN(2 * PI * input_freq * j * 16./40.e6) + NINT((2.31/2) * electrons) + sq_wave_array (array_pointer)
!                END DO
!                
!                ! Non-integer portion of each cycle length is added onto frequency shift to create a steadily changing frequency
!                remainder = square_total_real - INT(square_total_real)
!                square_total_real = (old_square_total_real * (1 - ((freq_shift/100)*i/((mult-1)*signal_end/(2500*length_multiplier))))) + remainder
!                
!            END DO   
                               
            
            ! The rest of the array is filled with zeroes
            DO  i = (mult * INT(square_total_real)), NPTS
                sq_wave_array (i) = 0
            END DO        
        
            WRITE (46,778) sq_wave_array
            778 FORMAT (F25.10)
        END IF
        
    
        IF  (noise_channel < 3) THEN
            
!!            number_of_sections = FLOOR(400./length_multiplier)
!!            n_main = CEILING(rand_4 * number_of_sections)
!            ! For simplicity, just take the first trapping event of each file to add test signal
!            n_main = 1
!!            n_main = CEILING(rand_4 * 12)
!!            IF  (n_main == 0)   THEN
!!                n_main = 1
!!            ELSE IF (n_main == 13)  THEN
!!                n_main = 12
!!            END IF
!               
!            ! Randomly select a noise file
!            j = FLOOR(rand_3 * number_of_files) 
!            IF  (j ==  number_of_files) THEN
!                j = number_of_files
!            END IF
                       
            ! Open input noise file
            unit = 1  
            PRINT *, "Enter noise file name (chA#####.txt or chA#####.dat)..."
            READ *, filename        
            WRITE (97,*) "Enter noise file name (chA#####.txt or chA#####.dat)..."
            WRITE (97,*) filename
            
            ! Select file section        
            PRINT *, "Enter section # on which to perform DFT..."
            READ *, n_main 
            WRITE (97,*) "Enter section # on which to perform DFT..."
            WRITE (97,*) n_main                  
                  
            IF  (data_type == 1)    THEN
            
!                IF  (noise_channel == 1)    THEN                
!                    WRITE (file_number, fmt) j+10000
!                    filename = 'chA'//TRIM(file_number)//'.dat' 
!                ELSE IF (noise_channel == 2)    THEN                
!                    WRITE (file_number, fmt) j+10000
!                    filename = 'chB'//TRIM(file_number)//'.dat'            
!                END IF          
            
                PRINT *, "Now reading from ", TRIM(filename)
                OPEN (UNIT=unit, FILE=filename, STATUS="old", ACTION="read", IOSTAT=stat, FORM="binary")             
            
                input_array = 0
                y_input = 0
                ! Read input file into array "y_input"
                READ (unit, END=10) input_array     ! Read filename values and store in array "y"
   10           CONTINUE   
            
                DO  k = 1, NPTS10
                    y_input(k) = input_array(k)
                END DO
            
            ELSE IF (data_type == 2)   THEN
                   
!                IF  (noise_channel == 1)    THEN                
!                    WRITE (file_number, fmt) j+10000
!                    filename = 'chA'//TRIM(file_number)//'.txt' 
!                ELSE IF (noise_channel == 2)    THEN                
!                    WRITE (file_number, fmt) j+10000
!                    filename = 'chB'//TRIM(file_number)//'.txt'            
!                END IF    
             
                PRINT *, "Now reading from ", TRIM(filename)
                OPEN (UNIT=unit, FILE=filename, STATUS="old", ACTION="read", IOSTAT=stat)           
                      
                y_input = 0
                
                ! Read input file into array "y_input"
                READ (unit, *, END=81) y_input    ! Read filename values and store in array "y"
   81           CONTINUE 

            END IF

            CLOSE (unit)               

            !Switch each pair of points since each pair gets switched in the data file.
            DO  j = 2, NPTS10, 2
               ytmp = y_input(j-1)
               y_input(j-1) = y_input(j)
               y_input(j) = ytmp
            END DO

            !Add test signal onto noise file for troubleshooting/optimization purposes    
            DO  i = 1, signal_end
               !y_input(i) = y_input(i) + sq_wave_array (NINT(rand_5 * old_square_total_real) + i)
                y_input(i) = y_input(i) + sq_wave_array (i)
            END DO

            WRITE (47,779) y_input
  779       FORMAT (F25.10)

            ! Calculate the derivative of y_input
            DO  j = 1, NPTS
                IF  (j == 1)    THEN
                    derivative (j) = (y_input(j+1) - y_input(j))/((j+1)-j)
                ELSE IF (j == NPTS)  THEN
                    derivative (j) = (y_input(j) - y_input(j-1))/(j-(j-1))
                ELSE           
                    derivative (j) = 0.5 * (((y_input(j+1) - y_input(j))/((j+1)-j)) + ((y_input(j) - y_input(j-1))/(j-(j-1)))) 
                END IF
            END DO        
            
              
            derivative_peaks = 0 
            pos_peaks = 0
            neg_peaks = 0
            
            DO  j = 1, NPTS
                IF  (derivative(j) > 5000)  THEN
                    IF  (pos_peaks == 0)    THEN
                        pos_peaks = pos_peaks + 1
                        derivative_peaks(1,pos_peaks) = j
                    ELSE IF (j > (derivative_peaks(1,pos_peaks) + 50))  THEN
                        pos_peaks = pos_peaks + 1
                        derivative_peaks(1,pos_peaks) = j
                    END IF
                    
                    IF  (pos_peaks == 400)  THEN
                        GOTO 16
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
                        GOTO 16
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
                                        
               !CALL fft_filter(input, NPTS, fft_filter_output)
                CALL fft_filter(input, fft_filter_output)
                input = fft_filter_output
                
            ELSE IF (((derivative_peaks(2,(n_main)) - derivative_peaks(1,n_main)) > ((2500 * REAL(max_trap_time + (load_time / 2.))) - 50)) &
                & .AND. ((derivative_peaks(2,(n_main)) - derivative_peaks(1,n_main)) < ((2500 * REAL(max_trap_time + (load_time / 2.))) + 50))) THEN
                
                derivative_peak = derivative_peaks(1,n_main) + (2500 * load_time / 2.)

                input = 0
                DO  i = 0, ((2500 * max_trap_time) - 1)
                    input(i+1) = y_input(derivative_peak + i)
                END DO
               !CALL fft_filter(input, NPTS, fft_filter_output)
                CALL fft_filter(input, fft_filter_output)
                input = fft_filter_output
                                 
            END IF                 
               
           !! Write name of file currently being analyzed to output files
           !WRITE (26,3000) filename, "Section_#:", n_main
           !3000 FORMAT (A10, TR5, A10, TR5, I2)
            
        ELSE
            n_main = 1
            y_input = 0
        END IF
        
        
        ! Choose a random starting point in the cycle to realistically approximate unknown aspect of real signal
        start = FLOOR(rand_5 * (old_square_total_real + 1))
        IF  (start == (old_square_total_real + 1))  THEN
           start = old_square_total_real
        END IF
        
        section = n_main
        n_main = 1     
            
        ! Add test signal onto noise file for troubleshooting/optimization purposes
        IF (noise_channel == 3) THEN    
            DO  i = 1, signal_end
              !input(i) = input(i) + sq_wave_array (NINT(rand_5 * old_square_total_real) + i)
               input(i) = input(i) + sq_wave_array (i)
            END DO
            WRITE (47,780) input
  780       FORMAT (F25.10)
        END IF
        
        
       !section = n_main
        ! Analyze test signal (added to noise file) to determine accuracy and limits of program
        CALL peakfinder_VI(filename, input, n_main, channel, charge, length_multiplier, max_trap_time, dead_time, single_event, multiple_event, low_mass_cutoff, high_mass_cutoff, sig_cutoff, section)
        
        ! FFT outputs for both windowed and non-windowed data
        CALL fftransform(input, n_main, signal_end, fft_magnitude, length_multiplier, dead_time)
        WRITE (31,999) fft_magnitude
  999   FORMAT (F25.10)  
        
        CALL fft_no_window(input, n_main, signal_end, fft_magnitude, length_multiplier, dead_time) 
        WRITE (32,999) fft_magnitude
      
   16   CONTINUE
      
        ! Call end clock time to determine program running time
        ! total_time is in seconds
        CALL CLOCKX(end_run_time)
        total_time = (end_run_time - start_run_time) / 1.E6
    
        ! Close all output files and keep them (i.e., don't delete them)
        CLOSE (26, STATUS = "KEEP")         
        CLOSE (31, STATUS = "KEEP") 
        CLOSE (32, STATUS = "KEEP") 
        CLOSE (46, STATUS = "KEEP")
        CLOSE (47, STATUS = "KEEP")   
            
                               
!===============================================================================
!Actions to take for creating simulated signal, adding onto a noise file, and analyzing
    ELSE IF (channel == 4)  THEN

        ! Call initial clock time (to measure program running time)
        CALL CLOCKX(start_run_time)

        expt_num = INT(400 / length_multiplier) - 1

        CALL simulations(channel, data_type, fmt, length_multiplier, max_trap_time, load_time, dead_time, expt_num)

        ! Call end clock time to determine program running time
        ! total_time is in seconds
        CALL CLOCKX(end_run_time)
        total_time = (end_run_time - start_run_time) / 1.E6


!===============================================================================
!Actions to take to look at data from a specific file & section
    ELSE IF (channel == 5)  THEN
               
        ! Open file to store info from FFT scan across data section
        unit = 21
        filename = "peak_info_FFT_section.txt"
        CALL open_new_file (unit, filename)            
        
        ! Open file to store info from FFT scan across data section
        unit = 26
        filename = "full_info_FFT_section.txt"
        CALL open_new_file (unit, filename)            
        
        ! Headers for verbose output files
        WRITE (21,221) "File", "Section_#", "Signal_Length", "Max_Freq", "Magnitude"
        221 FORMAT (A4, TR10, A10, TR5, A13, T49, A9, T71, A10)
        
        WRITE (26,241) "Section_#", "End_of_Signal", "Max_Freq_Win", "Mag_Win", "m/z_Win", "Charge_Win", "Mass_Win"
        241 FORMAT (T15, A10, TR5, A16, TR5, A12, TR10, A7, TR10, A7, TR10, A10, TR10, A8)
    
        CALL CLOCKX(start_run_time)     
        
        ! Perform FFT on selected file & section
        CALL fft_of_selected_section(channel, length_multiplier, load_time, data_type)
        
        ! Call end clock time to determine program running time
        ! total_time is in microseconds
        CALL CLOCKX(end_run_time)
        total_time = (end_run_time - start_run_time) / 1.E6
    
        ! Close all output files and keep them (i.e., don't delete them)
        CLOSE (21, STATUS = "KEEP") 
        CLOSE (26, STATUS = "KEEP") 
    
     
!===============================================================================
!Investigate frequency shift      
    ELSE IF (channel == 6 .OR. channel == 7)  THEN            
 
        IF  (channel == 7)  THEN

            ! Open file to store info from FFT scan across data section
            unit = 27
            filename = "freq_shift_data.txt"
            CALL open_new_file (unit, filename)            

            ! Headers for verbose output files
            WRITE (27,232) "File", "Section_#", "Freq_Avg_1", "Freq_Avg_1_Stdev", "Freq_Avg_2", "Freq_Avg_2_Stdev", "Inflection_Pt", "Freq_Shift", &
                           & "Mass_Avg_1", "Mass_Avg_1_Stdev", "Mass_Avg_2", "Mass_Avg_2_Stdev", "Mass_Shift", &
                           & "m/z_Avg_1", "m/z_Avg_1_Stdev", "m/z_Avg_2", "m/z_Avg_2_Stdev", "m/z_Shift", &
                           & "Charge_Avg_1", "Charge_Avg_1_Stdev", "Charge_Avg_2", "Charge_Avg_2_Stdev", "Charge_Shift"
            232 FORMAT (A4, TR10, A10, TR5, A10, TR5, A16, TR5, A10, TR5, A16, TR5, A13, TR5, A10, TR8, &
                        & A10, TR7, A16, TR7, A10, TR7, A16, TR5, A10, TR7, &
                        & A9, TR5, A15, TR5, A9, TR5, A15, TR5, A9, TR5, &
                        & A12, TR5, A18, TR5, A12, TR5, A18, TR5, A12)
        
        END IF

        CALL CLOCKX(start_run_time)

        ! Perform FFT on selected file & section
        CALL find_frequency_shift(channel, length_multiplier, load_time, dead_time)
 
        ! Call end clock time to determine program running time
        ! total_time is in microseconds
        CALL CLOCKX(end_run_time)
        total_time = (end_run_time - start_run_time) / 1.E6
 
        IF  (channel == 7)  THEN
            CLOSE (27, STATUS = "KEEP")
        END IF


!===============================================================================
    ELSE IF (channel == 8)  THEN
      
        PRINT *, "Enter last file number..."
        READ *, number_of_files
        WRITE (97,*) "Enter last file number..."
        WRITE (97,*) number_of_files
        number_of_files = number_of_files + 1 - 10000 

        IF  (derivative_peaks(1,1) <= (2500 * (400 - (INT(400 / length_multiplier) * length_multiplier)))) THEN        
            expt_num = INT(400 / length_multiplier)        
        ELSE IF (derivative_peaks(1,1) > (2500 * (400 - (INT(400 / length_multiplier) * length_multiplier)))) THEN  
            expt_num = INT(400 / length_multiplier) - 1
        END IF        

        CALL CLOCKX(start_run_time)     

        ! Perform FFT on selected file & section
        CALL frequency_data(channel, length_multiplier, load_time, dead_time, total_single_events, total_multiple_events, number_of_files, expt_num, run_number)
        
        ! Call end clock time to determine program running time
        ! total_time is in microseconds
        CALL CLOCKX(end_run_time)
        total_time = (end_run_time - start_run_time) / 1.E6
        
        ! Write program running time to "complete peak data" file
        WRITE(21,222) "End: Program running time =", total_time, "seconds"
        WRITE(21,223) "Single ion trapping events =", total_single_events, "/", (expt_num * number_of_files)
        WRITE(21,224) "Multiple ion trapping events =", total_multiple_events, "/", (expt_num * number_of_files)
        WRITE(21,225) "No ion trapping events =", ((expt_num * number_of_files) - total_multiple_events - total_single_events), "/", (expt_num * number_of_files)
           
        ! Close all output files and keep them (i.e., don't delete them)
        CLOSE (21, STATUS = "KEEP")
        CLOSE (22, STATUS = "KEEP")  
        
                
!===============================================================================
!Actions to take for calibration of trap        
    ELSE IF (channel == 9)  THEN
    
        CALL CLOCKX(start_run_time)     
        
        expt_num = INT(400 / length_multiplier) - 1
        
        ! Calibration options subroutine
        CALL trap_charge_calibration(length_multiplier, load_time, expt_num, data_type)
    
        ! Call end clock time to determine program running time
        ! total_time is in seconds
        CALL CLOCKX(end_run_time)
        total_time = (end_run_time - start_run_time) / 1.E6
    
    
!===============================================================================
!Noise analysis
    ELSE IF (channel == 10)  THEN
        
        CALL CLOCKX(start_run_time)
        
        ! Calculates voltage and standard deviation
        CALL noise_analysis(length_multiplier, load_time, dead_time, data_type)
           
        ! Call end clock time to determine program running time
        ! total_time is in seconds
        CALL CLOCKX(end_run_time)
        total_time = (end_run_time - start_run_time) / 1.E6
    
    
!===============================================================================
!Find voltage from data obtained on lab computer (for calibration purposes)
    ELSE IF (channel == 11)  THEN
        
        CALL CLOCKX(start_run_time)
        
        ! Calculates voltage and standard deviation
        CALL voltage_finder
           
        ! Call end clock time to determine program running time
        ! total_time is in seconds
        CALL CLOCKX(end_run_time)
        total_time = (end_run_time - start_run_time) / 1.E6
    
    END IF 
    
!=========================================================================================
!=========================================================================================


!Write data to event types file
    IF (channel == 1 .OR. channel == 2 .OR. channel == 12) THEN
        WRITE (single_ion_section_char, "(I2)") single_ion_section
        WRITE (multiple_ion_section_char, "(I2)") multiple_ion_section
        WRITE (no_ion_section_char, "(I2)") no_ion_section
        DO i = 1, number_of_events
            IF (single_ion_section(i) == 0 .AND. multiple_ion_section(i) == 0 .AND. no_ion_section(i) == 0) EXIT
            IF (single_ion_filename(i) == '0') THEN
                single_ion_filename(i) = ''; single_ion_section_char(i) = ''
            END IF
            IF (multiple_ion_filename(i) == '0') THEN
                multiple_ion_filename(i) = ''; multiple_ion_section_char(i) = ''
            END IF
            IF (no_ion_filename(i) == '0') THEN
                no_ion_filename(i) = ''; no_ion_section_char = ''
            END IF
            WRITE (20, 693) single_ion_filename(i), single_ion_section_char(i), multiple_ion_filename(i), multiple_ion_section_char(i), no_ion_filename(i), no_ion_section_char(i)
            693 FORMAT (A12, TR6, A4, TR8, A12, TR10, A4, TR14, A12, TR4, A4)
        END DO
        CLOSE (UNIT=20, STATUS='KEEP')
    END IF
    
 !Print program running time to screen
    PRINT *, "End: Program running time=", total_time, " seconds"          
       
    !PAUSE "Press Enter to continue..."
    
    
    END PROGRAM FFT
    
    
    
!*****************************************************************************************
!Subroutine for opening and old file
    SUBROUTINE open_old_file (open_unit, open_name, stat)
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: open_unit          !Unit number for r/w files
    CHARACTER (*), INTENT(IN) :: open_name    !Dummy file name for r/w files
    INTEGER :: stat                           !Error message when opening files 

    INTEGER :: status                         !Error message when opening files 
    
    PRINT *, "Now reading from ", TRIM(open_name)
    OPEN (UNIT=open_unit, FILE=open_name, STATUS="old", ACTION="read", IOSTAT=status, FORM="binary")
    stat = status    
    END SUBROUTINE open_old_file
    
    
!*****************************************************************************************
!Subroutine for creating new file
    SUBROUTINE open_new_file (new_unit, new_name)
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: new_unit           !Unit number for r/w files
    CHARACTER (*), INTENT(IN) :: new_name     !Dummy file name for r/w files

    INTEGER :: status                         !Error message when opening files 
    
    OPEN (UNIT=new_unit, FILE=new_name, STATUS="replace", ACTION="write", IOSTAT=status)
           ! What do STATUS="replace", ACTION="write", and IOSTAT=status mean?
    PRINT *, "Error creating file? ", new_name, "-0- means success! ", status
    END SUBROUTINE open_new_file
    
       
!*****************************************************************************************
    SUBROUTINE avg_stdev (array, n_start, n_end, avg, stdev)
    IMPLICIT NONE

    REAL, DIMENSION (*) :: array
    INTEGER :: n_start, n_end 
    REAL    :: avg, stdev

    REAL    :: sum , diffsq
    INTEGER :: i

    sum = 0
    DO  i = n_start, n_end
      sum = sum + array(i)
    END DO            
    avg = sum / (n_end - n_start + 1)  
    diffsq = 0 
    DO  i = n_start, n_end
      diffsq = diffsq + ((array(i) - avg)**2)
    END DO                    
    stdev = SQRT(diffsq/(n_end - n_start))  
    END SUBROUTINE avg_stdev
    
    
!*****************************************************************************************
!Subroutine for rolling average high-pass filter (here 9.375 kHz)
    SUBROUTINE filter (filter_input, filter_output)
    USE trap_mod
    IMPLICIT NONE

   !INTEGER, PARAMETER :: NPTS=1048576
    REAL(4), DIMENSION(1:NPTS), INTENT(IN)  :: filter_input
    REAL(4), DIMENSION(1:NPTS), INTENT(OUT) :: filter_output

    REAL(4) :: filter_summation
    INTEGER :: i, j

!Take 100 point rolling average and subtract from each point.
!Acts as a high pass filter ... gets rid of everyting below 9.375 kHz.
    DO  i = 101, NPTS-100
       filter_summation = 0
       DO  j = i-100, i+100
          filter_summation = filter_summation + filter_input(j)
       END DO
       filter_output(i) = filter_input(i) - (filter_summation/201)
    END DO

!Make sure everything before rolling average start is zero
    DO  i = 1, 100
       filter_output(i) = 0
    END DO 

!Make sure everything before rolling average start is zero
    DO  i = NPTS-99, NPTS
       filter_output(i) = 0
    END DO 

    END SUBROUTINE filter
     
     
!*****************************************************************************************
! Subroutine for FFT high-pass filter. This routine applies a high-pass filter to time-
! domain values in fft_input, and returns the filtered output in fft_filter_output. The
! routine applies an FFT to the input, zeros out the low frequencies of the transformed
! data, and transforms back to the time-domain. Currently, the frequency cut-off is
! 1 kHz.

    SUBROUTINE fft_filter(fft_input, fft_filter_output)
    USE mkl_dfti
    USE trap_mod
    USE omp_lib     !DKB-debug
    IMPLICIT NONE

   !INTEGER, PARAMETER :: NPTS = 1048576                  !Number of datapoints
   !INTEGER, INTENT(IN)  :: npts                          !Number of datapoints
    REAL(4), INTENT(IN)  :: fft_input(NPTS)               !Input values
    REAL(4), INTENT(OUT) :: fft_filter_output(NPTS)       !Output values

    TYPE(DFTI_DESCRIPTOR), POINTER :: My_Desc1_Handle
    INTEGER :: Status                                     !Sets DFT Parameters
    INTEGER :: istat                                      !Allocate status
    REAL(4), ALLOCATABLE :: fft_forward_output(:)         !DFT output values

    INTEGER  :: mythrd = 0            !DKB-debug
    CHARACTER(12), PARAMETER :: procname="fft_filter"  !DKB-debug

!#ifdef _OPENMP                        !DKB-debug
!    mythrd = omp_get_thread_num()     !DKB-debug
!#endif                                !DKB-debug
    !write(6,"(a,' : 1693 : thread=',i4,'  npts=',i8)") procname,mythrd,npts   !DKB-debug
    ALLOCATE(fft_forward_output(NPTS))

    Status = DftiCreateDescriptor(My_Desc1_Handle, DFTI_SINGLE, DFTI_REAL, 1, NPTS)
        !Format of DftiCreateDescriptor: DftiCreateDescriptor(descriptor handle, precision,
        !domain of forward (time to frequency) transform, dimension, length)
        !"Precision" may be DFTI_SINGLE or DFTI_DOUBLE. "Domain" may be DFTI_COMPLEX or DFTI_REAL.
    Status = DftiSetValue(My_Desc1_Handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
        !Format of DftiSetValue: DftiSetValue(descriptor handle, configuration parameter, configuration value)
        !See Tables 11-3 and 11-5 of Intel Math Kernel Library Reference Manual.
        !DFTI_PLACEMENT places where the output goes. DFTI_NOT_INPLACE: "Output does not overwrite input"
    Status = DftiSetValue(My_Desc1_Handle, DFTI_NUMBER_OF_USER_THREADS, 1)
    Status = DftiSetValue(My_Desc1_Handle, DFTI_FORWARD_SCALE, 1.0)
        !DFTI_FORWARD_SCALE is "scale factor for forward transform."
    Status = DftiSetValue( My_Desc1_Handle, DFTI_BACKWARD_SCALE, 1.0/REAL(NPTS))
        !According to p2478 of the MKL user's guide, setting the forward scale to 1 and the backward scale to 1/NPTS makes the backward
        !transform the inverse of the forward transform.
    Status = DftiCommitDescriptor(My_Desc1_Handle)
        !This initializes the descriptor described by DftiCreateDescriptor.
    Status = DftiComputeForward(My_Desc1_Handle, fft_input, fft_forward_output)
        !This actually computes the forward FT using the factor e^(-i*2*Pi/n).
        !Is the n above NPTS? Thus, the greater NPTS, the better the resolution?

    !write(6,"(a,' : 1716 : thread=',i4)") procname,mythrd        !DKB-debug
    fft_forward_output(1:840) = 0      ! 1 kHz filter
        !From where did the value 840 come? How does each element of the input correspond to the output?
   !fft_forward_output(1:8390) = 0    ! 10 kHz filter
    !write(6,"(a,' : 1720 : thread=',i4)") procname,mythrd        !DKB-debug
 
    Status = DftiComputeBackward(My_Desc1_Handle, fft_forward_output, fft_filter_output)
        !This computes the backward FT using the factor e^(i*2*Pi/n).
        !After the filter has been applied in frequency space, this transforms back to the time domain.
    Status = DftiFreeDescriptor(My_Desc1_Handle)
        !This frees the memory allocated for the descriptor.
    DEALLOCATE(fft_forward_output)
    !write(6,"(a,' : 1730 : thread=',i4)") procname,mythrd        !DKB-debug

    END SUBROUTINE fft_filter
 

!*****************************************************************************************
! This is subroutine is similar to fft_filter, except it applies a 10 kHz high pass
! filter.

    SUBROUTINE ten_khz_fft_filter(fft_input, fft_filter_output)
    USE mkl_dfti
    USE trap_mod
    USE omp_lib     !DKB-debug
    IMPLICIT NONE

   !INTEGER, PARAMETER :: NPTS = 1048576                !Number of datapoints
    REAL(4), INTENT(IN)  :: fft_input(NPTS)             !Input values
    REAL(4), INTENT(OUT) :: fft_filter_output(NPTS)     !Output values

    TYPE(DFTI_DESCRIPTOR), POINTER :: My_Desc1_Handle
    INTEGER :: Status                                   !Sets DFT Parameters
    INTEGER :: istat                                    !Allocate status
    REAL(4), ALLOCATABLE :: fft_forward_output(:)       !DFT output values

    INTEGER  :: mythrd = 0            !DKB-debug
    CHARACTER(12), PARAMETER :: procname="ten_khz_fft*"  !DKB-debug
!#ifdef _OPENMP                        !DKB-debug
!    mythrd = omp_get_thread_num()     !DKB-debug
!#endif                                !DKB-debug

    !write(6,"(a,' : 1751 : thread=',i4,'  npts=',i8)") procname,mythrd,npts   !DKB-debug
    ALLOCATE(fft_forward_output(NPTS))

    Status = DftiCreateDescriptor(My_Desc1_Handle, DFTI_SINGLE, DFTI_REAL, 1, NPTS)
    Status = DftiSetValue(My_Desc1_Handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
    Status = DftiSetValue(My_Desc1_Handle, DFTI_NUMBER_OF_USER_THREADS, 1)
    Status = DftiSetValue(My_Desc1_Handle, DFTI_FORWARD_SCALE, 1.0)
    Status = DftiSetValue( My_Desc1_Handle, DFTI_BACKWARD_SCALE, 1.0/REAL(NPTS))
    Status = DftiCommitDescriptor(My_Desc1_Handle)
    Status = DftiComputeForward(My_Desc1_Handle, fft_input, fft_forward_output)

   !fft_forward_output(1:840)  = 0    !  1 kHz filter
    fft_forward_output(1:8390) = 0    ! 10 kHz filter

    Status = DftiComputeBackward(My_Desc1_Handle, fft_forward_output, fft_filter_output)
    Status = DftiFreeDescriptor(My_Desc1_Handle)
    DEALLOCATE(fft_forward_output)
    !write(6,"(a,' : 1768 : thread=',i4)") procname,mythrd        !DKB-debug

    END SUBROUTINE ten_khz_fft_filter


!*****************************************************************************************
!Band pass FFT filter
    SUBROUTINE band_pass_fft_filter(fft_input, freq_input, fft_filter_output)
    USE mkl_dfti
    USE trap_mod
    IMPLICIT NONE

   !INTEGER, PARAMETER   :: NPTS = 1048576           !Number of datapoints
    REAL(4), INTENT(IN)  :: fft_input(NPTS)          !Input values
    REAL(4), INTENT(IN)  :: freq_input               !Frequency of signal
    REAL(4), INTENT(OUT) :: fft_filter_output(NPTS)  !Output values

    TYPE (DFTI_DESCRIPTOR), POINTER :: My_Desc1_Handle
    INTEGER :: Status                                !Sets DFT Parameters
    INTEGER :: istat                                 !Allocate status
    INTEGER :: low_cutoff                            !Low frequency cutoff for band pass filter
    INTEGER :: high_cutoff                           !High frequency cutoff for band pass filter
    REAL(4), ALLOCATABLE :: fft_forward_output(:)    !DFT output values

    ALLOCATE(fft_forward_output(NPTS))
    Status = DftiCreateDescriptor(My_Desc1_Handle, DFTI_SINGLE, DFTI_REAL, 1, NPTS)
    Status = DftiSetValue(My_Desc1_Handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
    Status = DftiSetValue(My_Desc1_Handle, DFTI_NUMBER_OF_USER_THREADS, 1)
    Status = DftiSetValue(My_Desc1_Handle, DFTI_FORWARD_SCALE, 1.0)
    Status = DftiSetValue(My_Desc1_Handle, DFTI_BACKWARD_SCALE, 1.0/REAL(NPTS))
    Status = DftiCommitDescriptor(My_Desc1_Handle)
    Status = DftiComputeForward(My_Desc1_Handle, fft_input, fft_forward_output)

    low_cutoff  = INT((freq_input-250) / 40.E6 * 16 * NPTS)
    high_cutoff = INT((freq_input+250) / 40.E6 * 16 * NPTS)
    fft_forward_output(1:(low_cutoff*2)) = 0
    fft_forward_output((high_cutoff*2):NPTS) = 0

    Status = DftiComputeBackward(My_Desc1_Handle, fft_forward_output, fft_filter_output)
    Status = DftiFreeDescriptor(My_Desc1_Handle)
    DEALLOCATE(fft_forward_output)

    END SUBROUTINE band_pass_fft_filter
    
    
!*****************************************************************************************
!Subroutine for FFT with a window function

    SUBROUTINE fftransform(data_in, n_fft, signal_end, fft_mag, len_mult, dead_time)
    USE mkl_dfti
    USE trap_mod
    IMPLICIT NONE
    
   !INTEGER, PARAMETER :: NPTS0 = 65536            !Base number of datapoints  
   !INTEGER, PARAMETER :: FREQ_MULT = 16           !Multiplier for length of FFT     
   !INTEGER, PARAMETER :: NPTS = FREQ_MULT*NPTS0   !Number of datapoints  
   !REAL(4), PARAMETER :: PI = 3.14159265                      

    REAL(4), INTENT(IN)  :: data_in(NPTS)    !Input values 
    INTEGER, INTENT(IN)  :: n_fft            !Pointer for file section
    INTEGER, INTENT(IN)  :: signal_end       !Length of signal on which to perform DFT
    REAL(4), INTENT(OUT) :: fft_mag(NPTS/2)  !DFT magnitude values                          
   !REAL(4), INTENT(OUT) :: phase(NPTS/2)    !DFT phase values   
    INTEGER, INTENT(IN)  :: len_mult         !Multiplier for length of each trapping
                                             !   event (multiple of 1 ms or 2,500 datapoints)
    INTEGER, INTENT(IN)  :: dead_time        !Number of points of dead time due to A250
                                             !   recovering from impulse noise of trap gate pulse
    
    TYPE(DFTI_DESCRIPTOR), POINTER :: My_Desc1_Handle  
    INTEGER :: Status              !Sets DFT Parameters
    INTEGER :: istat               !Allocate status               
    INTEGER :: i                   !Loop variables
    INTEGER :: nDKB
    REAL(4) :: sigma_fft           !Sigma value for Gaussian window function 
   !REAL(4) :: alpha_fft           !Alpha value for Hann-Poisson window function

    REAL(4), ALLOCATABLE :: fft_input(:)       !Input values after applying rolling average filter
    REAL(4), ALLOCATABLE :: fft_output(:)      !DFT output values

    ALLOCATE(fft_input(NPTS))
    ALLOCATE(fft_output(NPTS))
    
!Select section of data file and cut off first "dead_time" points to account for pulse
!generator impulse noise and possible packet of untrapped ions passing through the
!detection tube.
    DO  i = 1+(n_fft-1)*(2500*len_mult)+dead_time, (n_fft-1)*(2500*len_mult)+signal_end
       fft_input(i-(n_fft-1)*(2500*len_mult)-dead_time) = data_in(i)
    END DO

    nDKB = signal_end-(dead_time+1)
    
!Make sure everything after signal of interest is zero
    DO  i = signal_end-(dead_time-1), NPTS
       fft_input(i) = 0
    END DO                      
    
!Gaussian Window
    sigma_fft = 0.45
    DO  i = 0, nDKB
       fft_input(i+1) = EXP(-0.5*(((i-REAL(nDKB)/2)/(sigma_fft*(REAL(nDKB)/2)))**2)) &
                        * fft_input(i+1) 
    END DO  
                   
!!Hann Window
!    DO  i = 0, nDKB
!       IF  (i == 0)   THEN
!          fft_input(i+1) = (0.5 + 0.5*cos((2*PI*i)/REAL(nDKB))) * fft_input(i+1)  
!       ELSE
!          fft_input(i+1) = (0.5 - 0.5*cos((2*PI*i)/REAL(nDKB))) * fft_input(i+1)
!       END IF
!    END DO

!!Hamming Window
!    DO  i = 0, nDKB
!       IF  (i == 0)   THEN
!          fft_input(i+1) = (0.54 + 0.46*cos((2*PI*i)/REAL(nDKB))) * fft_input(i+1)  
!       ELSE
!          fft_input(i+1) = (0.54 - 0.46*cos((2*PI*i)/REAL(nDKB))) * fft_input(i+1)
!       END IF
!    END DO  
    
!!Hann-Poisson Window
!    alpha_fft = 3
!    DO  i = 0, nDKB
!       IF  (i == 0)   THEN
!          fft_input(i+1) = (0.5 + 0.5*cos((2*PI*i)/(REAL(nDKB)/2)))  &
!                    * exp(-alpha_fft*(i/(REAL(nDKB)/2))) * fft_input(i+1)  
!       ELSE
!          fft_input(i+1) = (0.5 - 0.5*cos((2*PI*i)/(REAL(nDKB)/2)))   &
!                    * exp(-alpha_fft*(i/(REAL(nDKB)/2))) * fft_input(i+1)
!       END IF
!    END DO
    
    Status = DftiCreateDescriptor(My_Desc1_Handle, DFTI_SINGLE, DFTI_REAL, 1, NPTS)
    Status = DftiSetValue(My_Desc1_Handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
    Status = DftiSetValue(My_Desc1_Handle, DFTI_NUMBER_OF_USER_THREADS, 1)
    Status = DftiSetValue(My_Desc1_Handle, DFTI_FORWARD_SCALE, 1.0)
    Status = DftiCommitDescriptor(My_Desc1_Handle)
    Status = DftiComputeForward(My_Desc1_Handle, fft_input, fft_output)
    Status = DftiFreeDescriptor(My_Desc1_Handle)
    
!Calculate magnitude of DFT from real and imaginary parts.
!Real parts of DFT are in odd rows [(2*i)-1], Imaginary parts are in even rows (2*i) {for i = 1, NPTS}.
!Real parts of DFT are in even rows [(2*i)-1], Imaginary parts are in odd rows (2*i) {for i = 0, (NPTS-1)}.
    fft_mag = 0
    DO  i = 1, NPTS/2
       fft_mag(i) = sqrt( fft_output(2*i-1)**2 + fft_output(2*i)**2 )
      !phase (i) = atan( fft_output(2*i)/fft_output(2*i-1) )    ! arctan (imag/real)
    END DO
    
    DEALLOCATE(fft_input,fft_output)
    END SUBROUTINE fftransform
    
    
!*****************************************************************************************
!Subroutine for FFT without a window function

    SUBROUTINE fft_no_window(data_in, n_fft, signal_end, fft_mag, len_mult, dead_time)
    USE mkl_dfti
    USE trap_mod
    IMPLICIT NONE

   !INTEGER, PARAMETER :: NPTS0 = 65536            !Base number of datapoints  
   !INTEGER, PARAMETER :: FREQ_MULT = 16           !Multiplier for length of FFT     
   !INTEGER, PARAMETER :: NPTS = FREQ_MULT*NPTS0   !Number of datapoints  
   !REAL(4), PARAMETER :: PI = 3.14159265                      

    REAL(4), INTENT(IN)  :: data_in(NPTS)    !Input values 
    INTEGER, INTENT(IN)  :: n_fft            !Pointer for file section
    INTEGER, INTENT(IN)  :: signal_end       !Length of signal on which to perform DFT
    REAL(4), INTENT(OUT) :: fft_mag(NPTS/2)  !DFT magnitude values
   !REAL(4), INTENT(OUT) :: phase(NPTS/2)    !DFT phase values
    INTEGER, INTENT(IN)  :: len_mult         !Multiplier for length of each trapping
                                             !   event (multiple of 1 ms or 2,500 datapoints)
    INTEGER, INTENT(IN)  :: dead_time        !Number of points of dead time due to A250
                                             !   recovering from impulse noise of trap gate pulse

    TYPE(DFTI_DESCRIPTOR), POINTER :: My_Desc1_Handle
    INTEGER :: Status              !Sets DFT Parameters
    INTEGER :: istat               !Allocate status
    INTEGER :: i                   !Loop variables

    REAL(4), ALLOCATABLE :: fft_input(:)       !Input values after applying rolling average filter
    REAL(4), ALLOCATABLE :: fft_output(:)      !DFT output values

    ALLOCATE(fft_input(NPTS))
    ALLOCATE(fft_output(NPTS))

!Select section of data file and cut off first "dead_time" points to account for pulse
!generator impulse noise and possible packet of untrapped ions passing through the
!detection tube.
    DO  i = 1+(n_fft-1)*(2500*len_mult)+dead_time, (n_fft-1)*(2500*len_mult)+signal_end
       fft_input(i-(n_fft-1)*(2500*len_mult)-dead_time) = data_in(i)
    END DO

!Make sure everything after signal of interest is zero.
    DO  i = signal_end-(dead_time-1), NPTS
       fft_input(i) = 0
    END DO

    Status = DftiCreateDescriptor(My_Desc1_Handle, DFTI_SINGLE, DFTI_REAL, 1, NPTS)
    Status = DftiSetValue(My_Desc1_Handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
    Status = DftiSetValue(My_Desc1_Handle, DFTI_NUMBER_OF_USER_THREADS, 1)
    Status = DftiSetValue(My_Desc1_Handle, DFTI_FORWARD_SCALE, 1.0)
    Status = DftiCommitDescriptor(My_Desc1_Handle)
    Status = DftiComputeForward(My_Desc1_Handle, fft_input, fft_output)
    Status = DftiFreeDescriptor(My_Desc1_Handle)

!Calculate magnitude of DFT from real and imaginary parts.
!Real parts of DFT are in odd rows [(2*i)-1], Imaginary parts are in even rows (2*i) {for i = 1, NPTS}.
!Real parts of DFT are in even rows [(2*i)-1], Imaginary parts are in odd rows (2*i) {for i = 0, (NPTS-1)}.
    fft_mag = 0
    DO  i = 1, NPTS/2
       fft_mag(i) = sqrt( fft_output(2*i-1)**2 + fft_output(2*i)**2 )
      !phase (i) = atan( fft_output(2*i)/fft_output(2*i-1) )    ! arctan (imag/real)
    END DO

    DEALLOCATE(fft_input,fft_output)
    END SUBROUTINE fft_no_window


!*****************************************************************************************
!Subroutine for scanning across the section of the data file with a window of length "increment"
    SUBROUTINE fft_scan_start (data_in, n_fft, signal_start, increment, fft_mag, len_mult, freq_mult)
    USE mkl_dfti
    IMPLICIT NONE
    
    INTEGER, PARAMETER :: NPTS0 = 65536                 !Base number of datapoints  
    REAL(4), PARAMETER :: PI = 3.14159265

    REAL(4), INTENT(IN)  :: data_in(*)        !Raw input values (nominally 1048576)
    INTEGER, INTENT(IN)  :: n_fft             !Pointer for file section
    INTEGER, INTENT(IN)  :: signal_start      !Pointer for start of moving DFT
    INTEGER, INTENT(IN)  :: increment         !Length of DFT window 
    REAL(4), INTENT(OUT) :: fft_mag(*)        !DFT magnitude values (nominally 524288)
   !REAL(4), INTENT(OUT) :: phase(524288)     !DFT phase values (nominally 524288)
    INTEGER, INTENT(IN)  :: len_mult          !Multiplier for length of each trapping
                                              !   event (multiple of 1 ms or 2,500 datapoints)
    INTEGER, INTENT(IN)  :: freq_mult         ! Multiplier for length of FFT (base is 2^16)

    TYPE (DFTI_DESCRIPTOR), POINTER :: My_Desc1_Handle  
    INTEGER :: npts                                     !Number of datapoints 
    INTEGER :: Status                                   !Sets DFT Parameters
    INTEGER :: istat                                    !Allocate status               
    INTEGER :: i                                        !Loop variables
    REAL(4) :: sigma_fft                                !Sigma value for Gaussian window function 
   !REAL(4) :: alpha_fft                                !Alpha value for Hann-Poisson window function
    
    REAL(4), ALLOCATABLE :: fft_input (:)               ! Input values after applying rolling average filter 
    REAL(4), ALLOCATABLE :: fft_output (:)              ! DFT output values  
    
    npts = (freq_mult * NPTS0)
    ALLOCATE (fft_input(1:npts), fft_output(1:npts))
               
    DO  i = (1 + ((n_fft-1)*(2500*len_mult)) + signal_start), (((n_fft-1)*(2500*len_mult)) + signal_start + increment)
       fft_input(i - ((n_fft-1)*(2500*len_mult)) - signal_start) = data_in(i)
    END DO     
         
!Make sure everything after signal of interest is zero
    DO  i = (increment + 1), npts
       fft_input(i) = 0
    END DO                   
         
!Gaussian Window
    sigma_fft = 0.45
    DO  i = 0, increment - 1
       fft_input(i+1) = ( EXP( -0.5 * ((  (i - ((REAL(increment-1))/2)) / (sigma_fft * ((REAL(increment-1))/2))  )**2)   )  ) * fft_input(i+1)
    END DO
                  
!    !   Hann Window
!    DO  i = 0, (signal_end - (dead_time+1))
!        IF  (i == 0)   THEN
!            fft_input(i+1) = ( 0.5 + (0.5 * cos ((2*PI*i)/REAL(increment - 1))) ) * fft_input(i+1)  
!        ELSE
!            fft_input(i+1) = ( 0.5 - (0.5 * cos ((2*PI*i)/REAL(increment - 1))) ) * fft_input(i+1)
!        END IF
!    END DO

!    !   Hamming Window
!    DO  i = 0, increment - 1
!        IF  (i == 0)   THEN
!            fft_input(i+1) = ( 0.54 + (0.46 * cos ((2*PI*i)/REAL(increment - 1))) ) * fft_input(i+1)  
!        ELSE
!            fft_input(i+1) = ( 0.54 - (0.46 * cos ((2*PI*i)/REAL(increment - 1))) ) * fft_input(i+1)
!        END IF
!    END DO  
    
!    !   Hann-Poisson Window
!    alpha_fft = 3
!    DO  i = 0, increment - 1
!        IF  (i == 0)   THEN
!            fft_input(i+1) = (( 0.5 + (0.5 * cos ((2*PI*i)/(REAL(increment - 1)/2))) )) * ( exp(-alpha_fft * (i/(REAL(increment - 1)/2))) ) * fft_input(i+1)  
!        ELSE
!            fft_input(i+1) = ( 0.5 - (0.5 * cos ((2*PI*i)/(REAL(increment - 1)/2))) ) * ( exp(-alpha_fft * (i/(REAL(increment - 1)/2))) ) * fft_input(i+1)
!        END IF
!    END DO
    
    Status = DftiCreateDescriptor(My_Desc1_Handle, DFTI_SINGLE, DFTI_REAL, 1, npts)   
    Status = DftiSetValue(My_Desc1_Handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
    Status = DftiSetValue(My_Desc1_Handle, DFTI_NUMBER_OF_USER_THREADS, 1)
    Status = DftiSetValue(My_Desc1_Handle, DFTI_FORWARD_SCALE, 1.0)
    Status = DftiCommitDescriptor(My_Desc1_Handle)
    Status = DftiComputeForward(My_Desc1_Handle, fft_input, fft_output)
    Status = DftiFreeDescriptor(My_Desc1_Handle)
        
    fft_mag(1:npts/2) = 0
        
 !Calculate magnitude of DFT from real and imaginary parts
 !Real parts of DFT are in odd rows [(2*i)-1], Imaginary parts are in even rows (2*i) {for i = 1, npts}
 !Real parts of DFT are in even rows [(2*i)-1], Imaginary parts are in odd rows (2*i) {for i = 0, (npts-1)}
    DO  i = 1, (freq_mult * 32768)
       fft_mag(i) = sqrt( ((fft_output ((2*i)-1))**2) + ((fft_output (2*i))**2) )
      !phase (i) = atan( (fft_output (2*i)) / (fft_output ((2*i)-1)) )    ! arctan (imag/real)
    END DO
        
    DEALLOCATE  (fft_input, fft_output)
  
    END SUBROUTINE fft_scan_start
    
    
!*****************************************************************************************
    SUBROUTINE parabolic_interpolation (f2int, I2, I1, I3, Iv, fv)
    IMPLICIT NONE
    
    INTEGER :: f2int                             ! Integer read in local_max
    REAL(4) :: I1,I2,I3                          ! DFT peak maximum (I2) and adjacent points
    REAL(4) :: Iv                                ! Peak magnitude determined by parabolic interpolation
    REAL(4) :: fv                                ! "Pointer" for peak maximum determined by parabolic interpolation
    
    REAL(8) :: A, B, C                           ! Quadratic coefficients: I = A*f**2 + B*f + C. Real(8) substantially improves fit
    REAL(8) :: f1,f2,f3                          ! Real form of f2int as well as adjacent points
    INTEGER :: i                                 ! DO loop index
    REAL(4) :: freq_spacing                      ! Spacing between points used to fill in parabola
    REAL(4), DIMENSION(1:10001) :: freq_point    ! Frequencies of points to fill in parabola
    REAL(4), DIMENSION(1:10001) :: mag_point     ! Magnitudes of points to fill in parabola
      
!Three points are (local_max_less_one,DFT_mag_less_one), (local_max, DFT_mag),(local_max_plus_one,DFT_mag_plus_one),
!abbreviated to (f1,I1),(f2,I2),(f3,I3). Vertex of parabola at (interp_max,interp_mag), abbreviated to (fm,Im).
    f2 = REAL(f2int); f1 = f2 - 1.; f3 = f2 + 1.
    
    A = ((I1 - I3)*(f1 - f2) - (I1 - I2)*(f1 - f3))/((f1**2 - f3**2)*(f1 - f2) - (f1**2 - f2**2)*(f1 - f3))

    B = (I1 - I2 - A*(f1**2 - f2**2))/(f1 - f2)

    C = I1 - A*f1**2 - B*f1
    

 !Define points to fill in parabola
    freq_spacing = (f3 - f1)/10000.
    DO i = 1, 10001
       freq_point(i) = f1 + freq_spacing*(i - 1)
       mag_point(i) = A*freq_point(i)**2 + B*freq_point(i) + C
    END DO
    
!Scan through points to find vertex
    Iv = 0.
    DO i = 1, 10001
       IF (mag_point(i) > Iv) THEN
          Iv = mag_point(i)
          fv = freq_point(i)
       END IF
    END DO
    
!This analytical method of calculating fv and Iv may suffer from rounding errors (huge
!numbers subtracted to get small numbers).
   !fv = -B/(2*A)
   !Iv = -B**2/(4*A) + C
    
    END SUBROUTINE parabolic_interpolation
