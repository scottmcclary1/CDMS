    SUBROUTINE voltage_finder
    
    IMPLICIT NONE        
    CHARACTER (len = 50) :: filename                    ! Filename for r/w files 
    INTEGER :: nvals = 0                                ! Number of datapoints              
    INTEGER :: i, j                                     ! Loop variables   
    INTEGER :: unit                                     ! Unit number for r/w files
    INTEGER :: status                                   ! Error message when opening files
    INTEGER :: n_high                                   ! Number of datapoints with voltage > midpoint*1.05
    INTEGER :: n_low                                    ! Number of datapoints with voltage < midpoint*0.95
    
    REAL(4) :: max                                      ! Maximum value in array
    REAL(4) :: min                                      ! Minimum value in array
    REAL(4) :: midpoint                                 ! Midpoint between max and min values in array
    REAL(4) :: sum_high                                 ! Sum of all points > midpoint*1.05
    REAL(4) :: sum_low                                  ! Sum of all points < midpoint*0.95
    REAL(4) :: sum_of_squares_high                      ! Sum of (data point - avg_high)^2 
    REAL(4) :: sum_of_squares_low                       ! Sum of (data point - avg_low)^2 
    REAL(4) :: avg_high_sum                             ! Sum of avg_high for each file
    REAL(4) :: avg_low_sum                              ! Sum of avg_low for each file
    REAL(4) :: stdev_high_sq_sum                        ! Sum of square of std_dev_high for each file
    REAL(4) :: stdev_low_sq_sum                         ! Sum of square of std_dev_low for each file
    REAL(4) :: final_avg_high                           ! Avg of all high values for all 8 files
    REAL(4) :: final_avg_low                            ! Avg of all low values for all 8 files
    REAL(4) :: final_stdev_high                         ! Square root of stdev_high_sq_sum (cf. propagation of error)
    REAL(4) :: final_stdev_low                          ! Square root of stdev_low_sq_sum (cf. propagation of error)
    REAL(4) :: final_stdev                              ! Final standard deviation
    REAL(4) :: volts                                    ! Difference between final_avg_high and final_avg_low
    REAL(4) :: charge                                   ! Number of electrons corresponding to input voltage across 9 x 6.8 pF capacitors
    REAL(4) :: final_error_95                           ! 95% Confidence Interval
    REAL(4) :: final_error_99                           ! 99% Confidence Interval
    REAL(4) :: final_error_999                          ! 99.90% Confidence Interval
    
    REAL(4), DIMENSION (1:1048560) :: array_1           ! Input array 
    REAL(4), DIMENSION (1:8) :: avg_high_array          ! Array for avg value of datapoints > midpoint*1.05 for each of 8 files
    REAL(4), DIMENSION (1:8) :: avg_low_array           ! Array for avg value of low datapoints < midpoint*0.95 for each of 8 files         
    REAL(4), DIMENSION (1:8) :: stdev_high_array        ! Array for std dev of high datapoints for each of 8 files 
    REAL(4), DIMENSION (1:8) :: stdev_low_array         ! Array for std dev of low datapoints for each of 8 files         
    
    
    
        
    PRINT *, "You must have 8 files, numbered 1-8..."      
    PAUSE "Press Enter to continue..."
       
    unit = 42
    filename = "input_voltage.txt"
    CALL open_new_file (unit, filename)
    
    ! Select and read in file
    unit = 1     
    DO  i = 1, 8
        filename = CHAR(i+48)//'.txt'   
        
        
        PRINT *, "Now reading from ", filename
        OPEN (UNIT=unit, FILE=filename, STATUS="old", ACTION="read", IOSTAT=status) 
            
        nvals = 0
        DO
            READ (1, *, IOSTAT=status) array_1      ! Read file values and store in dummy variable "temp"
            IF (status /= 0) EXIT                   ! Stop reading values when end of file is reached
            nvals = nvals + 1                       ! Add count to number of datapoints
        END DO
        
        CLOSE (unit)    
      
        
        ! Assign midpoint between max value and min value
        ! Within +/- 5% of that midpoint value is a "dead zone"
        ! Assign all values > midpoint*1.05 to high value
        ! Assign all values < midpoint*0.95 to low value
        midpoint = (MAXVAL(array_1) + MINVAL(array_1))/2    
        sum_high = 0
        sum_low = 0
        n_high = 0
        n_low = 0
        DO  j = 1, 1048560
            IF  (array_1(j) > (midpoint*1.05))  THEN
                n_high = n_high + 1
                sum_high = sum_high + array_1(j)
            ELSE IF (array_1(j) < (midpoint*0.95))  THEN
                n_low = n_low + 1
                sum_low = sum_low + array_1(j)
            END IF
        END DO
        
        ! Find average of both high and low values
        avg_high_array(i) = sum_high/n_high
        avg_low_array(i) = sum_low/n_low       
        
        ! Find standard deviation of both high and low values    
        sum_of_squares_high = 0
        sum_of_squares_low = 0
        DO  j = 1, 1048560
            IF  (array_1(j) > (midpoint*1.05))  THEN
                sum_of_squares_high = sum_of_squares_high + ( (array_1(j) - avg_high_array(i))**2 )
            ELSE IF (array_1(j) < (midpoint*0.95))  THEN
                sum_of_squares_low = sum_of_squares_low + ( (array_1(j) - avg_low_array(i))**2 )
            END IF
        END DO    
        
        stdev_high_array(i) = SQRT( sum_of_squares_high / (n_high - 1) )
        stdev_low_array(i) = SQRT( sum_of_squares_low / (n_low - 1) )  
        
    END DO     
    
    ! Find average and standard deviation of all 8 files
    avg_high_sum = 0
    avg_low_sum = 0
    stdev_high_sq_sum = 0
    stdev_low_sq_sum = 0
    DO  i = 1, 8
        avg_high_sum = avg_high_sum + avg_high_array(i)
        avg_low_sum = avg_low_sum + avg_low_array(i)
        stdev_high_sq_sum = stdev_high_sq_sum + (stdev_high_array(i))**2
        stdev_low_sq_sum = stdev_low_sq_sum + (stdev_low_array(i))**2
    END DO
    
    final_avg_high = avg_high_sum/8
    final_avg_low = avg_low_sum/8
    final_stdev_high = SQRT(stdev_high_sq_sum)
    final_stdev_low = SQRT(stdev_low_sq_sum)    
        
    volts = final_avg_high - final_avg_low 
    
    charge = volts * 6.8E-12 / 9. / 1.6022E-19
    
    ! Calculate standard deviation and confidence intervals for all 8 files
    final_stdev = SQRT((final_stdev_high)**2 + (final_stdev_low)**2)
    final_error_95 = SQRT((1.96*final_stdev_high/SQRT(REAL(n_high)))**2 + (1.96*final_stdev_low/SQRT(REAL(n_low)))**2)
    final_error_99 = SQRT((2.576*final_stdev_high/SQRT(REAL(n_high)))**2 + (2.576*final_stdev_low/SQRT(REAL(n_low)))**2)
    final_error_999 = SQRT((3.291*final_stdev_high/SQRT(REAL(n_high)))**2 + (3.291*final_stdev_low/SQRT(REAL(n_low)))**2)
    
    PRINT *, "Charge (electrons):        ", charge
    PRINT *, "Input Voltage:             ", volts*1000, "mV"
    PRINT *, "95% Confidence Interval:   ", final_error_95*1000, "mV"
    PRINT *, "99% Confidence Interval:   ", final_error_99*1000, "mV"
    PRINT *, "99.90% Confidence Interval:", final_error_999*1000, "mV"
    PRINT *, "Standard Deviation:        ", final_stdev*1000, "mV"
    
    WRITE (42,888)  "Charge_(e)", "Input Voltage_(mV)", "95%_Conf._Int.", "99%_Conf._Int.", "99.90%_Conf._Int.", "Std._Dev."
    WRITE (42,889)  charge, volts*1000, final_error_95*1000, final_error_99*1000, final_error_999*1000, final_stdev*1000
    888 FORMAT (A10, TR5, A18, TR5, A14, TR5, A14, TR5, A17, TR5, A9)
    889 FORMAT (F9.2, TR7, F7.5, TR16, E12.5, TR8, E12.5, TR8, E12.5, TR8, E12.5)
               
    ! Close all output files and keep them (i.e., don't delete them)
    CLOSE (42, STATUS = "KEEP")         
            
    END SUBROUTINE voltage_finder