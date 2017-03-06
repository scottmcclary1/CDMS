    SUBROUTINE harmonic_analysis(channel, freq_mult, fft_magnitude, fund_freq_win_1, fund_mag_win_1, temp_signal_end, harm_freq, harm_mag, t)
    ! See versions of code (Mar_10_2014 to Apr_01_2014) for different way of performing this analysis. Current method more flexible and streamlined
    
    IMPLICIT NONE
    INTEGER :: channel                  ! Channel read in so output only written for analyzing single trapping events
    INTEGER :: j, i
    INTEGER :: t                        ! Flag used in user-selected options for channels 3 and 5
    INTEGER :: freq_mult                ! Multiplier for length of FFT (base is 2^16)
    REAL(4) :: mag_cutoff    ! FFT magnitude analysis cutoff value
    REAL(4), DIMENSION (1:524288) :: fft_magnitude      ! FFT magnitude output
    REAL(4) :: fund_freq_win_1          ! Fundamental frequency as determined by the first window
    REAL(4) :: fund_mag_win_1           ! Magnitude of fundamental frequency as determined by the first window
    REAL(4) :: mag_max                  ! Maximum value of fft_magnitude in a desired range
    INTEGER :: harm1_flag               ! Used to terminate loop after REAL(harm_freq_index(1)) found
    INTEGER, DIMENSION(1:10) :: harm_freq_index
                                        ! Location in fft_magnitude of first 10 harmonics
    REAL(4), DIMENSION(1:10) :: harm_freq
                                        ! Frequency of first 10 harmonics
    REAL(4), DIMENSION(1:10) :: harm_mag
                                        ! Magnitude of first 10 harmonics
    REAL(4), DIMENSION(1:10) :: harm_index_low, harm_index_high
                                        ! Low and high ranges of locations in fft_magnitude to look for harmonics
    REAL(4) :: freq_ratio               ! Ratio of jth peak to fundamental
    INTEGER :: temp_signal_end          ! Temporary pointer for end of signal
    
    ! Variables for parabolic_interpolation
    INTEGER :: local_max                ! Pointer for local maximum
    REAL(4) :: DFT_mag, DFT_mag_less_one, DFT_mag_plus_one 
                                        ! Values of fft_magnitude corresponding to max of peak and adjacent peaks
    REAL(4) :: interp_max               ! Position of peak max determined by parabolic interpolation
    REAL(4) :: interp_mag               ! Magnitude of peak max determined by parabolic interpolation
    
    

    ! Find fundamental frequency
    harm_index_low(1) = FLOOR(0.95*fund_freq_win_1*2**20*freq_mult/40.E6 + 1)
    harm_index_high(1) = CEILING(1.05*fund_freq_win_1*2**20*freq_mult/40.E6 + 1)
    mag_max = 0
    DO j = harm_index_low(1) + 1, harm_index_high(1)
        IF (fft_magnitude(j) > mag_max) THEN
            mag_max = fft_magnitude(j)
            harm_freq_index(1) = j
        END IF
    END DO
    local_max = harm_freq_index(1)      ! Assign value in array as scalar so it can be used with parabolic_interpolation
    DFT_mag = fft_magnitude(local_max)
    DFT_mag_less_one = fft_magnitude(local_max - 1)
    DFT_mag_plus_one = fft_magnitude(local_max + 1)
    
    CALL parabolic_interpolation (local_max,DFT_mag,DFT_mag_less_one,DFT_mag_plus_one,interp_mag,interp_max)
    
    harm_freq(1) = (interp_max - 1.)*40.E6/(2**20*freq_mult)
    harm_mag(1) = interp_mag
!    harm_freq(1) = (REAL(harm_freq_index(1)) - 1.)*40.E6/(2**20*freq_mult)
!    harm_mag(1) = fft_magnitude(harm_freq_index(1))
    
    ! Search within 5% of where higher-order harmonics should be; define maximum in those ranges as value of the harmonic
    harm_index_low(2) = FLOOR(1.95*REAL(harm_freq_index(1))); harm_index_high(2) = CEILING(2.05*REAL(harm_freq_index(1)))
    harm_index_low(3) = FLOOR(2.95*REAL(harm_freq_index(1))); harm_index_high(3) = CEILING(3.05*REAL(harm_freq_index(1)))
    harm_index_low(4) = FLOOR(3.95*REAL(harm_freq_index(1))); harm_index_high(4) = CEILING(4.05*REAL(harm_freq_index(1)))
    harm_index_low(5) = FLOOR(4.95*REAL(harm_freq_index(1))); harm_index_high(5) = CEILING(5.05*REAL(harm_freq_index(1)))
    harm_index_low(6) = FLOOR(5.95*REAL(harm_freq_index(1))); harm_index_high(6) = CEILING(6.05*REAL(harm_freq_index(1)))
    harm_index_low(7) = FLOOR(6.95*REAL(harm_freq_index(1))); harm_index_high(7) = CEILING(7.05*REAL(harm_freq_index(1)))
    harm_index_low(8) = FLOOR(7.95*REAL(harm_freq_index(1))); harm_index_high(8) = CEILING(8.05*REAL(harm_freq_index(1)))
    harm_index_low(9) = FLOOR(8.95*REAL(harm_freq_index(1))); harm_index_high(9) = CEILING(9.05*REAL(harm_freq_index(1)))
    harm_index_low(10) = FLOOR(9.95*REAL(harm_freq_index(1))); harm_index_high(10) = CEILING(10.05*REAL(harm_freq_index(1)))
    
    DO j = 2, 10 
        mag_max = 0
        DO i = harm_index_low(j) + 1, harm_index_high(j)
            IF (fft_magnitude(i) > mag_max) THEN
                mag_max = fft_magnitude(i)
                harm_freq_index(j) = i
            END IF
        END DO
        local_max = harm_freq_index(j)      ! Assign value in array as scalar so it can be used with parabolic_interpolation
        DFT_mag = fft_magnitude(local_max)
        DFT_mag_less_one = fft_magnitude(local_max - 1)
        DFT_mag_plus_one = fft_magnitude(local_max + 1)
        
        CALL parabolic_interpolation (local_max,DFT_mag,DFT_mag_less_one,DFT_mag_plus_one,interp_mag,interp_max)
        
        harm_freq(j) = (interp_max - 1.)*40.E6/(2**20*freq_mult)
        harm_mag(j) = interp_mag
!        harm_freq(j) = (REAL(harm_freq_index(j)) - 1)*40.E6/(2**20*freq_mult)
!        harm_mag(j) = MAXVAL(fft_magnitude(harm_index_low(j):harm_index_high(j)))
    END DO
    
    IF ((channel == 3 .OR. channel == 5) .AND. t == 1) THEN
        WRITE (108, '(I16,20F15.1)') &
            temp_signal_end, harm_freq(1), harm_mag(1), harm_freq(2), harm_mag(2), harm_freq(3), harm_mag(3), harm_freq(4), harm_mag(4), &
            harm_freq(5), harm_mag(5), harm_freq(6), harm_mag(6), harm_freq(7), harm_mag(7), harm_freq(8), harm_mag(8), &
            harm_freq(9), harm_mag(9), harm_freq(10), harm_mag(10)
    END IF

    END SUBROUTINE harmonic_analysis