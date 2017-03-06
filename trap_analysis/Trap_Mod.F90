      MODULE  trap_mod

      REAL(4), PARAMETER :: PI = 3.14159265                            

      INTEGER, PARAMETER :: NPTS0 =  65536            !base number of datapoints (2^16)     
      INTEGER, PARAMETER :: FREQ_MULT0 = 16           !multiplier for length of FFT
      INTEGER, PARAMETER :: NPTS = FREQ_MULT0*NPTS0   !number of datapoints for FFT
      INTEGER, PARAMETER :: NPTS10 = 1000000          !acutal number of datapoints

      END MODULE  trap_mod
