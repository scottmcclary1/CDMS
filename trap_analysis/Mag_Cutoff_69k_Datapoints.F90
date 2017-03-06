    SUBROUTINE mag_cutoff_69k (window_length, frequency, mag_cutoff)
    
    IMPLICIT NONE
    
    REAL(4), INTENT(OUT) :: mag_cutoff                         ! FFT magnitude analysis cutoff value
    REAL(4), INTENT(IN) :: frequency
    REAL(4) :: a
    REAL(4) :: b
    REAL(4) :: c

    INTEGER, INTENT(IN) :: window_length 


    IF  (window_length == 1000)    THEN    
        a = 16376.49992

        b =	20130.15946

        c = 7560.51439

    ELSE IF (window_length == 2000)    THEN    
        a = 21733.70096

        b =	27463.22642

        c = 13674.96471

    ELSE IF (window_length == 3000)    THEN    
        a = 26521.41687

        b =	27496.86829

        c = 14021.6337

    ELSE IF (window_length == 4000)    THEN    
        a = 30542.70291

        b =	28015.56381

        c = 14551.80623

    ELSE IF (window_length == 5000)    THEN    
        a = 34215.30689

        b =	27899.50731

        c = 14483.62682

    ELSE IF (window_length == 6000)    THEN    
        a = 37610.61404

        b =	27430.8803

        c = 14114.38722

    ELSE IF (window_length == 7000)    THEN    
        a = 40692.56123

        b =	27203.11112

        c = 13970.63969

    ELSE IF (window_length == 8000)    THEN    
        a = 43470.05296

        b =	27285.53572

        c = 14082.85569

    ELSE IF (window_length == 9000)    THEN    
        a = 46061.2553

        b =	27394.40009

        c = 14184.74633

    ELSE IF (window_length == 10000)    THEN    
        a = 48520.55432

        b =	27468.60865

        c = 14248.96764

    ELSE IF (window_length == 11000)    THEN    
        a = 50921.42518

        b =	27405.7946

        c = 14208.78871

    ELSE IF (window_length == 12000)    THEN    
        a = 53226.6241

        b =	27317.46407

        c = 14148.01335

    ELSE IF (window_length == 13000)    THEN    
        a = 55412.34924

        b =	27299.97355

        c = 14151.24152

    ELSE IF (window_length == 14000)    THEN    
        a = 57527.71685

        b =	27251.87302

        c = 14122.84288

    ELSE IF (window_length == 15000)    THEN    
        a = 59551.95838

        b =	27241.1884

        c = 14115.24916

    ELSE IF (window_length == 16000)    THEN    
        a = 61490.60705

        b =	27272.87833

        c = 14136.09181

    ELSE IF (window_length == 17000)    THEN    
        a = 63364.99987

        b =	27311.4448

        c = 14159.88022

    ELSE IF (window_length == 18000)    THEN    
        a = 65194.69612

        b =	27333.58914

        c = 14173.81083

    ELSE IF (window_length == 19000)    THEN    
        a = 66966.5353

        b =	27370.76503

        c = 14199.06521

    ELSE IF (window_length == 20000)    THEN    
        a = 68702.38385

        b =	27384.69151

        c = 14203.94927

    ELSE IF (window_length == 21000)    THEN    
        a = 70400.45748

        b =	27384.28898

        c = 14196.70437

    ELSE IF (window_length == 22000)    THEN    
        a = 72054.27949

        b =	27392.19602

        c = 14200.12373

    ELSE IF (window_length == 23000)    THEN    
        a = 73673.24823

        b =	27397.97042

        c = 14203.89011

    ELSE IF (window_length == 24000)    THEN    
        a = 75252.92886

        b =	27416.10918

        c = 14219.55581

    ELSE IF (window_length == 25000)    THEN    
        a = 76796.0328

        b =	27442.12192

        c = 14242.28171

    ELSE IF (window_length == 26000)    THEN    
        a = 78306.76038

        b =	27469.90422

        c = 14265.86834

    ELSE IF (window_length == 27000)    THEN    
        a = 79790.48977

        b =	27490.21444

        c = 14284.43823

    ELSE IF (window_length == 28000)    THEN    
        a = 81244.4152

        b =	27508.44185

        c = 14301.46653

    ELSE IF (window_length == 29000)    THEN    
        a = 82665.86277

        b =	27534.67255

        c = 14325.56821

    ELSE IF (window_length == 30000)    THEN    
        a = 84072.03155

        b =	27543.34041

        c = 14333.43549

    ELSE IF (window_length == 31000)    THEN    
        a = 85453.51049

        b =	27549.82482

        c = 14339.13146

    ELSE IF (window_length == 32000)    THEN    
        a = 86819.53089

        b =	27545.37526

        c = 14336.33788

    ELSE IF (window_length == 33000)    THEN    
        a = 88163.44427

        b =	27542.38408

        c = 14334.92792

    ELSE IF (window_length == 34000)    THEN    
        a = 89495.91833

        b =	27527.4332

        c = 14323.63239

    ELSE IF (window_length == 35000)    THEN    
        a = 90809.30035

        b =	27512.91157

        c = 14311.18188

    ELSE IF (window_length == 36000)    THEN    
        a = 92105.85041

        b =	27496.42743

        c = 14296.22917

    ELSE IF (window_length == 37000)    THEN    
        a = 93376.22565

        b =	27492.91082

        c = 14291.73405

    ELSE IF (window_length == 38000)    THEN    
        a = 94631.72444

        b =	27487.27296

        c = 14286.61562

    ELSE IF (window_length == 39000)    THEN    
        a = 95863.88825

        b =	27491.64436

        c = 14288.51426

    ELSE IF (window_length == 40000)    THEN    
        a = 97079.27935

        b =	27496.844

        c = 14290.10323

    ELSE IF (window_length == 41000)    THEN    
        a = 98278.9193

        b =	27501.55684

        c = 14290.47143

    ELSE IF (window_length == 42000)    THEN    
        a = 99457.45963

        b =	27515.55318

        c = 14299.6273

    ELSE IF (window_length == 43000)    THEN    
        a = 100617.4198

        b =	27535.61371

        c = 14314.6942

    ELSE IF (window_length == 44000)    THEN    
        a = 101770.4541

        b =	27548.58511

        c = 14325.31082

    ELSE IF (window_length == 45000)    THEN    
        a = 102911.0686

        b =	27559.37389

        c = 14332.03625

    ELSE IF (window_length == 46000)    THEN    
        a = 104046.7424

        b =	27561.76885

        c = 14333.15087

    ELSE IF (window_length == 47000)    THEN    
        a = 105161.3113

        b =	27575.95343

        c = 14344.7729

    ELSE IF (window_length == 48000)    THEN    
        a = 106269.8259

        b =	27582.74492

        c = 14352.35235

    ELSE IF (window_length == 49000)    THEN    
        a = 107366.3703

        b =	27593.2337

        c = 14363.77354

    ELSE IF (window_length == 50000)    THEN    
        a = 108455.2586

        b =	27598.43204

        c = 14368.87153

    ELSE IF (window_length == 51000)    THEN    
        a = 109540.3693

        b =	27597.69385

        c = 14368.74179

    ELSE IF (window_length == 52000)    THEN    
        a = 110611.747

        b =	27600.33165

        c = 14372.69515

    ELSE IF (window_length == 53000)    THEN    
        a = 111674.3125

        b =	27603.34875

        c = 14379.22615

    ELSE IF (window_length == 54000)    THEN    
        a = 112728.9079

        b =	27603.91287

        c = 14383.79174

    ELSE IF (window_length == 55000)    THEN    
        a = 113778.4676

        b =	27601.29063

        c = 14382.7484

    ELSE IF (window_length == 56000)    THEN    
        a = 114821.627

        b =	27594.77878

        c = 14376.8947

    ELSE IF (window_length == 57000)    THEN    
        a = 115846.3685

        b =	27600.80883

        c = 14384.27458

    ELSE IF (window_length == 58000)    THEN    
        a = 116869.0753

        b =	27598.81714

        c = 14388.47827

    ELSE IF (window_length == 59000)    THEN    
        a = 117880.3806

        b =	27600.94958

        c = 14396.73044

    ELSE IF (window_length == 60000)    THEN    
        a = 118888.6165

        b =	27597.3646

        c = 14393.76593

    ELSE IF (window_length == 61000)    THEN    
        a = 119899.3254

        b =	27583.74315

        c = 14377.14249      
        
    ELSE IF (window_length == 62000)    THEN    
        a = 120888.5199

        b =	27584.4658

        c = 14377.65121

    ELSE IF (window_length == 63000)    THEN    
        a = 121869.1579

        b =	27584.72945

        c = 14386.1281

    ELSE IF (window_length == 64000)    THEN    
        a = 122835.5931

        b =	27592.30895

        c = 14404.96353

    ELSE IF (window_length == 65000)    THEN    
        a = 123815.0971

        b =	27579.41499

        c = 14383.9582

    ELSE IF (window_length == 66000)    THEN    
        a = 124848.0302

        b =	27506.0329

        c = 14269.47669

    ELSE IF (window_length == 67000)    THEN    
        a = 125923.9135

        b =	27385.02866

        c = 14056.80322

    ELSE IF (window_length == 68000)    THEN    
        a = 127209.9025

        b =	27082.10206

        c = 13497.65244

    ELSE IF (window_length == 69000)    THEN    
        a = 124248.6133

        b =	34877.67072

        c = 13305.36363

    END IF
          
    
    mag_cutoff = (a * EXP(b / (frequency + c))) + 20000
    
    
    END SUBROUTINE  mag_cutoff_69k