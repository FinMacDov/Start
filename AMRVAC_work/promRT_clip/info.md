-----------------------------------------------------------------------------
 -----------------------------------------------------------------------------
 |         __  __ ____ ___        _    __  __ ______     ___    ____         |
 |        |  \/  |  _ \_ _|      / \  |  \/  |  _ \ \   / / \  / ___|        |
 |        | |\/| | |_) | |_____ / _ \ | |\/| | |_) \ \ / / _ \| |            |
 |        | |  | |  __/| |_____/ ___ \| |  | |  _ < \ V / ___ \ |___         |
 |        |_|  |_|_|  |___|   /_/   \_\_|  |_|_| \_\ \_/_/   \_\____|        |
 -----------------------------------------------------------------------------
 -----------------------------------------------------------------------------
 Warning in ReadParameters: No save condition for file            3
 Warning in ReadParameters: No save condition for file            4
 Warning in ReadParameters: No save condition for file            5
 typelimited to predictor for RK
 Using          150  cells in dimension            1
 level one dx(           1 )=   6.6666666666666671E-003
 Using          150  cells in dimension            2
 level one dx(           2 )=   6.6666666666666671E-003
 Error estimation is Lohner's scheme
 Reading from inifile: promRT25D.par
 snapshotini         :           -1
 slicenext           :            0
 collapsenext        :            0
 Filenameini         : data
 Converting?         :  F
                                                                 
 number of modes=   50.000000000000000     
 2.5D or 3D MHD Rayleigh-Taylor instability in prominence
 -------------------------------------------------------------------------------
 Startup phase took :             1.036 sec
 -------------------------------------------------------------------------------
 BCs before Advance took :        0.002 sec
 -----------------------------------------------------------------------------
 Saving visual data. Coordinate directions and variable names are:
           1 X         
           2 Y         
           3 rho       rho       
           4 v1        v1        
           5 v2        v2        
           6 v3        v3        
           7 p         p         
           8 b1        b1        
           9 b2        b2        
          10 b3        b3        
          11 psi       psi       
          12 tr        tr        
          13 T         T         
          14 beta      beta      
 time =   0.0000000000000000     
 -----------------------------------------------------------------------------
 Total timeloop took        :     2425.240 sec
 Time spent on Regrid+Update:      121.843 sec
                  Percentage:         5.02 %
 Time spent on IO in loop   :       31.812 sec
                  Percentage:         1.31 %
 Time spent on BC           :      527.253 sec
                  Percentage:        21.74 %
 Time spent on run          :     2393.428 sec
 Total time spent on IO     :       31.812 sec
 Total timeintegration took :     2425.243 sec
 -------------------------------------------------------------------------------
 Finished AMRVAC in :          2426.279 sec
 -------------------------------------------------------------------------------

This simulation is ran in CGS units. It has a slightly tilted magnetic feild through z-x axis.

&savelist
        itsave(1,1)=0
        itsave(1,2)=0
        dtsave(1)     = 0.01
        dtsave(2)     = 0.10 




