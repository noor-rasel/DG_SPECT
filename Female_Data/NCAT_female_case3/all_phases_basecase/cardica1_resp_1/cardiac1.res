


              SIMIND Monte Carlo Simulation Program    V6.2.1
------------------------------------------------------------------------------
 Phantom(S): h2o       Crystal...: nai       InputFile.: cardiac1          
 Phantom(B): bone      BackScatt.: pmt       OutputFile: cardiac1          
 Collimator: pb_sb2    SourceRout: smap      SourceFile: cardiac1_act_1    
 Cover.....: al        ScoreRout.: scattwin  DensityMap: cardiac1_atn_1    
------------------------------------------------------------------------------
 PhotonEnergy.......: 140          pr-lehr  PhotonsPerProj....: 5309830        
 EnergyResolution...: 10           SPECT    Activity..........: 1              
 MaxScatterOrder....: 3            Random   DetectorLenght....: 25             
 DetectorWidth......: 25.6         Phantom  DetectorHeight....: 0.9525         
 UpperEneWindowTresh: 154          Resolut  Distance to det...: 25             
 LowerEneWindowTresh: 126          IntFile  ShiftSource X.....: 0              
 PixelSize  I.......: 0.8                   ShiftSource Y.....: 0              
 PixelSize  J.......: 0.8                   ShiftSource Z.....: 0              
 HalfLength S.......: 25.6                  HalfLength P......: 25.6           
 HalfWidth  S.......: 0.1                   HalfWidth  P......: 0.1            
 HalfHeight S.......: 0.1                   HalfHeight P......: 0.1            
 SourceType.........: XcatBinMap            PhantomType.......: XcatBinMap   
------------------------------------------------------------------------------
 GENERAL DATA
 keV/channel........: 1                     CutoffEnergy......: 0              
 Photons/Bq.........: 0.879                 StartingAngle.....: 0              
 CameraOffset X.....: 0                     CoverThickness....: 0              
 CameraOffset Y.....: 0                     BackscatterThickn.: 0              
 MatrixSize I.......: 64                    IntrinsicResolut..: 0.36           
 MatrixSize J.......: 64                    AcceptanceAngle...: 3.31774        
 Emission type......: 2                     Initial Weight....: 0.16554        
 NN ScalingFactor...: 5                     Energy Channels...: 512            
------------------------------------------------------------------------------
 SPECT DATA
 RotationMode.......: -360                  Nr of Projections.: 64             
 RotationAngle......: 5.625                 Projection.[start]: 1              
 Orbital fraction...: 1                     Projection...[end]: 64             
------------------------------------------------------------------------------
 COLLIMATOR DATA FOR ROUTINE: Analytical          
 CollimatorCode.....: pr-lehr               CollimatorType....: Parallel 
 HoleSize X.........: 0.14                  Distance X........: 0.018          
 HoleSize Y.........: 0.16166               Distance Y........: 0.09642        
 CenterShift X......: 0.079                 X-Ray flag........: F              
 CenterShift Y......: 0.13683               CollimThickness...: 2.7            
 HoleShape..........: Hexagonal             Space Coll2Det....: 0              
 CollDepValue [57]..: 0                     CollDepValue [58].: 0              
 CollDepValue [59]..: 0                     CollDepValue [60].: 0              
------------------------------------------------------------------------------
 IMAGE-BASED PHANTOM DATA
 RotationCentre.....:  33, 33               Bone definition...: 1170           
 CT-Pixel size......: 0.8                   Slice thickness...: 0.8            
 StartImage.........: 1                     No of CT-Images...: 64             
 MatrixSize I.......: 64                    CTmapOrientation..: 0              
 MatrixSize J.......: 64                    StepSize..........: 0.6            
 CenterPoint I......: 33                    ShiftPhantom X....: 0              
 CenterPoint J......: 33                    ShiftPhantom Y....: 0              
 CenterPoint K......: 33                    ShiftPhantom Z....: 0              
------------------------------------------------------------------------------
  Scattwin results: Window file: scattwin.win        
  
  Win  WinAdded  Range(keV)   ScaleFactor
   1       0    126.0 - 154.0    1.00
  
  Win    Total    Scatter   Primary  S/P-Ratio S/T Ratio  Cps/MBq
   1   0.699E+04 0.000E+00 0.699E+04 0.000E+00 0.000E+00 0.109E+03
  
  Win  Geo(Air)  Pen(Air)  Sca(Air)  Geo(Tot)  Pen(Tot)  Sca(Tot)
   1   100.00%     0.00%     0.00%   100.00%     0.00%     0.00%
  
 Header file........: cardiac1.h00                            
 Image file.........: cardiac1.a00                            
------------------------------------------------------------------------------
 INTERACTIONS IN THE CRYSTAL
 MaxValue spectrum..: 450.1          
 MaxValue projection: 0.2324         
 CountRate spectrum.: 110.6          
 CountRate E-Window.: 109.3          
------------------------------------------------------------------------------
 PHOTONS AFTER COLLIMATOR AND WITHIN ENER-WIN
 Geometric..........:   0.00 %         100.00 %
 Penetration........:   0.00 %           0.00 %
 Scatter in collim..:   0.00 %           0.00 %
 X-rays in collim...:   0.00 %           0.00 %
------------------------------------------------------------------------------
 CALCULATED DETECTOR PARAMETERS
 Efficiency E-Window: 0.8652          3.26 % (1SD)
 Efficiency spectrum: 0.8755         
 Sensitivity Cps/MBq: 109.2543       
 Sensitivity Cpm/uCi: 242.5445       
------------------------------------------------------------------------------
 Simulation started.: 2024:02:13 22:41:36
 Simulation stopped.: 2024:02:13 22:45:03
 Elapsed time.......: 0 h, 3 m and 27 s
 DetectorHits.......: 4926492        
 DetectorHits/CPUsec: 28199          
------------------------------------------------------------------------------
 SIMIND built 2021:07:31 with INTEL Linux  compiler
 Random number generator: Intel RAN
 Comment:EMISSION
 Energy resolution as function of 1/sqrt(E)
 Inifile: simind.ini
 Command: cardiac1 /SF:1 /NN:5 /fa:1
