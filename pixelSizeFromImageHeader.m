


% to extract pixel size form the header of raw image ~mh973

%% objective selection
% uncomment the conversionFactor of your objective
%     if (matches(objective,"Fetcho 4x")){
%         conversionFactor=392.9777;
%         } else if (matches(objective,"Fetcho 25x")){
%         conversionFactor=70.1628;
%         } else if (matches(objective,"Setup 1: 20x Olympus")){
%         conversionFactor=242.6515;
%         }  else if (matches(objective,"Setup 1: 20x Zeiss")){
%         conversionFactor=230.1013;
%         }  else if (matches(objective,"Setup 1: 20x Zeiss Coverslip Corrected")){
%         conversionFactor=229.7129;
%         } else if (matches(objective,"Setup 2: 20x Olympus")){
%         conversionFactor=248.5030;
%         }  else if (matches(objective,"Setup 2: 20x Zeiss")){
%         conversionFactor=234.3450;
%         }  else if (matches(objective,"Setup 2: 20x Zeiss Coverslip Corrected")){
%         conversionFactor=233.4116;
%         }  else if (matches(objective,"Setup 3: 20x Zeiss Coverslip Corrected")){
%         conversionFactor=146.180065;
%         }  else if (matches(objective,"Setup 3: 20x Zeiss")){
%          conversionFactor=150.263401;
%             } else if (matches(objective,"Setup 3: 40x Zeiss Oil Immersion")){
%         conversionFactor=36.8060;
%         } else if (matches(objective,"Setup 3: 4x Zeiss")){
%          conversionFactor=651.3437;
%         } else if (matches(objective,"Setup 3: 20x Olympus")){
%         conversionFactor=143.1174;
%         }  else if (matches(objective,"Setup 3: 10x Olympus")){
%         conversionFactor=157.0023;
%         }  else if (matches(objective,"Setup 3: 25x Olympus")){
         conversionFactor=125.2536;
%         }

%% Calculation 

[FileName,PathName] = uigetfile('*.*','Select the image file');
imInfo=imfinfo([PathName,'/',FileName]);
eval(imInfo(1).ImageDescription);
SI_version=state.software.version;


if (SI_version == 3.6)
    safast=2*state.acq.scanAmplitudeX;
    saslow=2*state.acq.scanAmplitudeY;
    voltsPerDegree=state.init.opticalDegreesConversion;
elseif (SI_version == 3.7)
    safast=state.acq.scanAngularRangeFast;
    saslow=state.acq.scanAngularRangeSlow;
    voltsPerDegree=state.init.voltsPerOpticalDegree;
    fastX=state.acq.fastScanningX;
else 
    safast=state.init.scanAngularRangeReferenceFast;
    saslow=state.init.scanAngularRangeReferenceSlow;
    voltsPerDegree=state.init.voltsPerOpticalDegree;
    fastX=state.acq.fastScanningX;
end

zoom=state.acq.zoomFactor;

voltsf=safast/zoom*voltsPerDegree/2;
FOVf=voltsf*conversionFactor;
voltss=saslow/zoom*voltsPerDegree/2;
FOVs=voltss*conversionFactor;

getHeight=imInfo(1).Height;
getWidth=imInfo(1).Width;

if (fastX==1) 
	ysize=FOVs/getHeight
	xsize=FOVf/getWidth
else 
	ysize=FOVf/getHeight
	xsize=FOVs/getWidth
end

