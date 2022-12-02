function varargout = arbitraryLinescanPreprocess(varargin)
% FUNCTION arbitraryLinescanPreprocess()
% FUNCTION [results, tifFilename] = arbitraryLinescanPreprocess()
% FUNCTION [results, tifFilename] = arbitraryLinescanPreprocess(setup,objective,vesselChannel,getVesselWidthSnap,getVesselWidthLine,nearLineMicrons,nearDiamMicrons)
%
% This functions pre-processes raw data from a scanimage aribtrary 
% linescan (with one or more ROIs) and outputs the linescan data to a 
% tif stack that can be fed into scripts for measuring blood flow speed 
% (i.e. the sn-lab/Blood-flow repo). Aditionally, this script pulls the 
% ms/line from the metadata, calculates the um/pixel conversion factor, 
% and can estimate the vessel width using either a single-frame snapshot 
% taken alongside the linescan or using a line ROI drawn sideways across
% the vessel.
% 
% DEPENDENCIES:
% saveastiff() (included in Blood-Flow repo: https://github.com/sn-lab/Blood-flow) 
%              (originally from normcorre: https://github.com/flatironinstitute/NoRMCorre)
% align_data() (included in Blood-Flow repo: https://github.com/sn-lab/Blood-flow)
% Pixel-to-Micron (repo: https://github.com/sn-lab/Pixel-to-Micron)
% Image Processing Toolbox
% Statistics and Machine Learning Toolbox
% Curve Fitting Toolbox
%
% OPTIONAL INPUTS:
% setup: the imaging setup used (e.g. 3)
% objective: the objective used, which must be an exact match to one listed
%             in the Pixel-to-Micron repo (e.g. 'Olympus 20x')
% vesselChannel: imaging channel (among the saved channels only) with 
%                vasculature label (e.g. 4)
% getVesselWidthSnap: whether to analyze a snapshot image of the vessel to
%                     calculate the vessel width (e.g. 'yes' or 'no')
% getVesselWidthLine: whether to analyze a line ROI across the vessel to
%                     calculate the vessel width (e.g. 'yes' or 'no')
% nearLineMicrons: maximum orthogonal distance the scan position can be
%                  from the blood flow ROI to be included
% nearDiamMicrons: maximum orthogonal distance the scan position can be
%                  from the diameter ROI to be included
%(if these 7 inputs aren't given, a dialog box will ask for them)
%
% ADDITIONAL REQUIRED FILES:
% When calling this function, you will be asked to navigate to arbitrary
% linescan data (a .pmt.dat file). An associated meta.txt file and 
% scnnr.dat file should be in the same folder, with the same base filename 
% as the .pmt.dat file. Additionally, if you are opting to calculate the 
% vessel width using a snaphot, before or after taking the linescan you
% should have acquired a single frame snapshot of the field of view used
% to take the linescan: this is what is used by this script to calculate
% the diameter of the vessel being scanned. If this function finds a lone
% .tif file in the same folder as the linescan file, it will assume it to
% be this snapshot; otherwise, you will be asked to navigate to the
% associated snapshot.
%   TLDR: 
%   Required files to run this function:
%   ...pmt.dat (from linescan)
%   ...meta.txt (from linescan)
%   ...scnnr.dat (from linescan)
%   ....tif (from framescan snapshot) -- optional, for getVesselWidthSnap
%
% OPTIONAL OUTPUTS:
%   results: a matlab struct containing the following values:
%                   1. bloodFlowROIs
%                   2. diameterROIs
%                   3. vesselDiameterInMicronsSnap
%                   4. vesselDiameterInMicronsLine
%                   5. linescanMicronsPerPixel
%                   6. linescanMsPerLine
%   tifFilename: filename of the .tif file containing pre-processed
%                linescan data, to be fed into the Blood-flow pipeline
% 
% In addition, 3-5 types of files are saved by this script:
%   1. A .tif stack of the blood flow linescan data, compatible with flow
%       measurement scripts extractVelTiffShared.m and velocity_from_tif.m
%       -- the filename of this tif stack is returned from the function
%   2  A .mat file of the results (e.g. microns/pixel, vessel width, 
%       ms/line) -- these results are also shown in a dialogbox when
%       finished and are returned from the function in the "results" struct
%   3. A .tif of plots of the ROIs and overlayed scan positions (to verify
%       that this function calculated the lines/positions correctly)
%   4. (optional) A .tif image of the linescan vessel diameter estimation 
%   5. (optional) A .tif image of the snapshot vessel diameter estimation 
%
% 9/19/22, Matt Isaacson, mdi22@cornell.edu


%% settings
%default parameter settings
setup = 3;
objective = 'Olympus 25x (1300)';
vesselChannel = 4;
getVesselWidthLine = 'yes';
getVesselWidthSnap = 'no';
nearLineMicrons = 4; %how close the scan position must be to a linescan ROI to be included
nearDiamMicrons = 8; %how close the scan position must be to a diameter ROI to be included

%other settings
incorporateFillfraction = true; %for debugging/comparing to other scripts that may not consider fillfraction
tifoptions.color = false; %tif files are in greyscale
tifoptions.compress = 'no'; %tif files are uncompressed
tifoptions.message = false; %don't show warning messages
tifoptions.append = false; %don't append onto existing files
tifoptions.overwrite = true; %overwrite existing files
tifoptions.big = false; %.tif files will be <4 GB

%check inputs
if nargin==7
    setup = varargin{1};
    objective = varargin{2};
    vesselChannel = varargin{3};
    getVesselWidthSnap = varargin{4};
    getVesselWidthLine = varargin{5};
    nearLineMicrons = varargin{6};
    nearDiamMicrons = varargin{7};
else %ask for inputs
    prompt = {'setup:','objective:','vessel channel:','get vessel width from linescan?','get vessel width from snapshot?','max distance (in microns) to blood flow ROI:','max distance (in microns) to diameter ROI:',};
    dlgtitle = 'Inputs';
    dims = [1 50];
    definput = {num2str(setup),objective,num2str(vesselChannel),getVesselWidthLine,getVesselWidthSnap,num2str(nearLineMicrons),num2str(nearDiamMicrons)};
    answers = inputdlg(prompt,dlgtitle,dims,definput);
    if isempty(answers)
        return;
    end
    setup = str2double(answers{1});
    objective = answers{2};
    vesselChannel = str2double(answers{3});
    getVesselWidthLine = answers{4};
    getVesselWidthSnap = answers{5};
    nearLineMicrons = str2double(answers{6});
    nearDiamMicrons = str2double(answers{7});
end

%validate inputs
if ~ischar(setup)
    setup = num2str(setup);
end
assert(any(str2double(setup)==[1 2 3 4]),'setup must be 1, 2, 3, or 4');
assert(strncmpi(getVesselWidthSnap,'y',1)|strncmpi(getVesselWidthSnap,'n',1),'answer to "get vessel width?" must be "yes" or "no"');
assert(strncmpi(getVesselWidthLine,'y',1)|strncmpi(getVesselWidthLine,'n',1),'answer to "get vessel width?" must be "yes" or "no"');
assert(exist('saveastiff','file')==2,'The function "saveastiff()" is not found in Matlab''s paths')


%% check for necessary files
%get linescan data filename/path
[pmtFilename, path] = uigetfile('*.pmt.dat','select the .pmt.dat file');
if pmtFilename==0
    return;
end

extInd = strfind(pmtFilename,'.pmt.dat'); %index to separate filename from extension
assert(~isempty(extInd)&&length(extInd)==1,'selected file is not a .pmt.dat file');

%check that metadata .txt file exists in the same folder, with the same base name
metaFilename = [pmtFilename(1:extInd-1) '.meta.txt'];
assert(isfile(fullfile(path,metaFilename)),['cannot find metadata txt file ' metaFilename]);

%check that scanner .dat file exists in the same folder, with the same base name
scnnrFilename = [pmtFilename(1:extInd-1) '.scnnr.dat'];
assert(isfile(fullfile(path,scnnrFilename)),['cannot find scanner .dat file ' scnnrFilename]);

%get snapshot .tif filename
if strncmpi(getVesselWidthSnap,'y',1)
    tifFiles = ls(fullfile(path,'*.tif'));
    if ~isempty(tifFiles) && size(tifFiles,1)==1
        snapName = tifFiles;
    else
        snapName = uigetfile(fullfile(path,'*.tif'),'select the snapshot .tif file');
        if snapName==0
            return;
        end
    end
end


%% get various metadata from linescan metadata text file
%get date created
metaFileDir = dir(fullfile(path,metaFilename));
date = metaFileDir.date;

%load metadata text
fileID = fopen(fullfile(path,metaFilename),'r');
metadataTxt = fscanf(fileID,'%c');

%get imaging parameters
channels = eval(char(extractBetween(metadataTxt, 'SI.hChannels.channelSave = ', newline)));
numChannels = length(channels);
assert(vesselChannel<=numChannels,['vesselchannel ' num2str(vesselChannel) ' not found in imaging data']);
% framesPerSlice = eval(char(extractBetween(metadataTxt,'SI.hStackManager.framesPerSlice = ', newline)));
samplesPerFrame = eval(char(extractBetween(metadataTxt,'SI.hScan2D.lineScanSamplesPerFrame = ', newline)));
scannerSamplesPerFrame = eval(char(extractBetween(metadataTxt,'SI.hScan2D.lineScanFdbkSamplesPerFrame = ', newline)));
frameRate = eval(char(extractBetween(metadataTxt,'SI.hRoiManager.scanFrameRate = ', newline)));
sampleRateLinescan = eval(char(extractBetween(metadataTxt,'SI.hScan2D.sampleRate = ', newline)));
linescanMsPerLine = 1000/frameRate;

%get metadata for "line" ROIs, size and location (in scan angle degrees)
sizeText = extractBetween(metadataTxt, 'sizeXY": [' , ']');
centerText = extractBetween(metadataTxt, 'centerXY": [' , ']');
roiTypeText = extractBetween(metadataTxt, 'mroi.stimulusfunctions.','"');
lineRois = find(cellfun(@strcmp,roiTypeText,repmat({'line'},size(roiTypeText))));
numRois = length(lineRois);
sizeXY = nan(numRois,2);
centerXY = nan(numRois,2);
for r = 1:numRois
    sizeXY(r,:) = str2num(sizeText{lineRois(r)});
    centerXY(r,:) = str2num(centerText{lineRois(r)}); 
end
plotColors = hsv(2+numRois);

%get metadata for the imaging field of view size
nominalFov = eval(char(extractBetween(metadataTxt,'SI.hScan2D.nominalFovCornerPoints = ', newline)));
nominalFov = reshape(nominalFov,[4 2]);
zoom = eval(char(extractBetween(metadataTxt,'SI.hRoiManager.scanZoomFactor = ',newline)));
fillFraction = eval(char(extractBetween(metadataTxt,'SI.hScan2D.fillFractionSpatial = ', newline)));
actualFov = nominalFov/zoom;
if incorporateFillfraction
    actualFov = fillFraction*actualFov;
end
fovRanges = max(abs(diff(actualFov))); %size of fov [x y]
fovMins = min(actualFov); %min value of fov [x y]
    
%get metadata for fov size in pixels
imageW = eval(char(extractBetween(metadataTxt,'SI.hRoiManager.pixelsPerLine = ', newline)));
imageH = eval(char(extractBetween(metadataTxt,'SI.hRoiManager.linesPerFrame = ', newline)));

%get metadata for scanner feedback
recFeedback = eval(char(extractBetween(metadataTxt,'SI.hScan2D.recordScannerFeedback = ', newline)));
assert(recFeedback==1,'Scanner position data was not logged; use the "NoPosition" script instead');
sampleRateFeedback = eval(char(extractBetween(metadataTxt,'SI.hScan2D.sampleRateFdbk = ', newline)));
numFeedbackChannels = eval(char(extractBetween(metadataTxt,'SI.hScan2D.lineScanNumFdbkChannels = ', newline)));

%if snapshot exists, use header info instead
if strncmpi(getVesselWidthSnap,'y',1)
    snapshot = imread(fullfile(path,snapName),vesselChannel);
    imgdescr = '';
    % Tiff doesn't always work on macs?
%     T = Tiff(fullfile(path,snapName));
%     setDirectory(T,vesselChannel)
%     snapshot = read(T);
%     [imageH,imageW] = size(snapshot);
%     imgdescr = getTag(T,'ImageDescription');
%     date = datetime(eval(char(extractBetween(imgdescr, 'epoch = ', newline))));
end

%get microns/pix measureent
framescanMicronsPerPixel = getPixelSize(setup, date, objective, zoom); %microns per pixel
assert(~isempty(framescanMicronsPerPixel),'could not get the pixel size for this imaging configuration with Pixel-to-Micron');
if diff(framescanMicronsPerPixel)==0 %This measurement should be the same in x and y, but maybe not always...
    framescanMicronsPerPixel = framescanMicronsPerPixel(1); ...so this removes the redundancy
else
    error('um/pixel conversion factors are not equal in x and y dimensions with this imaging configuration; this case isn''t supported yet');
end
switch setup
    case '3'
        calibrationImageSize = [1024 1024];
    case '4'
        calibrationImageSize = [512 512];
    otherwise
        error(['Image resolution during pixel-to-micron calibration is unknown for setup ' setup ' (maybe Nicole knows?)']);
end
framescanMicronsPerPixel = framescanMicronsPerPixel*(calibrationImageSize(1)/imageH);
degreesPerPixel = fovRanges(1)/imageW;
micronsPerDegree = framescanMicronsPerPixel/degreesPerPixel;

%calculate length of linescan line segments (in degrees, pixels, and microns)
lineXDeg = nan(numRois,2);
lineYDeg = nan(numRois,2);
lineXPix = nan(numRois,2);
lineYPix = nan(numRois,2);
centerXPix = nan(numRois,1);
centerYPix = nan(numRois,1);
lineDistPix = nan(numRois,1);
lineDistMicrons = nan(numRois,1);
lineAngle = nan(numRois,1);
for r = 1:numRois
    lineXDeg(r,:) = [centerXY(r,1)-0.5*sizeXY(r,1) centerXY(r,1)+0.5*sizeXY(r,1)]; %x [start stop] of line, in degrees?
    lineYDeg(r,:) = [centerXY(r,2)+0.5*sizeXY(r,2) centerXY(r,2)-0.5*sizeXY(r,2)]; %y [start stop] of line, in degrees?
    lineXPix(r,:) = 1+(imageW-1)*((lineXDeg(r,:)-fovMins(1))/fovRanges(1)); %x [start stop] of line, in pixels
    lineYPix(r,:) = 1+(imageH-1)*((lineYDeg(r,:)-fovMins(2))/fovRanges(2)); %y [start stop] of line, in pixels
    centerXPix(r) = mean(lineXPix(r,:)); %center x of line, in pixels
    centerYPix(r) = mean(lineYPix(r,:)); %center y of line, in pixels
    lineDistPix(r) = sqrt(diff(lineXPix(r,:))^2 + diff(lineYPix(r,:))^2); %line distance, in pixels of snapshot
    lineDistMicrons(r) = lineDistPix(r)*framescanMicronsPerPixel; %line distance, in microns
    lineDistPix(r) = round(lineDistPix(r));
    lineAngle(r) = atan2(sizeXY(r,2),sizeXY(r,1));
end


%% identify the ROIs
if numRois==1
    flowRois = 1;
    diamRois = [];
else %have user verify that the line is drawn correctly
    flowRois = 1:2:numRois;
    diamRois = 2:2:numRois;
    figure();
    ax = axes();
    if strncmpi(getVesselWidthSnap,'y',1)
        imshow(snapshot)
        hold on
    end
    legendNames = cell(1,numRois);
    for r = 1:numRois
        plot(lineXPix(r,:),lineYPix(r,:),'Color',plotColors(r,:),'LineWidth',3); %plot the linescan line over the snapshot
        hold on
        legendNames{r} = ['ROI ' num2str(r)];
    end
    set(ax,'YDir','reverse')
    xlim([1 imageW]);
    ylim([1 imageH]);
    legend(legendNames)
    title('Line position')
    
    prompt={'ROIs for blood flow measurement:','ROIs for diameter measurement:'};
    name='ROI Identification';
    numlines=2;
    defaultanswer={num2str(flowRois),num2str(diamRois)};
    answer = inputdlg(prompt,name,numlines,defaultanswer);
    
    flowRois = str2num(answer{1});
    diamRois = str2num(answer{2});
    assert(any(flowRois'==(1:numRois),2),'Invalid blood flow ROI number entered')
    assert(any(diamRois'==(1:numRois),2),'Invalid diameter ROI number entered')
    close(gcf)
end
roiNames = [];
for r = 1:numRois
    if any(r==flowRois)
        roiNames = [roiNames {['blood flow ROI ' num2str(r)]}];
    elseif any(r==diamRois)
        roiNames = [roiNames {['diameter ROI ' num2str(r)]}];
    else
        roiNames = [roiNames {['excluded ROI ' num2str(r)]}];
    end
end


%% pre-process linescan data
waitfig = waitbar(0,'Loading arbitrary linescan data...');

% get linescan position data from scanner
fileID = fopen(fullfile(path,scnnrFilename));
data = fread(fileID,'single');
scannerPos = reshape(data,numFeedbackChannels,scannerSamplesPerFrame,[]);
scannerPos = permute(scannerPos,[2 3 1]);
medScannerPos = squeeze(median(scannerPos,2)); %[x y]
% numScannerSamples = size(scannerPos,2);
scannerPosTime = (0:scannerSamplesPerFrame-1)/sampleRateFeedback;
% fullScannerPosTime = (0:numScannerSamples-1)/sampleRateFeedback;

%read linescan data
fileID = fopen(fullfile(path,pmtFilename));
data = fread(fileID,'int16');
linescanData = reshape(data,numChannels,samplesPerFrame,[]);
linescanData = linescanData(vesselChannel,:,:);
linescanData = permute(linescanData,[2 3 1]);
numLines = size(linescanData,2);
linescanTime = (0:samplesPerFrame-1)/sampleRateLinescan;
% fullLineTime = (0:framesPerSlice-1)/frameRate;

%make sure that the scan and data samples line up in time
if length(scannerPosTime)~=length(linescanTime) || any(scannerPosTime~=linescanTime)
    if (scannerPosTime(end)/linescanTime(end))>0.95
        %resample scanner pos to match linescan sample rate
        medScannerPosFull = align_data(medScannerPos(:,1),scannerPosTime',linescanTime');
        medScannerPosFull = fillmissing(medScannerPosFull,'makima');
        medScannerPosFull(:,2) = align_data(medScannerPos(:,2),scannerPosTime',linescanTime');
        medScannerPosFull(:,2) = fillmissing(medScannerPosFull(:,2),'makima');
        fitInds = false(size(linescanTime)); %no fit needed
    else 
        fprintf(['Only ' num2str(100*round(scannerPosTime(end)/linescanTime(end))) ' percent of scan position was logged -- missing positions will be estimated.\n']);
        %fit fourier series to position data to replace missing points and interpolate values
        period = (samplesPerFrame/sampleRateLinescan);
        tripleTime = [scannerPosTime'; period+scannerPosTime'; (2*period)+scannerPosTime'];
        triplePos = [medScannerPos; medScannerPos; medScannerPos];
        goodnessOfFit = 0;
        numFourier = numRois*3-1;
        while goodnessOfFit<0.95
            numFourier = numFourier+1;
            assert(numFourier<11,'scanner position curve could not be fit well with less than 11 terms of a fourier series');
            [fX, gX] = fit(tripleTime,triplePos(:,1),['fourier' num2str(numFourier)]);
            [fY, gY] = fit(tripleTime,triplePos(:,2),['fourier' num2str(numFourier)]);
            goodnessOfFit = gX.adjrsquare*gY.adjrsquare;
        end
        medScannerPosFull = [fX(linescanTime) fY(linescanTime)];
        fitInds = linescanTime>scannerPosTime(end);
    end
end

%process linescan data for each ROI
roiInds = false(numRois,samplesPerFrame);
linearizedLinescans = cell(numRois,1);
linearizedMicronsPerPixel = nan(numRois,1);
for r = 1:numRois
    waitfig = waitbar(0.1 + (0.9/numRois)*(r-1),waitfig,['Linearizing data for ROI ' num2str(r) '...']);
    
    %rotate ROI so that it's horizontal
    lineXY = ([lineXDeg(r,:)' lineYDeg(r,:)']) - repmat(centerXY(r,:),[2 1]);
    lineXYAngle = atan2(lineXY(:,2),lineXY(:,1));
    lineXYRadius = sqrt(lineXY(:,1).^2 + lineXY(:,2).^2);
    rotLineXY = lineXY;
    rotLineXYAngle = lineXYAngle + lineAngle(r);
    rotLineXY(:,1) = lineXYRadius.*cos(rotLineXYAngle); %after rotating, y should be 0
    rotLineXY(:,2) = lineXYRadius.*sin(rotLineXYAngle);
    assert(sum(abs(rotLineXY(:,2)))/sum(abs(rotLineXY(:,1)))<0.01,'bug in code: roi not rotated correctly');

    %rotate scanner position around ROI center to match the rotated ROI
    rotScanPos = medScannerPosFull - repmat(centerXY(r,:),[samplesPerFrame 1]); %centered on diam ROI
    scanAngle = atan2(rotScanPos(:,2),rotScanPos(:,1));
    scanRadius = sqrt(rotScanPos(:,1).^2 + rotScanPos(:,2).^2);
    rotScanAngle = scanAngle+lineAngle(r);
    rotScanPos(:,1) = scanRadius.*cos(rotScanAngle); %rotate scanner pos to match diam ROI
    rotScanPos(:,2) = scanRadius.*sin(rotScanAngle);

    %find linescan positions that are near the ROI(s)
    if any(r==flowRois)
        nearRoiIdx = abs(rotScanPos(:,2))<=(nearLineMicrons/micronsPerDegree); %restrict to +/- nearRoiMicrons
    else
        nearRoiIdx = abs(rotScanPos(:,2))<=(nearDiamMicrons/micronsPerDegree); %restrict to +/- nearRoiMicrons
    end
    CC = bwconncomp(nearRoiIdx);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [~,biggestComp] = max(numPixels);
    nearRoiIdx = false(1,samplesPerFrame);
    nearRoiIdx(CC.PixelIdxList{biggestComp}) = true;
    
    %restrict roiInds to a single direction along the ROI (ignore turnaround during "pauses" that are still near the ROI)
    scanVelocityAlongRoi = [0 diff(rotScanPos(:,1)')];
    scanVelocityAlongNearRoi = zeros(1,samplesPerFrame);
    scanVelocityAlongNearRoi(nearRoiIdx) = scanVelocityAlongRoi(nearRoiIdx);
    medianVelNearRoi = median(scanVelocityAlongNearRoi(nearRoiIdx));
    scanDirectionNearRoiCenter = medianVelNearRoi/abs(medianVelNearRoi);
    roiInds(r,:) = nearRoiIdx & (scanDirectionNearRoiCenter*scanVelocityAlongNearRoi)>0;
    roiRows = find(roiInds(r,:));
    
    %calculate the largest MicronsPerPixel of this linescan
    maxMicronsPerPixel = diff(rotScanPos(roiRows,1)')*micronsPerDegree; %(D/P)(M/D) = M/P
    actualMPPsign = maxMicronsPerPixel(round(end/2))/abs(maxMicronsPerPixel(round(end/2)));
    maxMicronsPerPixel = max(maxMicronsPerPixel*actualMPPsign);
    
    %bin the linescan points linearly across the ROI
    linearSamplesPerFrame = round(lineDistMicrons(r)/maxMicronsPerPixel);
    alignedX = linspace(rotLineXY(1,1),rotLineXY(2,1),linearSamplesPerFrame);
    unalignedX = rotScanPos(roiRows,1)';
    unalignedY = linescanData(roiRows,:);
    tmp = align_data(unalignedY(:,1)',unalignedX,alignedX); %linearize the 1st one to see where it starts/stops
    startInd = find(~isnan(tmp),1);
    stopInd = find(~isnan(tmp),1,'last');
    tmpLinearizedLinescans = zeros([linearSamplesPerFrame numLines]);
    for i = 1:numLines
        if mod(i,1000)==0
            waitfig = waitbar(0.1 + (0.9/numRois)*(r-1) + (0.9/numRois)*(i/numLines),waitfig);
        end
        tmp = align_data(unalignedY(:,i)',unalignedX,alignedX);
        tmp2 = tmp(startInd:stopInd); %fill in any missing data, but only inside the linescan
        tmp2 = fillmissing(tmp2,'linear');
        tmp(startInd:stopInd) = tmp2;
        tmp(isnan(tmp)) = 0;
        tmpLinearizedLinescans(:,i) = tmp;
    end
    linearizedLinescans(r) = {tmpLinearizedLinescans};
    linearizedMicronsPerPixel(r) = lineDistMicrons(r)/linearSamplesPerFrame;
end
close(waitfig)


%% create figure of roi lines and full scannerpos
marker_size = 8;
line_width = 5;
fig1 = figure('Position',[100 100 1500 400]);
ax = subplot(1,3,1);
for r = 1:numRois
    plot(lineXDeg(r,:),lineYDeg(r,:),'Color',plotColors(r,:),'LineWidth',line_width);
    hold on
end
loggedRows = find(~fitInds);
interpRows = find(fitInds);
scatter(medScannerPosFull(loggedRows,1),medScannerPosFull(loggedRows,2),marker_size,'filled','MarkerEdgeColor','none','MarkerFaceColor',plotColors(numRois+1,:));
scatter(medScannerPosFull(interpRows,1),medScannerPosFull(interpRows,2),marker_size,'filled','MarkerEdgeColor','none','MarkerFaceColor',plotColors(numRois+2,:));
legend([roiNames {'logged scanner positions','interpolated scanner positions'}]);
axis equal
xlimits = ax.XLim;
xlimits = [xlimits(1)-0.5 xlimits(2)+0.5];
xlim(xlimits);
ylimits = ax.YLim;
ylimits = [ylimits(1)-1 ylimits(2)];
ylim(ylimits);
ylabel('y scan angle (deg)')
xlabel('x scan angle (deg)')
set(gca,'YDir','reverse')

subplot(1,3,2)
for r = 1:numRois
    plot(lineXDeg(r,:),lineYDeg(r,:),'Color',plotColors(r,:),'LineWidth',line_width);
    hold on
end
legendNames = roiNames;
for r = 1:numRois
    scatter(medScannerPosFull(find(roiInds(r,:)),1),medScannerPosFull(find(roiInds(r,:)),2),marker_size,'filled','MarkerEdgeColor','none','MarkerFaceColor',plotColors(1+numRois,:));
    legendNames = [legendNames {['ROI ' num2str(r) ' positions']}];
end
scatter(medScannerPosFull(find(~any(roiInds,1)),1),medScannerPosFull(find(~any(roiInds,1)),2),marker_size,'filled','MarkerEdgeColor','none','MarkerFaceColor',[0.4 0.4 0.4]);
legend([legendNames {'excluded positions'}]);
axis equal
xlim(xlimits);
ylim(ylimits);
ylabel('y scan angle (deg)')
xlabel('x scan angle (deg)')
set(gca,'YDir','reverse')

subplot(1,3,3)
for r = 1:numRois
    plot(lineXDeg(r,:),lineYDeg(r,:),'Color',plotColors(r,:),'LineWidth',line_width);
    hold on
end
colorVec = repmat(plotColors(numRois+1,:),[samplesPerFrame 1]);
brightnessVec = repmat(mean(linescanData,2)/max(mean(linescanData,2)),[1 3]);
scatter(medScannerPosFull(:,1),medScannerPosFull(:,2),marker_size,colorVec.*brightnessVec,'filled');
legend([roiNames {'mean pixel value'}]);
axis equal
xlim(xlimits);
ylim(ylimits);
ylabel('y scan angle (deg)')
xlabel('x scan angle (deg)')
set(gca,'YDir','reverse')

%save figure as tiff
saveas(fig1,fullfile(path,['AL3_AllROIs_' pmtFilename(1:extInd-1) '.tif']));
  

%% create tif file of linescan data (for line ROI(s) only)
for r = flowRois
    linearData = linearizedLinescans{r};
    samplesPerLine = size(linearData,1);
    
    %reshape data to create a multiframe tif
    possibleLinesPerFrame = [500:100:1000 501:1024]; %my subjective preference for image sizes
    validDivisors = rem(numLines,possibleLinesPerFrame)==0;
    if any(validDivisors)
        linesPerFrame = possibleLinesPerFrame(find(validDivisors,1));
        numFrames = ceil(numLines/linesPerFrame);
        linescanTif = reshape(linearData,[samplesPerLine linesPerFrame numFrames]);
    else
        linesPerFrame = samplesPerLine;
        numFrames = ceil(numLines/samplesPerLine);
        linescanTif = nan([samplesPerLine linesPerFrame numFrames]);
        linescanTif(1:numel(linearData)) = linearData;
    end
    linescanTif = permute(linescanTif,[2 1 3]); %change lines to horizontal

    %save linescan data to tif file
    if numel(linescanTif)>=40000000 %this much data may create a .tif file close to or over 4 GB
        tifoptions.big = true; % "big" option uses 64-bit addressing
    end
    linescanSavename = ['AL1_ROI' num2str(r) '_LinescanCh' num2str(vesselChannel) '_' pmtFilename(1:extInd-1) '.tif'];
    tifFilename = fullfile(path,linescanSavename);
    saveastiff(linescanTif,tifFilename,tifoptions,metadataTxt);
end


%% (optional) get vessel diameter from linescan diameter ROI
vesselDiameterInMicronsLine = [];
if strncmpi(getVesselWidthLine,'y',1)
    for r = diamRois
        diamImage = linearizedLinescans{r};
        diamImageH = size(diamImage,1);
        diamImageW = 512;
        diamImage = imresize(diamImage,[diamImageH diamImageW]); %create image of these inds,
        diamImage = diamImage/prctile(diamImage,99,'all');

        %have the user calculate the edges of the vessel
        fig4 = figure(); %have user move lines for diameter
        ax4 = axes();
        imshow(diamImage)
        topPos = diamImageH*0.25;
        botPos = diamImageH*0.75;
        top = drawline('Position',[0 topPos; diamImageW topPos],'LineWidth',0.5,'label','drag to vessel top','LabelAlpha',0.2);
        bot = drawline('Position',[0 botPos; diamImageW botPos],'LineWidth',0.5,'label','drag to vessel bottom','LabelAlpha',0.2);
        title('estimate vessel diameter')
        fig4.Position = [100 100 diamImageW diamImageH];
        ax4.Position = [0 0 1 1];
        c = uicontrol(fig4,'String','set vessel edges','Callback','uiresume(gcbf)','FontSize',12);
        c.Position(3:4) = [150 30];
        uiwait(gcf)
        clear c
        diamRows = find(roiInds(r,:));
        topInd = diamRows(round(mean(top.Position(:,2)))); %index of scanner position on top side of vessel
        botInd = diamRows(round(mean(bot.Position(:,2)))); %index of scanner position on bottom side of vessel

        %find vessel width
        diamInPixels = abs(topInd-botInd);
        diamInMicrons = diamInPixels*linearizedMicronsPerPixel(r); %median microns (in y) of the vessel
        vesselDiameterInMicronsLine = [vesselDiameterInMicronsLine diamInMicrons];
        
        %save figure as tif
        fontColor = [1 1 1];
        backgroundColor = [0 0 0];
        bot.Label = '';
        top.Label = '';
        text(0.01*diamImageW,6,['microns/pixel: ' num2str(linearizedMicronsPerPixel(r))],'FontSize',8,'Color',fontColor,'BackgroundColor',backgroundColor,'Margin',1);
        text(0.01*diamImageW,22,['median vessel width (pixels): ' num2str(diamInPixels)],'FontSize',8,'Color',fontColor,'BackgroundColor',backgroundColor,'Margin',1);
        text(0.01*diamImageW,38,['median vessel width (microns): ' num2str(diamInMicrons)],'FontSize',8,'Color',fontColor,'BackgroundColor',backgroundColor,'Margin',1);

        vesselWidthImage = getframe(fig4);
        vesselWidthImage = vesselWidthImage.cdata;
        vesselWidthImage = imresize(vesselWidthImage,[diamImageH diamImageW]);
        saveastiff(vesselWidthImage,fullfile(path,['AL4_ROI' num2str(r) '_LinescanVesselWidth_' pmtFilename(1:extInd-1) '.tif']),tifoptions);
        close(gcf)
    end
end


%% (optional) get vessel width from snapshot
vesselDiameterInMicronsSnap = [];
if strncmpi(getVesselWidthSnap,'y',1)
    for r = flowRois
        %crop the snapshot around the line
        fig = figure();
        ax = axes();
        imshow(snapshot);
        fig.Position = [0 0 imageW imageH];
        ax.Position = [0 0 1 1];
        cropRange = 200; %# of pixels around the line center to crop
        xlim([centerXPix(r)-cropRange centerXPix(r)+cropRange]); %crop in x
        ylim([centerYPix(r)-cropRange centerYPix(r)+cropRange]); %crop in y
        snapWithLineCropped = getframe(fig);
        snapWithLineCropped = mean(snapWithLineCropped.cdata,3);
        snapWithLineCropped = imresize(snapWithLineCropped,[2*cropRange 2*cropRange]);
        close(gcf);

        %rotate and crop the image so that the drawn vessel line is straight
        lineAngle = atan2(diff(lineYPix(r,:)),diff(lineXPix(r,:)));
        snapRotated = imrotate(snapWithLineCropped,rad2deg(lineAngle),'bilinear','crop');
        snapRotated = snapRotated/prctile(snapRotated,99,'all'); %scale brightness to more visible range
        fig = figure();
        ax = axes();
        imshow(snapRotated)
        xlim([cropRange-lineDistPix(r)/2 cropRange+lineDistPix(r)/2]); %crop image in x so it only includes the line
        fig.Position = [0 0 lineDistPix(r) 2*cropRange];
        ax.Position = [0 0 1 1];

        %get the image of the straightened/cropped vessel
        snapVessel = getframe(fig);
        snapVessel = mean(snapVessel.cdata,3);
        snapVessel = imresize(snapVessel,[2*cropRange lineDistPix(r)]);
        snapVessel = snapVessel/prctile(snapVessel,99,'all');
        close(fig)

        %draw lines around the vessel and let the user move them
        fig3 = figure();
        fig3.Position = [0 0 2*lineDistPix(r) 4*cropRange];
        ax3 = axes();
        imshow(snapVessel)
        topPos = 0.5*cropRange;
        botPos = 1.5*cropRange;
        top = drawline('Position',[0 topPos; lineDistPix(r) topPos],'LineWidth',0.5,'label','drag to vessel top','LabelAlpha',0.2);
        bot = drawline('Position',[0 botPos; lineDistPix(r) botPos],'LineWidth',0.5,'label','drag to vessel bottom','LabelAlpha',0.2);
        title('estimate vessel diameter')
        c = uicontrol(fig3,'String','set vessel edges','Callback','uiresume(gcbf)','FontSize',12);
        c.Position(3:4) = [150 30];
        uiwait(gcf)
        diamInPixels = abs(mean(top.Position(:,2)) - mean(bot.Position(:,2)));
        diamInMicrons = diamInPixels*framescanMicronsPerPixel; %median microns (in y) of the vessel
        vesselDiameterInMicronsSnap = [vesselDiameterInMicronsSnap diamInMicrons];
        
        %create vessel width estimation summary figure
        fontColor = [1 1 1];
        backgroundColor = [0 0 0];
        fig3.Position = [0 0 lineDistPix(r) 2*cropRange];
        ax3.Position = [0 0 1 1];
        bot.Label = '';
        top.Label = '';
        clear c
        text(0.01*lineDistPix(r),6,['microns/pixel: ' num2str(framescanMicronsPerPixel)],'FontSize',8,'Color',fontColor,'BackgroundColor',backgroundColor,'Margin',1);
        text(0.01*lineDistPix(r),22,['median vessel width (pixels): ' num2str(diamInPixels)],'FontSize',8,'Color',fontColor,'BackgroundColor',backgroundColor,'Margin',1);
        text(0.01*lineDistPix(r),38,['median vessel width (microns): ' num2str(diamInMicrons)],'FontSize',8,'Color',fontColor,'BackgroundColor',backgroundColor,'Margin',1);

        %save figure as tif
        vesselWidthImage = getframe(fig3);
        vesselWidthImage = vesselWidthImage.cdata;
        vesselWidthImage = imresize(vesselWidthImage,[2*cropRange lineDistPix(r)]);
        saveastiff(vesselWidthImage,fullfile(path,['AL5_ROI' num2str(r) '_SnapshotVesselWidth_' pmtFilename(1:extInd-1) '.tif']),tifoptions,imgdescr);
        close(gcf)
    end
end


%% display/save final results
%save results to a .mat file
resultsSavename = ['AL2_Results_'  pmtFilename(1:extInd-1) '.mat'];
results.bloodFlowROIs = flowRois;
results.diameterROIs = diamRois;
results.vesselDiameterInMicronsSnap = vesselDiameterInMicronsSnap;
results.vesselDiameterInMicronsLine = vesselDiameterInMicronsLine;
results.bloodflowMicronsPerPixel = linearizedMicronsPerPixel(flowRois);
results.linescanMsPerLine = linescanMsPerLine;
save(fullfile(path,resultsSavename),'results');

%print results in a messagebox
% txt1 = strcat("framescan microns/pixel: ",num2str(micronsPerPixel));
txt1 = strcat("blood flow ROIs: ", num2str(flowRois));
txt2 = strcat("diameter ROIs: ", num2str(diamRois));
txt3 = strcat("estimated vessel diameter from snapshot (microns): ", num2str(vesselDiameterInMicronsSnap));
txt4 = strcat("estimated vessel diameter from linescan (microns): ", num2str(vesselDiameterInMicronsLine));
txt5 = strcat("blood flow microns/pixel: ", num2str(linearizedMicronsPerPixel(flowRois)));
txt6 = strcat("linescan ms/line: ", num2str(linescanMsPerLine));
txt7 = strcat("results saved in '", path, "'");

msgText = ["Pre-processing complete";"";txt1;txt2];
if strncmpi(getVesselWidthSnap,'y',1)
    msgText = [msgText;"";txt3];
end
if strncmpi(getVesselWidthLine,'y',1)
    msgText = [msgText;"";txt4];
end
msgText = [msgText;"";txt5;txt6;"";txt7];
msgbox(msgText,'Results','help');

if nargout>0
    varargout{1} = results;
end
if nargout>1
    varargout{2} = tifFilename;
end