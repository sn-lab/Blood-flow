%% Calculate linescan velocity (Tiff)
linescan.calcLinescanVelTiff()

%% Clean (threshold) linescan velocity
linescan.cleanLinescanVel()


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate linescan velocity batch from spreadsheet (Tiff)
[file, path] = uigetfile({'*.csv;*.xls;*.xlsb;*.xlsm;*.xlsx;*.xltm;*.xltx;*.ods'}, 'Select Linescan Settings File');
T = readtable(fullfile(path,file));
% TODO: check if this is a double first
T.MaxLines = str2double(T.MaxLines);
linescan.calcLinescanVelTiff(T);

%% Plot linescan against image (Tiff)
imagefile = '*.tif';
linescanfile = '*.mat';

% Load linescan file
Vel = load(linescanfile);

% Read Tiff
hTiffReader = util.io.readTiffStack(filepath);
I = permute(hTiffReader.data(), [2,1,3]);
delete(hTiffReader);

% Reshape linescan to 2D
nPix = size(I, 2);
I = permute(I, [1,3,2]);
I = reshape(I, [], nPix);
nLines = size(I, 1);
I = I(1:Vel.Settings.MaxLines, Vel.Settings.MaskLeft:Vel.Settings.MaskRight);

% Plot
linescan.plotLinescan(I, mean(Vel.Result{:,1:2}, 2), vel.Result.Velocity);
