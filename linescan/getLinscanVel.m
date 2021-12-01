function Result = getLinescanVel(varargin)
% Process linescan files and calculate velocity
% Data file should be in greyscale, '.tif' format. This program assumes the data
% is stored as sequential images and that there is no gap in time between
% the last line of a frame and the first line of the next frame. 
% Part 1 allows user to choose parameters and a region of interest for each line
% scan file. 
% Saves filenames and parameters in a "bfdata.mat" file and also a ".csv" file that can
% be opened in Excel or a text program. 
% Part 2 uses parameters and filenames saved in "bfdata.mat" file to
% then calculate velocites in each file.
% OUTPUT: variable called 'Result' is saved in separate files for each .tif file       
%   1) line number 
%   2) time (ms)
%   3) Velocity (mm/s), positive veloctiy indicates RBCs going from left to
%   right
%   4) Sep (Seperability)
%   5) Angle (angle of stripes in radians)
%   6) blank
%   7) blank
%   8) blank
%   9) blank
%   10) blank
%
% Please reference the following publications:
% C. B. Schaffer, B. Friedman, N. Nishimura, L. F. Schroeder, P. S. Tsai,
% F. F. Ebner, P. D. Lyden, and D. Kleinfeld, 
% Two-photon imaging of cortical surface microvessels reveals 
% a robust redistribution in blood flow after vascular occlusion.  
% Public Library of Science Biology 4, e22 (2006).
% 
% D. Kleinfeld, P.P. Mitra, F. Helmchen, W. Denk, Fluctuations and 
% stimulus-induced changes in blood flow observed in individual capillaries 
% in layers 2 through 4 of rat neocortex. Proc Natl Acad Sci U S A 95, 15741- 
% 15746 (1998). 
%
% Questions, bugs, etc. please contact:
% Nozomi Nishimura
% nn62@cornell.edu
% last mod 02-07-09

% INPUTS
p = inputParser();
p.addOptional('Openfile','',@ischar);
% TODO: these should probably be required
p.addOptional('msPerLine',0.6,@(x) isnumeric(x)&&isscalar(x));
p.addOptional('umPerPx',0.1,@(x) isnumeric(x)&&isscalar(x));
% TODO: more strict validation functions
p.addOptional('WinSize',75,@(x) isnumeric(x)&&isscalar(x));
p.addOptional('WinPixelsDown',50,@(x) isnumeric(x)&&isscalar(x));
% TODO: these should probably be moved up to calling function
p.addOptional('WinLeft',[],@(x) isnumeric(x)&&isscalar(x));
p.addOptional('WinRight',[],@(x) isnumeric(x)&&isscalar(x));
p.addOptional('Maxlines',inf);
p.addOptional('UseAvg',false,@islogical);
p.addOptional('errorcheck',false,@islogical);
p.addParameter('Method','Radon',@(x) any(strcmp(x, {'Radon', 'SVD'})))
% TODO: allow user to flip image for opposite velocity?
% if I_sign==0
%     I=fliplr(I);
% end

p.parse(varargin{:});

% TODO: Probably don't need to reassign most of these
Tfactor = 1/p.Results.msPerLine; % ypixel per ms
Xfactor = p.Results.umPerPx; % microns per xpixel
WinSize = p.Results.WinSize;
WinPixelsDown = p.Results.WinPixelsDown;
WinLeft = p.Results.WinLeft;
WinRight = p.Results.WinRight;
Maxlines = p.Results.Maxlines;
UseAvg = p.Results.UseAvg;
errorcheck = p.Results.errorcheck;

%    % actual data used is only center circle ~70% of area (square window)
%  % number of pixels between top of last window and next window

% TODO: allow output as CSV?

% TODO: make sure all output files are being saved in input directory;
% Get input folder and file
if isempty(p.Results.Openfile)
    % Get file to open
    [fname,pname] = uigetfile('*.*');
    Openfile = fullfile(pname, fname);
    disp(Openfile);
else
    Openfile = p.Results.Openfile;
    [pname,fname] = fileparts(Openfile);
end



%   fileinfo = imfinfo(Openfile);

    % get time file was created
    info = dir(Openfile);
    filetime = info.date;
        
    %% Select bounding box



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate velocities or diameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\

    
% TODO: should use this default name or ask user?
    %Datafile = [char(strrep(OpenName(i),'.tif',['--wpd', num2str(WinPixelsDown)])), date, '.mat'];
    Datafile = [char(strrep(Openfile,'.tif',[' rawVel', num2str(WinPixelsDown),num2str(WinSize)])),'.mat'];


    % read in image data and reshape
    I = f_loadTiff(Openfile);
    nPix = size(I, 2);
    I = permute(I, [1,3,2]);
    I = reshape(I, [], nPix);
    nLines = size(I, 1);
    
    if isempty(WinLeft)
        WinLeft = 1;
    end
    
    if isempty(WinRight)
        WinRight = nPix;
    end
    
    Maxlines = min(Maxlines, nLines);
    
    % subtract average value of each column from image to take out vertical stripes
    % TODO: should this happen frame-by-frame, block-by-block or over
    % entire image??
    % Before this was happening block-by-block which seems a bit sus
    % Maybe should move out depening on which
    if UseAvg
        I = I-mean(I);
    end
%     [r c] = size(I);
%     avgPixValues = ones(r,1)*mean(I);
% 
%     if UseAvg
%         I = I-avgPixValues;
%     end
%     clear avgPixValues
%     % TODO: move this out into calling function
%     % Take out vertical stripes
%     blocksize= size(block);
%     avg = mean(block);
%     avgs = ones([blocksize(1), 1])*avg;
%     if useaverage
%         block = block-avgs;
%     end
%     clear avgs avg;
    
    % Loop through lines
    
    % Calculate block indices
    last = WinSize:WinPixelsDown:Maxlines;
    first = last - WinSize + 1;
    nWins = length(first);

    % Initialize outputs
    Result = zeros(nWins,10);
    
    % Create waitbar
    hWait = waitbar(0);
    
    % Pick function to calculate velocity
    switch p.Results.Method
        case 'Radon'
            getLinescanVelFcn = @(block) method.getLinescanVelRadon(block, Tfactor, Xfactor);
        case 'SVD'
            getLinescanVelFcn = @(block) method.getLinescanVelSVD(block, Tfactor, Xfactor);
    end
    
    for iWin = 1:1:nWins
        % TODO: use im2double instead?
        block = double(I(first(iWin):last(iWin), WinLeft: WinRight));
        veldata = getLinescanVelFcn(block);

        veldata(1) = first(iWin);
        veldata(2) = iWin*WinPixelsDown/Tfactor;
    % ---------------------------------------------
        % For Debugging
        if (errorcheck ==1) && (npoints< 20)
            subplot(2,1,1);imagesc(lines); f_niceplot;
            title(Openfile)
            subplot(2,1,2); imagesc(block); f_niceplot;title('block')
            angle = acot(veldata(3)/Xfactor/Tfactor)*180/pi;
            title([num2str(veldata(1)), ' vel:', num2str(veldata(3)), ' angle: ', num2str(angle)]);
            xlabel('press a key to continue');
            pause;
        end
        % ---------------------------------------------
        veldata(5) = -1*acot(veldata(3)/Xfactor/Tfactor);

        Result(iWin, :) = veldata;
    end

    waitbar(1, hWait, 'Done!')

    % TODO: this should be moved out into calling function 
    save(Datafile,'Result', 'Tfactor', 'WinPixelsDown','Openfile', 'filetime', 'WinLeft', 'WinRight', 'Tfactor', 'Xfactor');

%     clear data data1 cropped Result Rotdata time Data;

    close(hWait);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = f_niceplot

axis image; colormap gray;
set(gca, 'XTickLabel', [])
set(gca, 'YTickLabel', [])

function I = f_loadTiff(filename)
    T = Tiff(filename, 'r');
    
    % Extracting Matlab Type Name and Bits per Sample
    typeName = T.getTag('SampleFormat');
    if typeName == 1
       typeName = 'uint';
    elseif typeName  == 2
       typeName = 'int';
    else
       error('Tiff Sample Format is not supported. It must be UInt or Int');
    end
    bits = T.getTag('BitsPerSample');
    type = [typeName, num2str(bits)];

    % Get Image Height and Width
    width = T.getTag('ImageWidth');
    height = T.getTag('ImageLength');
            
    % Get number of images in Tiff Stack
    T.setDirectory(1);
    numImages = 1;
    while ~T.lastDirectory
       numImages = numImages + 1;
       T.nextDirectory();
    end

    I = zeros(height, width, numImages, type);
    % Setting data for each Image File Directory
    for IFD = 1:1:numImages
       T.setDirectory(IFD);
       I(:,:,IFD) = T.read();
    end
        
    close(T);    



