function Result = calcLinescanVel(varargin)
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
% [ADD NEW PUBLICATION HERE]
% 
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
% Nash Allan-Rahill
% na428@cornell.edu
% last modified Feb 23, 2021

% Parse and validate inputs
p = inputParser();
p.addRequired('I',@ismatrix);
p.addRequired('msPerLine',@(x) isnumeric(x)&&isscalar(x));
p.addRequired('umPerPx',@(x) isnumeric(x)&&isscalar(x));
% TODO: more strict validation functions
p.addOptional('WinSize',75,@(x) isnumeric(x)&&isscalar(x));
p.addOptional('WinStep',50,@(x) isnumeric(x)&&isscalar(x));
p.addOptional('errorcheck',false,@islogical);
% TODO: add validiation function? @(x) any(strcmp(x, {'fminsearch', 'legacy'})))
p.addParameter('Transform','Radon',@(x) any(strcmp(x, {'Radon', 'Rotate'})));
p.addParameter('Metric','Var',@(x) any(strcmp(x, {'Sep', 'Var'})));
p.addParameter('Optimizer','fminbnd');
p.addParameter('FilterVar', 0, @isfinite);
p.parse(varargin{:});

% TODO: Probably don't need to reassign most of these
I = p.Results.I;
Tfactor = 1/p.Results.msPerLine; % ypixel per ms
Xfactor = p.Results.umPerPx; % microns per xpixel
WinSize = p.Results.WinSize;
WinPixelsDown = p.Results.WinStep;
errorcheck = p.Results.errorcheck;
a = p.Results.FilterVar;
% TODO: allow user to flip image for opposite velocity?
% if I_sign==0
%     I=fliplr(I);
% end

% TODO: should use this default name or ask user?
    %Datafile = [char(strrep(OpenName(i),'.tif',['--wpd', num2str(WinPixelsDown)])), date, '.mat'];
% TODO: consider making this the default for all radon methods
% Force filtering at a = 25 for legacy Radon
if strcmp(p.Results.Transform, 'Radon') && strcmp(p.Results.Optimizer, 'radonlegacy')
    a = 25;
end
    


%% Filter Image (if requested)
% TODO: move this into subfunction?
% This filtering step seems to be very important for radon method in
% general
% TODO: should image be cast to double first?
if a
    %Create high pass flter using isotropic Gaussian
    I_size = size(I);
    y = (1:I_size(1))';
    x = 1:I_size(2);
    gaus = 1 - exp(-.5 * ( ((y-0.5*I_size(1)).^2 / a) + ((x-0.5*I_size(2)).^2 / a) ) );

    %Multiply frequency components of image by Gaussian filter
    I = fftshift(fft2(I));
    I = I.*gaus;
    I = ifft2(ifftshift(I));
    I = real(I);
end

%% Loop through lines
    
    % Calculate block indices
    last = WinSize:WinPixelsDown:size(I,1);
    first = last - WinSize + 1;
    nWins = length(first);

    % Initialize outputs
    Result = zeros(nWins,6);
    
    % Create waitbar
    % TODO: move this waitbar out to calcLinescanVelTiff?
    hWait = waitbar(0, 'Calculating linescan velocity', 'Name', 'Linescan');
    
    
    % Function handle for calculating linescan slope
    calcLinescanSlopeFcn = @(block) linescan.calcLinescanSlope(...
        block, 'Optimizer', p.Results.Optimizer, 'Transform',...
        p.Results.Transform, 'Metric', p.Results.Metric);

%     % Pick function to calculate velocity
%     switch p.Results.Method
%         case 'Radon'
% %             calcLinescanVelFcn = @(block) linescan.method.calcLinescanVelRadon(block, Tfactor, Xfactor);
%             calcLinescanSlopeFcn = @(block) linescan.method.calcLinescanSlopeRadon(block, p.Results.Optimizer);
%         case 'SVD'
% %             calcLinescanVelFcn = @(block) linescan.method.calcLinescanVelSVD(block, Tfactor, Xfactor);
%             calcLinescanSlopeFcn = @(block) linescan.method.calcLinescanSlopeSVD(block, p.Results.Optimizer);
%     end
    
    for iWin = 1:1:nWins
        % TODO: use im2double instead?
        block = double(I(first(iWin):last(iWin), :));
%         block = im2double(I(first(iWin):last(iWin), :));
        [dYdt, metric] = calcLinescanSlopeFcn(block);
        
        % TODO: change this to first, last? Or this is supposed to be time?
        Result(iWin,1) = first(iWin);
        % TODO: should time be average time or start time of block?
        Result(iWin,2) = iWin*WinPixelsDown/Tfactor;
        % Velocity
        Result(iWin,3) = dYdt*Xfactor*Tfactor;
        Result(iWin,4) = metric;
        
        
%         veldata = [0, 0, veldata, 0];
%         veldata(1) = first(iWin);
%         veldata(2) = iWin*WinPixelsDown/Tfactor;
    % ---------------------------------------------
        % For Debugging
        if errorcheck && npoints < 20
            subplot(2,1,1);imagesc(lines); f_niceplot;
            title(Openfile)
            subplot(2,1,2); imagesc(block); f_niceplot;title('block')
            angle = acot(veldata(3)/Xfactor/Tfactor)*180/pi;
            title([num2str(veldata(1)), ' vel:', num2str(veldata(3)), ' angle: ', num2str(angle)]);
            xlabel('press a key to continue');
            pause;
        end
        % ---------------------------------------------
%         veldata(4) = -1*acot(veldata(3)/Xfactor/Tfactor);

%         Result(iWin, :) = veldata;
        
        waitbar(iWin/nWins, hWait);
    end

    waitbar(1, hWait, 'Done!')
    
%     clear data data1 cropped Result Rotdata time Data;

    close(hWait);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function f_niceplot
    axis image; colormap gray;
    set(gca, 'XTickLabel', [])
    set(gca, 'YTickLabel', [])
end
