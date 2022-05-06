function Result = calcLinescanVel(varargin)
% calcLinescanVel Calculate velocity from linescan kymograph
%   Detailed explanation goes here
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


%% Parse and validate inputs
p = inputParser();
p.addRequired('I',@ismatrix);
p.addRequired('msPerLine',@(x) isnumeric(x)&&isscalar(x));
p.addRequired('umPerPx',@(x) isnumeric(x)&&isscalar(x));
% TODO: more strict validation functions
p.addOptional('WinSize',75,@(x) isnumeric(x)&&isscalar(x));
p.addOptional('WinStep',50,@(x) isnumeric(x)&&isscalar(x));
% TODO: add validiation function? @(x) any(strcmp(x, {'fminsearch', 'legacy'})))
p.addParameter('Invert',false,@islogical);
p.addParameter('Transform','Radon',@(x) any(strcmp(x, {'Radon', 'Rotate'})));
p.addParameter('Metric','Var',@(x) any(strcmp(x, {'Sep', 'Var'})));
p.addParameter('Optimizer','globalsearch');
p.addParameter('FilterVar', 0, @isfinite);
p.parse(varargin{:});

% TODO: Probably don't need to reassign most of these
I = p.Results.I;
Tfactor = 1/p.Results.msPerLine; % ypixel per ms
Xfactor = p.Results.umPerPx; % microns per xpixel
WinSize = p.Results.WinSize;
WinPixelsDown = p.Results.WinStep;
a = p.Results.FilterVar;

% Force filtering at a = 25 for legacy Radon
if strcmp(p.Results.Transform, 'Radon') && strcmp(p.Results.Optimizer, 'radonlegacy')
    a = 25;
end
    
%% Filter Image (if requested)
% TODO: move this into subfunction?
% This filtering step seems to be very important for radon method in
% general
% TODO: should image be cast to double first?
% TODO: should filter variance be a function of the window size?  
% Create filter kernel once
if a
    %Create high pass flter using isotropic Gaussian
    I_size = size(I);
    y = (1:WinSize)';
    x = 1:size(I,2);
    gaus = 1 - exp(-.5 * ( ((y-0.5*I_size(1)).^2 / a) + ((x-0.5*I_size(2)).^2 / a) ) );
end

%% Loop through lines
% TODO: could be more efficient about using the table to hold these
% variables so not duplicating/calculating every loop
    % Calculate block indices
    last = WinSize:WinPixelsDown:size(I,1);
    first = last - WinSize + 1;
    nWins = length(first);

    % Initialize outputs
    Result = table('Size',[nWins,6],'VariableTypes',{'double','double','duration','duration','double','double'},...
        'VariableNames', {'StartLine','EndLine','StartTime','EndTime','Velocity','Metric'});
    Result.Properties.VariableUnits = {'','','ms','ms','mm/s',''};
    
    % Create waitbar
    % TODO: move this waitbar out to calcLinescanVelTiff?--Create input
    % DisplayWaitbar argument?
    hWait = waitbar(0, 'Calculating linescan velocity', 'Name', 'Linescan');
    
    % Function handle for calculating linescan slope
    calcLinescanAngleFcn = @(block) linescan.calcLinescanAngle(...
        block,'Transform',p.Results.Transform, 'Metric', p.Results.Metric,...
        'Optimizer', p.Results.Optimizer);
    
    % Loop over windows
    for iWin = 1:1:nWins
        % TODO: use im2double instead?
        block = double(I(first(iWin):last(iWin), :));
        
        % Filter section
        if a; block = applyFilt(block, gaus); end
        
%         block = im2double(I(first(iWin):last(iWin), :));
        [angle, metric] = calcLinescanAngleFcn(block);
        
        % Slope is actually opposite tand(thetaMax) because tand provides slope
        % in cartesian coordinates as opposed to image coordinates where Y-axis
        % increases going down, rather than up.
        dXdt = -1/tand(angle);

        % TODO: change this to first, last? Or this is supposed to be time?
        Result.StartLine(iWin) = first(iWin);
        Result.EndLine(iWin) = last(iWin);
        % TODO: should time be average time or start time of block?
        Result.StartTime(iWin) = milliseconds((first(iWin)-1)/Tfactor);
        Result.EndTime(iWin) = milliseconds((last(iWin)-1)/Tfactor);
        
        % Calculate velocity from slope
        vel = dXdt*Xfactor*Tfactor;
        if p.Results.Invert; vel = -vel; end    % Reverse if(invert)
        Result.Velocity(iWin) = vel;
        
        Result.Metric(iWin) = metric;
        
        % Update waitbar
        waitbar(iWin/nWins, hWait);
    end

    % Update waitbar and close
    waitbar(1, hWait, 'Done!')
    close(hWait);
end

function I = applyFilt(I, gaus)
    %Multiply frequency components of image by Gaussian filter
    I = fftshift(fft2(I));
    I = I.*gaus;
    I = ifft2(ifftshift(I));
    I = real(I);
end

% function I = isoGaussFilt(I, a)
%     %Create high pass flter using isotropic Gaussian
%     I_size = size(I);
%     y = (1:I_size(1))';
%     x = 1:I_size(2);
%     gaus = 1 - exp(-.5 * ( ((y-0.5*I_size(1)).^2 / a) + ((x-0.5*I_size(2)).^2 / a) ) );
% 
%     %Multiply frequency components of image by Gaussian filter
%     I = fftshift(fft2(I));
%     I = I.*gaus;
%     I = ifft2(ifftshift(I));
%     I = real(I);
% end