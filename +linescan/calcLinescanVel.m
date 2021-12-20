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
p.addRequired('I',@ismatrix);
% TODO: these should probably be required
p.addRequired('msPerLine',@(x) isnumeric(x)&&isscalar(x));
p.addRequired('umPerPx',@(x) isnumeric(x)&&isscalar(x));
% TODO: more strict validation functions
p.addOptional('WinSize',75,@(x) isnumeric(x)&&isscalar(x));
p.addOptional('WinStep',50,@(x) isnumeric(x)&&isscalar(x));
p.addOptional('errorcheck',false,@islogical);
p.addParameter('Method','Radon',@(x) any(strcmp(x, {'Radon', 'SVD'})))
% TODO: allow user to flip image for opposite velocity?
% if I_sign==0
%     I=fliplr(I);
% end

p.parse(varargin{:});

% TODO: Probably don't need to reassign most of these
I = p.Results.I;
Tfactor = 1/p.Results.msPerLine; % ypixel per ms
Xfactor = p.Results.umPerPx; % microns per xpixel
WinSize = p.Results.WinSize;
WinPixelsDown = p.Results.WinStep;
errorcheck = p.Results.errorcheck;

%    % actual data used is only center circle ~70% of area (square window)
%  % number of pixels between top of last window and next window


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate velocities or diameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\

    
% TODO: should use this default name or ask user?
    %Datafile = [char(strrep(OpenName(i),'.tif',['--wpd', num2str(WinPixelsDown)])), date, '.mat'];

   
    
    

    % Loop through lines
    
    % Calculate block indices
    last = WinSize:WinPixelsDown:size(I,1);
    first = last - WinSize + 1;
    nWins = length(first);

    % Initialize outputs
    Result = zeros(nWins,6);
    
    % Create waitbar
    hWait = waitbar(0);
    
    % Pick function to calculate velocity
    switch p.Results.Method
        case 'Radon'
            getLinescanVelFcn = @(block) linescan.method.calcLinescanVelRadon(block, Tfactor, Xfactor);
        case 'SVD'
            getLinescanVelFcn = @(block) linescan.method.calcLinescanVelSVD(block, Tfactor, Xfactor);
    end
    
    for iWin = 1:1:nWins
        % TODO: use im2double instead?
        block = double(I(first(iWin):last(iWin), :));
%         block = im2double(I(first(iWin):last(iWin), :));
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
        veldata(6) = -1*acot(veldata(3)/Xfactor/Tfactor);

        Result(iWin, :) = veldata;
        
        waitbar(iWin/nWins, hWait);
    end

    waitbar(1, hWait, 'Done!')
    
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



