function [] = velocity_from_tif
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


% Ask to go through new files or not
button = questdlg('New files to look at individually?',...
    'New files','Yes','No', 'Yes');
if strcmp(button,'Yes')
    newfiles = 1;
elseif strcmp(button,'No')
    newfiles = 0;
    keepgoing = 1;
end

% let user choose files and select region of interest
if newfiles == 1; % User wants to enter new files
    errorcheck = 0 ;


    OpenNameTemp = [];
    FnameTemp = [];
    FileTimeTemp = [];
    WinLefts = [];
    WinRights = [];
    Npics = [];
    Tfactors = [];
    Xfactors =[];
    UseAvgs = [];

    Npic = 0;



    % Running parameters
    prompt={'Always subtract average? (N/Y)'};
    def={'N'};
    dlgTitle='Processing parameters';
    lineNo=1;
    answer=inputdlg(prompt,dlgTitle,lineNo,def,'on');
    if strcmp(answer{1},'N')
        alwaysuseavg = 0;
    else
        alwaysuseavg = 1;
        useavg = 1;
    end


    morefiles = 1;
    while morefiles
        % Get file to open         
        pause (0.1);
        [fname,pname] = uigetfile('*.*');
        Openfile = [pname, fname]
        cd(pname);

        fileinfo = imfinfo(Openfile);
        [maxframes, scrap] = size(fileinfo);

        % get time file was created
        info= dir(pname);
        fnames = {};
        times = {};
        nfiles = length(info);
        for i = 1:nfiles
            fnames = strvcat(fnames,char(info(i).name));
            times = strvcat(times,char(info(i).date));
            if strcmp(char(info(i).name), fname)
                filetime = char(info(i).date);
                continue
            end
        end

        
        % get info for each file

        % Calibration factors
        % Calculate Tfactor (number of pixels per ms)
        prompt={'ms per line', 'microns per pixel'};
        def={'2', '1'};
        dlgTitle='Conversion factors';
        lineNo=1;
        answer2=inputdlg(prompt,dlgTitle,lineNo,def,'on'); 
        Tfactor = 1/str2double(cell2mat(answer2(1))); % ypixel per ms
        Xfactor = str2double(cell2mat(answer2(2))); % microns per xpixel



        
        figdisp = figure;

        % USER CHOOSES RELEVANT AREA FOR ANALYSIS
        % show 1 frame at a time
        done = 0;
        framenumber = 1;

        showlines = imread(Openfile,1);
        [numlines, nx] = size(showlines);
        
        if numlines>500;
            showlines = showlines(1:500, :);
        end;
        
        imagesc(showlines); f_niceplot;
        title({[fname];['frame:', num2str(framenumber)]});

        % get coordinates of rrbox
        fignum = gcf;
        xlabel('Select region of interest');
        Roi = round(getrect);
        WinLeft = Roi(1);
        width = Roi(3);
        WinRight = Roi(1) + Roi(3);
        rectangle('Position', [WinLeft, 1, width, numlines],'EdgeColor', 'r');
        xlabel ('Space - keep this box, f-forward, b-back, s-skip forward,  n-new box');
        while not(done);
            waitforbuttonpress;
            fignum = gcf;
            pressed = get(fignum, 'CurrentCharacter');
            if pressed == ' ' % space for keep this box
                done = 1;
            elseif pressed == 'f' % forward 1 frame
                if framenumber < maxframes;
                    framenumber = framenumber +1;
                    showlines = imread(Openfile, framenumber);
                else
                    beep;
                    display ('no more frames')
                end

                imagesc(showlines); f_niceplot;
                title({[fname];['frame:', num2str(framenumber)]});

                rectangle('Position', [WinLeft, 1, width, numlines],'EdgeColor', 'r');
                xlabel ('Space - keep this box, f-forward, b-back, s-skip forward,  n-new box');

            elseif pressed == 'b' % back 1 frame
                if framenumber > 1
                    framenumber = framenumber - 1;
                    showlines = imread(Openfile, framenumber);
                else
                    beep;

                end

                imagesc(showlines); f_niceplot;
                title({[fname];['frame:', num2str(framenumber)]});

                rectangle('Position', [WinLeft, 1, width, numlines],'EdgeColor', 'r');
                xlabel ('Space - keep this box, f-forward, b-back, s-skip forward,  n-new box');


            elseif pressed == 's' % skip 10 frames forward
                if framenumber + 10 <= maxframes;
                    framenumber = framenumber + 10;
                    showlines = imread(Openfile, framenumber);
                else
                    beep;
                end

                imagesc(showlines); f_niceplot;
                title({[fname];['frame:', num2str(framenumber)]});

                rectangle('Position', [WinLeft, 1, width, numlines],'EdgeColor', 'r');
                xlabel ('Space - keep this box, f-forward, b-back, s-skip forward,  n-new box');

            elseif pressed == 'n'
                clf;
                imagesc(showlines); f_niceplot;
                title({[fname];['frame:', num2str(framenumber)]});
                xlabel('Select region of interest');
                Roi = round(getrect)
                WinLeft = Roi(1);
                width = Roi(3);
                WinRight = Roi(1) + Roi(3);
                rectangle('Position', [WinLeft, 1, width, numlines],'EdgeColor', 'r');
                xlabel ('Space - keep this box, f-forward, b-back, s-1skip forward,  n-new box');
            else
                beep;
                display (' not a good key')
            end; %if
        end; %loop while not done


        slope = 2;
        
        % Ask user if subtract average of linescans across from each block of
        % data
        if alwaysuseavg == 0
            button = questdlg('Subtract average across linescans?',...
                'Use average?','Yes','No', 'Yes');
            if strcmp(button,'Yes')
                useavg = 1;
            elseif strcmp(button,'No')
                useavg = 0;
            end
        end;
                
        % ask if more files?
        button = questdlg('More files?',...
            'Continue','Yes','No','Yes');
        if strcmp(button,'Yes')
            morefiles = 1;
        elseif strcmp(button,'No')
            morefiles = 0;
        end
        
        close (figdisp);

        
        % saveinfo
        OpenNameTemp = strvcat(OpenNameTemp, Openfile);
        FnameTemp = strvcat(FnameTemp, fname);
        FileTimeTemp = strvcat(FileTimeTemp,filetime); 
        WinLefts = vertcat(WinLefts, WinLeft);
        WinRights = vertcat(WinRights, WinRight);
        Tfactors = vertcat(Tfactors, Tfactor);
        Xfactors = vertcat(Xfactors, Xfactor);
        UseAvgs = vertcat(UseAvgs, useavg);

    end; % morefiles
    
    
    OpenName = cellstr(OpenNameTemp);
    Fname = cellstr(FnameTemp);
    FileTime = cellstr(FileTimeTemp);
    
    
    %  as a comma delimited text file
    pause(0.1);
    [filename, pathname] = uiputfile('*.csv', 'Comma delimited file save As');
    Datafile = [pathname, filename];
    
    Datafile2 = [pathname, filename, 'bfdata.mat'];
    save(Datafile2, 'OpenName','WinLefts', 'WinRights', 'Tfactors', 'Xfactors','UseAvgs'); 
    
    OpenName = strvcat('OpenName', OpenNameTemp);
    Fname = strvcat('filename', FnameTemp);
    FileTime = strvcat('Time',FileTimeTemp); 
    WinLefts = strvcat('WinLefts',num2str(WinLefts));
    WinRights = strvcat('WinRights', num2str(WinRights)); 
    Tfactors = strvcat('Tfactors', num2str(Tfactors));
    Xfactors = strvcat('Xfactors', num2str(Xfactors));
    UseAvgs = strvcat('UseAvgs', num2str(UseAvgs));

    commas = [];
    [lines, col] = size(UseAvgs);
    for i = 1:lines
        commas(i,1) = ',';
    end
    tosave = horzcat((OpenName), commas, (WinLefts),commas, (WinRights), ...
        commas,  (Tfactors), commas, (Xfactors), commas, (Fname), commas, ...
        (FileTime), commas, (UseAvgs));


    diary(Datafile)
    tosave
    diary off

    button = questdlg('Calculate velocites now?',...
        'Continue','Yes','No','Yes');
    if strcmp(button,'Yes')
        keepgoing = 1;
    elseif strcmp(button,'No')
        keepgoing = 0;
    end
    
end; % if new files



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate velocities or diameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if keepgoing 
    
    % % GET SAVED DATA (optional)
    if exist('Datafile2') == 1
        load (Datafile2);
    else
        pause (0.1);
        [filename3, pathname3] = uigetfile('*.*', 'parameter file (...bfdata.mat)');
        Datafile3 = [pathname3, filename3];
        load(Datafile3);
    end
    [nfiles,z] = size(OpenName);
    
    % For running continuously from setup 
    clear showlines info
    
    
    % Running parameters from user
    prompt={'Number of time pixels per data point', 'time pixels between data points',...
        'number of lines to process', 'Start with file #:', 'turn on display?'};
    def={'75', '50', '50000', '1', 'no'};
    dlgTitle='Processing parameters';
    lineNo=1;
    answer=inputdlg(prompt,dlgTitle,lineNo,def,'on');
    WinSize =  str2double(answer(1));   % actual data used is only center circle ~70% of area (square window)
    WinPixelsDown = str2double(answer(2)); % number of pixels between top of last window and next window
    Maxlines =  str2double(answer(3)); %  total number of lines
    startfilenumber = str2double(answer(4));
    if strcmp(answer(5), 'no')
        errorcheck = 0;
    else
        errorcheck = 1;
    end;

    
    for i =startfilenumber:nfiles %%%%%%%%% Loop through all files
        Openfile2 = char(OpenName(i));
        %Datafile = [char(strrep(OpenName(i),'.tif',['--wpd', num2str(WinPixelsDown)])), date, '.mat'];
        Datafile = [char(strrep(OpenName(i),'.tif',['SVDrawVel', num2str(WinPixelsDown)])),'.mat'];
        
        WinLeft = WinLefts(i,1); % leftmost pixel
        WinRight = WinRights(i,1); % rightmost pixel
        Tfactor = Tfactors(i, 1);
        Xfactor = Xfactors(i, 1);
        UseAvg = UseAvgs(i,1);
        
        
        % read in data
        
        FR1 = 1;
        FRLast = 1; 
        datachunk = [];

        % Loop through lines
        npoints = 0;
        first = 1; 
        last = first+WinSize;
        Result = [];
        Slope = 2;
        
        while last<Maxlines % loop thorugh lines
            lines = f_get_lines_from_tiff(Openfile2, first, last);
            if isempty(lines)
                break;
            end
            [tny, tnx] = size(lines);
            if tny < WinSize
                break
            end
            
            block = lines(:, WinLeft: WinRight);
            if errorcheck == 1;
                veldata = f_find_vel(block, Tfactor, Xfactor, Slope, UseAvg, 1);

            else
                veldata = f_find_vel(block, Tfactor, Xfactor, Slope, UseAvg);
            end
       
        veldata(1) = first;
        veldata(2) = npoints *WinPixelsDown/Tfactor;
        % ---------------------------------------------
            % For Debugging
            if (errorcheck ==1) & (npoints< 20)
                subplot(2,1,1);imagesc(lines); f_niceplot;
                title(Openfile2)
                subplot(2,1,2); imagesc(block); f_niceplot;title('block')
                angle = acot(veldata(3)/Xfactor/Tfactor)*180/pi;
                title([num2str(veldata(1)), ' vel:', num2str(veldata(3)), ' angle: ', num2str(angle)]);
                xlabel('press a key to continue');
                pause;
            end
            % ---------------------------------------------
            veldata(5) = -1*acot(veldata(3)/Xfactor/Tfactor);
          
            Result = vertcat(Result, veldata);
            first = first + WinPixelsDown
            last = first+WinSize;
            npoints = npoints+1;
        end
        
        save(Datafile,'Result', 'Tfactor', 'WinPixelsDown');
        
        clear data data1 cropped Result Rotdata Small ecog shutter ekg BioRad time Stimulus Data;

    end; % Loop for each file
    
    message = 'done'
    beep
    beep
    
end % if keepgoing 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = f_niceplot

axis image; colormap gray;
set(gca, 'XTickLabel', [])
set(gca, 'YTickLabel', [])

% --------------------------------------------------------------
function data = f_get_lines_from_tiff(filename, startline, endline)
% filename - path
% startline - 1st line to get
% endline - last line to get (inclusive)
% OUTPUT:
% data - matrix of values, [] if something is not valid


% temp
% startline = 204800;
% endline = 204900;
% [fname,pname] = uigetfile('*.*');
% filename = [pname, fname];


data = [];

% error check to see if really a file
fid = fopen(filename);
if fid ~= -1
    fclose(fid);
    first = imread(filename);
    [ny, nx] = size(first);
    
    % get first and last frame and line index
    startframe = floor(startline/ny) + 1;
    startinframeline = startline - (startframe-1)*ny;
    
    endframe = floor(endline/ny) + 1;
    endinframeline = endline - (endframe-1)*ny;
    
    

    for nframe = startframe:endframe
        try 
            tempframe = imread(filename, startframe);
        catch
            disp('invalid linenumber');
            break;
            
        end;
        
        if startframe == endframe
            data = tempframe(startinframeline:endinframeline, :);
            
        elseif nframe == startframe
            data = vertcat(data, tempframe(startinframeline:end,:));
        elseif nframe == endframe;
            data = vertcat(data, tempframe(1:endinframeline, :));
        else
            data = vertcat(data, tempframe);
        end;
        
        
    end;
else
    data = [];
    disp('invalid file');
end

data = double(data);


% ---------------------------------------------------------
function Result = f_find_vel(small, varargin)
% function Result = f_find_vel(small, (Tfactor), (Xfactor), (slope), (useaverage), (debug))
% based on RotMatOpt19_rand_time
% IN: small = 1 frame
%               Xfactor (microns/pixel)
%               Tfactor (pixels/ms)
%               slope= 1 for positive, 0 for neg, 2 for automatically pick
%               slope
%               useaverage = 1 to subtract our average across lines, 0 for
%               slow capilliaries; if not used, sets = 1;
%               debug: if exists, will show each frame
% OUT: Result: (preserved data structure)
%       unneeded numbers are = 0
%       column     3) Velocity (mm/s) + veloctiy is in x-dir (RBC's from
%                             left to right)
%                       4) Sep
%                       5) Angle (true angle of data unmodified by this
%                       function, positive is RBC move left to right)
%                       6) Flux
% 07-11-03: Make an option to not use average across the frame. Useful for
% slow capillaries.
%12-20-04: Finds Flux based on method developed empirically: Thresholds
% image of rotated block data by average; Takes an average projection, and
% thresholds; Finds derivative and zerocrossings to find RBC edges. Uses a
% ratio of standard deviation of intensities across time and space to
% reject some data points. 

%%%%TEMP
do_debug = 0;

% Get Tfactor and Xfactor from input parameters
if length(varargin) == 0 % no user parameters
    Xfactor = 205/500*250/512; % microns/pixel
    Tfactor = 1;
    slopeset = 0;
    useaverage = 1
elseif length(varargin) == 1 % user gave at least  Tfactor
    Tfactor= cell2mat(varargin(1));
    Xfactor = 205/500*250/512;
    slopeset = 0;
    useaverage = 1
elseif length(varargin) == 2 % gives Xfactor
    Tfactor= cell2mat(varargin(1));
    Xfactor = cell2mat(varargin(2));
    slopeset = 0;
    useaverage = 1
elseif length(varargin) == 3 % gives a slope
    Tfactor= cell2mat(varargin(1));
    Xfactor = cell2mat(varargin(2));
    slopenum = cell2mat(varargin(3));
    if slopenum ==2
        slopeset = 0;
    else % = 1 or =0;
        slopeset = 1;
    end
    useaverage = 1
elseif length(varargin) ==4
    Tfactor= cell2mat(varargin(1));
    Xfactor = cell2mat(varargin(2));
    slopenum = cell2mat(varargin(3));
    if slopenum ==2
        slopeset = 0;
    else % = 1 or =0;
        slopeset = 1;
    end
    useaverage = cell2mat(varargin(4));
elseif length(varargin) ==5 % use the debugging option
    Tfactor= cell2mat(varargin(1));
    Xfactor = cell2mat(varargin(2));
    slopenum = cell2mat(varargin(3));
    if slopenum ==2
        slopeset = 0;
    else % = 1 or =0;
        slopeset = 1;
    end
    useaverage = cell2mat(varargin(4));
    do_debug = 1;
end

block = small;

% Take out vertical stripes
blocksize= size(block);
avg = mean(block);
avgs = ones([blocksize(1), 1])*avg;
if useaverage
    block = block-avgs;
end
clear avgs, avg;

% Make data square
oldY= blocksize(1);
oldX = blocksize(2);

oldXs = ones(oldY,1)*[1:oldX];
oldYs = transpose(ones(oldX,1)*[1:oldY]);
if oldY > oldX;
    newX = oldY; newY = oldY;
    step = (oldX - 1)/(newX-1);
    Xs = ones(newY,1)*([1:step:oldX]);
    Ys = transpose(ones(newY,1)*[1:newY]);
    small = interp2(oldXs, oldYs,block, Xs, Ys);
    TfactorUse = Tfactor;
    XfactorUse = Xfactor/newX*oldX;
elseif oldY < oldX
    newX = oldX; newY = oldX;
    Xs = ones(newX,1)*[1:newY];
    step = (oldY - 1)/(newY-1);
    Ys = transpose(ones(newY,1)*([1:step:oldY]));
    small= interp2(oldXs, oldYs,block, Xs, Ys);
    TfactorUse = Tfactor*newY/oldY;
    XfactorUse = Xfactor;
    
else
    TfactorUse = Tfactor;
    XfactorUse = Xfactor;
end; % resize block


[WinSize, WinSizeX] = size(small);
% Pre-calculated numbers for RotateFindSVD, etc
MaxXRot = floor(WinSize/sqrt(2));
HalfMaxX = round(MaxXRot/2);
MidSmall = round(WinSize/2);

% PARAMETERS

Steps = 50;
XRAMP = ones(WinSize, 1)*[1:WinSize] - MidSmall;
YRAMP = MidSmall-[1:WinSize]'*ones(1, WinSize);
X = ones(MaxXRot, 1)*[1:MaxXRot] - HalfMaxX;
Y = HalfMaxX-[1:MaxXRot]'*ones(1, MaxXRot);
method = '*linear'; % method for interpolate in rotating image
SepTol = 0.01;


size(small);
% Left over variables from origianal program are set = 0
WinNumber = 0; Nframes = 0; WinPerFrame = 0; WinTop = 0; Period = 0;  WinPixelsDown = 0; 

Slope = 1; % SLOPE IS NOW FOUND AUTOMATICALLY

FoundMax = 0;       
if Slope==1;            
    MinTheta = 0;           % Starting negative value for angles of rotation
    MaxTheta = pi/2;       % Starting positive limit for angles of rotation
elseif Slope == 0 ;      
    MinTheta = -pi/2;        % Starting negative value for angles of rotation
    MaxTheta = 0;   % Starting positive limit for angles of rotation
end;

loops = 1;
OldSep = 0;
while (not(FoundMax))
    dTheta = (MaxTheta - MinTheta)/Steps;
    Sep = zeros(1, Steps+1);
    Angles = zeros(1, Steps+1);
    
    % loop for each value of dTheta
    for count = 1:Steps+1;        
        Angles(count) = MaxTheta - (count-1) * dTheta;
        [Sep(count), Rotdata] = RotateFindSVD(XRAMP, YRAMP, X, Y, small,Angles(count),method);
    end;
    
    [MaxSep, Index] = max(Sep);
    if Index==1;   % rotation is too large and positive
        if MaxTheta >= pi/2
            Result = [Nframes,WinTop, 50, 0, 0,0, (Nframes-1)*Period + 1/1000/TfactorUse*(WinNumber-1)*WinPixelsDown, Nframes,0,0];
            FoundMax = 1;
        else
            MaxTheta = MaxTheta + 3*dTheta;
            MinTheta = Angles(Index+1);
        end
    elseif Index == Steps +1; % rotation is too large and negative  
        if MinTheta <= -pi/2
            Result = [Nframes,WinTop,50, 0, 0,0,(Nframes-1) *Period + 1/1000/TfactorUse*(WinNumber-1)*WinPixelsDown, Nframes,0,0];
            FoundMax = 1;
        else
            MinTheta = MinTheta - 3*dTheta;
            MaxTheta = Angles(Index-1);
      end
else % found a good rotation
      if abs(MaxSep - OldSep)<SepTol;
            [scrap, Rotdata] = RotateFindSVD(XRAMP, YRAMP, X, Y, small,Angles(Index),method);
            
            % check orientation of rotated matrix
            vertavg = mean(Rotdata,1);
            horzavg = mean(Rotdata, 2);
            vertstd = std(vertavg);
            horzstd = std(horzavg);
            if horzstd> vertstd %lines are horizontal
                  vel = 1*TfactorUse*XfactorUse*abs(cot(Angles(Index)));
                  angle = Angles(Index);
                  angletrue =  acot(cot(Angles(Index))*TfactorUse*XfactorUse);  
                  % Threshold Rotdata before getting lineout
                  AvgRotdata = mean(mean(Rotdata));
                  ThreshRotdata = Rotdata > AvgRotdata;
                  LineOut1 = mean(ThreshRotdata,2);
                  xvar = mean(std(Rotdata,0,1));
                  tvar = mean(std(Rotdata,0,2));
            else % lines are vertical
                  vel = -1*TfactorUse*XfactorUse*abs(tan(Angles(Index)));
                  angle = atan(cot(Angles(Index)));
                  angletrue =  -acot(cot(Angles(Index))*TfactorUse*XfactorUse);
                  AvgRotdata = mean(mean(Rotdata));
                  ThreshRotdata = Rotdata > AvgRotdata;
                  LineOut1 = mean(ThreshRotdata,1);
                  xvar = mean(std(Rotdata,0,2));
                  tvar = mean(std(Rotdata,0,1));
            end
            
            % FLUX
            Flux = NaN;
            Flux2 = NaN;
            Flux3 = NaN;


            Result = [Nframes,WinTop, vel, MaxSep, angletrue, Flux,(Nframes-1)*Period + 1/1000/TfactorUse*(WinNumber-1)*WinPixelsDown, WinNumber, Flux2, Flux3];
            FoundMax = 1; %set flag for exiting loop for window
            
        else % new anlge range
            MaxTheta = Angles(Index)+2*dTheta;
            MinTheta = Angles(Index)-2*dTheta;
            OldSep = MaxSep;
        end
    end %if index
    
    loops = 1+ loops;
    if loops > 100
        Result = [Nframes,WinTop, 50, 0, 0, 0,(Nframes-1)*Period + 1/1000/TfactorUse*(WinNumber-1)*WinPixelsDown, Nframes,0,0];
        FoundMax = 1;
        vel =0;
    end;
end % while loop for thetas

% ------------------------------------------------------
function [seperability, Rotdata] = RotateFindSVD(XRAMP, YRAMP, X, Y,small,Theta,method)
%RotateFindSVD - rotates the center square matrix of small, returns seperability
% 090406 changed isnan
warpx = X*cos(Theta) +Y*sin(Theta) ;
                warpy = (-X*sin(Theta)+ Y*cos(Theta)) ;
                Rotdata = interp2(XRAMP, YRAMP, small, warpx, warpy, method);
                Rotdata(isnan(Rotdata))= mean(Rotdata(~isnan(Rotdata)));
                S = svd(Rotdata);
                seperability = S(1)^2/sum(S.^2);



