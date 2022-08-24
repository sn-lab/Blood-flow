
% extractVelTiffShared.m - Processes linescan tiff files to obtain RBC flow speed. User
% inputs microscope parameters (e.g. objective, magnification, frame rate,
% etc.), selects ROI, then inputs analysis parameters (e.g. window size)

% Jan 7th, 2012
% By Nozomi Nishimura and Thom P. Santisakultarm
% Chris Schaffer Lab, Cornell University, Ithaca, NY
% ts386@cornell.edu

% Output:  Result, Tfactor, WinPixelsDown
% below are the columns of Result
% 1)  starting line number (first)
% 2)  time (ms)
% 3)  velocity (mm/s). + veloctiy is in x-dir (RBC's from left to right)
% 4)  angle (true angle of stripes in degrees)

function [Result Tfactor WinPixelsDown] = extractVelTiffShared

clear
close all

disp('')
disp('extractVelTiffShared.m')
disp('')

dontaskslope = 0;  % user specifies slope

%% set up analysis parameters

% Ask to go through new files or not
button = questdlg('New files to look at individually?',...
    'New Files','Yes','No', 'Yes');
if strcmp(button,'Yes')
    newfiles = 1;
elseif strcmp(button,'No')
    newfiles = 0;
    keepgoing = 1;
end

if newfiles == 1; % User wants to enter new files
    errorcheck = 0;
    
    OpenNameTemp = [];
    FnameTemp = [];
    FileTimeTemp = [];
    WinLefts = [];
    WinRights = [];
    Slopes = [];
    NXs = [];
    Tfactors = [];
    Xfactors =[];
    UseAvgs = [];
    maxFrames = [];
    msPreAnalogPixel = [];
    
    % ask if subtract average for all files
    button = questdlg('Always subtract average?','Subtract Average','Yes','No','No');
    if strcmp(button,'No')
        alwaysuseavg = 0;
    else
        alwaysuseavg = 1;
        useavg = 1;
    end
    
    % ask if same acquisition rate and magnification for all files
    button = questdlg('Same acquisition rate and magnification for all files?','Acquisition Rate and Magnification','Yes','No','No');
    if strcmp(button,'No')
        sameAcqMag = 0;
    else
        sameAcqMag = 1;
        
        prompt = {'acquisition rate (ms/line)','magnification (um/pixel)'};
        def = {'0.6','0.1055'};%default values correspond to 3.39 (512x512) frames/s and 20x objective at 2 V scan mirror at 10x in MPScan program
        dlgTitle = 'Processing parameters';
        lineNo = 1;
        answer = inputdlg(prompt,dlgTitle,lineNo,def,'on');
        Tfactor = 1/str2double(cell2mat(answer(1))); % ypixel per ms
        Xfactor = str2double(cell2mat(answer(2))); % microns per xpixel
    end
    
    UseAna = 0;%no analog signals for this version of extractVelTiffShared.m
    
    morefiles = 1;
    while morefiles
        % Get file to open
        [fname,pname] = uigetfile('*.*');
        Openfile = [pname, fname];
        fprintf('showing:  %s\n',Openfile)
        cd(pname);
        
        % get time file was created
        fileinfo = imfinfo(Openfile);
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
        
        showlines = imread(Openfile,1);
        NX = fileinfo(1).Width;
        NY = fileinfo(1).Height;
        %% USER CHOOSES RELEVANT AREA FOR ANALYSIS
        % show 1 frame at a time
        done = 0;
        framenumber = 1;
        numlines = NY;
        maxframes = numel(fileinfo);
        imagesc(showlines); f_niceplot;
        title({fname;['frame:', num2str(framenumber),'/', num2str(maxframes)]});
        
        % get coordinates of rrbox
        fignum = gcf;
        Roi = round(getrect);
        WinLeft = Roi(1);
        width = Roi(3);
        WinRight = Roi(1) + Roi(3);
        rectangle('Position', [WinLeft, 1, width, numlines],'EdgeColor', 'r');
        xlabel ('Space - keep this box, f-forward, b-back, s-skip forward,  n-newbox');
        
        while not(done);
            waitforbuttonpress;
            fignum = gcf;
            pressed = get(fignum, 'CurrentCharacter');
            
            if pressed == ' ' % space for keep this box
                done = 1;
            elseif pressed == 'f' % forward 1 frame
                if framenumber < maxframes;
                    framenumber = framenumber +1;
                    showlines = imread(Openfile,framenumber);
                else
                    beep
                    display ('no more frames')
                end
                imagesc(showlines); f_niceplot;
                title({fname;['frame:', num2str(framenumber),'/', num2str(maxframes)]});
                
                rectangle('Position', [WinLeft, 1, width, numlines],'EdgeColor', 'r');
                xlabel ('Space - keep this box, f-forward, b-back, s-skip forward,  n-newbox');
                
            elseif pressed == 'b' % back 1 frame
                if framenumber > 1
                    framenumber = framenumber - 1;
                    showlines = imread(Openfile,framenumber);
                else
                    beep
                end
                
                imagesc(showlines); f_niceplot;
                title({fname;['frame:', num2str(framenumber),'/', num2str(maxframes)]});
                
                rectangle('Position', [WinLeft, 1, width, numlines],'EdgeColor', 'r');
                xlabel ('Space - keep this box, f-forward, b-back, s-skip forward,  n-newbox');
                
            elseif pressed == 's' % skip 10 frames forward
                if framenumber + 10 <= maxframes;
                    framenumber = framenumber + 10;
                    showlines = imread(Openfile,framenumber);
                else
                    beep;
                end
                
                imagesc(showlines); f_niceplot;
                title({fname;['frame:', num2str(framenumber),'/', num2str(maxframes)]});
                
                rectangle('Position', [WinLeft, 1, width, numlines],'EdgeColor', 'r');
                xlabel ('Space - keep this box, f-forward, b-back, s-skip forward,  n-newbox');
                
            elseif pressed == 'n'
                clf;
                imagesc(showlines); f_niceplot;
                title({fname;['frame:', num2str(framenumber),'/', num2str(maxframes)]});
                
                Roi = round(getrect);
                WinLeft = Roi(1);
                width = Roi(3);
                WinRight = Roi(1) + Roi(3);
                rectangle('Position', [WinLeft, 1, width, numlines],'EdgeColor', 'r');
                xlabel ('Space - keep this box, f-forward, b-back, s-1skip forward,  n-newbox');
            else
                beep;
                display (' not a good key')
            end
        end %loop while not done
        
        % Ask user for slope
        if dontaskslope == 0
            button = questdlg('slope', 'Slope of lines', ...
                'positive', 'negative', 'both', 'both');
            if strcmp(button, 'positive')
                slope = 1;
            elseif strcmp(button, 'negative')
                slope = 0;
            else
                slope = 2;
            end
        end
        
        % Ask user if subtract average of linescans across from each block of data
        if alwaysuseavg == 0
            button = questdlg('Subtract average across linescans?',...
                'Use average','Yes','No', 'Yes');
            if strcmp(button,'Yes')
                useavg = 1;
            elseif strcmp(button,'No')
                useavg = 0;
            end
        end;
        
        % Ask user for acquisition rate and magnification
        if ~sameAcqMag
            prompt = {'acquisition rate (ms/line)','magnification (um/pixel)'};
            def = {'0.6','0.1'};
            dlgTitle = 'Processing parameters';
            lineNo = 1;
            answer = inputdlg(prompt,dlgTitle,lineNo,def,'on');
            Tfactor = 1/str2double(cell2mat(answer(1))); % ypixel per ms
            Xfactor = str2double(cell2mat(answer(2))); % microns per xpixel
        end
        
        % ask if more files?
        button = questdlg('More files?',...
            'Continue','Yes','No','Yes');
        if strcmp(button,'Yes')
            morefiles = 1;
        elseif strcmp(button,'No')
            morefiles = 0;
        end
        
        %% save parameter file
        OpenNameTemp = strvcat(OpenNameTemp, Openfile);
        FnameTemp = strvcat(FnameTemp, fname);
        FileTimeTemp = strvcat(FileTimeTemp,filetime);
        WinLefts = vertcat(WinLefts, WinLeft);
        WinRights = vertcat(WinRights, WinRight);
        Slopes = vertcat(Slopes, slope);
        NXs = vertcat(NXs, NX);
        Tfactors = vertcat(Tfactors, Tfactor);
        Xfactors = vertcat(Xfactors, Xfactor);
        UseAvgs = vertcat(UseAvgs, useavg);
        maxFrames = vertcat(maxFrames,maxframes);
        
    end % morefiles
    
    OpenName = cellstr(OpenNameTemp);
    Fname = cellstr(FnameTemp);
    FileTime = cellstr(FileTimeTemp);
    
    % save as a comma delimited text file
    [filename, pathname] = uiputfile('*.csv', 'Comma delimited file save As');
    Datafile = [pathname, filename];
    
    Datafile2 = strrep(Datafile, '.csv', 'Parameters');
    save(Datafile2, 'OpenName', 'Slopes','WinLefts','WinRights','NXs', 'Tfactors', 'Xfactors','UseAvgs', 'maxFrames','UseAna');
    
    OpenName = strvcat('OpenName', OpenNameTemp);
    Fname = strvcat('filename', FnameTemp);
    FileTime = strvcat('Time',FileTimeTemp);
    WinLefts = strvcat('WinLefts',num2str(WinLefts));
    WinRights = strvcat('WinRights', num2str(WinRights));
    Slopes = strvcat('Slope',num2str(Slopes));
    NXs = strvcat('NX', num2str(NXs));
    Tfactors = strvcat('Tfactors', num2str(Tfactors));
    Xfactors = strvcat('Xfactors', num2str(Xfactors));
    UseAvgs = strvcat('UseAvgs', num2str(UseAvgs));
    maxFrames = strvcat('maxFrames', num2str(maxFrames));
    
    [lines, col] = size(Slopes);
    commas = char(44*ones(lines,1));
    
    tosave = horzcat((OpenName), commas, (NXs), commas, (WinLefts),commas, (WinRights), commas, (Slopes),commas, (Tfactors), commas, (Xfactors) , commas, (Fname), commas, (FileTime), commas, (UseAvgs), commas, (maxFrames));
    
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
    
end % if new files

%% Calculate velocities
if keepgoing
    errorcheck = 0;
    
    % GET SAVED DATA (optional)
    if exist('Datafile2','var') == 1
        load (Datafile2);
    else
        [filename3, pathname3] = uigetfile('*.*', 'parameter file');
        Datafile3 = [pathname3, filename3];
        load(Datafile3);
    end
    
    [nfiles,z] = size(OpenName);
    
    % For running continuously from setup
    clear showlines info
    
    % Running parameters from user
    prompt = {'WinPixelsDown', 'WinSize', 'Maxlines', 'Start with file #:'};
    def = {'50', '100', 'inf', '1'};
    dlgTitle = 'Processing parameters';
    lineNo = 1;
    answer = inputdlg(prompt,dlgTitle,lineNo,def,'on');
    WinPixelsDown = str2double(answer(1)); % number of pixels between top of last window and next window
    WinSize =  str2double(answer(2));   % EVEN NUMBER please! actual data used is only center circle ~70% of area (square window)
    Maxlines =  str2double(answer(3)); % TEMP total number of lines
    startfilenumber = str2double(answer(4));
    
    for i =startfilenumber:nfiles %%%%%%%%% Loop through all files
        Openfile2 = char(OpenName(i));
        
        addonName = [' rawVel ', num2str(WinPixelsDown), num2str(WinSize), '.mat'];
        Datafile = char(strrep(strrep(strrep(strrep(OpenName(i),'.tiff',addonName),'.TIFF',addonName),'.tif',addonName),'.TIF',addonName));
        
        fprintf('processing:  %s \n',Datafile)
        
        Slope = Slopes(i,1);
        WinLeft = WinLefts(i,1); % leftmost pixel
        WinRight = WinRights(i,1); % rightmost pixel
        NX = NXs(i,1);   % Pixels per frame in x
        Tfactor = Tfactors(i, 1);
        Xfactor = Xfactors(i, 1);
        UseAvg = UseAvgs(i,1);
        maxframes = maxFrames(i,1);
        
        FR1 = 1;
        FRLast = 1;
        datachunk = [];
        
        % Loop through lines
        npoints = 0;
        first = 1;
        last = first+WinSize;
        Result = [];
        
        prevstr = '';
        while last<Maxlines % loop thorugh lines
            lines = f_get_lines_from_tiff(Openfile2, first, last);
            
            if isempty(lines)
                break
            end
            
            [tny, tnx] = size(lines);
            
            if tny < WinSize
                break
            end
            
            block = lines(:, WinLeft: WinRight);
            veldata = f_find_vel_radon_shared1(block, Tfactor, Xfactor, Slope, UseAvg, 1); %corresponds to 2V mirror
            
            veldata(1) = first;
            veldata(2) = npoints *WinPixelsDown/Tfactor;
            
            %             ---------------------------------------------
            %             For Debugging
            %
            %             if (errorcheck ==1) && (npoints< 20)
            %                 subplot(4,1,1);imagesc(lines); f_niceplot;
            %                 title(Openfile)
            %                 subplot(4,1,2); imagesc(block); f_niceplot;title('block')
            %                 angle = acot(veldata(3)/Xfactor/Tfactor)*180/pi;
            %                 title([num2str(veldata(1)), ' vel:', num2str(veldata(3)), ' angle: ', num2str(angle)]);
            %                 Slope
            %                 pause
            %             end
            %
            %             ---------------------------------------------
            
            Result = vertcat(Result, veldata);
            first = first + WinPixelsDown;
            last = first+WinSize;
            npoints = npoints+1;
            str = ['Analyzed line ' num2str(first) ];
            refreshdisp(str,prevstr,first);
            prevstr = str;
        end

        save(Datafile,'Result', 'Tfactor', 'WinPixelsDown');
        
        clear Result
    end; % Loop for each file
    
    disp('done')
    beep
    beep
    
end % if keepgoing
end

%% f_find_vel_radon

% This function performs a Radon transform on an image prefiltered by
% Gaussian isotropic filter to determine angle. The angle is then converted
% to velocity in units of um/ms.

% Spring 2008
% By Nathan R. Cornelius
% Chris Schaffer Lab, Cornell University, Ithaca, NY
% ts386@cornell.edu

% I = Input image
% a = Gaussian filter variance (defaults to 25)
% I_sign = 1 for positive slope, 0 for negative, 2 for both
% Debug = 0 enables debug mode in which various images are displayed for
% each velocity calculation

% Modification History:
% 09/30/2008 Puifai changed to mirror voltage = 2 V by default
% 01/10/2010 Puifai added user option to subtract average

function Result = f_find_vel_radon_shared1(I, Tfactor, Xfactor, I_sign, UseAvg, Debug);

%tic
%error(nargchk(1,6,nargin));
if nargin<6, Debug=1; end
if nargin<5, UseAvg=0; end
if nargin<4, I_sign=1; end
if nargin<3, Tfactor=1; end
if nargin<3, Xfactor=1; end

if I_sign==0,
    I=fliplr(I);
end

I_size = size(I);
I = im2double(I);

%% subtract average value from image to take out vertical stripes

[r c] = size(I);
avgPixValues = ones(r,1)*mean(I);

if UseAvg
    I = I-avgPixValues;
end
clear avgPixValues

%% extract velocity by radon transform method


%figure(1), imshow(I);
%I=makeNoise(I);
%figure(2), imshow(I);

filtVar=25;     %Filter variance set to 25
a=filtVar;

%Create high pass flter using isotropic Gaussian
for x=1:I_size(1)
    for y=1:I_size(2)
        gaus(x,y) = 1 - exp(-.5 * ( ((x-0.5*I_size(1))^2 / a) + ((y-0.5*I_size(2))^2 / a) ) );
    end
end

%Multiply frequency components of image by Gaussian filter
I_fft = fftshift(fft2(I));
filtered = I_fft.*gaus;
final = ifft2(ifftshift(filtered));
final = real(final);

%Perform Radon transform
if I_sign~=2,
    theta = 90+radonAngles(181,4);
else    
    theta1 = 90 - radonAngles(181,4);
    theta1 = fliplr(theta1);
    theta2 = 90 + radonAngles(181,4);
    theta = [theta1 theta2];
end

%if I_sign==1,
    %theta = 90+linspace(0,90,181);
    %theta = 90+radonAngles(181);
%else
    %theta=linspace(0,90,181);
    %theta = radonAngles(181);
%end
[R, xp] = radon(final, theta);

%Determine maximum variance along columns
Variance=var(R);
[Y,n]=max(Variance);


%Debug mode
Debug = 1 ;
if Debug==0
    %Display the filtered image
    figure(2)
    imagesc(final);
    
    %Display the Radon transform image
    figure(3)
    iptsetpref('ImshowAxesVisible','on')
    imshow(R, [], 'Xdata',theta,'Ydata',xp,'InitialMagnification','fit')
    colormap(hot), colorbar
    ylabel('x'''); xlabel('\theta');
    
    %Show the plot of the variance
    figure(4)
    plot(Variance);
    
%     figure(5)
%     IR=iradon(R,theta);
%     imshow(IR,[]);
    %Wait for user
    pause;  
end
 

if I_sign==1
    thetaMax = theta(n)-90;
    velocity = -1*(1/tan(deg2rad(thetaMax)))*Tfactor*Xfactor;
elseif I_sign==0
    thetaMax = (theta(n)-90)*-1;
    velocity = (1/tan(deg2rad(-1*thetaMax)))*Tfactor*Xfactor;
else
    if theta(n)<90
        if n-1<1,
            n=2;
        end
        thetaMax = theta(n-1)-90;
    else
        thetaMax = theta(n)-90;
    end
    velocity = (1/tan(deg2rad(-1*thetaMax)))*Tfactor*Xfactor;
end

%toc

Result = [0, 0, velocity, thetaMax];

end

%% f_niceplot
function [] =f_niceplot

axis image; colormap gray;
set(gca, 'XTickLabel', [])
set(gca, 'YTickLabel', [])

end

%% f_get_lines_from_tiff
function data = f_get_lines_from_tiff(filename, startline, endline)
% filename - path
% startline - 1st line to get
% endline - last line to get (inclusive)
% OUTPUT:
% data - matrix of values, [] if something is not valid


% % temp
% startline = 4096;
% endline = startline+30;
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
    startframe = ceil(startline/ny);
    startinframeline = rem(startline-1, ny)+1;
    
    endframe = ceil(endline/ny);
    endinframeline = rem(endline-1, ny)+1;
    
    

    for nframe = startframe:endframe
        try 
%             tempframe = imread(filename, startframe);
            tempframe = imread(filename, nframe);
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
end

%% radonAngles
function theta = radonAngles(numInc,x);

numInc=x*numInc;
theta(1)=deg2rad(1);
thetaLast=deg2rad(89);
velocityHigh=1/tan(theta(1));
velocityLow=1/tan(thetaLast);

velInc = (velocityHigh - velocityLow)/numInc;

for i=1:x*173
    thetaNew=atan(1 / ((1/tan(theta(i))) - velInc));
    theta(i+1)=thetaNew;
end

theta = rad2deg(theta);

n=1;
for i=22:1/x:89,
    thetaLin(n)=i;
    n=n+1;
end

theta = [theta, thetaLin];
end


function refreshdisp(str,prevstr,iteration)

if ~exist('iteration','var')
    iteration=2;
end

if iteration==1
    fprintf(str)
else
    fprintf(char(8*ones(1,length(prevstr))));
    fprintf(str);
end

end