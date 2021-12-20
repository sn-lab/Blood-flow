function Result = calcLinescanVelSVD(varargin)
% function Result = f_find_vel(small, (Tfactor), (Xfactor), (slope), (useaverage), (debug))
% based on RotMatOpt19_rand_time
% IN: small = 1 frame
%               Xfactor (microns/pixel)
%               Tfactor (pixels/ms)
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

% TODO: setup input parser to handle optional inputs.
p = inputParser();
p.addRequired('small')
% TODO: these probably should be addRequired (not guessed)
p.addOptional('Tfactor', 1, @(x) isnumeric(x)&&isscalar(x)); % microns/pixel
p.addOptional('Xfactor', 205/500*250/512, @(x) isnumeric(x)&&isscalar(x)); % microns/pixel
p.parse(varargin{:});

small = p.Results.small;
Xfactor = p.Results.Xfactor;
Tfactor = p.Results.Tfactor;

block = small;

% Make data square
blocksize = size(block);
oldY = blocksize(1);
oldX = blocksize(2);

% TODO: Try imresize instead??
[oldXs,oldYs] = meshgrid(1:oldX,1:oldY);
if oldY > oldX
    newX = oldY; newY = oldY;
    step = (oldX - 1)/(newX-1);
    
    [Xs, Ys] = meshgrid(1:step:oldX, 1:newY);
    small = interp2(oldXs, oldYs, block, Xs, Ys);
    TfactorUse = Tfactor;
    XfactorUse = Xfactor/newX*oldX;
elseif oldY < oldX
    newX = oldX; newY = oldX;
    step = (oldY - 1)/(newY-1);
    
    [Xs, Ys] = meshgrid(1:newY, 1:step:oldY);
    small = interp2(oldXs, oldYs, block, Xs, Ys);
    TfactorUse = Tfactor*newY/oldY;
    XfactorUse = Xfactor;
    
else
    TfactorUse = Tfactor;
    XfactorUse = Xfactor;
end % resize block


[WinSize, ~] = size(small);
% Pre-calculated numbers for RotateFindSVD, etc
MaxXRot = floor(WinSize/sqrt(2));
HalfMaxX = round(MaxXRot/2);
MidSmall = round(WinSize/2);

% PARAMETERS
[XRAMP, YRAMP] = meshgrid((1:WinSize) - MidSmall, MidSmall-(1:WinSize));
[X, Y] = meshgrid((1:MaxXRot) - HalfMaxX, HalfMaxX-(1:MaxXRot));
method = '*linear'; % method for interpolate in rotating image
MinTheta = 0;           % Starting negative value for angles of rotation
MaxTheta = pi/2;       % Starting positive limit for angles of rotation
SepTol = 0.01;
Steps = 50;

% Find angle of maximum separability
fun = @(Theta) 1-(RotateFindSVD(XRAMP, YRAMP, X, Y,small,Theta,method));
[angle, MaxSep] = optimizeWithToolbox(fun, MinTheta, MaxTheta, SepTol, Steps);
% fun = @(Theta) RotateFindSVD(XRAMP, YRAMP, X, Y,small,Theta,method);
% [angle, MaxSep] = optimizeWithoutToolbox(fun, MinTheta, MaxTheta, SepTol, Steps);

% Determine whether maximum yields horizontal or vertical stripes
Rotdata = Rotate(XRAMP, YRAMP, X, Y, small, angle, method);

% Calculate velocity using angle and direction
if isLinescanHorizontal(Rotdata)
    vel = 1*TfactorUse*XfactorUse*abs(cot(angle));
    angletrue =  acot(cot(angle)*TfactorUse*XfactorUse);
else % lines are vertical
    vel = -1*TfactorUse*XfactorUse*abs(tan(angle));
    angletrue =  -acot(cot(angle)*TfactorUse*XfactorUse);
end

% TODO: what is this number?: (Nframes-1)*Period + 1/1000/TfactorUse*(WinNumber-1)*WinPixelsDown
% Result = [Nframes,WinTop, vel, MaxSep, angletrue, Flux,(Nframes-1)*Period + 1/1000/TfactorUse*(WinNumber-1)*WinPixelsDown, WinNumber, Flux2, Flux3];
% Left over variables from origianal program are set = 0
WinNumber = 0; Nframes = 0; WinPerFrame = 0; WinTop = 0; Period = 0;  WinPixelsDown = 0;
Flux = NaN; Flux2 = NaN; Flux3 = NaN;
Result = [Nframes, WinTop, vel, angletrue, MaxSep];


function [angle, MaxSep] = optimizeWithToolbox(fun, MinTheta, MaxTheta, SepTol, Steps)
%     fun = @(Theta) 1-(RotateFindSVD(XRAMP, YRAMP, X, Y,small,Theta,method));
    % TODO: set MaxFunEvals and MaxIter?
    % TODO: change the standard here in optimizeWithoutToolbox
    Steps = Steps*100;
    options = optimset('MaxIter', Steps, 'TolX', SepTol);
    [angle,MinSep,exitflag,output] = fminbnd(fun, MinTheta, MaxTheta, options);
    MaxSep = 1-MinSep;
    % TODO: Need to check exitflag, etc?
%     if ~exitflag
%         angle = NaN; %OR 50?
%         MaxSep = NaN; %OR 0?
%     end
        
    
function [MaxTheta, MaxSep] = optimizeWithoutToolbox(fun, MinTheta, MaxTheta, SepTol, Steps)
FoundMax = 0;

loops = 1;
OldSep = 0;
while (not(FoundMax))
    dTheta = (MaxTheta - MinTheta)/Steps;
    % TODO: below two lines could be moved out of while loop for efficiency
    Sep = zeros(1, Steps+1);
    Angles = zeros(1, Steps+1);
    
    % loop for each value of dTheta
    for count = 1:Steps+1       
        Angles(count) = MaxTheta - (count-1) * dTheta;
        Sep(count) = fun(Angles(count));
    end
    
    [MaxSep, Index] = max(Sep);
    if Index==1   % rotation is too large and positive
        if MaxTheta >= pi/2
%             MaxTheta = NaN;
%             MaxSep = 0;
%             Result = [Nframes,WinTop, 50, 0, 0,0, (Nframes-1)*Period + 1/1000/TfactorUse*(WinNumber-1)*WinPixelsDown, Nframes,0,0];
            FoundMax = 1;
        else
            MaxTheta = MaxTheta + 3*dTheta;
            MinTheta = Angles(Index+1);
        end
    elseif Index == Steps +1 % rotation is too large and negative  
        if MinTheta <= -pi/2
%             MaxTheta = NaN;
%             MaxSep = 0;
%             Result = [Nframes,WinTop,50, 0, 0,0,(Nframes-1) *Period + 1/1000/TfactorUse*(WinNumber-1)*WinPixelsDown, Nframes,0,0];
            FoundMax = 1;
        else
            MinTheta = MinTheta - 3*dTheta;
            MaxTheta = Angles(Index-1);
        end
else % found a good rotation
        if abs(MaxSep - OldSep)<SepTol
            FoundMax = 1; %set flag for exiting loop for window
        else % new angle range
            MaxTheta = Angles(Index)+2*dTheta;
            MinTheta = Angles(Index)-2*dTheta;
            OldSep = MaxSep;
        end
    end %if index
    
    loops = 1+ loops;
    if loops > 100
%             MaxTheta = NaN;
%             MaxSep = 0;
%         Result = [Nframes,WinTop, 50, 0, 0, 0,(Nframes-1)*Period + 1/1000/TfactorUse*(WinNumber-1)*WinPixelsDown, Nframes,0,0];
        FoundMax = 1;
    end
end % while loop for thetas
        

function tf = isLinescanHorizontal(Rotdata)
    % check orientation of rotated matrix
    vertavg = mean(Rotdata,1);
    horzavg = mean(Rotdata, 2);
    vertstd = std(vertavg);
    horzstd = std(horzavg);
    tf = horzstd> vertstd; %lines are horizontal
    
% ------------------------------------------------------
function [seperability, Rotdata] = RotateFindSVD(XRAMP, YRAMP, X, Y,small,Theta,method)
%RotateFindSVD - rotates the center square matrix of small, returns seperability
% 090406 changed isnan
    Rotdata = Rotate(XRAMP, YRAMP, X, Y,small,Theta,method);
    S = svd(Rotdata);
    seperability = S(1)^2/sum(S.^2);

function Rotdata = Rotate(XRAMP, YRAMP, X, Y,small,Theta,method)
    warpx = X*cos(Theta) +Y*sin(Theta) ;
    warpy = (-X*sin(Theta)+ Y*cos(Theta)) ;
    Rotdata = interp2(XRAMP, YRAMP, small, warpx, warpy, method);
    Rotdata(isnan(Rotdata))= mean(Rotdata(~isnan(Rotdata)));




