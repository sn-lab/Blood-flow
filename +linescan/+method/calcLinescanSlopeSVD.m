function [dYdt, MaxSep] = calcLinescanSlopeSVD(varargin)
% What does this do?

p = inputParser();
p.addRequired('I')
% TODO: are these necessary?
% TODO: these probably should be addRequired (not guessed)
% p.addRequired('Tfactor', @(x) isnumeric(x)&&isscalar(x)); % ms/line
% p.addRequired('Xfactor', @(x) isnumeric(x)&&isscalar(x)); % microns/pixel
p.addOptional('Optimizer', 'fminbnd', @ischar);
p.parse(varargin{:});

I = p.Results.I;
% Xfactor = p.Results.Xfactor;
% Tfactor = p.Results.Tfactor;

% Make data square?
% TODO: Try imresize instead??
% TODO: does it need to be square??
% sz = size(I);
% [Xq, Yq] = meshgrid(linspace(1,sz(2),max(sz)), linspace(1,sz(1),max(sz)));
% I = interp2(I, Xq, Yq);
% 
% % TODO: won't this return the same number??? doesn't matter?
% [WinSize, ~] = size(I);
% % Pre-calculated numbers for RotateFindSVD, etc
% MaxXRot = floor(WinSize/sqrt(2));
% HalfMaxX = round(MaxXRot/2);
% MidSmall = round(WinSize/2);

% OPTIMIZATION PARAMETERS
% [XRAMP, YRAMP] = meshgrid((1:WinSize) - MidSmall, MidSmall-(1:WinSize));
% [X, Y] = meshgrid((1:MaxXRot) - HalfMaxX, HalfMaxX-(1:MaxXRot));

sz = size(I);
[XRAMP, YRAMP] = meshgrid(1:sz(2), 1:sz(1));
X = XRAMP;
Y = YRAMP;

% method = '*linear'; % method for interpolate in rotating image
method = 'bilinear';
% TODO: does this need to be 0:180? or -90:90?
MinTheta = 0;           % Starting negative value for angles of rotation
MaxTheta = 90;       % Starting positive limit for angles of rotation
% TODO: tolerance should probably be for velocity, not for angle, that way
% it's linear
SepTol = 0.01;
Steps = 50;


% Find angle of maximum separability
switch p.Results.Optimizer
    case 'fminbnd'
        % TODO: should this just be -RotateFindSVD?? Why 1??
        fun = @(Theta) 1-(RotateFindSVD(XRAMP, YRAMP, X, Y,I,Theta,method));
        [angle, MaxSep] = optimizeSVDAngleFminbnd(fun, MinTheta, MaxTheta, SepTol, Steps);
    case 'legacy'
        fun = @(Theta) RotateFindSVD(XRAMP, YRAMP, X, Y,I,Theta,method);
        [angle, MaxSep] = optimizeSVDAngleLegacy(fun, MinTheta, MaxTheta, SepTol, Steps);
end

% Calculate velocity using angle and direction
% TODO: why isn't this using tan for both?
% TODO: can at least shrink this down to just if clause for flipping signs
% TODO: what's the difference betwen vel and angletrue?
dYdt = tand(angle);

% Determine whether angle of maximized SVD yields horizontal or vertical stripes
% Rotdata = RotateWithToolbox(XRAMP, YRAMP, X, Y, I, angle, method);
Rotdata = RotateWithToolbox(I, angle, method);
if ~isLinescanHorizontal(Rotdata)
    % TODO: is there more complicated math involved here?
    dYdt = -dYdt;
end

end
% TODO: what is this number?: (Nframes-1)*Period + 1/1000/TfactorUse*(WinNumber-1)*WinPixelsDown
% Result = [Nframes,WinTop, vel, MaxSep, angletrue, Flux,(Nframes-1)*Period + 1/1000/TfactorUse*(WinNumber-1)*WinPixelsDown, WinNumber, Flux2, Flux3];
% Left over variables from origianal program are set = 0
% WinNumber = 0; Nframes = 0; WinPerFrame = 0; WinTop = 0; Period = 0;  WinPixelsDown = 0;
% Flux = NaN; Flux2 = NaN; Flux3 = NaN;

% TODO: necessary to return Nframes and WinTop?
% Nframes = 0;
% WinTop = 0;
% Result = [Nframes, WinTop, vel, angletrue, MaxSep];

% TODO: should vel be renamed dYdt?
% Result = [vel, angletrue, MaxSep];


% if isLinescanHorizontal(Rotdata)
% %     vel = 1*TfactorUse*XfactorUse*abs(cot(angle));
% %     angletrue =  acot(cot(angle)*TfactorUse*XfactorUse);
%     vel = abs(cot(angle));
%     angletrue =  acot(cot(angle));
% else % lines are vertical
% %     vel = -1*TfactorUse*XfactorUse*abs(tan(angle));
% %     angletrue =  -acot(cot(angle)*TfactorUse*XfactorUse);
%     vel = -abs(tan(angle));
%     angletrue =  -acot(cot(angle));
% end

%% Optimization Functions
function [angle, MaxSep] = optimizeSVDAngleFminbnd(fun, MinTheta, MaxTheta, SepTol, Steps)
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
end


% TODO: rename this to legacy
function [MaxTheta, MaxSep, flag] = optimizeSVDAngleLegacy(fun, MinTheta, MaxTheta, SepTol, Steps)
FoundMax = 0;
flag = false;

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
        if MaxTheta >= 90
%             MaxTheta = NaN;
%             MaxSep = 0;
%             Result = [Nframes,WinTop, 50, 0, 0,0, (Nframes-1)*Period + 1/1000/TfactorUse*(WinNumber-1)*WinPixelsDown, Nframes,0,0];
            FoundMax = 1;
        else
            MaxTheta = MaxTheta + 3*dTheta;
            MinTheta = Angles(Index+1);
        end
    elseif Index == Steps +1 % rotation is too large and negative  
        if MinTheta <= -90
%             MaxTheta = NaN;
%             MaxSep = 0;
%             Result = [Nframes,WinTop,50, 0, 0,0,(Nframes-1) *Period + 1/1000/TfactorUse*(WinNumber-1)*WinPixelsDown, Nframes,0,0];
            FoundMax = 1;
        else
            MinTheta = MinTheta - 3*dTheta;
            MaxTheta = Angles(Index-1);
        end
else % found a good rotation
    flag = true;
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
end
        
%% Helper functions
function tf = isLinescanHorizontal(Rotdata)
    % check orientation of rotated matrix
    vertavg = mean(Rotdata,1);
    horzavg = mean(Rotdata, 2);
    vertstd = std(vertavg);
    horzstd = std(horzavg);
    tf = horzstd> vertstd; %lines are horizontal
end
% ------------------------------------------------------
function [seperability, Rotdata] = RotateFindSVD(XRAMP, YRAMP, X, Y,I,theta,method)
%RotateFindSVD - rotates the center square matrix of small, returns seperability
% 090406 changed isnan
%     Rotdata = RotateWithoutToolbox(XRAMP, YRAMP, X, Y,I,theta,method);
    Rotdata = RotateWithToolbox(I, theta, method);
    S = svd(Rotdata);
    seperability = S(1)^2/sum(S.^2);
end

function Rotdata = RotateWithoutToolbox(I, theta, method)
    % TODO: this seems to be adding an extra pixel of padding


    T = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
    sz = size(I);
    
    %     outputSize = getOutputBound(Rz,sz);
    xLimitsIn = [0.5 sz(2)+0.5];
    yLimitsIn = [0.5 sz(1)+0.5];
    % TODO: add/subtract 0.5????
%     xLimitsIn = [1, sz(2)];
%     yLimitsIn = [1, sz(1)];
    xCentered = xLimitsIn - mean(xLimitsIn);
    yCentered = yLimitsIn - mean(yLimitsIn);
    



    % TODO: can probably eliminate middle column
%     u = [xLimitsIn(1), mean(xLimitsIn), xLimitsIn(2)];
%     v = [yLimitsIn(1), mean(yLimitsIn), yLimitsIn(2)];
    
    % Form grid of boundary points and internal points used by
    % findbounds algorithm.
%     [U,V] = meshgrid(u,v);
    [U,V] = meshgrid(xCentered, yCentered);
    
    % TODO: is there a way to do this with array multiplication
    % Transform corner points to find boundaries
    X = T(1,1).*U + T(2,1).*V;
    Y = T(1,2).*U + T(2,2).*V;

    % XLimitsOut/YLimitsOut are formed from min and max of transformed points.
    XWorldLimitsOut = [min(X(:)), max(X(:))];
    YWorldLimitsOut = [min(Y(:)), max(Y(:))];
    
            
    % Use ceil to provide grid that will accomodate world limits at roughly
    % the target resolution. Assumes that the resolution is 1.
    numCols = ceil(diff(XWorldLimitsOut));
    numRows = ceil(diff(YWorldLimitsOut));

    outputSize = [numRows numCols];

    
    % Pad image for smooth edges
    I = [zeros(sz(1),1), I, zeros(sz(1),1)];
    I = [zeros(1,sz(2)+2); I; zeros(1,sz(2)+2)];
    
    
    sz = size(I);
    x = [1, sz(2)];
    y = [1, sz(1)];
    xCentered = x - mean(x);
    yCentered = y - mean(y);
    [X, Y] = meshgrid(xCentered(1):xCentered(2), yCentered(1):yCentered(2));
    
    
    [Xq, Yq] = meshgrid(XWorldLimitsOut(1):XWorldLimitsOut(2), YWorldLimitsOut(1):YWorldLimitsOut(2));
    % TODO: what is maximum size?
    
    % TODO: can probably clean this up by doing element-wise multiplication
    % rather than array multiplication
%     theta = -theta;
%     T = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
%     Xq = T(1,1).*Xq + T(2,1).*Yq;
%     Yq = T(1,2).*Xq + T(2,2).*Yq;
    
    Rq = T*[reshape(Xq, 1, []); reshape(Yq, 1, [])];
    Xq = reshape(Rq(1,:), size(Xq));
    Yq = reshape(Rq(2,:), size(Yq));
    
    % TODO: the input image needs to be padded with zeros for same behavior
    % as imrotate
    Rotdata = interp2(X, Y, I, Xq, Yq, method, 0);

    % TODO: this neeeds to rotate around center i.e. need to have X and Y query
    % and index vectors neg to pos.
    %     warpx = X*cosd(theta) +Y*sind(theta) ;
    %     warpy = (-X*sind(theta)+ Y*cosd(theta)) ;
    %     Rotdata = interp2(XRAMP, YRAMP, small, warpx, warpy, method);
    %     Rotdata(isnan(Rotdata))= mean(Rotdata(~isnan(Rotdata)));
    
%         Rotdata = imrotate(I, theta);
end


function Rotdata = RotateWithToolbox(I, theta, method)
% Use 'bilinear' for best results?
    Rotdata = imrotate(I, theta, method, 'crop');
end


