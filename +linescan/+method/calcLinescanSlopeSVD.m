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
method = 'bilinear';    % interpolation method for rotating image
bbox = 'crop';

% TODO: does this need to be 0:180? or -90:90?
% MinTheta = 0;           % Starting negative value for angles of rotation
MinTheta = 0;           % Starting negative value for angles of rotation
MaxTheta = 90;       % Starting positive limit for angles of rotation
Theta0 = 45;

% TODO: tolerance should probably be for velocity, not for angle, that way
% it's linear
% TODO: this is used as the tolerance for differences in the separabilty in
% the legacy function which is inconsistent with how it's use in other
% optimizers. Should redefine to be ThetaTol...

ThetaTol = 0.01;
Steps = 50;


% Find angle of maximum separability
switch p.Results.Optimizer
    case 'fminbnd'
        % TODO: should this just be -RotateFindSVD?? Why 1??
        fun = @(Theta) 1-(RotateSVDSep(I,Theta,method,bbox));
        [angle, MaxSep] = optimizeSVDAngleFminbnd(fun, MinTheta, MaxTheta, ThetaTol, Steps);
    case 'legacy'
        fun = @(Theta) RotateSVDSep(I,Theta,method,bbox);
        [angle, MaxSep] = optimizeSVDAngleLegacy(fun, MinTheta, MaxTheta, ThetaTol, Steps);
    case 'globalsearch'
        fun = @(Theta) 1-(RotateSVDSep(I,Theta,method,bbox));
        [angle, MaxSep] = optimizeSVDAngleGlobalSearch(fun, Theta0, MinTheta, MaxTheta, ThetaTol, Steps);
    case 'multistart'
        fun = @(Theta) 1-(RotateSVDSep(I,Theta,method,bbox));
        [angle, MaxSep] = optimizeSVDAngleMultiStart(fun, Theta0, MinTheta, MaxTheta, ThetaTol, Steps);
end

% TODO: Need to check exitflag, etc?
%     if ~exitflag
%         angle = NaN; %OR 50?
%         MaxSep = NaN; %OR 0?
%     end

% Calculate velocity using angle and direction
angle = -angle;

% TODO: return IRot?
% TODO: should bbox be same here?
IRot = imrotate(I, angle, method, bbox);
if isLinescanVertical(IRot)
    % Linescan has a negative slope
    % TODO: is there more complicated math involved here?
    angle = angle + 90;
end
dYdt = tand(angle);

end


%% Optimization Functions
function [angle, MaxSep, exitflag] = optimizeSVDAngleGlobalSearch(fun, Theta0, MinTheta, MaxTheta, ThetaTol)
    % TODO: these parameters can probably be tuned for faster and more robust optimization
    rng default % For reproducibility
%     opts = optimoptions(@fmincon,'Algorithm','sqp');
    problem = createOptimProblem('fmincon','objective',...
    fun,'x0',Theta0,'lb',MinTheta,'ub',MaxTheta);
    gs = GlobalSearch('XTolerance', ThetaTol, 'NumTrialPoints', 45, 'NumStageOnePoints', 9);
    [angle,MinSep,exitflag,~] = run(gs,problem);
    MaxSep = 1-MinSep;
end

function [angle, MaxSep, exitflag] = optimizeSVDAngleMultiStart(fun, Theta0, MinTheta, MaxTheta, ThetaTol)
    rng default % For reproducibility
%     opts = optimoptions(@fmincon,'Algorithm','sqp');
    problem = createOptimProblem('fmincon','objective',...
    fun,'x0',Theta0,'lb',MinTheta,'ub',MaxTheta);
    ms = MultiStart('XTolerance', ThetaTol);
    [angle,MinSep,exitflag] = run(ms,problem, 45);
    MaxSep = 1-MinSep;
end

function [angle, MaxSep, exitflag] = optimizeSVDAngleFminbnd(fun, MinTheta, MaxTheta, ThetaTol, Steps)
%     fun = @(Theta) 1-(RotateFindSVD(XRAMP, YRAMP, X, Y,small,Theta,method));
    % TODO: set MaxFunEvals and MaxIter?
    % TODO: change the standard here in optimizeWithoutToolbox
    Steps = Steps*100;
    options = optimset('MaxIter', Steps, 'TolX', ThetaTol);
    [angle,MinSep,exitflag] = fminbnd(fun, MinTheta, MaxTheta, options);
    MaxSep = 1-MinSep;
end



function [OptimalTheta, MaxSep, exitflag] = optimizeSVDAngleLegacy(fun, MinTheta, MaxTheta, ThetaTol, Steps)

% Initializations
exitflag = false;
loops = 1;
% OldSep = 0;
dTheta = (MaxTheta - MinTheta)/Steps;

while not(exitflag) && loops <= 100
    % TODO: below two lines could be moved out of while loop for efficiency
    Angles = MinTheta:dTheta:MaxTheta;
    
    nAngles = length(Angles);
    Sep = zeros(1, nAngles);
    % loop for each value of Angles
    for iAngle = 1:nAngles
        Sep(iAngle) = fun(Angles(iAngle));
    end
    
    [MaxSep, Index] = max(Sep);
    if Index==1   % rotation is too large and positive
        if MaxTheta >= 90
            % TODO: what is this case handling and why is it setting the
            % FoundMax flag to true?
%             MaxTheta = NaN;
%             MaxSep = 0;
%             Result = [Nframes,WinTop, 50, 0, 0,0, (Nframes-1)*Period + 1/1000/TfactorUse*(WinNumber-1)*WinPixelsDown, Nframes,0,0];
            exitflag = true;
        else
            MaxTheta = MaxTheta + 3*dTheta;
            % TODO: why not increase this as well?
            MinTheta = Angles(Index+1);
        end
    elseif Index == nAngles % rotation is too large and negative  
        if MinTheta <= -90
            % TODO: what is this case handling and why is it setting the
%             MaxTheta = NaN;
%             MaxSep = 0;
%             Result = [Nframes,WinTop,50, 0, 0,0,(Nframes-1) *Period + 1/1000/TfactorUse*(WinNumber-1)*WinPixelsDown, Nframes,0,0];
            exitflag = true;
        else
            MinTheta = MinTheta - 3*dTheta;
            MaxTheta = Angles(Index-1);
        end
else % found a good rotation
%         if abs(MaxSep - OldSep)<ThetaTol
        [~,inds] = maxk(Sep, 2);
        if abs(diff(Angles(inds))) < ThetaTol
            % TODO: this makes no fucking sense because if you run the same
            % angle twice, you'll get a separation of 0 even though you
            % neveer checked other points...
            % TODO: matybe this should be assigned at the beginning of the
            % while loop
            OptimalTheta = Angles(Index);
            exitflag = true; %set flag for exiting loop for window
            % TODO: this should be setting the MaxTheta to Angles(Index)
        else % new angle range
            % TODO: this is going to end up repeating the same angles which
            % is an unfortunate waste of time...
            MaxTheta = Angles(Index)+2*dTheta;
            MinTheta = Angles(Index)-2*dTheta;
            dTheta = dTheta/2;  % NAR Addition
%             OldSep = MaxSep;
        end
    end %if index
    
    loops = loops + 1;
end % while loop for thetas
end
        
%% Objective Function
function [seperability, IRot] = RotateSVDSep(I,theta,method,bbox)
% function [seperability, Rotdata] = RotateFindSVD(XRAMP, YRAMP, X, Y,I,theta,method)
%RotateFindSVD - rotates I, calculates SVD, and returns seperability
%     Rotdata = RotateWithoutToolbox(XRAMP, YRAMP, X, Y,I,theta,method);
    IRot = imrotate(I, theta, method, bbox);
    S = svd(IRot);
    % TODO: is there a better weighting scheme here to deal with values
    % outside of original image e.g. set to NaN? These seem to be affecting
    % the separability. i.e. separability decreases with more 0 values
    % (near 45 degrees)
    seperability = S(1)^2/sum(S.^2);
end

%% Helper functions
function tf = isLinescanVertical(Rotdata)
    % Check if linescan is vertical
    vertstd = std(mean(Rotdata,1));
    horzstd = std(mean(Rotdata, 2));
    tf = vertstd > horzstd; %lines are horizontal
end
