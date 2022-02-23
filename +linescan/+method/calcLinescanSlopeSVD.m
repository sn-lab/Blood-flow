function [dYdt, MaxSep, IRot] = calcLinescanSlopeSVD(varargin)
%CALCLINESCANSLOPESVD Summary of this function goes here
%   Detailed explanation goes here

% Parse inputs
p = inputParser();
p.addRequired('I')
p.addOptional('Optimizer', 'fminbnd', @ischar);
p.parse(varargin{:});
I = p.Results.I;

% TODO: does input data need to be square??


% Transformation parameters
method = 'bilinear';    % interpolation method for rotating image
bbox = 'crop';

% Optimization parameters
% TODO: does this need to be 0:180? or -90:90?
MinTheta = 0;           % Minimum rotation angle
MaxTheta = 90;          % Maximum rotation angle
Theta0 = 45;            % Starting rotation angle
% TODO: tolerance should probably be for velocity, not for angle, that way
% it's linear
ThetaTol = 0.01;
Steps = 50;


% Construct objective function



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


% TODO: Need to modify output based on exitflag or just return the
% exitflag?
%     if ~exitflag
%         angle = NaN; %OR 50?
%         MaxSep = NaN; %OR 0?
%     end


% Calculate velocity using angle and direction
angle = -angle;
% TODO: should bbox be same here?
IRot = imrotate(I, angle, method, bbox);
if isLinescanVertical(IRot) % Linescan has a negative slope
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


% TODO: probably need to restore this to true legacy
function [OptimalTheta, MaxSep, exitflag] = optimizeSVDAngleLegacy(fun, MinTheta, MaxTheta, ThetaTol, Steps)

% Initializations
exitflag = false;
outOfRange = false;
loops = 1;
MinThetaIter = MinTheta;
MaxThetaIter = MaxTheta;
dTheta = (MaxThetaIter - MinThetaIter)/Steps;

while not(exitflag) && loops <= 100 && not(outOfRange)
    % Compute objective function for angle values from MinThetaIter to
    % MaxThetaIter in steps of dTheta
    Angles = MinThetaIter:dTheta:MaxThetaIter;
    nAngles = length(Angles);
    Sep = zeros(1, nAngles);
    % loop for each value of Angles
    for iAngle = 1:nAngles
        Sep(iAngle) = fun(Angles(iAngle));
    end
    
    [MaxSep, Index] = max(Sep);
    OptimalTheta = Angles(Index);
    if Index==1   % rotation is too large and positive
        if MaxThetaIter >= MaxTheta
            outOfRange = true;
        else
            MaxThetaIter = MaxThetaIter + 3*dTheta;
            % TODO: why not increase this as well?
            MinThetaIter = Angles(Index+1);
        end
    elseif Index == nAngles % rotation is too large and negative  
        if MinThetaIter <= MinTheta
            outOfRange = true;
        else
            MinThetaIter = MinThetaIter - 3*dTheta;
            MaxThetaIter = Angles(Index-1);
        end
    else % found a good rotation
        [~,inds] = maxk(Sep, 2);
        if abs(diff(Angles(inds))) < ThetaTol
            exitflag = true; %set flag for exiting loop for window
        else % new angle range
            % TODO: this is going to end up repeating the same angles which
            % is an unfortunate waste of time...
            MaxThetaIter = Angles(Index)+2*dTheta;
            MinThetaIter = Angles(Index)-2*dTheta;
            dTheta = dTheta/2;  % NAR Addition
%             OldSep = MaxSep;
        end
    end
    
    loops = loops + 1;
end % while loop for thetas
end

%% Transformation Functions
function [IRot] = RotateImage(I,theta,method,bbox)
%     Rotdata = RotateWithoutToolbox(XRAMP, YRAMP, X, Y,I,theta,method);
    IRot = imrotate(I, theta, method, bbox);
end

% TODO: radon transform function would go here
% function R = RadonImage(I,theta)
% end

%% Metric Functions
function sep = separability(J)
    S = svd(J);
    % TODO: is there a better weighting scheme here to deal with values
    % outside of original image e.g. set to NaN? These seem to be affecting
    % the separability. i.e. separability decreases with more 0 values
    % (near 45 degrees)
    sep = S(1)^2/sum(S.^2);
end

%% Helper functions
function tf = isLinescanVertical(Rotdata)
    % Check if linescan is vertical
    vertstd = std(mean(Rotdata,1));
    horzstd = std(mean(Rotdata, 2));
    tf = vertstd > horzstd; %lines are horizontal
end
