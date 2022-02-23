function [dXdt, fval, IRot] = calcLinescanSlope(varargin)
% TODO: consider renaming this to calcLinescanAngle?
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here


% I = Input image
% options.I_sign = 1 for positive slope, 0 for negative, 2 for both

p = inputParser();
p.addRequired('I')
% TODO: maybe this should be passed as a function handle to allow user to
% pass other arguments such as options.I_sign and customize optimizer

% TODO: add validiation function? @(x) any(strcmp(x, {'fminsearch', 'legacy'})))
% TODO: should this be a function handle instead? Better validation?
p.addParameter('Transform','Radon',@(x) any(strcmp(x, {'Radon', 'Rotate'})));
p.addParameter('Metric','Var',@(x) any(strcmp(x, {'Sep', 'Var'})));
p.addParameter('Optimizer','fminbnd');
p.addOptional('options', struct(), @isstruct);
p.parse(varargin{:});
% TODO: allow optimset - optimizer settings e.g. options.I_sign
I = p.Results.I;
options = p.Results.options;

% Transformation parameters
% TODO: reconsider min/max theta
method = 'bilinear';    % interpolation method for rotating image
bbox = 'crop';
MaxTheta = 90;          % Maximum rotation angle
switch p.Results.Transform
    case 'Radon'
        TransformFun = @(angle) RadonImage(I,angle);
        MinTheta = -90;
    case 'Rotate'
        TransformFun = @(angle) RotateImage(I,angle,method,bbox);
        % TODO: does input data need to be square??
        MinTheta = 0;           % Minimum rotation angle
end

% TODO: change the function to allow handle inputs so don't have to do this
% silly switch case
switch p.Results.Metric
    case 'Var'
        MetricFun = @variability;
    % TODO: should this be renamed to separability?
    case 'Sep'
        MetricFun = @separability;
end


% Optimization parameters
% TODO: does this need to be 0:180? or -90:90?
Theta0 = 45;            % Starting rotation angle
% TODO: tolerance should probably be for velocity, not for angle, that way
% it's linear
ThetaTol = 0.01;
Steps = 50;

% Objective function is inverted because these otpimizers are minimizers
ObjectiveFun = @(angle) -MetricFun(TransformFun(angle));

% Optimization parameters
% TODO: construct objective function

% TODO: can we just use dX/dt or are we also worried about tolerance due to
% different Tfactor/Xfactor i.e. will these cause nonlinear????
% TODO: maybe can switch back to just output angle
% TODO: pull slope calculation out of all objective functions???
switch p.Results.Optimizer
    case 'fminbnd'
        [angle, fval] = optimizeSVDAngleFminbnd(ObjectiveFun, MinTheta, MaxTheta, ThetaTol, Steps);
    case 'svdlegacy'
        [angle, fval] = optimizeSVDAngleLegacy(ObjectiveFun, MinTheta, MaxTheta, ThetaTol, Steps);
    case 'binarysearch'
        [angle, fval] = optimizeRadonAngleLoop(ObjectiveFun,I, [MinTheta,MaxTheta], ThetaTol);
    case 'exhaustive'
        [angle, fval] = optimizeRadonAngleExhaustive(ObjectiveFun,I, options);
    case 'radonlegacy'
        options.Thetas = 'orig';    % Or 'simple'
        [angle, fval] = optimizeRadonAngleExhaustive(ObjectiveFun,I, options);
    case 'globalsearch'
        nTrialPts = 45;
        nStage1pts = 9;
        [angle, fval] = optimizeSVDAngleGlobalSearch(ObjectiveFun, Theta0, MinTheta, MaxTheta, ThetaTol, nTrialPts, nStage1pts);
    case 'multistart'
        k = 15;
        [angle, fval] = optimizeSVDAngleMultiStart(ObjectiveFun, Theta0, MinTheta, MaxTheta, ThetaTol, k);
    otherwise
        error([p.Results.Optimizer ' is not a valid optimizer']);
end

% Revert function value to positive
fval = -fval;
    
% TODO: Need to modify output based on exitflag or just return the
% exitflag?
%     if ~exitflag
%         angle = NaN; %OR 50?
%         MaxSep = NaN; %OR 0?
%     end


% TODO: should bbox be same here?
% TODO: should image be rotated even if SVD wasn't used? Should this
% function really output the rotated image?
IRot = imrotate(I, angle, method, bbox);

% Deal with "Rotate" Tranform which only requires angles (0,90) but needs a
% check to see if output is vertical or horizontal
if strcmp(p.Results.Transform, 'Rotate')
    angle = -angle;
    if isLinescanVertical(IRot) % Linescan has a negative slope
        angle = angle + 90;
    end
else
    angle = angle - 90;
end

% Slope is actually opposite tand(thetaMax) because tand provides slope
% in cartesian coordinates as opposed to image coordinates where Y-axis
% increases going down, rather than up.
dXdt = -tand(angle);

end








%% Optimization Functions
% For setting an appropriate tolerance, it makes most sense to perform
% optimization on slope, since the slope is linearly related to velocity,
% however, this is a difficult optimization problem since velocity is
% unconstrained. In reality, we have to set limits on velocity and...
% (-inf,inf)convergence of optimization function is more consistent when restricted
% peformed on

function [angle, Y, exitflag] = optimizeSVDAngleFminbnd(fun, MinTheta, MaxTheta, ThetaTol, Steps)
%     fun = @(Theta) 1-(RotateFindSVD(XRAMP, YRAMP, X, Y,small,Theta,method));
    % TODO: set MaxFunEvals and MaxIter?
    % TODO: change the standard here in optimizeWithoutToolbox
    % TODO: move "Steps" out?
    Steps = Steps*100;
    options = optimset('MaxIter', Steps, 'TolX', ThetaTol);
    [angle,Y,exitflag] = fminbnd(fun, MinTheta, MaxTheta, options);
end


% (SVD) Optimization Functions
function [angle, MaxSep, exitflag] = optimizeSVDAngleGlobalSearch(fun, Theta0, MinTheta, MaxTheta, ThetaTol, nTrialPts, nStage1pts)
    % TODO: these parameters can probably be tuned for faster and more robust optimization
    rng default % For reproducibility
%     opts = optimoptions(@fmincon,'Algorithm','sqp');
    problem = createOptimProblem('fmincon','objective',...
    fun,'x0',Theta0,'lb',MinTheta,'ub',MaxTheta);
    gs = GlobalSearch('XTolerance', ThetaTol, 'NumTrialPoints', nTrialPts, 'NumStageOnePoints', nStage1pts);
    [angle,MinSep,exitflag,~] = run(gs,problem);
    % TODO: move this out?
    MaxSep = 1-MinSep;
end

function [angle, fval, exitflag] = optimizeSVDAngleMultiStart(fun, Theta0, MinTheta, MaxTheta, ThetaTol, k)
    rng default % For reproducibility
%     opts = optimoptions(@fmincon,'Algorithm','sqp');
    problem = createOptimProblem('fmincon','objective',...
    fun,'x0',Theta0,'lb',MinTheta,'ub',MaxTheta);
    ms = MultiStart('XTolerance', ThetaTol);
    [angle,fval,exitflag] = run(ms,problem, k);
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
    
    [MaxSep, Index] = min(Sep);
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
        [~,inds] = mink(Sep, 2);
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


function [dXdt, fval] = optimizeRadonAngleExhaustive(I, options)
% This function performs the radon transform on a precalculated list of
% angles (approximately evenly spaced in velocity units), and only afer
% picks the angle from that list that that produced the radon transform
% with maximum variance.

    % I_sign improves efficiency (2x) at the expense of enforcing a
    % unidrectional flow constraint. Image is only radon transformed over
    % one quadrant of angles e.g. [0,90) for I_sign = 1 (positive slope) or
    % (-90, 0] for I_sign = -1 (negative slope).
    
    % I_sign = -1   % slope is negative
    % I_sign = 0    % slope is unknown/unspecified
    % I_sign = 1    % slope is positive
    
    
    % Initialize theta vector for radon transform
    switch options.thetas
        case 'orig'
            thetaHalf = radonAngles(181,4);
        case 'simple'
            % TODO: should this range be extended? 
            thetaHalf = radonAnglesSimple(tand([0, 89]), 0.06);
    end
    
    
    % NOTE: it is not critical for theta to be sorted... not worth running
    % the fliplr?
    switch options.I_sign
        case -1
            theta = 90-fliplr(thetaHalf);
        % TODO: should this be switched to NaN?
        case 0
            theta = [90-fliplr(thetaHalf), 90+thetaHalf];
        case 1
            theta = 90+thetaHalf;
    end
    
    
    % TODO: move this into separate "objective function"
    % For radon transform, theta is angle from vertical?
    [R, xp] = radon(I, theta);

    % TODO: this assumes that there is some variance--but if the window is
    % small and there aren't many (or any) clear stripes then that doesn't
    % work. But then again maybe nothing will work in that case.
    % Determine maximum variance along columns
    Variance = var(R);
    [fval,iTheta] = max(Variance);


    %Debug mode
    if options.Debug
        displayDebugger(R, theta, xp, Variance)
    end
    
    % For cartesian convention. Should be in range (-90, 90)
    thetaMax = theta(iTheta)-90;
    % Slope is actually opposite tand(thetaMax) because tand provides slope
    % in cartesian coordinates as opposed to image coordinates where Y-axis
    % increases going down, rather than up.
    dXdt = -tand(thetaMax);

    
    
    % TODO: these should return angles prepared properly for radon
    % transform i.e. range (0, 180)
    % TODO: is it worth setting a speed limit?
    % TODO: what is appropriate thetaTol?
    
    % TODO: there should be a faster way to do this at the specified
    % precision i.e. linearly spaced along velocity
    % TODO: make this persistent so only has to be run once?
    % numIng
    % TODO: what are these function inputs??
    function theta = radonAngles(numInc,x)
        numInc=x*numInc;
        velocityHigh = tand(89);
        velocityLow = tand(1);
        
        % TODO: is this close enough to original or should go back to true
        % legacy compatibility???
        velInc = (velocityHigh - velocityLow)/numInc;
        theta = atand(1./(velocityHigh:-velInc:(velocityHigh-velInc*x*173)));

        thetaLin = 22:1/x:89;

        theta = [theta, thetaLin];
    end

    % thetas = radonAnglesSimple(tand([1, 89]), 0.06);
    % thetas = radonAnglesSimple([57.290, 0.0175]),0.06);
    % TODO: velRange should be divisible by velTol so endpoints are
    % included
    function thetas = radonAnglesSimple(velRange, velTol)        
        % TODO: should this be corrected for radon? i.e. atand(-vel)+90
%         vels = velRange(1):velTol:velRange(2);
        vels = velRange(2):-velTol:velRange(1);
        thetas = atand(-vels)+90;
    end

end


%% (Expermental) Optimization Functions
% ```radon``` is O(theta) i.e. linear relationship between runtime and
% length of theta vector. So, try to minimize the number of thetas that
% radon is evaluated for by doing some sort of optimization. This function
% essentially uses a binary search to optimize theta.
% TODO: tolerance should be in velocity units, not theta units....
function [theta, fval] = optimizeRadonAngleRecurse(I, thetaRange, velTol)
    % |----'----!----'----|
    % Check quarters (') of range to narrow down range to one half of
    % original, either [|, !], or [!, |]. If first quarter has a greater
    % variance than the third quarter, new range is [|, !], otherwise it is
    % [!, |].
    midTheta = mean(thetaRange);
%     iTol = (midTheta-thetaRange(1))/2;
    thetas = [mean([thetaRange(1), midTheta]), mean([midTheta, thetaRange(2)])];
    [R, xp] = radon(I, thetas);

    % TODO: this assumes that there is some variance--but if the window is
    % small and there aren't many (or any) clear stripes then that doesn't
    % work. But then again maybe nothing will work in that case.
    %Determine maximum variance along columns
    Variance=var(R);
    [fval,n]=max(Variance);

    % Determine newThetaRange, which is half of the original thetaRange
    if n == 1
        newThetaRange = [thetaRange(n), midTheta];
    else
        newThetaRange = [midTheta, thetaRange(n)];
    end

    % Calculate difference in velocity for 
    % iTol = abs(diff(tand(newThetaRange)));
    iTol = min( tand(newThetaRange)-tand(mean(newThetaRange)) );
    % TODO: does this need to be if/else or just if?
    if abs(iTol) < velTol
        theta = thetas(n);
    else
        [theta, fval] = optimizeRadonAngleRecurse(I, newThetaRange, velTol);
    end
end

% Same function as above but using a while loop instead of recursion
function [theta, fval, flag] = optimizeRadonAngleLoop(I, thetaRange, velTol)
    % |----'----!----'----|
    % Check quarters (') of range to narrow down range to one half of
    % original, either [|, !], or [!, |]. If first quarter has a greater
    % variance than the third quarter, new range is [|, !], otherwise it is
    % [!, |].
    i = 1;
    maxIter = 5000;
    iTol = Inf;
    % Loop until tolerance is less than requested, or maximum number of
    % iterations is exceeded.
    while iTol > velTol && i <= maxIter
        midTheta = mean(thetaRange);
        thetas = [mean([thetaRange(1), midTheta]), mean([midTheta, thetaRange(2)])];
        [R, xp] = radon(I, thetas);

        % TODO: this assumes that there is some variance--but if the window is
        % small and there aren't many (or any) clear stripes then that doesn't
        % work. But then again maybe nothing will work in that case.
        %Determine maximum variance along columns
        Variance=var(R);
        [fval,n]=max(Variance);

        % Determine newThetaRange, which is half of the original thetaRange
        theta = thetas(n);
        if n == 1
            thetaRange = [thetaRange(n), midTheta];
        else
            thetaRange = [midTheta, thetaRange(n)];
        end

        % Calculate difference in velocity for 
%         iTol = min(tand(thetaRange)-tand(theta));
        iTol = abs(diff(tand(thetaRange)));
        % TODO: does this need to be if/else or just if?
        % TODO: this may be less readable than a while loop
        i = i+1;
    end
    
    flag = i <= maxIter;
    % TODO: if i > maxIter, set the flag;
end



%% Transformation Functions
function IRot = RotateImage(I,theta,method,bbox)
    IRot = imrotate(I, theta, method, bbox);
end

% TODO: radon transform function would go here
function R = RadonImage(I,theta)
    R = radon(I, theta);
end

%% Metric Functions
function sep = separability(J)
    S = svd(J);
    % TODO: is there a better weighting scheme here to deal with values
    % outside of original image e.g. set to NaN? These seem to be affecting
    % the separability. i.e. separability decreases with more 0 values
    % (near 45 degrees)
    sep = S(1)^2/sum(S.^2);
end

% TODO: define this so peforms nanmean and nanvar across both dimensions,
% checks which is greater (or less) and returns index-1 wheree 0/1 would
% signify vertical/horizontal
function vari = variability(J)
    vari = -var(J);
end

%% Helper functions
function tf = isLinescanVertical(Rotdata)
    % Check if linescan is vertical
    vertstd = std(mean(Rotdata,1));
    horzstd = std(mean(Rotdata, 2));
    tf = vertstd > horzstd; %lines are horizontal
end

% TODO: generalize this to allow other optimization functions to use
function displayDebugger(I, R, theta, xp, Variance)
    %Display the filtered image
    subolot(1,3,1)
    imagesc(I);

    %Display the Radon transform image
    subplot(1,3,2)
    iptsetpref('ImshowAxesVisible','on')
    imshow(R, [], 'Xdata',theta,'Ydata',xp,'InitialMagnification','fit')
    colormap(hot), colorbar
    ylabel('x'''); xlabel('\theta');

    %Show the plot of the variance
    subplot(1,3,3)
    plot(Variance);

%     figure(5)
%     IR=iradon(R,theta);
%     imshow(IR,[]);
    %Wait for user
    pause; 
end
