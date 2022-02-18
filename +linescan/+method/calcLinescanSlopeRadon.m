function [dXdt, Y] = calcLinescanSlopeRadon(varargin)
% TODO: consider renaming this to calcLinescanAngleRadon?

% This function performs a Radon transform on an image prefiltered by
% Gaussian isotropic filter to determine angle. The angle is then converted
% to velocity in units of um/ms.

% I = Input image
% options.I_sign = 1 for positive slope, 0 for negative, 2 for both

p = inputParser();
p.addRequired('I')
% TODO: are these necessary?
% TODO: these probably should be addRequired (not guessed)
% p.addRequired('Tfactor', @(x) isnumeric(x)&&isscalar(x)); % ms/line
% p.addRequired('Xfactor', @(x) isnumeric(x)&&isscalar(x)); % microns/pixel
% TODO: maybe this should be passed as a function handle to allow user to
% pass other arguments such as options.I_sign and customize optimizer

% TODO: should this be a function handle instead?
p.addOptional('Optimizer', 'fminsearch', @ischar); % microns/pixel
p.addOptional('options', struct(), @isstruct);
p.parse(varargin{:});
% TODO: allow optimset - optimizer settings e.g. options.I_sign
I = p.Results.I;
options = p.Results.options;
% Xfactor = p.Results.Xfactor;
% Tfactor = p.Results.Tfactor;


%Perform Radon transform

% TODO: can we just use dX/dt or are we also worried about tolerance due to
% different Tfactor/Xfactor i.e. will these cause nonlinear????
% TODO: maybe can switch back to just output angle
switch p.Results.Optimizer
    case 'fminbnd'
        [dXdt, Y] = optimizeRadonAngleFminbnd(I, 0.01);
    case 'binarysearch'
        [thetaMax, Y] = optimizeRadonAngleLoop(I, [0,180], 0.01);
        dXdt = -tand(thetaMax-90);
    case 'exhaustive'
%         [thetaMax, Y] = optimizeRadonAngleLegacy(final, [0,180], 0.01);
        % TODO: switch to passing in angle range rather than I_sign
        % TODO: maybe switch this to 'exhaustive' and just use exhaustive
        % with appropriate settings for legacy compatibility
        [dXdt, Y] = optimizeRadonAngleExhaustive(I, options);
    case 'legacy'
%         [thetaMax, Y] = optimizeRadonAngleLegacy(final, [0,180], 0.01);
        % TODO: switch to passing in angle range rather than I_sign
        % TODO: maybe switch this to 'exhaustive' and just use exhaustive
        % with appropriate settings for legacy compatibility
        options.Thetas = 'orig';    % Or 'simple'
        [dXdt, Y] = optimizeRadonAngleExhaustive(I, options);
    otherwise
        error([p.Results.Optimizer ' is not a valid optimizer']);
end
 
%     % Slope is actually opposite tand(thetaMax) because tand provides slope
%     % in cartesian coordinates as opposed to image coordinates where Y-axis
%     % increases going down, rather than up.
%     dXdt = -tand(thetaMax);

end

%% Optimization Functions
% For setting an appropriate tolerance, it makes most sense to perform
% optimization on slope, since the slope is linearly related to velocity,
% however, this is a difficult optimization problem since velocity is
% unconstrained. In reality, we have to set limits on velocity and 
% (-inf,inf)convergence of optimization function is more consistent when restricted
% peformed on

function [dXdt, Y, flag, output] = optimizeRadonAngleFminbnd(I, thetaTol)
    fun = @(theta) -var(radon(I, theta));
    options = optimset('MaxIter', 5000, 'TolX', thetaTol);
    [theta,Y,flag,output] = fminbnd(fun, 0, 180, options);
    Y = -Y;
    % Slope is actually opposite tand(thetaMax) because tand provides slope
    % in cartesian coordinates as opposed to image coordinates where Y-axis
    % increases going down, rather than up.
    dXdt = -tand(theta-90);
end


function [dXdt, Y] = optimizeRadonAngleExhaustive(I, options)
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
    [Y,iTheta] = max(Variance);


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
function [theta, Y] = optimizeRadonAngleRecurse(I, thetaRange, velTol)
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
    [Y,n]=max(Variance);

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
        [theta, Y] = optimizeRadonAngleRecurse(I, newThetaRange, velTol);
    end
end

% Same function as above but using a while loop instead of recursion
function [theta, Y, flag] = optimizeRadonAngleLoop(I, thetaRange, velTol)
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
        [Y,n]=max(Variance);

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

%% Objective Functions
function [theta, Y, R] = radonMaxVar(I, theta)
% TODO: move this into separate "objective function"
    % For radon transform, theta is angle from vertical?
    [R, xp] = radon(I, theta);

    % TODO: this assumes that there is some variance--but if the window is
    % small and there aren't many (or any) clear stripes then that doesn't
    % work. But then again maybe nothing will work in that case.
    % Determine maximum variance along columns
    Variance = var(R);
    [Y,iTheta] = max(Variance);
    

    %Debug mode
    if options.Debug
        displayDebugger(I, R, theta, xp, Variance)
    end
    
    % For cartesian convention. Should be in range (-90, 90)
    thetaMax = theta(iTheta)-90;
end

%% Debugging
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