function [dYdt, Y] = calcLinescanVelRadon(varargin)
%% TODO: change this to calcLinescanAngleRadon
%% TODO: should this be renamed dXdt?

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

% Modification History:
% 09/30/2008 Puifai changed to mirror voltage = 2 V by default
% 01/10/2010 Puifai added user option to subtract average

%tic
%error(nargchk(1,6,nargin));
p = inputParser();
p.addRequired('I')
% TODO: these probably should be addRequired (not guessed)
% p.addOptional('Tfactor', 1, @(x) isnumeric(x)&&isscalar(x)); % microns/pixel
% p.addOptional('Xfactor', 205/500*250/512, @(x) isnumeric(x)&&isscalar(x)); % microns/pixel
% TODO: maybe this should be passed as a function handle to allow user to
% pass other arguments such as I_sign and customize optimizer

% TODO: should this be a function handle instead?
p.addOptional('Optimizer', 'fminsearch', @ischar); % microns/pixel
p.parse(varargin{:});
% TODO: allow optimset - optimizer settings e.g. I_sign
I = p.Results.I;
% Xfactor = p.Results.Xfactor;
% Tfactor = p.Results.Tfactor;
I_sign = 0;


% if I_sign==0
%     I=fliplr(I);
% end

%% extract velocity by radon transform method


%Perform Radon transform
% TODO: this should be moved out to filter full image all at once so don't
% have to recalculate kernel on each call of this function
% TODO: this should probably be done in a recursive way. Start with large
% steps e.g. 1 or 2 degree steps, find max, then reduce range and increment
% size
% Can't the sign be determined by whether the angle of maximum variance of
% the radon transform is more than 90 or less than 90? (negative if more
% than 90 and postive if less?
% Okay, so I_sign allows you to force all the input data to have a positive
% slope by flipping the image left to right if it has a negative slope,
% that way, you know the max will be in the range (90, 180] and therefore,
% only have to do the transform over a smaller range (one quadrant rather
% than two)
% However, this forces an assumption of uni-directional flow which is
% likely the case but it's possible that flow slightly reverses in some
% cases like stalls

% TODO: can we just use dY/dt or are we also worried about tolerance due to
% different Tfactor/Xfactor i.e. will these cause nonlinear????
switch p.Results.Optimizer
    case 'fminsearch'
        % TODO: change this to dY/dt
        [dYdt, Y] = optimizeRadonAngleFminsearch(I, 0.01);
    case 'fminbnd'
        [dYdt, Y] = optimizeRadonAngleFminbnd(I, 0.01);
    case 'binarysearch'
        [thetaMax, Y] = optimizeRadonAngleRecurse(I, [0,180], 0.01);
        dYdt = tand(thetaMax);
   case 'fminunc'
        [dYdt, Y] = optimizeRadonAngleFminunc(I, 0.01);
    case 'legacy'
%         [thetaMax, Y] = optimizeRadonAngleLegacy(final, [0,180], 0.01);
        % TODO: switch to passing in angle range rather than I_sign
        [dYdt, Y] = optimizeRadonAngleLegacy(I, 2);
    otherwise
        error([p.Results.Optimizer ' is not a valid optimizer']);
end


% %Debug mode
% Debug = 0;
% if Debug
%     showRadon();
% end
 

% TODO: why is this necessary? Can't you just do theta over range of [0,
% 180) and then take tand(maxTheta). On interval [0,90) -> velocity is
% positive, on interval (90,180) -> velocity is negative
% thetaMax = theta(n);
% velocity = tand(thetaMax)*Tfactor*Xfactor;
% velocity = velocity/Tfactor/Xfactor;
% thetaMax = NaN;
% TODO: always return dYdt? instead of theta?

% Result = [dYdt, thetaMax, Y];
% Result = [dYdt, Y];
end


%% radonAngles
% maxVel = 60;
% minVel = -60;
% tol = 0.1;
% thetas = [atand(0:tol:maxVel), atand(minVel:tol:0)+180];
% % OR
% maxSpeed = 60;
% thetaHalf = atand(0:tol:maxSpeed);
% thetas = [thetaHalf, -fliplr(thetaHalf)+180];

function showRadon()
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

%% Optimization Functions
% TODO: why isn't this using fminbnd?
% TODO: figure out why this isn't working well

% For setting an appropriate tolerance, it makes most sense to perform
% optimization on slope, since the slope is linearly related to velocity,
% however, this is a difficult optimization problem since velocity is
% unconstrained. In reality, we have to set limits on velocity and 
% (-inf,inf)convergence of optimization function is more consistent when restricted
% peformed on 
function [dYdt, Y, flag] = optimizeRadonAngleFminsearch(I, slopeTol)
    fun = @(dYdt) -RadonVarFromSlope(I, dYdt);
    options = optimset('MaxIter', 5000, 'TolX', slopeTol);
    [dYdt,Y,flag,output] = fminsearch(fun, 0, options);
    Y = -Y;
end
    
function [dYdt, Y, flag] = optimizeRadonAngleFminunc(I, slopeTol)
    fun = @(dYdt) -RadonVarFromSlope(I, dYdt);
%     dYdt = optimvar('dYdt');
%     prob = optimproblem('Objective', fun(dYdt));
%     x0.dYdt = 0;
% TODO: what is appropriate starting point for fminunc? 0 seems like would
% never occur in practice
    [dYdt,Y,flag,output] = fminunc(fun,0);
    Y = -Y;
end

% TODO: what is appropriate thetaTol?
function [dYdt, Y, flag] = optimizeRadonAngleFminbnd(I, thetaTol)
    fun = @(theta) -RadonVarFromAngle(I, theta);
    options = optimset('MaxIter', 5000, 'TolX', thetaTol);
    [theta,Y,flag,output] = fminbnd(fun, 0, 180, options);
    Y = -Y;
    % Slope is actually opposite tand(thetaMax) because tand provides slope
    % in cartesian coordinates as opposed to image coordinates where Y-axis
    % increases going down, rather than up.
    dYdt = -tand(theta);
end


%% Objective Function
% TODO: move this inside fminsearch function for better organization?
function Y = RadonVarFromSlope(I, dYdt)
    % Calculate variance of image at angle specified by slope (dYdt)
    theta = atand(dYdt);
    Y = RadonVarFromAngle(I, theta);
%     R = radon(I, theta);
%     Y = var(R);
end

function Y = RadonVarFromAngle(I, theta)
    % Calculate variance of image at specified angle (theta)
    R = radon(I, theta);
    Y = var(R);
end

%%
function [dYdt, Y] = optimizeRadonAngleLegacy(final, I_sign)
    % Okay, so I_sign allows you to force all the input data to have a positive
    % slope by flipping the image left to right if it has a negative slope,
    % that way, you know the max will be in the range (90, 180] and therefore,
    % only have to do the transform over a smaller range (one quadrant rather
    % than two)

    if I_sign~=2
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
    
    % For radon transform, theta appears to be angle from vertical
    [R, xp] = radon(final, theta);

    % TODO: this assumes that there is some variance--but if the window is
    % small and there aren't many (or any) clear stripes then that doesn't
    % work. But then again maybe nothing will work in that case.
    % Determine maximum variance along columns
    Variance=var(R);
    [Y,iTheta]=max(Variance);


    %Debug mode
    Debug = 0;
    if Debug
        displayDebugger(R, theta, xp, Variance)
    end


    % TODO: why is this necessary? Can't you just do theta over range of [0,
    % 180) and then take tand(maxTheta). On interval [0,90) -> velocity is
    % positive, on interval (90,180) -> velocity is negative
    % TODO: need angle convention
    % TODO; theta-90?
    
    % For cartesian convention. Should be in range (-90, 90)
    thetaMax = theta(iTheta)-90;
    % Slope is actually opposite tand(thetaMax) because tand provides slope
    % in cartesian coordinates as opposed to image coordinates where Y-axis
    % increases going down, rather than up.
    dYdt = -tand(thetaMax);

    % if I_sign==1
    %     thetaMax = theta(n)-90;
    %     velocity = -1*(1/tand(thetaMax))*Tfactor*Xfactor;
    % elseif I_sign==0
    %     thetaMax = (theta(n)-90)*-1;
    %     velocity = (1/tand(-1*thetaMax))*Tfactor*Xfactor;
    % else
    %     if theta(n)<90
    %         if n-1<1
    %             n=2;
    %         end
    %         thetaMax = theta(n-1)-90;
    %     else
    %         thetaMax = theta(n)-90;
    %     end
    %     velocity = (1/tand(-1*thetaMax))*Tfactor*Xfactor;
    % end

    %toc

    
    % maxVel = 60;
    % minVel = -60;
    % tol = 0.1;
    % thetas = [atand(0:tol:maxVel), atand(minVel:tol:0)+180];
    % % OR
    % maxSpeed = 60;
    % thetaHalf = atand(0:tol:maxSpeed);
    % thetas = [thetaHalf, -fliplr(thetaHalf)+180];
    
    function displayDebugger(final, R, theta, xp, Variance)
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
    
    % TODO: there should be a faster way to do this at the specified
    % precision i.e. linearly spaced along velocity
    % TODO: make this persistent so only has to be run once?
    function theta = radonAngles(numInc,x)
        numInc=x*numInc;
        thetaFirst = deg2rad(1);
        thetaLast=deg2rad(89);
        velocityHigh=1/tan(thetaFirst);
        velocityLow=1/tan(thetaLast);

        velInc = (velocityHigh - velocityLow)/numInc;

        theta = zeros(1, x*173+1);
        theta(1) = thetaFirst;
        for i=1:x*173
            thetaNew=atan(1 / ((1/tan(theta(i))) - velInc));
            theta(i+1)=thetaNew;
        end

        theta = rad2deg(theta);
        
        % TODO: this can be vectorized
        thetaLin = zeros(1,x*(89-22)+1);
        n=1;
        for i=22:1/x:89
            thetaLin(n)=i;
            n=n+1;
        end

        theta = [theta, thetaLin];
    end

end


%% (Expermental) Optimization Functions
% % ```radon``` is O(theta) i.e. linear relationship between runtime and
% % length of theta vector. So, try to minimize the number of thetas that
% % radon is evaluated for by doing some sort of optimization. This function
% % essentially uses a binary search to optimize theta.
% % TODO: also write this as a while loop to see if more readable etc.
% % TODO: tolerance should be in velocity units, not theta units....
% % TODO: this can probably be optimized faster with the optimization
% % toolbox, this is just a binary search
% % TODO: this function is currently non-functional
% function [theta, Y] = optimizeRadonAngleRecurse(I, thetaRange, velTol)
%     % |----'----!----'----|
%     % Check quarters (') of range to narrow down range to one half of
%     % original, either [|, !], or [!, |]. If first quarter has a greater
%     % variance than the third quarter, new range is [|, !], otherwise it is
%     % [!, |].
%     midTheta = mean(thetaRange);
% %     iTol = (midTheta-thetaRange(1))/2;
%     thetas = [mean([thetaRange(1), midTheta]), mean([midTheta, thetaRange(2)])];
%     [R, xp] = radon(I, thetas);
% 
%     % TODO: this assumes that there is some variance--but if the window is
%     % small and there aren't many (or any) clear stripes then that doesn't
%     % work. But then again maybe nothing will work in that case.
%     %Determine maximum variance along columns
%     Variance=var(R);
%     [Y,n]=max(Variance);
% 
%     % Determine newThetaRange, which is half of the original thetaRange
%     if n == 1
%         newThetaRange = [thetaRange(n), midTheta];
%     else
%         newThetaRange = [midTheta, thetaRange(n)];
%     end
% 
%     % Calculate difference in velocity for 
%     % iTol = abs(diff(tand(newThetaRange)));
%     iTol = min( tand(newThetaRange)-tand(mean(newThetaRange)) );
%     % TODO: does this need to be if/else or just if?
%     if abs(iTol) < velTol
%         theta = thetas(n);
%     else
%         [theta, Y] = optimizeRadonAngle(I, newThetaRange, velTol);
%     end
% end
% 
% function [theta, Y, flag] = optimizeRadonAngleLoop(I, thetaRange, velTol)
%     % |----'----!----'----|
%     % Check quarters (') of range to narrow down range to one half of
%     % original, either [|, !], or [!, |]. If first quarter has a greater
%     % variance than the third quarter, new range is [|, !], otherwise it is
%     % [!, |].
%     i = 1;
%     maxIter = 5000;
%     iTol = Inf;
%     % Loop until tolerance is less than requested, or maximum number of
%     % iterations is exceeded.
%     while iTol > velTol && i <= maxIter
%         midTheta = mean(thetaRange);
%         thetas = [mean([thetaRange(1), midTheta]), mean([midTheta, thetaRange(2)])];
%         [R, xp] = radon(I, thetas);
% 
%         % TODO: this assumes that there is some variance--but if the window is
%         % small and there aren't many (or any) clear stripes then that doesn't
%         % work. But then again maybe nothing will work in that case.
%         %Determine maximum variance along columns
%         Variance=var(R);
%         [Y,n]=max(Variance);
% 
%         % Determine newThetaRange, which is half of the original thetaRange
%         theta = thetas(n);
%         if n == 1
%             thetaRange = [thetaRange(n), midTheta];
%         else
%             thetaRange = [midTheta, thetaRange(n)];
%         end
% 
%         % Calculate difference in velocity for 
% %         iTol = min(tand(thetaRange)-tand(theta));
%         iTol = abs(diff(tand(thetaRange)));
%         % TODO: does this need to be if/else or just if?
%         % TODO: this may be less readable than a while loop
%         i = i+1;
%     end
%     
%     flag = i <= maxIter;
%     % TODO: if i > maxIter, set the flag;
% end