function Result = calcLinescanVelRadon(varargin)
%%

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
p.addOptional('Tfactor', 1, @(x) isnumeric(x)&&isscalar(x)); % microns/pixel
p.addOptional('Xfactor', 205/500*250/512, @(x) isnumeric(x)&&isscalar(x)); % microns/pixel
% TODO: maybe this should be passed as a function handle to allow user to
% pass other arguments such as I_sign and customize optimizer
p.addOptional('Optimizer', @ischar); % microns/pixel
p.parse(varargin{:});
I_sign = 0;

I = p.Results.I;
Xfactor = p.Results.Xfactor;
Tfactor = p.Results.Tfactor;

% if I_sign==0
%     I=fliplr(I);
% end

I_size = size(I);

%% extract velocity by radon transform method
%figure(1), imshow(I);
%I=makeNoise(I);
%figure(2), imshow(I);

a=25;     %Filter variance set to 25

%Create high pass flter using isotropic Gaussian
% TODO: this should be moved out to filter full image all at once so don't
% have to recalculate kernel on each call of this function
gaus = zeros(I_size);
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

switch p.Results.Optimizer
    case 'fminbnd'
        [velocity, Y] = optimizeRadonAngleFminsearch(final, 0.01);
    case 'optimizer'
        [thetaMax, Y] = optimizeRadonAngleFminbnd(final, [0,180], 0.01);
    case 'binarysearch'
        [theta, Y] = optimizeRadonAngleRecurse(final, [0,180], 0.01);
    case 'legacy'
        
end


%Debug mode
Debug = 0;
if Debug
    showRadon();
end
 

% TODO: why is this necessary? Can't you just do theta over range of [0,
% 180) and then take tand(maxTheta). On interval [0,90) -> velocity is
% positive, on interval (90,180) -> velocity is negative
% thetaMax = theta(n);
% velocity = tand(thetaMax)*Tfactor*Xfactor;
velocity = velocity/Tfactor/Xfactor;
thetaMax = NaN;

Result = [0, 0, velocity, thetaMax, Y];
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


% ```radon``` is O(theta) i.e. linear relationship between runtime and
% length of theta vector. So, try to minimize the number of thetas that
% radon is evaluated for by doing some sort of optimization. This function
% essentially uses a binary search to optimize theta.
% TODO: also write this as a while loop to see if more readable etc.
% TODO: tolerance should be in velocity units, not theta units....
% TODO: this can probably be optimized faster with the optimization
% toolbox, this is just a binary search
% TODO: this function is currently non-functional
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
        [theta, Y] = optimizeRadonAngle(I, newThetaRange, velTol);
    end
end

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


% TODO: why isn't this using fminbnd?
function [vel, Y, flag] = optimizeRadonAngleFminsearch(I, velTol)
    fun = @(vel) -VelAngleRadonVar(I, vel);
    options = optimset('MaxIter', 5000, 'TolX', velTol);
%     [theta,Y,flag,output] = fminbnd(fun, thetaRange(1), thetaRange(2), options);
    [vel,Y,flag,output] = fminsearch(fun, 0, options);
    Y = -Y;
end


function optimizeRadonAngleLegacy()
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
    [R, xp] = radon(final, theta);

    % TODO: this assumes that there is some variance--but if the window is
    % small and there aren't many (or any) clear stripes then that doesn't
    % work. But then again maybe nothing will work in that case.
    %Determine maximum variance along columns
    Variance=var(R);
    [Y,n]=max(Variance);


    %Debug mode
    Debug = 0;
    if Debug
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


    % TODO: why is this necessary? Can't you just do theta over range of [0,
    % 180) and then take tand(maxTheta). On interval [0,90) -> velocity is
    % positive, on interval (90,180) -> velocity is negative
    thetaMax = theta(n);
    velocity = tand(thetaMax)*Tfactor*Xfactor;

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

        thetaLin = zeros(1,x*(89-22)+1);
        n=1;
        for i=22:1/x:89
            thetaLin(n)=i;
            n=n+1;
        end

        theta = [theta, thetaLin];
    end

end



function Y = VelAngleRadonVar(I, vel)
    % Calculate angle from requested velocity, with wrapping
    theta = atand(vel);
    if vel < 0
        theta = theta+180;
    end
    
    R = radon(I, theta);
    Y = var(R);
end