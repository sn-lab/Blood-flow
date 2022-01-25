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
p.parse(varargin{:});
I_sign = 0;

I = p.Results.I;
Xfactor = p.Results.Xfactor;
Tfactor = p.Results.Tfactor;

if I_sign==0
    I=fliplr(I);
end

I_size = size(I);

%% extract velocity by radon transform method
%figure(1), imshow(I);
%I=makeNoise(I);
%figure(2), imshow(I);

filtVar=25;     %Filter variance set to 25
a=filtVar;

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
