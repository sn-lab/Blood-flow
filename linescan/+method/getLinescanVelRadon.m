function Result = getLinescanVelRadon(varargin)
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

I = p.Results.I;
Xfactor = p.Results.Xfactor;
Tfactor = p.Results.Tfactor;

I_size = size(I);

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
        if n-1<1
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


%% radonAngles
function theta = radonAngles(numInc,x)

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
for i=22:1/x:89
    thetaLin(n)=i;
    n=n+1;
end

theta = [theta, thetaLin];
end
