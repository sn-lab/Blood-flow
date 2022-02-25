% NOTE: this is not a perfect replacement for MATLAB's "imrotate" -- very
% slightly different results
% TODO: maybe allow user to turn on/off smooth edges
% TODO: varargin??
% Theta in degrees
function IRot = imrotate(I, theta, method, bbox, fillval, smooth)
    p = inputParser;
    p.addRequired('I')
    p.addRequired(theta, @isreal);
    p.addOptional('method', 'nearest')  % Let interp2 validate this input
    p.addOptional('bbox', 'loose', @(x) any(strcmp(x, {'loose', 'crop'})));
    p.addOptional('fillval', 0);
    p.addOptional('smooth', false);
    p.parse(varargin{:});
    
    
    if rem(theta, 90) == 0
        theta = rem(theta, 360);    % Remove multiples of 360
        IRot = rot90(I, theta/90);
    else

        T = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];

        sz = size(I);

        switch bbox
            case 'crop'
                warning('Not yet supported')
                xLimitsOutInt = [1, sz(2)];
                yLimitsOutInt = [1, sz(1)];
            case 'loose'
                % Calculate corner points
                xLimitsIn = [1 sz(2)];
                yLimitsIn = [1 sz(1)];
                xMid = mean(xLimitsIn);
                yMid = mean(yLimitsIn);
                xCentered = xLimitsIn - xMid;
                yCentered = yLimitsIn - yMid;

                [xCorners,yCorners] = meshgrid(xCentered, yCentered);

                % Transform corner points to find output boundaries
                xCornersOut = T(1,1).*xCorners + T(2,1).*yCorners;
                yCornersOut = T(1,2).*xCorners + T(2,2).*yCorners;

                % XLimitsOut/YLimitsOut are min and max of transformed points
                xLimitsOut = [min(xCornersOut(:)), max(xCornersOut(:))];
                yLimitsOut = [min(yCornersOut(:)), max(yCornersOut(:))];

                % Round away from zero
                xLimitsOutInt = ceil(abs(xLimitsOut)).*sign(xLimitsOut);
                yLimitsOutInt = ceil(abs(yLimitsOut)).*sign(yLimitsOut);
            otherwise
                error(['Invalid ''bbox'' parameter ''', bbox '''. Must be either ''loose'' or ''crop'''])
        end

        % Pad image for smooth edges
        if smooth
            I = [zeros(sz(1),1), I, zeros(sz(1),1)];
            I = [zeros(1,sz(2)+2); I; zeros(1,sz(2)+2)];
            xMid = xMid + 1;
            yMid = yMid + 1;
        end


        [XqOrig, YqOrig] = meshgrid(xLimitsOutInt(1):xLimitsOutInt(2), yLimitsOutInt(1):yLimitsOutInt(2));
        % Apply rotation matrix
        % Cleaner to do it vector-by-vector rather than constructing an xy
        % coordinate matrix for matrix multiplication
        Xq = T(1,1).*XqOrig + T(1,2).*YqOrig + xMid;
        Yq = T(2,1).*XqOrig + T(2,2).*YqOrig + yMid;

        IRot = interp2(I, Xq, Yq, method, fillval);
    end
end