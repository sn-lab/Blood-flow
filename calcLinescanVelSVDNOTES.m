% function Result = f_find_vel(small, (Tfactor), (Xfactor), (slope), (useaverage), (debug))
% based on RotMatOpt19_rand_time
% IN: small = 1 frame
%               Xfactor (microns/pixel)
%               Tfactor (pixels/ms)
% OUT: Result: (preserved data structure)
%       unneeded numbers are = 0
%       column     3) Velocity (mm/s) + veloctiy is in x-dir (RBC's from
%                             left to right)
%                       4) Sep
%                       5) Angle (true angle of data unmodified by this
%                       function, positive is RBC move left to right)
%                       6) Flux
% 07-11-03: Make an option to not use average across the frame. Useful for
% slow capillaries.
%12-20-04: Finds Flux based on method developed empirically: Thresholds
% image of rotated block data by average; Takes an average projection, and
% thresholds; Finds derivative and zerocrossings to find RBC edges. Uses a
% ratio of standard deviation of intensities across time and space to
% reject some data points.