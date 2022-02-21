# Blood-flow
This repository contains the most up-to-date version of the linescan code used for measuring blood flow from linescans

### In here you will find:

- 1. **scim_openTif.m** (only use if you forgot the ms/line you used to image. You can also use the pipeline to get this number)
	Input: Raw image stacks
	Output: number of ms/line

- 2. **pixelSizeFromImageHeader.m** (you should use the pipeline to get this number. THIS IS NOT UP-TO-DATE WITH THE CONVERSION FACTORS FOR ALL THE SETUPS. If you decide to use it, make sure you are using the appropriate conversion factors)
	Input: Raw image stacks
	Output: Conversion factor for um/pixels of your image.

- 3. **extractVelTiffShared.m** extracts the raw blood flow velocity from each file using a **radon transform** algorithm. requires um/pixel and ms/line, also requires Texas Red or blood plasma label channel only
	Input: Blood labelling channel
	Output: .mat file with the RAW velocity data

- 4. **velocity_from_tif.m** extracts the raw blood flow velocity from each file using a **SVD algorithm**. requires um/pixel and ms/line, also requires Texas Red or blood plasma label channel only
	Input: Blood labelling channel
	Output: .mat file with the RAW velocity data
	
- 5. **view_velocities_save_data.m** extracts the final velocity values by filtering noise and background from the initial measurements. This is done with input from the user
	Input: rawVel.mat files from extractVelTiffShared.m 
	Output: rawVel ExcludedPts.mat and excel file with average values of velocity and 	standard deviation 

- 6. **avgDiameterinROI.m** extracts the diameter of the vessel by selecting a rectangular area on the vessel manually. gives value in pixels and user needs to use um/pixel conversion factor to transform to microns

**UN THE FILES IN THE ORDER THEY APPEAR IN THIS README. YOU MAY SKIP #1 AND #2 IF YOU ALREADY HAVE MS/LINE AND PIXEL/UM conversion factors**

## From share_instructions_2012-03-08.pdf (Revised on 3/8/2012):
### extractVelTiffShared.m

This Matlab function calculates red blood cell velocity from linescan files. This code was originally developed by Nozomi Nishimura, Nathan Cornelius, and Puifai Santisakultarm in David Kleinfeld’s and Chris Schaffer’s labs.

This software is freely available to researchers in academia or all non-commercial, non-profit making organizations. Please reference the following publications:

T. P. Santisakultarm, N. R. Cornelius, N. Nishimura, A. I. Schafer, R. T. Silver, P. C. Doerschuk,
W. L. Olbricht, and C. B. Schaffer. In Vivo Two-photon Excited Fluorescent Microscopy Reveals Cardiac- and Respiration-dependenct Pulsatile Blood Flow in Cortical Blood Vessels in Mice. Am J Physiol Heart Circ Physiol (2012).

C. B. Schaffer, B. Friedman, N. Nishimura, L. F. Schroeder, P. S. Tsai, F. F. Ebner, P. D. Lyden, and D. Kleinfeld, Two-photon imaging of cortical surface microvessels reveals a robust redistribution in blood flow after vascular occlusion. Public Library of Science Biology 4, e22 (2006).

This is an ongoing project, and updates and new algorithms are anticipated. Please send us your email address if you would like to be informed of any changes. Questions, bugs, discussion, etc. please contact:

Puifai Santisakultarm ts386@cornell.edu

Step by step instructions

Installation – copy the files ‘extractVelTiffShared.m’ into the current directory in Matlab, or add the location of the files to the Matlab path.

1.	New Files
•	New files to look at individually?
Select ‘Yes’ to define the region of interest (ROI) for files. Or if ROIs are already saved, click ‘No’ and you will be asked for the parameter file with filenames and ROIs. Skip to step 9.

2.	Subtract Average
•	Always subtract average?
‘Yes’ means that at each time step the average intensity across the spatial dimension is calculated and subtracted. This removes vertical stripes in the data that result from variations in intensity across a blood vessel. If flow is very slow and the stripes are vertical, do not subtract the average. If ‘No,’ the program will ask the user whether the average should be subtracted for each file.

3.	Acquisition Rate and Magnification
•	Same acquisition rate and magnification for all files? If ‘Yes,’ Processing Parameters window will pop up.
o	Acquisition rate (ms/line)
o	Magnification (um/pixel)
If ‘No,’ the program will ask for acquisition rate and magnification for each file.

4.	Select file to open
Click on the ‘.tif or .tiff’ file with linescan data. This program assumes greyscale images. If the data is stored as multiple frames, the program assumes that there is no gap in time between the last line of a frame and the first line of the next frame. The frame can be any pixel size.
 

5.	Select region of interest
Use the mouse to draw a box that defines the left and right borders of the region of interest (ROI). The larger the ROI, the slower the analysis will run. Typically, we choose ROIs less than 100 pixels wide for processing speed. If the data is noisy, more pixels may be useful. You can now scroll through the images with the following keys:
‘f’ – goes forward 1 frame ‘b’ – goes back 1 frame
‘s’ – skips forward 10 frames
Press ‘space’ to keep the ROI and go to the next step, or press ‘n’ to redefine the ROI.

6.	Slope of lines
•	Slope
Pick ‘both’ when unsure or if the slope of the linescan file changes

7.	Use Average?
•	Subtract average across linescans?

8.	Continue
•	More files?
Press ‘yes’ to define the ROI for more files, or ‘no’ to process the data.

9.	Comma delimited file save as
Enter a file name to save the filenames and the ROIs. This will save a ‘.csv’ file which you can open in Excel. It will also save a ‘Parameters.mat’ file that saves the same data in a format for Matlab.

10.	Continue
•	Calculate velocities now?

11.	Processing parameters
•	WinPixelsDown
Number of time pixels between velocity measurements. Often, this is smaller than the number of time pixels per data point so that there is overlap between data points.
•	WinSize
This is the number of pixel lines used to find a velocity for a given point in time. This will depend on the scan speed of your microscope and the your required time resolution. Typically, this is chosen to be short relative to the cardiac cycles
•	Maxlines
Total number of pixel lines to process. If the file has less than this number of lines, the analysis will end when the end of the file is reached. The program will then go on to the next file.
•	Start with file #:
Files will be processed in the order that the ROIs were selected. You do not have to start with the first file.

The command window will show the line number of the data being processed. When the processing is done, the message ‘done’ will appear and the program will beep.
For each ‘.tif’ file ,a new file will be generated with the file name of the ‘.tif’ file appended by ‘ rawVel (WinPixelsDown)(WinSize).mat’. Please note that if you run the same ‘.tif’ file more than one time, the file containing the results with be overwritten.

The data is saved in a matrix with 4 columns called ‘Result.’
i.	Starting line number
ii.	Time (ms)
iii.	Velocity (mm/s)
 
iv.	Angle of stripes
The algorithm does generate some incorrect velocities. Often, these are clustered near 0 mm/s or at very large values. We recommend inspecting each file visually and manually selecting the range of valid velocities. Usually, the cluster of invalid data is quite distinct from valid data.
