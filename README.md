# Blood-flow
This repository contains the most up-to-date version of the linescan code used for measuring blood flow from linescans

##In here you will find:

-scim_openTif.m (only use if you forgot the ms/line you used to image)
	Input: Raw image stacks
	Output: number of ms/line

-pixelSizeFromImageHeader.m (you should use the pipeline to get this number. THIS IS NOT UP-TO-DATE WITH THE CONVERSION FACTORS FOR ALL THE SETUPS. If you decide to use it, make sure you are using the appropriate conversion factors)
	Input: Raw image stacks
	Output: Conversion factor for um/pixels of your image.

-extractVelTiffShared.m extracts the raw blood flow velocity from each file using a radon transform algorithm. requires um/pixel and ms/line, also requires Texas Red or blood plasma label channel only
	Input: Blood labelling channel
	Output: .mat file with the RAW velocity data

-view_velocities_save_data.m extracts the final velocity values by filtering noise and background from the initial measurements. This is done with input from the user
	Input: rawVel.mat files from extractVelTiffShared.m 
	Output: rawVel ExcludedPts.mat and excel file with average values of velocity and 	standard deviation 

- avgDiameterinROI.m extracts the diameter of the vessel by selecting a rectangular area on the vessel manually. gives value in pixels and user needs to use um/pixel conversion factor to transform to pixel