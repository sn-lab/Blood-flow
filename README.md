# Blood-flow
This repository contains the most up-to-date version of the linescan code used for measuring blood flow from linescans

### In here you will find:

- 1. **scim_openTif.m** (only use if you forgot the ms/line you used to image. You can also use the pipeline to get this number)
	Input: Raw image stacks
	Output: number of ms/line

- 2. **pixelSizeFromImageHeader.m** (you should use the pipeline to get this number. THIS IS NOT UP-TO-DATE WITH THE CONVERSION FACTORS FOR ALL THE SETUPS. If you decide to use it, make sure you are using the appropriate conversion factors)
	Input: Raw image stacks
	Output: Conversion factor for um/pixels of your image.

- 3. **extractVelTiffShared.m** extracts the raw blood flow velocity from each file using a radon transform algorithm. requires um/pixel and ms/line, also requires Texas Red or blood plasma label channel only
	Input: Blood labelling channel
	Output: .mat file with the RAW velocity data

- 4. **view_velocities_save_data.m** extracts the final velocity values by filtering noise and background from the initial measurements. This is done with input from the user
	Input: rawVel.mat files from extractVelTiffShared.m 
	Output: rawVel ExcludedPts.mat and excel file with average values of velocity and 	standard deviation 

- 5. **avgDiameterinROI.m** extracts the diameter of the vessel by selecting a rectangular area on the vessel manually. gives value in pixels and user needs to use um/pixel conversion factor to transform to microns

**UN THE FILES IN THE ORDER THEY APPEAR IN THIS README. YOU MAY SKIP #1 AND #2 IF YOU ALREADY HAVE MS/LINE AND PIXEL/UM conversion factors**
