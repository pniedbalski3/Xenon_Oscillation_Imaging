# Xenon_Oscillation_Imaging
Work in Progress code for analyzing regional oscillations in the xenon RBC signal.

## Authorship
Author: Peter J. Niedbalski 

Contact: pniedbalski@kumc.edu


## Main Function: wiggle_imaging.m
This function is the overall controller for performing oscillation analysis. This is a work in progress, and I am actively looking at multiple possible ways to generate these images. 

First global oscillations in RBC k0 are calculated using the function Wiggle_Tools.bin_wiggles. This code performs a bandpass filter then finds the mean difference between the peaks and the troughs of the data.

As of 10/5/2021, 3 methods of calculating regional oscillations are performed:

- Bin "high" and "low" regions of k0, then create two key images from RBC data. The scaled difference between these key images is then the RBC oscillation amplitude (This is the method published in JAPPL)

The two additional methods that I am exploring start the same way:
Reconstruct key images in increments of 200 ms (i.e., create a key image for every 0.2 s of data). This results in >30 key images. The mean of these images shows an oscillation. Then, voxel-by-voxel:
- Perform a sine fitting to get amplitude and phase of oscillation (I'm not convinced that this is ideal, because the data is very noisy... can we really get an accurate phase measurement?)
- Find the mean difference between peaks and troughs

The jury is still out on which way is best. Way number 1 is published, but I think there are several fatal flaws - most notably how do you bin people with noisy data/small oscillations?
Ways 2 and 3 basically require using a fast Dixon acquisition approach - Using a TR of 8 ms (dissolved to dissolved), .2 s gives 25 projections per key - that is, we're really undersampled. With TR = 15 (standard), that number is down to 13. It's unclear whether this massive undersampling provides good enough images to get sensible values from.
However, ways 2 and 3 do mitigate concerns with binning k0 - It doesn't matter what the raw k0 oscillations look like, because we are leaving our analysis of oscillations until we get to image space. My hope is that this will help with repeatability.

I'm still working to optimize these and run some comparisons. 

### Inputs
- disfid: all dissolved FIDs, shaped as (NPts x NFIDs)
- gasfid: all gas FIDs, shaped as (NPts x NFIDs)
- traj: trajectories for reconstructing dissolved images
- H1_im: anatomic image (purely used to make pretty pictures. could pass a matrix of zeros instead)
- Gas_Image: low-resolution gas image from dissolved imaging
- H1_Mask: Mask from anatomic image
- Vent_Mask: H1_Mask with ventilation defect voxels removed
- RBC_Mask: H1_Mask with ventilation defect and RBC defect voxels removed
- RBC2Bar: RBC to barrier ratio
- TR: Repetition time (in ms)
- ImSize: desired image size
- scanDateStr: string with the scan date for reporting purposes
- write_path: path where you want images and reports to be saved 

### Outputs: 
No outputs. This function writes out several figures and .mat files to the path specified by write_path. I've recently made a lot of edits to this to generate fewer images write out fewer things until I'm able to figure out optimal methods. 
