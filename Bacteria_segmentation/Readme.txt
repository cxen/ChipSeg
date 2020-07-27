'segmentation_main_offline.m' quantifies the fluorescence of single bacterial cells; the path with raw images in the code needs to be changed by the user.
The script works together with the following Matlab files, which should be in the same folder:
'segmentation_GF.m' to segment single cells using the Otsu algorithm,
'fluorescence_eval_Init.m' to calculate the fluorescence to define a fluorescence threshold, and
'fluorescence_eval_GF.m' to calculate the fluorescence at every sampling time.

The user is required to run 'segmentation_main_offline' and it will be asked to crop the region where cells are located and the background zone.

The folder with the raw images needs to contain the phase contrast and fluorescence images.
