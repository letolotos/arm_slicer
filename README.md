# arm_slicer
A tool for slicing spiral structure of disk galaxies.

This program is to be used for slicing spiral arms of disk galaxies. With its help you can measure the widths of spiral arms and estimate the asymmetry of their profiles. The code is written in Python3.

For arm_slicer to work succesfully on your computer you should have DS9 software installed since you'll be working with it interactively. Also make sure that you have the following Python modules installed on your machine: Pyfits, SciPy, NumPy, and Matplotlib. 

So, let us have three images (corresponding to three spectral bands) for the galaxy of interest with names file_name_g.fits, file_name_i.fits, and file_name_r.fits.

To run the code just type:

python arm_slicer.py file_name bands

in Terminal. 

For the set of bands g, r, and i bands keyword is simply "gir":

python arm_slicer.py file_name gir

The folder you're running the program at should contain fits files for the set of bands for the galaxy of interest, arm_slicer.py, and the file galaxy_info.dat with input parameters for the code.
