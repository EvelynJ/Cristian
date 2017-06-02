Codes for Cristian Herrera for the kinematics project with NGC 1407 and IC 719

This code is written in IDL, and is designed to logarithmically bin the datacube in the wavelength direction (necessary step for the kinematics measurements), measure the S/N per spaxel, Voronoi bin the datacube and finally measure the kinematics with pPXF.

Getting started:
  - download the ppxf and voronoi binning codes in IDL
  http://www-astro.physics.ox.ac.uk/~mxc/software/
  - download the IDL astronomy library if it's not already installed
  https://idlastro.gsfc.nasa.gov/homepage.html
  - download the coyote codes. These are useful codes, especially for plotting etc
  http://www.idlcoyote.com/programs/zip_files/coyoteprograms.zip
  - download the contents of this repository
  - set the IDL paths to the codes in this repository and those you've downloded.
  e.g. export IDL_PATH=+/home/ejohnsto/GitHub/Cristian/:$IDL_PATH
  - in the "kinematics_input.txt" file included in this repositry, set up the necessary information. this will tell the code the directory paths to the data, the pixel values for the centre of the galaxy, the S/N limit and target etc. feel free to ask if you want me to explain any of the points in that file.
  - decide which steps of the code you want to run. in the text file, look at those options starting with D). You have four steps to this code- log-rebin the data, voronoi bin the data, measure the kinematics, and plot the kinematics. To run each step, choose 'y' (for yes), and to skip it, use 'n' (for no).
  - in the terminal window, open IDL by typing >  idl
  - once IDL has started, run the code (from the directory you kinematics_input.txt file is in) >  measure_kinematics,'kinematics_input.txt'
