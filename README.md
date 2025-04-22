# Relevant code for a study of precursor abundance and transformation rate as a function of distance from the edge and time since deposition

Relevant publication here:

Scripts

In MATLAB
Precursor_Phase_Distance.m
  ##developed by Isabelle Marie LeCloux and Zoë Rechav
  ##function: collect proportion of precursor in each pixel of a photoshop image of a coral skeleton area
  ##returns: .csv file containing the proportion of a specified precursor for each pixel in the image as a function of distance from the edge
  ##input
    all_masks: mask idenifying specified precursor, type: 8 bit greyscale , format: .png
    contiguous_mask: mask idenifying coral skeleton edge of area, type: 8 bit greyscale , format: .png
    pMap_image: proportion map of specified precursor, type: 16 bit grayscale, format: .tiff
    fov: field of view of image acquisitioned in PEEM, 45 micron or 56 micron
    length: depth to collect proportions into the coral skeleton, default 4 micron
    pix_wid: width of pMap_image, 1030 pixels
    csv_name: filepath of .csv file output

pp_tester_code.m
  ##developed by Isabelle Marie LeCloux and Zoë Rechav
  ##script to call function Precursor_Phase_Distance

in PYTHON
fitting_functions.py
  ##developed by Zoë Rechav
  ##library of functions used to perform fits
  ##basic function structures and fitting structures are provided
  ##can modify for your own functions as needed
