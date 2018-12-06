  This folder includes functions and data involved in SMDisplay.m, which are extracted from the GIFT toolbox (http://mialab.mrn.org/software/gift/index.html). 

  Functions:

  (1) icatb_getColormap.m: gets colormap based on number of components, whether the data is absoulte value and if there is a structural image. A matlab file specified by icatb_defaults contains the colorsmaps it uses.
  (2) icatb_overlayImages.m: puts components on top of structural image. 

  Datasets: 

  (1) 2DSMshow.nii: component file for volume information

  (2) niisetV: volume information, including
        V: spm - 3D normalized

  (3) parameters_returnResizedImage_ch2betr.mat: parameters setting for icatb_returnResizedImage, includuing
        anatomicalPlane: anatomical plane
        slices_in_mm: Slices in mm (range)
        dataType: data type ('real' or 'complex')
        file_numbers: number of files
        mask_file: full file path of the mask file

  See the GIFT toolbox for more information.
  

