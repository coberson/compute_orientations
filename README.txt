These mathematica files make it possible to compute and show the local orientation vector field of a 3D image, using the structure tensor.  
The issuing vector field can be seen as a description of the local orientation of a 3D image and can be used to describe the orientation of objects in such an image, for example in biomedical images.
The code was written in 2012, during my time at the Workgroup Imaging, University of Muenster, as part of a joint project with Christoph Brune and Hendrik Dirks.


The co_3D_abstract.pdf file is a short paper justifying the algorithm used in the compute_orientations_3D functions.
The main function is in compute_orientations_3D.m 
Are needed fspecial3.m and showVectorfield3D.m. 

ellipsoid_gf generates an ellipsoid for testing purposes.

