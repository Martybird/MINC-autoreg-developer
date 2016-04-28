# MINC-autoreg-developer
A light version of MINC-autoreg that provides a special version of minctracc that performs 2D multi-slice (e.g., tracked ultrasound slices) to 3D volume image registration.



*If you are using this software, please cite the following articles:

Yan CX, Goulet B, Chen SJ, Tampieri D, Collins DL. Validation of automated ultrasound-CT registration of vertebrae. Int J Comput Assist Radiol Surg. 2012 Jul;7(4):601–10.

Yan, C. X., Goulet, B., Tampieri, D., & Collins, D. L. (2012). Ultrasound-CT registration of vertebrae without reconstruction. International journal of computer assisted radiology and surgery, 7(6), 901–909.

C.X.B. Yan, B. Goulet, J. Pelletier, S.J.S. Chen, D. Tampieri and D.L. Collins, “Towards Accurate, Robust and Practical Ultrasound-CT Registration of Vertebrae for Image-Guided Spine Surgery,” International Journal of Computer Assisted Radiology and Surgery, 2011 Jul;6(4):523–37.

* How to use it?


./minctracc -clobber -debug -est_center -lsq6 -xcorr -slices_to_volume -nearest_neighbour -source_lattice -transformation initial_transform.xfm -simplex 4 -tol 0.001 -threshold 0 -1 -step 1 1 1 us_slices_00???.mnc ct_volume.mnc final_transform.xfm
