SE-Sync is a certifiably correct algorithm for synchronization over the special Euclidean group, a common problem arising in the context of 2D and 3D geometric estimation (e.g. pose-graph SLAM and camera motion estimation).

We are making this software freely available in the hope that it will be useful to others.  If you use SE-Sync in your own work, please cite our papers:

@inproceedings{Rosen2016Certifiably,
title = {A Certifiably Correct Algorithm for Synchronization over the Special {Euclidean} Group},
author = {Rosen, D.M. and Carlone, L. and Bandeira, A.S. and Leonard, J.J.},
booktitle = {Intl. Workshop on the Algorithmic Foundations of Robotics (WAFR)},
month = dec,
year = {2016},
address = {San Francisco, CA},
}

@techreport{Rosen2016SESync,
title = {{SE-Sync}: A Certifiably Correct Algorithm for Synchronization over the Special {Euclidean} Group},
author = {Rosen, D.M. and Carlone, L. and Bandeira, A.S. and Leonard, J.J.},
institution = {Computer Science and Artificial Intelligence Laboratory, Massachusetts Institute of Technology},
address = {Cambridge, MA},
number = {MIT-CSAIL-TR-2017-002},
year = {2017},
month = feb,
}

and the following paper of Absil et al., which describes the Riemannian trust-region (RTR) method that SE-Sync employs:

@article{Absil2007Trust,
title = {Trust-Region Methods on {Riemannian} Manifolds},
author = {Absil, P.-A. and Baker, C.G. and Gallivan, K.A.},
journal = {Found.\ Comput.\ Math.},
volume = {7},
number = {3},
pages = {303--330},
year = {2007},
month = jul,
}

If you use the MATLAB implementation of SE-Sync, please also cite the following reference for the Manopt toolbox, which provides the MATLAB implementation of RTR that the SE-Sync toolbox employs:

@article{Boumal2014Manopt,
  title={{Manopt}, a {MATLAB} Toolbox for Optimization on Manifolds.},
  author={Boumal, N. and Mishra, B. and Absil, P.-A. and Sepulchre, R.},
  journal={Journal of Machine Learning Research},
  volume={15},
  number={1},
  pages={1455--1459},
  year={2014}
}

If you use the C++ implementation of SE-Sync, please also cite the following reference for the ROPTLIB library, which provides the C++ implementation of RTR that the SE-Sync C++ library employs:

@techreport{Huang16ROPTLIB,
title = {{ROPTLIB}: An Object-Oriented {C++} Library for Optimization on {Riemannian} Manifolds},
author = {Huang, W. and Absil, P.-A. and Gallivan, K.A. and Hand, P.},
institution = {Florida State University},
number = {FSU16-14},
year = {2016},
}


==== Copyright and License ====

The C++ and MATLAB implementations of SE-Sync contained herein are copyright (C) 2016 - 2017 by David M. Rosen, and are distributed under the terms of the GNU General Public License (GPL) version 3 (or later).  Please see the files LICENSE.txt and COPYING.txt for more information.

Contact: dmrosen@mit.edu


==== Getting Started ====

MATLAB:

To use the MATLAB implementation of SE-Sync, simply place the 'MATLAB' folder in any convenient (permanent) location, and then run the script MATLAB/import_SE_Sync.m.  Congrats!  SE-Sync is now ready to go :-).  For a minimal working example, see MATLAB/examples/main.m

C++:

The C++ implementation of SE-Sync can be built and exported as a CMake project.  For a minimal working example, see C++/examples/main, which provides a simple command-line utility for processing .g2o files.
