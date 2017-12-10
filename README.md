# SE-Sync


**SE-Sync** is a *certifiably correct* algorithm for synchronization over the special Euclidean group, a common problem arising in the context of 2D and 3D geometric estimation (for example, [pose-graph SLAM](http://domino.informatik.uni-freiburg.de/teaching/ws11/robotics2/pdfs/ls-slam-tutorial.pdf),  camera motion estimation, and sensor network localization).  

A detailed description of the algorithm and its implementation can be found in our [technical report](https://github.com/david-m-rosen/SE-Sync/blob/master/SE-Sync%20-%20A%20Certifiably%20Correct%20Algorithm%20for%20Synchronization%20over%20the%20Special%20Euclidean%20Group.pdf).



## Getting Started

### MATLAB

To use the MATLAB implementation of SE-Sync, simply place the 'MATLAB' folder in any convenient (permanent) location, and then run the script MATLAB/import_SE_Sync.m.  Congrats!  SE-Sync is now ready to go :-).  For a minimal working example, see MATLAB/examples/main.m

### C++

The C++ implementation of SE-Sync can be built and exported as a CMake project.  For a minimal working example, see C++/examples/main, which provides a simple command-line utility for processing .g2o files.

## References

We are making this software freely available in the hope that it will be useful to others. If you use SE-Sync in your own work, please cite [our](https://github.com/david-m-rosen/SE-Sync/blob/master/A%20Certifiably%20Correct%20Algorithm%20for%20Synchronization%20over%20the%20Special%20Euclidean%20Group.pdf) [papers](https://github.com/david-m-rosen/SE-Sync/blob/master/SE-Sync%20-%20A%20Certifiably%20Correct%20Algorithm%20for%20Synchronization%20over%20the%20Special%20Euclidean%20Group.pdf):

```
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
```

and the following [paper](https://pdfs.semanticscholar.org/90b8/a3b089509dfea2fb83b2e49d77a443b2a3f7.pdf) of Absil et al., which describes the Riemannian trust-region (RTR) method that SE-Sync employs:

```
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
```

If you use the MATLAB implementation of SE-Sync, please also cite the following [reference](http://www.jmlr.org/papers/volume15/boumal14a/boumal14a.pdf) for the [Manopt toolbox](https://www.manopt.org/), which provides the MATLAB implementation of RTR that the SE-Sync toolbox employs:

```
@article{Boumal2014Manopt,
  title={{Manopt}, a {MATLAB} Toolbox for Optimization on Manifolds.},
  author={Boumal, N. and Mishra, B. and Absil, P.-A. and Sepulchre, R.},
  journal={Journal of Machine Learning Research},
  volume={15},
  number={1},
  pages={1455--1459},
  year={2014}
}
```
If you use the C++ implementation of SE-Sync, please also cite the following [reference](https://www.math.fsu.edu/~whuang2/papers/ROPTLIB.htm) for the [ROPTLIB library](https://github.com/whuang08/ROPTLIB), which provides the C++ implementation of RTR that the SE-Sync C++ library employs:

```
@techreport{Huang16ROPTLIB,
title = {{ROPTLIB}: An Object-Oriented {C++} Library for Optimization on {Riemannian} Manifolds},
author = {Huang, W. and Absil, P.-A. and Gallivan, K.A. and Hand, P.},
institution = {Florida State University},
number = {FSU16-14},
year = {2016},
}
```


## Copyright and License 

The C++ and MATLAB implementations of SE-Sync contained herein are copyright (C) 2016 - 2017 by David M. Rosen, and are distributed under the terms of the GNU General Public License (GPL) version 3 (or later).  Please see the files [LICENSE.txt](https://github.com/david-m-rosen/SE-Sync/blob/master/LICENSE.txt) and [COPYING.txt](https://github.com/david-m-rosen/SE-Sync/blob/master/COPYING.txt) for more information.

Contact: dmrosen@mit.edu
