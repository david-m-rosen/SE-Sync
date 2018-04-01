# SE-Sync


**SE-Sync** is a *certifiably correct* algorithm for performing *synchronization over the special Euclidean group*: estimate the values of a set of unknown *poses* (positions and orientations in Euclidean space) given noisy measurements of a subset of their pairwise relative transforms.  This problem frequently arises in the context of 2D and 3D geometric estimation; for example, the foundational problems of [pose-graph SLAM](http://domino.informatik.uni-freiburg.de/teaching/ws11/robotics2/pdfs/ls-slam-tutorial.pdf) (in robotics), [camera motion estimation](http://cmp.felk.cvut.cz/ftp/articles/pajdla/Martinec-Pajdla-CVPR-2007.pdf) (in computer vision), and [sensor network localization](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3741000/) (in distributed sensing) all require synchronization over the special Euclidean group.  SE-Sync improves upon prior methods by exploiting a novel (convex) *semidefinite relaxation* of the special Euclidean synchronization problem to directly search for *globally* optimal solutions, and is capable of producing a *computational certificate of correctness* (global optimality) in the (typical) case that a global minimizer is found.

A detailed description of the algorithm and its implementation can be found in our [technical report](https://github.com/david-m-rosen/SE-Sync/blob/master/SE-Sync%20-%20A%20Certifiably%20Correct%20Algorithm%20for%20Synchronization%20over%20the%20Special%20Euclidean%20Group.pdf).



## Getting Started

### MATLAB

To use the MATLAB implementation of SE-Sync, simply place the 'MATLAB' folder in any convenient (permanent) location, and then run the script MATLAB/import_SE_Sync.m.  Congrats!  SE-Sync is now ready to go :-).  For a minimal working example, see [MATLAB/examples/main.m](https://github.com/david-m-rosen/SE-Sync/blob/master/MATLAB/examples/main.m)

### C++

The C++ implementation of SE-Sync can be built and exported as a CMake project.  For a minimal working example, see [C++/examples/main](https://github.com/david-m-rosen/SE-Sync/blob/master/C%2B%2B/examples/main.cpp), which provides a simple command-line utility for processing .g2o files.

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


## Copyright and License 

The C++ and MATLAB implementations of SE-Sync contained herein are copyright (C) 2016 - 2018 by David M. Rosen, and are distributed under the terms of the GNU Lesser General Public License (LGPL) version 3 (or later).  Please see the [LICENSE](https://github.com/david-m-rosen/SE-Sync/blob/master/LICENSE) for more information.

Contact: drosen2000@gmail.com
