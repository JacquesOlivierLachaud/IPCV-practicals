# IPCV-practicals

Practicals for the IPCV lecture on digital geometry processing.

![DGtal logo](https://github.com/DGtal-team/DGtal/blob/main/doc/images/logoDGtal-small.png)

They are mainly based on a tutorial given at DGMM2022. Their original versions were written by

* [David Coeurjolly](https://perso.liris.cnrs.fr/david.coeurjolly)
* [Jacques-Olivier Lachaud](http://www.lama.univ-savoie.fr/pagesmembres/lachaud/People/LACHAUD-JO/person.html)
* [Tristan Roussillon](https://perso.liris.cnrs.fr/tristan.roussillon/)
* [Bertrand Kerautret](https://liris.cnrs.fr/page-membre/bertrand-kerautret)

## The slides

* Keynote: [part1](https://perso.liris.cnrs.fr/david.coeurjolly/talk/dgtal-tutorial/DGtalDGMM2022-part1.key), [part2](https://perso.liris.cnrs.fr/david.coeurjolly/talk/dgtal-tutorial/DGtalDGMM2022-part2.key)
* pdf export: [part1](https://perso.liris.cnrs.fr/david.coeurjolly/talk/dgtal-tutorial/DGtalDGMM2022-part1.pdf), [part2](https://perso.liris.cnrs.fr/david.coeurjolly/talk/dgtal-tutorial/DGtalDGMM2022-part2.pdf)
 
## Get the code

To get the code, the easiest way is to `clone` the github repository. E.g. with the command-line git client:

```
git clone https://github.com/JacquesOlivierLachaud/IPCV-practicals
```

## First build

To build DGtal related examples, you would need:
  - C++17 enabled compiler (most c++ compilers are C+17)
  - a [cmake](https://cmake.org) client (at least 3.25)
  

After cloning IPCV-practicals, set-up the compiling environment in a terminal with
```
cd IPCV-practicals
mkdir build
cd build
# configure the compilation in Release mode (for performance).
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_POLICY_VERSION_MINIMUM=3.5 ..
```

It will automatically download `boost`, `eigen`, `polyscope` and `dgtal`.
This should finish with the following lines:
```
-- Configuring done (1.9s)
-- Generating done (0.1s)
-- Build files have been written to: ${HOME}/.../IPCV-practicals/build
```

After that just type `make` to compile all the provided programs, some
are examples, others are the skeleton c++ files of the different
practicals.

## Data

A few 3D vol images are provided in `data` subdirectory. You may have a look at some binary or gray-level 3d image with `volViewer`. In `build` directory, you may try:

```
# A 3d shape (binary image)
./volViewer ../data/fertility-128.vol
# A thoracic CT scan, showing fat, lungs, heart, and the pulmonary vascular system.
./volViewer ../data/PA05-img-low.vol
# The segmentation of its pulmonary vascular system.
./volViewer ../data/PA05-lbl-low.vol
```

## The tutorials

After cloning the code, you may now have a look at each practical

* TP1 [Filtering and segmentation of 3D CT thoracic image](https://codimd.math.cnrs.fr/s/VSydn_j1G)
    > (ready) 3d image filtering, digital surfaces
* TP2 [Homotopic thinning](https://codimd.math.cnrs.fr/s/bVqUMFbM5) (wip)
    > (wip) 3d shapes, digital topology, cubical complexes

* TP3 [Geometric estimations on 3d digital surfaces](https://codimd.math.cnrs.fr/s/s2pNRQuga)
    > (ready) 3d shapes, geometric estimators, multigrid convergence


## Troubleshooting

Normally, `cmake` should install automatically `DGtal` and required
dependencies like `boost`, `zlib`, `polyscope`, `eigen`. However I
have not checked if this work on all platforms. Here are a few more
hints for specialized installations.

To build DGtal related examples, you would need:
  - C++17 enabled compiler (most c++ compilers are C+17)
  - a [cmake](https://cmake.org) client (at least 3.25)
  - [boost](http://boost.org) (>= 1.50).
  - the zlib package (already installed in many linux/macos distributions)

On linux, deps could be installed with:
```
sudo apt-get cmake git libboost-all-dev build-essential
```

For windows users, you can follow these [installation](https://github.com/DGtal-team/DGtal-Tutorials-DGMM2022/blob/main/windowsDGtalInstall.md) steps.

For polyscope based practicals, you may need X11/OpenGL headers (e.g. `sudo apt-get install xorg-dev libglu1-mesa-dev freeglut3-dev mesa-common-dev` on ubuntu). For more details, check [polyscope Building page](https://polyscope.run/building/).

From cmake, you can create the practicals, either from the commandline, or from a cmake GUI. E.g., from the commandline (using a Makefile target):

```
mkdir build
cd build
cmake ..
make
```

(`cmake .. -G Xcode` for an Xcode project, or `cmake .. -G "CodeBlocks - Unix Makefiles"` for a codeblocks one, or VisualStudio on MS Windows using the GUI)


By default, `cmake` will clone a copy of the DGtal repository, set up all the dependencies and build a first `helloworld` program.
