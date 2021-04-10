<p align="center"><img width=80% src="https://github.com/Waller-Lab/SpectralDiffuserCam/blob/gh-pages/resources/website_overview_figure.jpg"></p>

# [Miniscope3D](https://waller-lab.github.io/Miniscope3D/)



## Paper 
[Miniscope3D: optimized single-shot miniature 3D fluorescence microscopy](https://www.nature.com/articles/s41377-020-00403-7)

Please cite the following paper when using this code or data:

Yanny, Kyrollos, et al. "Miniscope3D: optimized single-shot miniature 3D fluorescence microscopy." Light: Science & Applications 9.1 (2020): 1-13.

```
@article{yanny2020miniscope3d,
  title={Miniscope3D: optimized single-shot miniature 3D fluorescence microscopy},
  author={Yanny, Kyrollos and Antipa, Nick and Liberti, William and Dehaeck, Sam and Monakhova, Kristina and Liu, Fanglin Linda and Shen, Konlin and Ng, Ren and Waller, Laura},
  journal={Light: Science \& Applications},
  volume={9},
  number={1},
  pages={1--13},
  year={2020},
  publisher={Nature Publishing Group}
}

```


## Contents

1. [Data](#Data)
2. [Setup](#Setup)
3. [Description](#Description)

## Data
Sample data (needed to run the code) can be found [here](https://drive.google.com/drive/folders/1dmfzkTLFZZFUYW8GC6Vn6SOuZiZq47SS?usp=sharing)

This includes the following files:
 * calibration.mat - includes the calibratated point spread function, filter function, and wavelength list
 * four sample raw measurements


## Setup
Clone this project using: 
```
git clone https://github.com/Waller-Lab/SpectralDiffuserCam.git
```

The dependencies can be installed by using:
```
conda env create -f environment.yml
source activate SpectralDiffuserCam
```

Please place the downloaded data in SampleData folder in the Python and/or Matlab folders.

[Reconstruction Demo.ipynb](https://github.com/Waller-Lab/SpectralDiffuserCam/blob/master/Python/Reconstruction%20Demo.ipynb) contains an example reconstruction in Python. 

[reconstruction_demo.m](https://github.com/Waller-Lab/SpectralDiffuserCam/blob/master/Matlab/reconstruction_demo.m) contains an example reconstruction in Matlab.

We recommend running this code on a GPU, but it can also be run on a CPU (much slower!). 

## Description 
This repository contains code in both Python and Matlab that is needed to process raw Spectral DiffuserCam images and reconstruct 3D hyperspectral volumes from the raw 2D measurements.  Four example raw images are provided, along with the calibrated point spread function and spectral filter function. Both the Python and Matlab versions support GPU acceleration. In Python, this is accomplished using cupy. We use FISTA for our reconstructions with a 3D total variation prior. 
