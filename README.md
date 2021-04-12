<p align="center"><img width=100% src="https://waller-lab.github.io/Miniscope3D/Miniscope3D_webpage_files/setup.png"></p>

# [Miniscope3D](https://waller-lab.github.io/Miniscope3D/)



## Paper 
[Miniscope3D: optimized single-shot miniature 3D fluorescence microscopy](https://www.nature.com/articles/s41377-020-00403-7)

Please cite the following paper when using this code or data:


```
Yanny, Kyrollos, et al. "Miniscope3D: optimized single-shot miniature 3D fluorescence microscopy." Light: Science & Applications 9.1 (2020): 1-13.

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
git clone https://github.com/Waller-Lab/Miniscope3D.git
```
### 1) Adaptive Stitching with Nanoscribe
This code performs adaptive stitching for randomly spaced microlenses that are multifocal and have aberrations added to them. The adaptive stitching is needed for the microlenses to be correctly 3D printed with Nanoscribe. The code places the stitching artifacts at the boundaries of the microlenses and thus maximizes optical quality. 

The demonstration here will be for the 36 multifocal microlenses with astigmatism and tilt added. This is the phase mask used in the paper. The code also provides the necessary .gwl files needed for the print-job on Nanoscribe

### 2) 3D Reconstruction
Reconstruction code will be available soon.

### 3) Microlens Phase Mask Optimization
Phase mask optimization code will be available soon.

## Description 
This repository contains code for the Miniscope3D paper. The code includes: 1) 3D Reconstruction with Shif-Varying model (matlab). This code reconstructs 3D volumes from a single 2D image. 2) Microlens phase mask optimization (python). This code optimizes a multi-focal randomly spaced mircolens array by addting aberrations and changing focal lengths and positions to imrpove 3D reconstructions (by minimizing PSF cross correlations at different depths). 3) Adaptive Stitching with Nanoscribe (python). This code performs the adaptive stitching needed to 3D print the phase mask with Nanoscribe.
