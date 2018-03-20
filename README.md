# MultiDIC: a MATLAB Toolbox for Multi-View 3D Digital Image Correlation

This is the offical repository for:

```
MultiDIC: a MATLAB Toolbox for Multi-View 3D Digital Image Correlation
Dana Solav, Kevin M. Moerman, Aaron M. Jaeger, Hugh M. Herr
DOI: TBD
```

Please cite this paper if you use this toolbox.

# Table of contents
- [Project Summary](#Summary)  
- [Installation](#Installation)  
- [Getting started](#Start)
- [License](#License)  
- [Citing](#Cite)
- [Contributing](#Contributing)  


# Project summary <a name="Summary"></a>
MultiDIC (Multi Digital Image Correlation) is an open-source MATLAB toolbox by [Dana Solav](media.mit.edu/people/danask/). Three-dimensional (stereo) Digital Image Correlation (3D-DIC) is an important technique for measuring the mechanical behavior of materials. MultiDIC was developed to allow fast calibration even for a large number of cameras, and be easily adaptable to different experimental requirements. It integrates the 2D-DIC subset-based software [Ncorr] (github.com/justinblaber/ncorr_2D_matlab) with several camera calibration algorithms to reconstruct 3D surfaces from multiple stereo image pairs. Moreover, it contains algorithms for merging multiple surfaces, and for computing and plotting displacement, deformation and strain measures. High-level scripts allow users to perform 3D-DIC analyses with minimal interaction with MATLAB syntax, while proficient MATLAB users can use stand-alone functions and data-structures to write custom scripts for specific experimental requirements. Comprehensive documentation, user guide, and sample data are included.

# Installation <a name="Installation"></a>  
To install Multi DIC simply follow these two steps:

### 1. Get a copy of MultiDIC
Use one of these two options:      
a. Download and unzip the latest [zip file](https://github.com/MultiDIC/MaDICT/archive/master.zip).   
b. Clone MultiDIC using: `git clone https://github.com/MultiDIC/MaDICT/.git`.

### 2. Install (or add to path)    
Use one of these two options:          
a. In Matlab, navigate to the unzipped MaDICT folder and type installMultiDIC in the command window, and hit enter.   
b. In Matlab, Add the MaDICT folder (with subfolders) to the path and save the path definitions.

# Getting started <a name="Start"></a>

# License <a name="License"></a>
MultiDIC is provided under:

# Citing <a name="Cite"></a>

# Contributing <a name="Contributing"></a>
