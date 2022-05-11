# MultiDIC: a MATLAB Toolbox for Multi-View 3D Digital Image Correlation   

- [Summary](#Summary)  
- [Installation](#Installation)  
- [Getting started](#Start)
- [Citing](#Cite)
- [Contributing](#Contributing)  
- [Application Highlights](#Applications)
- [Publications](#Publications)
- [License](#License)  

## Summary <a name="Summary"></a>
This is the official repository associated with the paper: [MultiDIC: an Open-Source Toolbox for Multi-View 3D Digital Image Correlation](https://ieeexplore.ieee.org/document/8371235/) (DOI:[ 10.1109/ACCESS.2018.2843725](https://ieeexplore.ieee.org/document/8371235/)).  
MultiDIC is an open-source MATLAB toolbox by [Dana Solav](https://www.solavlab.com/). Three-dimensional (stereo) Digital Image Correlation (3D-DIC) is an important technique for measuring the mechanical behavior of materials. MultiDIC was developed to allow fast calibration and data merging even for a large number of cameras, and to be easily adaptable to different experimental requirements. If you're interested in using only 2 cameras, you should check out [DuoDIC](https://github.com/SolavLab/DuoDIC/). MultiDIC integrates the 2D-DIC subset-based software [Ncorr](https://www.github.com/justinblaber/ncorr_2D_matlab) with several camera calibration algorithms to reconstruct 3D surfaces from multiple stereo image pairs. Moreover, it contains algorithms for merging multiple surfaces, and for computing and visualizing 3D displacement, deformation and strain measures. High-level scripts allow users to perform 3D-DIC analyses with minimal interaction with MATLAB syntax, while proficient MATLAB users can also use stand-alone functions and data-structures to write custom scripts for specific experimental requirements. Comprehensive documentation, [instruction manual](https://github.com/MultiDIC/MultiDIC/blob/master/docs/pdf/MultiDIC_v_1_1_0_instruction_manual.pdf), and [sample data](https://github.com/MultiDIC/MultiDIC/tree/master/sample_data) are included.  
Check out this [short video](https://www.youtube.com/watch?v=DC9ifDJ7lvo&t) for a quick demonstration of MultiDIC's capabilities:  

[![](https://img.youtube.com/vi/DC9ifDJ7lvo/0.jpg)](https://www.youtube.com/watch?v=DC9ifDJ7lvo&t)

## Installation <a name="Installation"></a>  
### System Requirements
MultiDIC was developed and tested on 64-bit Windows 10 and has not yet been tested on other platforms.        
#### MATLAB
MultiDIC was developed on MATLAB versions R2017a and R2017b, and has not yet been tested on prior versions.  

MATLAB toolbox dependencies:
* Image Processing Toolbox
* Computer Vision System Toolbox
* Statistics and Machine Learning Toolbox


### Installation Instructions
To install MultiDIC simply follow these two steps:
#### 1. Get a copy of MultiDIC
Use **one** of these two options:      
**a.** Download and unzip the latest [zip file](https://github.com/MultiDIC/MultiDIC/archive/master.zip).   
**b.** Clone MultiDIC using: `git clone https://github.com/MultiDIC/MultiDIC/.git`.

#### 2. Install (or add to path)    
Use **one** of these two options:          
**a.** In MATLAB, navigate to the unzipped MultiDIC folder, type `installMultiDIC` in the command window, and hit Enter.   
**b.** In MATLAB, Add the MultiDIC folder (with subfolders) to the path and save the path definitions.


## Getting started <a name="Start"></a>
Check out the [instruction manual](https://github.com/MultiDIC/MultiDIC/blob/master/docs/pdf/MultiDIC_v_1_1_0_instruction_manual.pdf). It should have all the information you need to get started.


## Citing <a name="Cite"></a>   
This is the official repository for:

[MultiDIC: an Open-Source Toolbox for Multi-View 3D Digital Image Correlation](https://ieeexplore.ieee.org/document/8371235/)   
Dana Solav, Kevin M. Moerman, Aaron M. Jaeger, Katia Genovese, Hugh M. Herr   
DOI: [10.1109/ACCESS.2018.2843725](https://ieeexplore.ieee.org/document/8371235/)

Please cite this paper if you use the toolbox.


## Contributing <a name="Contributing"></a>   
If you wish to contribute code/algorithms to this project, or to propose a collaboration study, please send an email to danask (at) mit (dot) edu.

## Application Highlights <a name="Applications"></a>
### These are some examples of figures obtained using MultiDIC:
#### Shape change and skin deformation induced by joint movement and muscle contraction
<img src="docs/img/Shank2D_corr_204_205.gif">     
<img src="docs/img/ShankFull_L1_L2.gif">       

#### Shape change and skin deformation induced by indentation    
<img src="docs/img/indentation_204_205_DispMgn_onImages.gif">
<img src="docs/img/indentation_3D_Lamda1_Lamda2_reducedLight.gif">

## Publications <a name="Publications"></a>
### MultiDIC has been used in the following publications:
1. [MultiDIC: An Open-Source Toolbox for Multi-View 3D Digital Image Correlation](https://ieeexplore.ieee.org/abstract/document/8371235)
2. [A framework for measuring the time-varying shape and full-field deformation of residual limbs using 3D digital image correlation](https://ieeexplore.ieee.org/document/8625546)
3. [Decoding of facial strains via conformable piezoelectric interfaces](https://www.nature.com/articles/s41551-020-00612-w)
4. [Seeing the Big Picture: Improving The Prosthetic Design Cycle Using 360° 3D Digital Image Correlation](https://www.techrxiv.org/articles/preprint/Seeing_the_Big_Picture_Improving_The_Prosthetic_Design_Cycle_Using_360_3D_Digital_Image_Correlation/12705722)
5. [Implementing 3D Digital Image Correlation to Study the Dynamics of Globally-Coupled Multistable Structures](https://kb.osu.edu/handle/1811/87248)
6. [Using Digital Image Correlation to Quantify Skin Deformation with Von Frey Monofilaments](https://ieeexplore.ieee.org/abstract/document/9663031/)
7. [A Comparative Study of 2D and 3D Digital Image Correlation Approaches for the Characterization and Numerical Analysis of Composite Materials](https://ieeexplore.ieee.org/abstract/document/9634011)
8. [Stereo Digital Image Correlation in MATLAB
](https://www.mdpi.com/2076-3417/11/11/4904/htm)
9. [Using Multi-Dimensional Digital Image Correlation to Measure Mechanical Properties of Additively Manufactured and Anisotropic Materials](https://www.proquest.com/docview/2550639051?pq-origsite=gscholar&fromopenview=true)
10. [Deformation Measurements of a Full Span Model with Adaptive Elasto-Flexible Membrane Wings](https://link.springer.com/chapter/10.1007/978-3-030-79561-0_51)



## License <a name="License"></a>
MultiDIC is provided under the [Apache-2.0 license](https://www.apache.org/licenses/). The [license file](https://www.github.com/MultiDIC/MultiDIC/blob/master/LICENSE) is found on the GitHub repository.
