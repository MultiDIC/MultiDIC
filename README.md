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
1. [Solav, Dana, et al. "MultiDIC: An open-source toolbox for multi-view 3D digital image correlation." _IEEE Access_ 6 (2018): 30520-30535.](https://ieeexplore.ieee.org/abstract/document/8371235)
2. [Solav, Dana, et al. "A framework for measuring the time-varying shape and full-field deformation of residual limbs using 3-D digital image correlation." _IEEE Transactions on Biomedical Engineering_ 66.10 (2019): 2740-2752.](https://ieeexplore.ieee.org/document/8625546)
3. [Sun, Tao, et al. "Decoding of facial strains via conformable piezoelectric interfaces." _Nature biomedical engineering_ 4.10 (2020): 954-972.](https://www.nature.com/articles/s41551-020-00612-w)
4. [Cabrera, Isaac A., et al. "Seeing the Big Picture: Improving the prosthetic design cycle using 360 3D digital image correlation." (2020).](https://www.techrxiv.org/articles/preprint/Seeing_the_Big_Picture_Improving_The_Prosthetic_Design_Cycle_Using_360_3D_Digital_Image_Correlation/12705722)
5. [DeNise, Adam. Implementing 3D Digital Image Correlation to Study the Dynamics of Globally-Coupled Multistable Structures. _Diss. The Ohio State University_, 2019.](https://kb.osu.edu/handle/1811/87248)
6. [Kao, Anika Rothlin, Chang Xu, and Gregory John Gerling. "Using Digital Image Correlation to Quantify Skin Deformation with Von Frey Monofilaments." _IEEE Transactions on Haptics_ (2021).](https://ieeexplore.ieee.org/abstract/document/9663031/)
7. [Pisonero, Javier, et al. "A Comparative Study of 2D and 3D Digital Image Correlation Approaches for the Characterization and Numerical Analysis of Composite Materials." _IEEE Access_ 9 (2021): 160675-160687.](https://ieeexplore.ieee.org/abstract/document/9634011)
8. [Atkinson, Devan, and Thorsten Hermann Becker. "Stereo Digital Image Correlation in MATLAB." Applied Sciences 11.11 (2021): 4904.](https://www.mdpi.com/2076-3417/11/11/4904/htm)
9. [Hebert, John J. "Using Multi-Dimensional Digital Image Correlation to Measure Mechanical Properties of Additively Manufactured and Anisotropic Materials."" _Diss. University of Louisiana at Lafayette_, 2020.](https://www.proquest.com/docview/2550639051?pq-origsite=gscholar&fromopenview=true)
10. [Pflüger, Jonathan, Yuting Chen, and Christian Breitsamter. "Deformation Measurements of a Full Span Model with Adaptive Elasto-Flexible Membrane Wings." _STAB/DGLR Symposium_. Springer, Cham, 2020.](https://link.springer.com/chapter/10.1007/978-3-030-79561-0_51)
11. [Gacek, Elizabeth, Arin M. Ellingson, and Victor H. Barocas. "In Situ Lumbar Facet Capsular Ligament Strains Due to Joint Pressure and Residual Strain." _Journal of biomechanical engineering_ 144.6 (2022): 061007.](https://asmedigitalcollection.asme.org/biomechanical/article/144/6/061007/1137926/In-Situ-Lumbar-Facet-Capsular-Ligament-Strains-Due?casa_token=-fX2MKVtvboAAAAA:Xg4p7p4Yngk9gVcIpgqUgXydaY_8OCdNjeAGut_qfLDG5MU9mLpEp5zWJsfZCrmm7CSqWeDLsQ)
12. [Li, Yuanzhe, et al. "Imaging dynamic three-dimensional traction stresses." _Science advances_ 8.11 (2022): eabm0984.](https://www.science.org/doi/full/10.1126/sciadv.abm0984)
13. [Jung, Yei Hwan, et al. "A wireless haptic interface for programmable patterns of touch across large areas of the skin." _Nature Electronics_ (2022): 1-12.](https://www.nature.com/articles/s41928-022-00765-3)
14. [Álvarez-Blanco, Mario, et al. "Crack morphology in lattice-core specimens made of biopolymer via fused deposition modelling." _Procedia Structural Integrity_ 39 (2022): 379-386.](https://www.sciencedirect.com/science/article/pii/S2452321622003195)
15. [Shiraishi, Yasuyuki, et al. "In Vitro Modelling for Bulging Sinus Effects of an Expanded Polytetrafluoroethylene Valved Conduit Based on High-Speed 3D Leaflet Evaluation." 2022 44th Annual International Conference of the IEEE Engineering in Medicine & Biology Society (EMBC). IEEE, 2022.](https://ieeexplore.ieee.org/abstract/document/9871676)
16. [Ni, Xinchen, et al. "Soft shape-programmable surfaces by fast electromagnetic actuation of liquid metal networks." Nature communications 13.1 (2022): 1-9.](https://www.nature.com/articles/s41467-022-31092-y)
17. [Kang, Youn J., et al. "Soft skin-interfaced mechano-acoustic sensors for real-time monitoring and patient feedback on respiratory and swallowing biomechanics." npj Digital Medicine 5.1 (2022): 1-13.](https://www.nature.com/articles/s41746-022-00691-w)
18. [Sun, Weifu, et al. "Digital image correlation-aided non-destructive buckling load prediction of cylindrical shells." International Journal of Solids and Structures 254 (2022): 111941.](https://www.sciencedirect.com/science/article/pii/S0020768322003961)
19. [Lee, Duncan RC. "Design and Clinical Evaluation of a Digital Transtibial Prosthetic Interface." Diss. Massachusetts Institute of Technology, 2022.](https://dspace.mit.edu/handle/1721.1/144501)
20. [Zhu, Xing, et al. "Identification of crack initiation and damage thresholds in sandstone using 3D digital image correlation." Theoretical and Applied Fracture Mechanics 122 (2022): 103653.](https://www.sciencedirect.com/science/article/pii/S0167844222003974)
21. [Kim, Jin-Tae, et al. "Dynamics of plosive consonants via imaging, computations, and soft electronics." Proceedings of the National Academy of Sciences 119.46 (2022): e2214164119.](https://www.pnas.org/doi/full/10.1073/pnas.2214164119)
22. [Kim, Jin-Tae, et al. "Mechanics of vibrotactile sensors for applications in skin-interfaced haptic systems." Extreme Mechanics Letters 58 (2023): 101940.](https://www.sciencedirect.com/science/article/pii/S2352431622002164)
23. [Ghadimi, Hamed, et al. "Effects of Printing Layer Orientation on the High-Frequency Bending-Fatigue Life and Tensile Strength of Additively Manufactured 17-4 PH Stainless Steel." Materials 16.2 (2023): 469.](https://www.mdpi.com/1996-1944/16/2/469)
24. [Álvarez-Blanco, Mario. "Influence of material extrusion parameters on fracture mechanisms of polylactic acid under three-point bending." Engineering Fracture Mechanics (2023): 109223.](https://www.sciencedirect.com/science/article/pii/S0013794423001819)
25. [Huang, Senbin, et al. "Experimental investigations of optimized 3D Printing Planar X-joints manufactured by stainless steel and high-strength steel." Engineering Structures 285 (2023): 116054.](https://www.sciencedirect.com/science/article/pii/S0141029623004686)


## License <a name="License"></a>
MultiDIC is provided under the [Apache-2.0 license](https://www.apache.org/licenses/). The [license file](https://www.github.com/MultiDIC/MultiDIC/blob/master/LICENSE) is found on the GitHub repository.
