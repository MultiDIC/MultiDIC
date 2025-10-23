# MultiDIC: a MATLAB Toolbox for Multi-View 3D Digital Image Correlation   

- [Summary](#summary)  
- [Installation](#installation)  
- [Getting started](#getting-started)
- [Citing](#citing)
- [Contributing](#contributing)  
- [License](#license)  
- [Application Highlights](#application-highlights)
- [Publications](#publications)

## Summary
This is the official repository associated with the paper: [MultiDIC: an Open-Source Toolbox for Multi-View 3D Digital Image Correlation](https://ieeexplore.ieee.org/document/8371235/) (DOI:[ 10.1109/ACCESS.2018.2843725](https://ieeexplore.ieee.org/document/8371235/)).  
MultiDIC is an open-source MATLAB toolbox by [Dana Solav](https://www.solavlab.com/). Three-dimensional (stereo) Digital Image Correlation (3D-DIC) is an important technique for measuring the mechanical behavior of materials. MultiDIC was developed to allow fast calibration and data merging even for a large number of cameras, and to be easily adaptable to different experimental requirements. If you're interested in using only 2 cameras, you should check out [DuoDIC](https://github.com/SolavLab/DuoDIC/). MultiDIC integrates the 2D-DIC subset-based software [Ncorr](https://www.github.com/justinblaber/ncorr_2D_matlab) with several camera calibration algorithms to reconstruct 3D surfaces from multiple stereo image pairs. Moreover, it contains algorithms for merging multiple surfaces, and for computing and visualizing 3D displacement, deformation and strain measures. High-level scripts allow users to perform 3D-DIC analyses with minimal interaction with MATLAB syntax, while proficient MATLAB users can also use stand-alone functions and data-structures to write custom scripts for specific experimental requirements. Comprehensive documentation, [instruction manual](https://github.com/MultiDIC/MultiDIC/blob/master/docs/pdf/MultiDIC_v_1_1_0_instruction_manual.pdf), and [sample data](https://github.com/MultiDIC/MultiDIC/tree/master/sample_data) are included.  
Check out this [short video](https://www.youtube.com/watch?v=DC9ifDJ7lvo&t) for a quick demonstration of MultiDIC's capabilities:  

[![](https://img.youtube.com/vi/DC9ifDJ7lvo/0.jpg)](https://www.youtube.com/watch?v=DC9ifDJ7lvo&t)

## Installation
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


## Getting started
Check out the [instruction manual](https://github.com/MultiDIC/MultiDIC/blob/master/docs/pdf/MultiDIC_v_1_1_0_instruction_manual.pdf). It should have all the information you need to get started.


## Citing
This is the official repository for:

[MultiDIC: an Open-Source Toolbox for Multi-View 3D Digital Image Correlation](https://ieeexplore.ieee.org/document/8371235/)   
Dana Solav, Kevin M. Moerman, Aaron M. Jaeger, Katia Genovese, Hugh M. Herr   
DOI: [10.1109/ACCESS.2018.2843725](https://ieeexplore.ieee.org/document/8371235/)

Please cite this paper if you use the toolbox.


## Contributing
If you wish to contribute code/algorithms to this project, or to propose a collaboration study, please send an email to danask (at) mit (dot) edu.

## License
MultiDIC is provided under the [Apache-2.0 license](https://www.apache.org/licenses/). The [license file](https://www.github.com/MultiDIC/MultiDIC/blob/master/LICENSE) is found on the GitHub repository.

## Application Highlights
### These are some examples of figures obtained using MultiDIC:
#### Shape change and skin deformation induced by joint movement and muscle contraction
<img src="docs/img/Shank2D_corr_204_205.gif">     
<img src="docs/img/ShankFull_L1_L2.gif">       

#### Shape change and skin deformation induced by indentation    
<img src="docs/img/indentation_204_205_DispMgn_onImages.gif">
<img src="docs/img/indentation_3D_Lamda1_Lamda2_reducedLight.gif">

## Publications
### MultiDIC has been used in the following publications:
1. Solav, Dana, et al. "MultiDIC: An open-source toolbox for multi-view 3D digital image correlation." _IEEE Access_ 6 (2018): 30520-30535. [https://doi.org/10.1109/ACCESS.2018.2843725](https://doi.org/10.1109/ACCESS.2018.2843725)
1. Solav, Dana, et al. "A framework for measuring the time-varying shape and full-field deformation of residual limbs using 3-D digital image correlation." _IEEE Transactions on Biomedical Engineering_ 66.10 (2019): 2740-2752. [https://doi.org/10.1109/TBME.2019.2895283](https://doi.org/10.1109/TBME.2019.2895283)
1. Sun, Tao, et al. "Decoding of facial strains via conformable piezoelectric interfaces." _Nature biomedical engineering_ 4.10 (2020): 954-972. [https://doi.org/10.1038/s41551-020-00612-w](https://doi.org/10.1038/s41551-020-00612-w)
1. Cabrera, Isaac A., et al. "Seeing the Big Picture: Improving the prosthetic design cycle using 360 3D digital image correlation." (2020). [https://doi.org/10.36227/techrxiv.12705722.v1](https://doi.org/10.36227/techrxiv.12705722.v1)
1. DeNise, Adam. Implementing 3D Digital Image Correlation to Study the Dynamics of Globally-Coupled Multistable Structures. _Diss. The Ohio State University_, 2019. [https://kb.osu.edu/handle/1811/87248](https://kb.osu.edu/handle/1811/87248)
1. Kao, Anika Rothlin, Chang Xu, and Gregory John Gerling. "Using Digital Image Correlation to Quantify Skin Deformation with Von Frey Monofilaments." _IEEE Transactions on Haptics_ (2021). [https://doi.org/10.1109/TOH.2021.3138350](https://doi.org/10.1109/TOH.2021.3138350)
1. Pisonero, Javier, et al. "A Comparative Study of 2D and 3D Digital Image Correlation Approaches for the Characterization and Numerical Analysis of Composite Materials." _IEEE Access_ 9 (2021): 160675-160687. [https://doi.org/10.1109/ACCESS.2021.3132393](https://doi.org/10.1109/ACCESS.2021.3132393)
1. Atkinson, Devan, and Thorsten Hermann Becker. "Stereo Digital Image Correlation in MATLAB." _Applied Sciences_ 11.11 (2021): 4904. [https://doi.org/10.3390/app11114904](https://doi.org/10.3390/app11114904)
1. Hebert, John J. "Using Multi-Dimensional Digital Image Correlation to Measure Mechanical Properties of Additively Manufactured and Anisotropic Materials."" _Diss. University of Louisiana at Lafayette_, 2020. [Link to thesis](https://www.proquest.com/docview/2550639051?pq-origsite=gscholar&fromopenview=true)
1. Pflüger, Jonathan, Yuting Chen, and Christian Breitsamter. "Deformation Measurements of a Full Span Model with Adaptive Elasto-Flexible Membrane Wings." _STAB/DGLR Symposium_. Springer, Cham, 2020. [https://doi.org/10.1007/978-3-030-79561-0_51](https://doi.org/10.1007/978-3-030-79561-0_51)
1. Gacek, Elizabeth, Arin M. Ellingson, and Victor H. Barocas. "In Situ Lumbar Facet Capsular Ligament Strains Due to Joint Pressure and Residual Strain." _Journal of biomechanical engineering_ 144.6 (2022): 061007. [https://doi.org/10.1115/1.4053993](https://doi.org/10.1115/1.4053993)
1. Li, Yuanzhe, et al. "Imaging dynamic three-dimensional traction stresses." _Science advances_ 8.11 (2022): eabm0984. [https://doi.org/10.1126/sciadv.abm0984](https://doi.org/10.1126/sciadv.abm0984)
1. Jung, Yei Hwan, et al. "A wireless haptic interface for programmable patterns of touch across large areas of the skin." _Nature Electronics_ (2022): 1-12. [https://doi.org/10.1038/s41928-022-00765-3](https://doi.org/10.1038/s41928-022-00765-3)
1. Álvarez-Blanco, Mario, et al. "Crack morphology in lattice-core specimens made of biopolymer via fused deposition modelling." _Procedia Structural Integrity_ 39 (2022): 379-386. [https://doi.org/10.1016/j.prostr.2022.03.106](https://doi.org/10.1016/j.prostr.2022.03.106)
1. Shiraishi, Yasuyuki, et al. "In Vitro Modelling for Bulging Sinus Effects of an Expanded Polytetrafluoroethylene Valved Conduit Based on High-Speed 3D Leaflet Evaluation." _2022 44th Annual International Conference of the IEEE Engineering in Medicine & Biology Society (EMBC). IEEE_, 2022. [https://doi.org/10.1109/EMBC48229.2022.9871676](https://doi.org/10.1109/EMBC48229.2022.9871676)
1. Ni, Xinchen, et al. "Soft shape-programmable surfaces by fast electromagnetic actuation of liquid metal networks." _Nature communications_ 13.1 (2022): 1-9. [https://doi.org/10.1038/s41467-022-31092-y](https://doi.org/10.1038/s41467-022-31092-y)
1. Kang, Youn J., et al. "Soft skin-interfaced mechano-acoustic sensors for real-time monitoring and patient feedback on respiratory and swallowing biomechanics." _npj Digital Medicine_ 5.1 (2022): 1-13. [https://doi.org/10.1038/s41746-022-00691-w](https://doi.org/10.1038/s41746-022-00691-w)
1. Sun, Weifu, et al. "Digital image correlation-aided non-destructive buckling load prediction of cylindrical shells." _International Journal of Solids and Structures_ 254 (2022): 111941.[https://doi.org/10.1016/j.ijsolstr.2022.111941](https://doi.org/10.1016/j.ijsolstr.2022.111941)
1. Lee, Duncan RC. "Design and Clinical Evaluation of a Digital Transtibial Prosthetic Interface." _Diss. Massachusetts Institute of Technology_, 2022. [Link to thesis](https://dspace.mit.edu/handle/1721.1/144501)
1. Zhu, Xing, et al. "Identification of crack initiation and damage thresholds in sandstone using 3D digital image correlation." _Theoretical and Applied Fracture Mechanics_ 122 (2022): 103653. [https://doi.org/10.1016/j.tafmec.2022.103653](https://doi.org/10.1016/j.tafmec.2022.103653)
1. Kim, Jin-Tae, et al. "Dynamics of plosive consonants via imaging, computations, and soft electronics." _Proceedings of the National Academy of Sciences_ 119.46 (2022): e2214164119. [https://doi.org/10.1073/pnas.2214164119](https://doi.org/10.1073/pnas.2214164119)
1. Kim, Jin-Tae, et al. "Mechanics of vibrotactile sensors for applications in skin-interfaced haptic systems." _Extreme Mechanics Letters_ 58 (2023): 101940. [https://doi.org/10.1016/j.eml.2022.101940](https://doi.org/10.1016/j.eml.2022.101940)
1. Ghadimi, Hamed, et al. "Effects of Printing Layer Orientation on the High-Frequency Bending-Fatigue Life and Tensile Strength of Additively Manufactured 17-4 PH Stainless Steel." _Materials_ 16.2 (2023): 469. [https://doi.org/10.1016/j.eml.2022.101940](https://doi.org/10.1016/j.eml.2022.101940)
1. Álvarez-Blanco, Mario, et al. "Influence of material extrusion parameters on fracture mechanisms of polylactic acid under three-point bending." _Engineering Fracture Mechanics_ 283 (2023): 109223. [https://doi.org/10.1016/j.engfracmech.2023.109223](https://doi.org/10.1016/j.engfracmech.2023.109223)
1. Huang, Senbin, et al. "Experimental investigations of optimized 3D Printing Planar X-joints manufactured by stainless steel and high-strength steel." _Engineering Structures_ 285 (2023): 116054. [https://doi.org/10.1016/j.engstruct.2023.116054](https://doi.org/10.1016/j.engstruct.2023.116054)
1. Arias-Blanco, Adrián, et al "Experimental and numerical analysis of the influence of intramedullary nail position on the cut-out phenomenon." _Computer Methods and Programs in Biomedicine_ (2023): 107734.[https://doi.org/10.1016/j.cmpb.2023.107734](https://doi.org/10.1016/j.cmpb.2023.107734)
1. Ghadimi, Hamed, et al. "Effects of printing layer orientation on the high-frequency bending-fatigue life and tensile strength of additively manufactured 17-4 PH stainless steel." Materials 16.2 (2023): 469. [https://doi.org/10.3390/ma16020469](https://doi.org/10.3390/ma16020469)
1. Kao, Anika R., M. Terry Loghmani, and Gregory J. Gerling. "Strain-based biomarkers at the skin surface differentiate asymmetries in soft tissue mobility associated with myofascial pain." _Journal of the Mechanical Behavior of Biomedical Materials_ 172, 107175 (2025). [https://doi.org/10.1016/j.jmbbm.2025.107175](https://doi.org/10.1016/j.jmbbm.2025.107175)
1. Pflueger, Jonathan, and Christian Breitsamter. "Gust response of an elasto-flexible morphing wing using fluid–structure interaction simulations." _Chinese Journal of Aeronautics_ 37.2 (2024): 45-57. [https://doi.org/10.1016/j.cja.2023.12.017](https://doi.org/10.1016/j.cja.2023.12.017)
1. Lin, Yu-Liang, et al. "Shaking table experiment on seismic response of a three-stage slope supported by anchoring lattice beam." _Soil Dynamics and Earthquake Engineering_ 187 (2024): 109003. [https://doi.org/10.1016/j.soildyn.2024.109003](https://doi.org/10.1016/j.soildyn.2024.109003)
1. Liu, Claire, et al. "Multifunctional materials strategies for enhanced safety of wireless, skin‐interfaced bioelectronic devices." _Advanced Functional Materials_ 33.34 (2023): 2302256. [https://doi.org/10.1002/adfm.202302256](https://doi.org/10.1002/adfm.202302256)
1. Cheng, Zhanqi, et al. "Crack development and fracture performance of recycled powder engineered cementitious composites (ECC): experimental investigation." _Construction and Building Materials_ 429 (2024): 136388. [https://doi.org/10.1016/j.conbuildmat.2024.136388](https://doi.org/10.1016/j.conbuildmat.2024.136388)
1. Kim, Jae-Hwan, et al. "A wirelessly programmable, skin-integrated thermo-haptic stimulator system for virtual reality." _Proceedings of the National Academy of Sciences_ 121.22 (2024): e2404007121. [https://doi.org/10.1073/pnas.2404007121](https://doi.org/10.1073/pnas.2404007121)
1. Lee, Duncan RC, et al. "A clinical comparison of a digital versus conventional design methodology for transtibial prosthetic interfaces." _Scientific Reports_ 14.1 (2024): 25833. [https://doi.org/10.1038/s41598-024-74504-3](https://doi.org/10.1038/s41598-024-74504-3)
1. Liu, Kaixin, et al. "Macro/meso-evolution of shear behavior of sand-steel interface in pipe jacking via digital image correlation (DIC) technology." _Tunnelling and Underground Space Technology_ 155 (2025): 106197. [https://doi.org/10.1016/j.tust.2024.106197](https://doi.org/10.1016/j.tust.2024.106197)
1. López-Rebollo, Jorge, et al. "Digital Image Correlation and Reliability-Based Methods for the Design of Structural Beams Made from Recycled Concrete Using Aggregates from Precast Rejects." _Applied Sciences_ 15.2 (2025): 656. [https://doi.org/10.3390/app15020656](https://doi.org/10.3390/app15020656)
1. Doumont, Donatien, et al. "3-D Reconstruction of Fingertip Deformation during Contact Initiation." _bioRxiv_ (2025): 2025-02. [https://doi.org/10.1101/2025.02.10.637366](https://doi.org/10.1101/2025.02.10.637366)
1. Yuan, Shaoke, et al. "Design, Modeling, and Control of a Soft Abdominal Compression Robot for Respiratory Assistance." _Soft Robotics_ (2025). [https://doi.org/10.1089/soro.2024.0088](https://doi.org/10.1089/soro.2024.0088)
1. Gautam, Kushagra, et al. "Using Digital Image Correlation to Analyze Deformation in Wood-based Liquid Deposition Modelling."(Proceedings of the 30th International Conference of the Association for Computer-Aided Architectural Design Research in Asia (CAADRIA) 2025, Volume 2, 49-58. [https://doi.org/10.52842/conf.caadria.2025.2.049](https://doi.org/10.52842/conf.caadria.2025.2.049)
1. Kirk, C. D., et al. "An Integrated High-Speed Microstructural Characterization Method Using Simultaneous XRD, Stereo-DIC, and PCI." _Journal of Dynamic Behavior of Materials_ (2025): 1-12.[ https://doi.org/10.1007/s40870-025-00476-8](https://doi.org/10.1007/s40870-025-00476-8)
1. Sato, Kai, and Mayuko Nishio. "Uncertainty quantification using digital-image-correlation data for load-carrying capacity analysis of structural members with local damages." _Structures_. Vol. 79. Elsevier, 2025. [https://doi.org/10.1016/j.istruc.2025.109431](https://doi.org/10.1016/j.istruc.2025.109431)
1. Bertsch, J. Michael, Anika R. Kao, and Gregory J. Gerling. "Capturing 3D Skin Deformation and Finger Kinematics in Active Touch with Compliant Materials." _2025 IEEE World Haptics Conference (WHC). IEEE_, 2025. [https://doi.org/10.1109/WHC64065.2025.11123409](https://doi.org/10.1109/WHC64065.2025.11123409)
1. Bai, Lu, et al. "Study on the uniaxial tensile fatigue performance of magnesium phosphate cement-based engineered cementitious composites." _Journal of Building Engineering_ (2025): 114168. [https://doi.org/10.1016/j.jobe.2025.114168](https://doi.org/10.1016/j.jobe.2025.114168)
1. Cheng, Zhanqi, et al. "Experimental study on crack development and fracture performance of 3D printed green engineered cementitious composites (ECC)." _Theoretical and Applied Fracture Mechanics_ 141 (2026): 105292. [https://doi.org/10.1016/j.tafmec.2025.105292](https://doi.org/10.1016/j.tafmec.2025.105292)
