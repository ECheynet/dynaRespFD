# dynaRspFD
The dynamic response of a suspension bridge to wind turbulence is computed in the frequency domain.

[![View Buffeting response of a suspension bridge (frequency domain) on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://se.mathworks.com/matlabcentral/fileexchange/51970-buffeting-response-of-a-suspension-bridge-frequency-domain)
[![DOI](https://zenodo.org/badge/249144990.svg)](https://zenodo.org/badge/latestdoi/249144990)


The estimation of the displacement response of a large civil engineering structure to wind turbulence is based on the buffeting theory [1, 2, 5]. Ref. [5] contains the theoretical background I have used for the function dynaRespFD3. In the present script, the structure in question is a suspension bridge modelled using the theory of continuous beams [3]. The buffeting response is computed in the frequency domain using the quasi-steady theory.  Modal coupling was assumed negligible, which is generally well verified for most of the wind velocities recorded in full scale [4]. The present script is a  simplified version of the one used in [6]. 

The present script computes the lateral, vertical and torsional displacement response. A multi-modes approach is used. Some knowledge in the field of random vibration analysis and wind loading on structures are advised for proper use of this script. 


The present submission contains 

•	dynaRespFD3.m : Function that calculates the displacement response spectrum of the bridge

•	Two example files Example_1.m and Example_2.m 

•	Two .mat files bridgeModalProperties.mat and DynamicDispl.mat that are used in the  2 examples. 

Any question, comment or suggestion to improve the submission is welcomed.


References

 [1] Davenport, A.G., The response of slender line-like structures to a gusty wind, Proceedings of the Institution of Civil Engineers, Vol. 23, 1962, pp. 389 – 408. 
 
[2] Scanlan, R. H. (1978). The action of flexible bridges under wind, II: Buffeting theory. Journal of Sound and vibration, 60(2), 201-211.

[3] http://www.mathworks.com/matlabcentral/fileexchange/51815-suspension-bridge--eigen-frequency-and-mode-shapes-benchmark-solutions 

[4] Thorbek, L. T., & Hansen, S. O. (1998). Coupled buffeting response of suspension bridges. Journal of Wind Engineering and Industrial Aerodynamics, 74, 839-847.

[5] Hjorth-Hansen, E. (1993). Fluctuating drag, lift and overturning moment for a line-like structure predicted (primarily) from static, mean loads. Wind Engineering, Lecture note no, 2.

[6] Cheynet, E., Jakobsen, J. B., & Snæbjörnsson, J. (2016). Buffeting response of a suspension bridge in complex terrain. Engineering Structures, 128, 474-487.   http://dx.doi.org/10.1016/j.engstruct.2016.09.060
