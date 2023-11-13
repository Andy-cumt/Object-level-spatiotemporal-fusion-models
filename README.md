# Object-level-spatiotemporal-fusion-models

A flexible object-level (OL) processing strategy is proposed to raise the efficiency of several popluar spatiotemporal fusion methods including STARFM, ESTARFM, Fit-FC, FSDAF and FSDAF 2.0. The OL fusion versions of STARFM, ESTARFM, and Fit-FC can better preserve the structural information, and were 102.89–113.71, 92.77–115.73, and 30.51–36.15 times faster than their original methods (the 1.1 version of OL-Fit-FC is about 50 times faster than Fit-FC).

In addition, an object-level hybrid spatiotemporal fusion method (OL-HSTFM) were proposed, which incorporates the efficiency advantage of object-level fusion strategy, spectral accuracy advantage of the three-step method (Fit-FC), and the spatial accuracy advantage of the spatial and temporal adaptive reflectance fusion model (STARFM).

Please double-click OLfusionv1.1.exe to install OL-fusionAPP and MATLAB Runtime version 9.9 (R2020b). Make sure you are connected to the network and anti-virus software (e.g., 360) is turned off during installation.

OL-STARFM can run on Google Earth engine:https://code.earthengine.google.com/6d9409377fe18e8a974e169ccfadadbb

D. Guo, W. Shi, H. Zhang and M. Hao, "A Flexible Object-Level Processing Strategy to Enhance the Weight Function-Based Spatiotemporal Fusion Method," in IEEE Transactions on Geoscience and Remote Sensing, vol. 60, pp. 1-11, 2022, Art no. 4414811, doi: 10.1109/TGRS.2022.3212474.

D. Guo and W. Shi, "Object-Level Hybrid Spatiotemporal Fusion: Reaching a Better Trade-Off Among Spectral Accuracy, Spatial Accuracy and Efficiency," in IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing, doi: 10.1109/JSTARS.2023.3310195.

The complete datasets, quantitative comparison data for OL—fusion experiments are available in https://github.com/Andy-cumt/Spatiotemporal-fusion-data

Paper Errata:
There is a clerical error in Eq.(17) in paper "A Flexible Object-Level Processing Strategy to Enhance the Weight Function-Based Spatiotemporal Fusion Method," where the fine image of the predicted phase should be involved in the calculation, rather than the base phase.
