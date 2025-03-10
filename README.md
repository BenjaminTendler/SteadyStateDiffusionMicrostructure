#  Investigating time-independent and time-dependent diffusion phenomena using steady-state diffusion MRI - Software
This repository contains software to characterise the diffusion-weighted steady-state free precession (DW-SSFP) signal based on the proposed frameworks in "Tendler, Investigating time-independent and time-dependent diffusion phenomena using steady-state diffusion MRI, Scientific Reports 2025, [https://doi.org/10.1038/s41598-025-87377-x](https://doi.org/10.1038/s41598-025-87377-x)".

Details are as follows:

**TimeIndependent/Dependent...Example.m**

These scripts provides example implementations of the time independent & dependent framework proposed in the manuscript using both conventional and oscillating DW-SSFP gradients.

**Figures**

This folder contains the scripts used for synthesising many of the Figures in the manuscript. These scripts can be used to replicate data in the manuscript, or identify different modes of investigation with the proposed frameworks. 

**SupportingInformationFigures**

This folder contains the scripts used for synthesising many of the Supporting Information Figures in the manuscript. These scripts can be used to replicate data in the manuscript, or identify different modes of investigation with the proposed frameworks. 

**AnalyticalModels**

This folder contains the scripts used to estimate the DW-SSFP signal for (1) diffusion gradients of a fixed duration, (2) oscillating diffusion gradients and (3) a diffusion tensor, alongside code to evaluate transverse period approximations of DW-SSFP. Implementations are described in Appendix [1-3] in the manuscript. They are based on "Freed et al., Steady-state free precession experiments and exact treatment of diffusion in a uniform gradient, J. Chem. Phys 2001, [10.1063/1.1389859](https://doi.org/10.1063/1.1389859)".

**MCOutputs**

This folder contains the Monte Carlo timeseries data for the DW-SSFP & DW-SE investigations displayed in Figure 8 of the manuscript.

**bin**

This folder contains the code database for the proposed frameworks. All code was generated in house, with the exception of Fourier-based convolution scripts (used in the time-independent framework) (Bruno Luong 2024; [FFT-based convolution](https://www.mathworks.com/matlabcentral/fileexchange/24504-fft-based-convolution), MATLAB Central File Exchange), and colormap scripts (used for generating the power-spectrum representations) (Matteo Courthoud 2021; [centered-colormap-master](https://github.com/matteocourthoud/centered-colormap?tab=readme-ov-file), GitHub).

**Notes**

You may need to compile the Fourier-based convolution script defined in bin/CONVNFFT. To achieve this, please run 'bin/CONVFFT/convnfft_install.m.

**Copyright**

Copyright, 2024, University of Oxford. All rights reserved

Any questions please contact benjamin.tendler@ndcn.ox.ac.uk

