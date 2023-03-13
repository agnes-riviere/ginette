![Ginette-2](Ginette-2.png)
==================================================================================
Ginette is a 2-D variably saturated groundwater flow model with integrated 1-D surface flow for the stream. It rigorously simulates water, energy, and solute transport fluxes in porous media using physically-based equations. Ginette was initially developed at METIS (Sorbonne University) to deal with interactions between streams and aquifers as they fluctuate from a connected to a disconnected status. Numerical simulations of experimental laboratory results reproducing such conditions provided the opportunity to test the coupled 1D surface water–2D variably saturated groundwater code (Rivière et al., 2014). Then, Ginette was extended to incorporate coupled heat transfer and water flow in saturated porous media.  The code was compared to experimental data acquired on a complex laboratory system to provide validation on the physical processes and mathematical formulations, in particular for the representation of density change between frozen and liquid water (Rivière et al., 2019). The coupling of fluid flow and heat transfer, accounting for freezing and thawing processes, is implemented in the code for a fully saturated medium and unsaturated zone.

The code is now jointly developed by MINES Paris (PSL University), METIS (Sorbonne University), and Cerege (Marseille University).

 Real-world cryo-hydrogeological paleo-applications, which have been presented in conferences (e.g. Jost, 2011; Jost et al., 2014), were also proposed using Ginette, requiring some additional adaptation to the specific needs of basin-scale calculations. 
The code is applied to estimate the hydraulic and thermal properties of the hyporheic zone (stream-aquifer interface) from the pressure differential and the temperature profle ( Cucchi et al. 2018 and 2021).

Ginette strives to provide the user with complete simulation control. The command files can be created using either a Python or a R script.

Different applications are available in the application directory. 
1) The test cases of the interfrost group benchmark.
 The Interfrost group proposes a benchmark exercise dealing with "subsurface thermal hydrologic processes" as presented by (Painter et al., 2012) or within the field of cryohydrogeology (McKenzie et al. 2020). In the first phase of the project, we first limit our efforts to the more simple set of equations involving Darcy flow (fully saturated porous medium) coupled with heat transfer with advection and phase change. Extensions of the benchmark to Richard equations or including the air phase are considered for later phases of the project. The benchmark consists of some test cases inspired by existing literature (e.g., Mc Kenzie et al., 2007) as well as new ones. Some experimental cases in the cold room will complement the validation approach. In view of a second phase, the benchmark project is open as well to new or alternative cases reflecting a numerical or process oriented interest or answering a more general concern among the cold region community. 

2) Different test cases in the unsaturated zone are available.

3) The chain of estimation of the hydraulic and thermal hyporheic zone properties from LOMOS-mini sensor data ( Cucchi et al. 2018 (⟨10.1016/j.jhydrol.2017.10.074⟩. ⟨hal-01656455⟩). 



## Authors:
- Riviere, Agnes, agnes.riviere@mines_paristech.fr
- Goncalves, Julio, goncalves@cerege.fr 
- Jost, Anne, anne.jost@upmc.fr


## Contributions by:
- Texier, Jérôme, texier@cerege.fr 


## Userguide
Userguide available in [User_guide.md](User_guide.md)


## References:
- Rivière, A., Gonçalvès, J., Jost, A., & Font, M. (2014). Experimental and numerical assessment of transient stream-aquifer exchange during disconnection. Journal of Hydrology, 517, 574–583. https://doi.org/10.1016/j.jhydrol.2014.05.040
- Rivière, A., Jost, A., Gonçalvès, J. & Font, M. (2018). Pore water pressure evolution below a freezing front under saturated conditions: Large-scale laboratory experiment and numerical investigation. Cold Regions Science and Technology, 158, 76-94. https://doi.org/10.1016/j.coldregions.2018.11.005
- Rivière, A., Gonçalvès, J., Jost, A., (2020). agnes-riviere/ginette: Ginette-2020-09 (Version 2020-09). Zenodo. http://doi.org/10.5281/zenodo.4058821
- Dangeard, M., Rivière, A., Bodet, L., Schneider, S., Guérin, R., Jugnot, D. et Maineult, A. River Corridor Model Constrained by Time‐Lapse Seismic Acquisition. Water Resources Research, American Geophysical Union, 2021, 57 (10), https://doi.org/10.1029/2020WR028911
- Cucchi, K., Flipo, N., Rivière, A. and Rubin, Y. Estimating Hydrothermal Properties and High-Frequency Fluxes From Geophysical Measurements in the Hyporheic Zone. Frontiers in Water, Frontiers, 2021, 3, https://doi.org/10.3389/frwa.2021.700274
- Cucchi, K.,  Rivière, A., Baudin, A., Berrhouma, A., Durand, V., Rejiba, F. et Robin, Y. LOMOS-mini: A coupled system quantifying transient water and heat exchanges in streambeds. Journal of Hydrology, Elsevier, 2018, 561, https:/10.1016/j.jhydrol.2017.10.074


## Message:
"If you use this software, please cite it as below."
Rivière, A., Gonçalvès, J., Jost, A., Ginette,   [![DOI](https://zenodo.org/badge/242535776.svg)](https://zenodo.org/badge/latestdoi/242535776)


authors:
  - Rivière Agnès
    orcid: https://orcid.org/0000-0002-6002-3189
  - Gonçalvès Julio
    orcid: https://orcid.org/0000-0003-0047-4233
  - Jost Anne
    orcid: https://orcid.org/0000-0002-0925-3376
    
    
title: agnes-riviere/ginette: Ginette-2020-09
version: 2020-09
date-released: 2017-12-18

