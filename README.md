![Ginette-2](Ginette-2.png)
==================================================================================
Ginette is a 2-D variably saturated groundwater flow model with integrated 1-D surface flow for the stream. It rigorously simulates water, energy, and solute transport fluxes in porous media using physically-based equations. The coupling of fluid flow and heat transfer, accounting for freezing and thawing processes, is implemented in the code for a fully saturated medium. Ginette was initially developed at METIS (Sorbonne University) to deal with interactions between streams and aquifers as they fluctuate from a connected to a disconnected status. Numerical simulations of experimental laboratory results reproducing such conditions provided the opportunity to test the coupled 1D surface water–2D variably saturated groundwater code (Rivière et al., 2014). Then, Ginette was extended to incorporate coupled heat transfer and water flow in saturated porous media.  The code was compared to experimental data acquired on a complex laboratory system to provide validation on the physical processes and mathematical formulations, in particular for the representation of density change between frozen and liquid water (Rivière et al., 2019).

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
  - family-names: Rivière
    given-names: Agnès
    orcid: https://orcid.org/0000-0002-6002-3189
  - family-names: Gonçalvès
    given-names: Julio
    orcid: https://orcid.org/0000-0003-0047-4233
  - family-names: Jost
    given-names: Anne
    orcid: https://orcid.org/0000-0002-0925-3376
    
    
title: agnes-riviere/ginette: Ginette-2020-09
version: 2020-09
date-released: 2017-12-18

## License
_Created in 2008_. Copyright ©
- Sorbonne Universite
- Universite Aix Marseille
- Mines ParisTech


-------------------------------------------------------------------------
This software is governed by the CeCILL license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILLerging pull requests, you can allow any combin
license as circulated by CEA, CNRS and INRIA at the following URL
http://www.cecill.info. 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.

-------------------------------------------------------------------------
Ce logiciel est régi par la licence CeCILL soumise au droit français et
respectant les principes de diffusion des logiciels libres. Vous pouvez
utiliser, modifier et/ou redistribuer ce programme sous les conditions
de la licence CeCILL telle que diffusée par le CEA, le CNRS et l'INRIA 
sur le site http://www.cecill.info.

En contrepartie de l'accessibilité au code source et des droits de copie,
de modification et de redistribution accordés par cette licence, il n'est
offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,
seule une responsabilité restreinte pèse sur l'auteur du programme,  le
titulaire des droits patrimoniaux et les concédants successifs.

A cet égard  l'attention de l'utilisateur est attirée sur les risques
associés au chargement,  à l'utilisation,  à la modification et/ou au
développement et à la reproduction du logiciel par l'utilisateur étant 
donné sa spécificité de logiciel libre, qui peut le rendre complexe à 
manipuler et qui le réserve donc à des développeurs et des professionnels
avertis possédant  des  connaissances  informatiques approfondies.  Les
utilisateurs sont donc invités à charger  et  tester  l'adéquation  du
logiciel à leurs besoins dans des conditions permettant d'assurer la
sécurité de leurs systèmes et ou de leurs données et, plus généralement, 
à l'utiliser et l'exploiter dans les mêmes conditions de sécurité. 

Il n'est pas permis de vendre ce logiciel ni de l'utiliser sous une prestation sans l'accord des développeurs.

Le fait que vous puissiez accéder à cet en-tête signifie que vous avez 
pris connaissance de la licence CeCILL, et que vous en avez accepté les
termes.


