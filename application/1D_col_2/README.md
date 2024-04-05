# Ginette application



# Utisation
compile ginette in the current folder
use compile_ginette.py to modify the setup and the parameters of the model


## Setup

### time step in s
dt=900
### duration of the simulation in days
nb_day=100

### state
##### 0 steady state
##### 1 transient state (dynamic state)
state=0

##### Geometry in meter
z_top=40
dz=0.04

#### apply unsaturated flow and thermal 
##### unsat =1 apply
##### unsat=0 cancel unsaturated zone
unsat=1

##### Number of zone 

nb_zone=1
##### user-defined  parameters zone 1 
##### intrinsic permeability [m2]  k=K*mu/(rho*g)
###### K hydraulic conductivity [m.s-1]
###### mu viscosity [Pa.s]
###### rho density [kg.m-3]
######g gravity  9.81 [m2.s-1]
val_k=3.33333333333333e-15
##### porosity
val_n=0.38 # \Phi
##### solid grain density rho_s=val_r  [kg.m-3]
val_r=2578
##### Van Genuchten parameters
val_a=2.00000 #m-1 alpha_vg
val_nVG= 1.23  # n_vg
val_swres=0.26 # S_wr



### Boundary conditions hydraulic head h=P/rho g + Z (in meter)
top=15 
bot=15

# Initial conditions







# License

Authors:
- ... Riviere , Agnes, agnes.riviere@mines_paristech.fr
- ... Goncalves , Julio, goncalves@cerege.fr 
- ... Jost , Anne, anne.jost@upmc.fr

_Created in 2008_. Copyright ©
- Sorbonne Universite
- Universite Aix Marseille
- Mines ParisTech

[describe functionalities and technical features of your software]

_**Update coming up later**_

-------------------------------------------------------------------------
This software is governed by the CeCILL license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
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

Le fait que vous puissiez accéder à cet en-tête signifie que vous avez 
pris connaissance de la licence CeCILL, et que vous en avez accepté les
termes.

