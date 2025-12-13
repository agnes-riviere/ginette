!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC GINETTE CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! AUCUNE MODIFICATION DU CODE N'IMPLIQUE LES DEVELLOPEURS AGNES RIVIERE JULIO GONCALVES ANNE JOST      C
! contacts agnes.riviere@mines_paristech.fr                                    C
!       goncalves@cerege.fr                                           C
!       anne.jost@upmc.fr                                         C
! Toute utilisation ou copie de ce code implique la citation des auteurs                    C
!                                                         C
! Toute vente commerciale en temps que service ou logiciel doit obtenir l'accord des auteurs         C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
program pression_ecoulement_transport_thermique
   implicit double precision(a - h, o - x, z), integer(I - N)
   implicit CHARACTER*5(y)
   double precision dx, dz, dt, al, az, reptop, repbot, rho1
   double precision g, amu, akrx, akx, akrz, akz, omp, sss, ans, asp,swo
   double precision swres, allg, alt, qre, hbot, xberg, unitsortie
   double precision rug, pent, qriva, hriv, aklit, aklitv, unitsim, sum
   double precision tempriv, elit, akdrain, edrain, crconvp, crconvc
   double precision pas
   integer :: compteur_div
   integer nc, nr, n, iec, irp, ith, nitt, ixy, ii
   integer imaille, itopo, ibuilt, icolone
   integer nci, nri, nclog, ilog, irho
   integer ipermh, ipermv, iom, iss, ia2, ivg
   integer iteration, itsortie, icalvit
   integer mcol, nmaille1, nmaille2, nmaille3, nmaille4
   integer nmaille5, nmaille6, nmaille7, nmaille8, nmaille9
   integer nmaille10
   integer iparo, iriv, iqriv, itr
   CHARACTER(5) ::  ytypakrice, ytypsice
   CHARACTER(3) :: ytest, yunconfined, ysolv
   integer, dimension(:), allocatable :: icol_ind, irow_ptr, ibas
   double precision, dimension(:), allocatable::  val, b, pr, x, z
   integer, dimension(:, :), allocatable :: icl, iclc, iclt, ivois
   integer, dimension(:), allocatable :: inum, inum2, izone, jzone, icol
   double precision, dimension(:), allocatable:: dswdp, om, rhold
   double precision, dimension(:), allocatable:: akr, ak, akrv, akv, rho
   double precision, dimension(:), allocatable:: bm, am, sw, pro
   double precision, dimension(:, :), allocatable:: zs, zso, zsoo, zl, zlo
   double precision, dimension(:, :), allocatable:: zloo
   double precision, dimension(:), allocatable:: temp, qtherm, asun
   double precision, dimension(:), allocatable:: ansun, conco, tempo
   double precision, dimension(:, :), allocatable:: valclt, valclc, valcl
   double precision, dimension(:, :), allocatable:: valclto
   double precision, dimension(:), allocatable:: rhos, alanda, cps
!      dimension icpiso(:)
   double precision, dimension(:), allocatable::dl, def, defo, dsidtempo
   double precision, dimension(:), allocatable::alph, dsidtempoo
   double precision, dimension(:), allocatable::vxm, vxp, vzp, vzm
   double precision, dimension(:), allocatable::prk, tempk, conck, conc
   double precision, dimension(:), allocatable::chg, alandas, ss, topo, bot
   integer, dimension(:), allocatable :: irow
   double precision, dimension(:), allocatable::sice, rhoi, siceo, siceoo
!     allocatable :: zbot(:),zaqui(:)
   double precision, dimension(:), allocatable::tempoo, swp, sicep, dswpdp
   double precision, dimension(:), allocatable::dsipdp, dsidp
   double precision, dimension(:), allocatable::dsidtemp, dsipdtemp
   double precision, dimension(:), allocatable::akzone, omzone, anszone
   double precision, dimension(:), allocatable::aspzone, dswdt
   double precision, dimension(:), allocatable::swreszone, alandazone
   double precision, dimension(:), allocatable::cpmzone, rhomzone
   double precision, dimension(:), allocatable::swresz
   double precision, dimension(:), allocatable::tempsol, qpluie
   double precision, dimension(:), allocatable::chgRD, chgRG, tempRD
   double precision, dimension(:), allocatable::tempRG
   double precision, dimension(:), allocatable::timeG, timeD, timeDTS
   double precision, dimension(:), allocatable::cRivG, cRivD, tempriver
   integer, dimension(:), allocatable :: id_river, id_rivert, id_ZH
   integer, dimension(:), allocatable :: id_RD, id_RG
   double precision, dimension(:), allocatable::chgbot, chgsurf, chgriver
   double precision, dimension(:), allocatable:: qbot, qsurf, xDTS
   double precision, dimension(:), allocatable::tempbot, tempsurf
   double precision, dimension(:, :), allocatable::  tempDTS, slopeRH
   double precision, dimension(:, :), allocatable::xpool, xriffle
   double precision, dimension(:), allocatable::qout_w, qout_wR
   double precision, dimension(:), allocatable::qin_w, qin_wR
   double precision, dimension(:), allocatable::qout_h, qout_hR
   double precision, dimension(:), allocatable::qin_h, qin_hR
   double precision, dimension(:), allocatable::qadvout_h, qadvout_hR
   double precision, dimension(:), allocatable::qadvin_h, qadvin_hR
   double precision, dimension(:), allocatable::qcondout_h, qcondout_hR
   double precision, dimension(:), allocatable::qcondin_h, qcondin_hR
   double precision, dimension(:), allocatable::qad, qcondu
   integer           :: nm
   LOGICAL           :: An_Error

!     integer(4) :: deg_max_gc
!     integer(4) :: sw_int(12)
!     real(8) :: sw_reel(10)

!   include 'dump.h'
!     common /entsor/ lec,imp
!     lec = 5
!     imp = 6

!     call GC_set_dump ('solgc',0)

!     nsys = 1
!     initb0 = 0
!     nivprop = 0
!     lmac5prop = 0
!     niter = 10
!     itrait = 1
!     istgrc = 1
!     i_cal_pc = 1
!     ieco = 0
!     epsgc = 1.e-9
!     deg_max_gc = 10
!     iappli = 1
!
!     sw_int(1) = nsys
!     sw_int(2) = initb0
!     sw_int(3) = nivprop
!     sw_int(4) = lmac5prop
!     sw_int(5) = niter
!     sw_int(6) = itrait
!     sw_int(7) = istgrc
!     sw_int(8) = i_cal_pc
!     sw_int(9) = ieco
!     sw_int(10) = deg_max_gc
!     sw_int(11) = message
!     sw_int(12) = iappli
!     sw_reel(1) = epsgc
!     sw_reel(2:10) = 0.

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!          MODIFICATIONS             C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC....Modification du format de swres f6.4 swres
!CC....Modification format porosite f6.4
!CC....Modification du format parametre chaleur latente format de chlat L (m2/s2) en D8.0
!CC....19-11-2014 TEST avec densite de la glace dans l'equation de la chaleur

!CC....13-07-2010 modif matp BILAN D EAU
!CC....13-07-2010 modif matt matc vitesse
!CC....11-2010 Loi AP
!CC....23-03-2011 variation de dz
!CC....modif vitesse, matp, matc, matt
!CC....modif vitesse x ajout de rho g * dz/dx
!CC....ajout condition cdt imposee vitesse
!CC....ajout permeabilite vertical
!CC....ajout coefficient d emmagasinement quand pas de zone non sature
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!      COMMENTAIRE+TODOLIST              C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC....si utilisation variation densite moyenne a changer dans matp!!!!
!cc au depart tsolidus et liquidus = GEl

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!      SIGNIFICATION                C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC....al longeur
!CC....az profondeur
!CC....nc nb de icolones
!CC....nr nb de ligne
!CC....qre debit pluie
!CC....rho1 densite eau
!CC....ia2 variation temporelle de la densite prise en compte dans matp
!CC....ans=n
!CC....as=alpha
!CC....ivois(ik,1)= voisin droite
!CC....ivois(ik,2)= voisin gauche
!CC....ivois(ik,3)= voisin haut
!CC....ivois(ik,4)= voisin ibas
!CC....nm nombre de mailles reelles
!CC....inum ancien numeros de maille avant recalcul
!CC....inum2 numero maille apres recalcul
!CC....am=largeur de la maille
! ECOULEMENT icl condition valcl valeur
!CC....ICL=-1 Flux impose sur une face
!CC....ICL=-2 potentiel impose sur une face
!CC....ICL=1 Mailles 'normale'
!CC.... ICL(ik,1)=-1, -2 ou 1 ik=1,...4: face droite,gauche,haute et BASSE
! TRANSPORT iclc condition valclc valeur
!CC....ICLT=-1 Flux impose sur une face
!CC....ICLT=-2 potentiel impose sur une face
!CC....ICLT=1 Mailles 'normale'
!CC....ICLT(ik,1)=-1, -2 ou 1 ik=1,...4: face droite,gauche,haute et BASSE
!CCCC  N.B.: LES FLUX IMPOSES (ICLC(i,j)=-1) SONT TRAITES VIA LE TERME ADVECTIF
! THERMIQUE iclt condition valclt valeur
!CC....ICLT=-1 Flux impose sur une face
!CC....ICLT=-2 potentiel impose sur une face
!CC....ICLT=1 Mailles 'normale'
!CC....ICLT(ik,1)=-1, -2 ou 1 ik=1,...4: face droite,gauche,haute et BASSE
!CCCCC  N.B.: LES FLUX IMPOSES (ICLC(i,j)=-1) SONT TRAITES VIA LE TERME ADVECTIF
!CC....pr(i)=pression
!CC....akr(i)=permeabilite relative
!CC....ak(i)=permeabilite intrinseque
!CC....sw(i)=saturation
!CC....om(i)=porosite
!CC....dswdp(i)=0.
!CC....rho(i)=densite
!CC....conc(i)=concentration
!CC....temp(i)=temperature
!CC....alanda(i)=conductivite thermauque
!CC....alandae conductivite thermique eau
!CC....alandam conductivite thermique milieu
!CC....igel=2 degel
!CC....igel=1 gel
!CC....igel=0 pas de variation de temperature
!CC....swo : old staturation time n-1
!CC....compteur_div option debug pour compter les iteration de temps pendant la convergence

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!         FICHIERS ENTREES                C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
   open (unit=11, file='E_coordonnee.dat', iostat=iocell)
   open (unit=12, file='E_def_maille.dat')
   open (unit=13, file='E_voisins.dat')
   open (unit=40, file='E_colonne.dat')
!   open(unit=21,file='E_row.dat')
!

!   open(unit=41,file='E_cl_drain_t.dat',iostat=io)
!   open(unit=42,file='E_cl_riv_pluie1_t.dat',iostat=iot)
!   open(unit=43,file='E_cl_riv_pluie2_t.dat',iostat=io2)
!   open(unit=44,file='E_cl_riv_pluie3_t.dat',iostat=io3)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!        lecture parameters           C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

   call lecture_parametre( &
      ilog, ipermh, ipermv, &
      iec, irp, ith, dt, nitt, unitsim, al, &
      az, igravit, imaille, itopo, reptop, repbot, &
      icolone, nmi, nci, nri, dx, dz, nclog, rho1, irho, g, amu, &
      akrx, akx, akrz, akz, iom, omp, iss, sss, ia2, yunconfined, ivg, &
      ans, asp, swres, itr, allg, alt, iriv, iqriv, qre, hbot, xberg, &
      rug, pent, qriva, hriv, aklit, aklitv, tempriv, elit, akdrain, &
      edrain, crconvp, crconvc, &
      iteration, itsortie, unitsortie, icalvit, mcol, nmaille1, nmaille2, &
      nmaille3, nmaille4, nmaille5, &
      nmaille6, nmaille7, nmaille8, nmaille9, nmaille10, iparo, ibuilt, &
      ysolv)

   ixy = igravit

   call lecture_parametre_thermique(irpth, crconvt, &
                                    ithec, idecouplage, ysupdp, &
                                    bll, blt, alandae, cpe, ilanda, &
                                    alandami, irhomi, rhosi, icpm, cpm, ymoycondtherm, &
                                    icycle, igel, igelzns, iomdegel, &
                                    ytypakrice, ytypsice, &
                                    alandai, rhoii, cpice, &
                                    chlat, rhog, alandag, cpg, dk, tsg, tlg, tsd, tld, cimp, swressi, &
                                    iaquitard, aktardx, aktardz, omptard, alandatard, sstard, &
                                    nrowtard, nmailleaqui, iexutoire, xexutoire, tempexutoire, &
                                    iinfil, itopomch, ibigridice, ytest, omega, cimpt)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!        LECTURE CDT INITALES                 C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
   call lecture_cdt_ini(chgi, conci, tempini, ichi, ichi2, iconci, &
                        itempi)

   if (iconci == 1) open (unit=22, file='E_conc_initiale.dat')
   if (itempi == 1) open (unit=23, file='E_temperature_initiale.dat')
!CC...charge imposee
   if (ichi == 1) open (unit=24, file='E_charge_initiale.dat')
!CC...VARIATION DE LA CHARGE INITIALE sur tout le model
   if (ichi2 == 1) open (unit=242, file='E_pression_initiale.dat')

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!        LECTURE CDT LIMITES                  C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
   call lecture_cdt_lt(icl_gauche, valcl_gauche, icl_droite, &
                       valcl_droite, icl_haut, valcl_haut, icl_bas, valcl_bas, iclc_gauche, &
                       valclc_gauche, iclc_droite, valclc_droite, iclc_haut, valclc_haut, &
                       iclc_bas, valclc_bas, &
                       iclt_gauche, valclt_gauche, iclt_droite, valclt_droite, iclt_haut, &
                       valclt_haut, iclt_bas, valclt_bas, iclect, icltherm, iclectchgt &
                       , iclchgt, icldrain, iclriviere, itlecture)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!      LECTURE FICHIER Variation spatiale     C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

   if (ipermh == 1) &
      open (unit=14, file='E_permeabilite_horizontale.dat')
   if (ipermv == 1) &
      open (unit=140, file='E_permeabilite_verticale.dat')
   if (iom == 1) &
      open (unit=15, file='E_porosite.dat')
   if (irho == 1) &
      open (unit=25, file='E_densite_fluide.dat')
   if (iss == 1) &
      open (unit=16, file='E_coefficient_emmagasinement.dat')
   if (irhomi == 1) &
      open (unit=26, file='E_densite_milieu.dat')
   if (ilanda == 1) &
      open (unit=17, file='E_conductivite_thermique.dat')
   if (icpm == 1) &
      open (unit=27, file='E_capacit_calorifique.dat')
!CC....Fichiers variations conditions aux limites
   if (iclect == 1) &
      open (unit=37, file='E_cl_ecoulement.dat')
   if (icltherm == 1) &
      open (unit=38, file='E_cl_thermique.dat')

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!      GEOMETRIE MODEL                    C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
   allocate(ibas(nci))
   allocate(topo(nci))
   allocate(bot(nci))

   select case (itopo)
   case (1)
      open (unit=10, file='E_geom.dat', form='formatted', &
            status='old', iostat=igeom)
      call count_file(10, igeom, nc)
      do i = 1, nc
         read (10, *) topo(i), bot(i)
      end do
   case (0)
      do i = 1, nci
         topo(i) = (reptop)
         bot(i) = (repbot)
      end do
!c endif itop
   end select

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!      LECTURE FICHIER IBUILT                 C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

   if (ytest == "VIT") then
!     open(unit=453,file='E_vitesse.dat')
!   open(unit=452,file='E_pres_t.dat',iostat=iop2)
      open (unit=45, file='E_charge_t.dat', iostat=iop)
      open (unit=68, file='E_temp_t.dat', iostat=ios)
      open (unit=23, file='E_temperature_initiale.dat')
   end if

   if (ytest == "TEX") then
      icolone = 1
      imaille = 0
      itopo = 1
      call open_file('E_zone.dat', An_Error, 32)
      call open_file('E_zone_parameter.dat', An_Error, 321)
      open (unit=222, file='E_BordRD.dat', iostat=iidRD)
      open (unit=223, file='E_BordRG.dat', iostat=iidRG)
      open (unit=224, file='E_Id_river.dat', iostat=iidrv)
      open (unit=2244, file='E_Id_river_max.dat', iostat=iidrvtest)
      open (unit=226, file='E_PluieR.dat', iostat=i226)
      open (unit=227, file='E_chargeT_RD.dat', iostat=iCRD)
      open (unit=228, file='E_chargeT_RG.dat')
      open (unit=232, file='E_chargeT_Riv.dat')

   end if

   if (ytest == "AVA") then
      icolone = 1
      imaille = 0
      itopo = 1
      call open_file('E_zone.dat', An_Error, 32)
      call open_file('E_zone_parameter.dat', An_Error, 321)
      open (unit=222, file='E_BordRD.dat', iostat=iidRD)
      open (unit=223, file='E_BordRG.dat', iostat=iidRG)
      open (unit=2242, file='E_bottomZH.dat', iostat=iidZH)
      open (unit=224, file='E_Id_river.dat', iostat=iidrv)
      open (unit=2244, file='E_Id_river_max.dat', iostat=iidrvtest)
      open (unit=225, file='E_TempSol.dat', iostat=iTsol)
      open (unit=226, file='E_PluieR.dat', iostat=i226)
      open (unit=227, file='E_chargeT_RD.dat', iostat=iCRD)
      open (unit=228, file='E_chargeT_RG.dat')
!     open(unit=229,file='E_tempT_RD.dat')
!     open(unit=230,file='E_tempT_RG.dat')
!     open(unit=231,file='E_tempT_Riv.dat')
      open (unit=232, file='E_chargeT_Riv.dat')

   end if

   if (ytest == "R2D") then
      icolone = 1
      imaille = 0
      itopo = 1

      call open_file('E_zone.dat', An_Error, 32)
      call open_file('E_zone_parameter.dat', An_Error, 321)
      open (unit=222, file='E_BordRD.dat', iostat=iidRD)
      open (unit=223, file='E_BordRG.dat', iostat=iidRG)
      open (unit=2242, file='E_bottomZH.dat', iostat=iidZH)
      open (unit=224, file='E_Id_river.dat', iostat=iidrv)
      open (unit=2244, file='E_Id_river_max.dat', iostat=iidrvtest)
 !     open (unit=225, file='E_TempSol.dat', iostat=iTsol)
 !     open (unit=226, file='E_PluieR.dat', iostat=i226)
      open (unit=227, file='E_chargeT_RD.dat', iostat=iCRD)
      open (unit=228, file='E_chargeT_RG.dat')
     open(unit=229,file='E_tempT_RD.dat')
     open(unit=230,file='E_tempT_RG.dat')
     open(unit=231,file='E_tempT_Riv.dat')
      open (unit=232, file='E_chargeT_Riv.dat')

   end if

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!    FICHIERS DTS                             C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

   if (ytest == "DTS") then
      icolone = 1
      imaille = 0
      ibuilt = 0
      itopo = 1

      call DTS_open(icolone, imaille, itopo, ligne, ligne1, &
                    ligne2, ligne3, ligne4, ligne6)
      allocate(cRivD(ligne))
      allocate(cRivG(ligne))

   end if

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!    FICHIERS Construction manuelle maillage  C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
   if (imaille == 0 .or. ibuilt == 1) then
      open (unit=11, file='E_coordonnee.dat')
      open (unit=12, file='E_def_maille.dat')
      open (unit=13, file='E_voisins.dat')
      open (unit=40, file='E_colonne.dat')
      open (unit=21, file='E_row.dat')
   end if

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!      LECTURE FICHIER ZHR                    C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
   if (ytest == "ZHR" .or. ytest == "ZHZ") then
      open (unit=45, file='E_charge_t.dat', iostat=iop)
      open (unit=68, file='E_temp_t.dat', iostat=ios)
      if (ios /= 0) then
         stop 'File E_temp_t.dat does not exist'
      end if
   end if
   if (ytest == "ZHZ") then
      call open_file('E_zone.dat', An_Error, 32)
      call open_file('E_zone_parameter.dat', An_Error, 321)
   end if

   if (ytest == "1DS") then
      open (unit=45, file='E_ec_bc_t.dat', iostat=iop)
      open (unit=68, file='E_temp_t.dat', iostat=ios)
      if (ios /= 0) then
         stop 'File E_temp_t.dat does not exist'
      end if
      call open_file('E_zone.dat', An_Error, 32)
      call open_file('E_zone_parameter.dat', An_Error, 321)
   end if

   if (ytest == "1DJ" .or. ytest == "ZND") then
      open (unit=45, file='E_ec_bc_t.dat', iostat=iop)
   end if

   if (ytest == "ZNS" .or. ytest == "1DJ" .or. ytest == "ZND") then
      call open_file('E_zone.dat', An_Error, 32)
      call open_file('E_zone_parameter.dat', An_Error, 321)
   end if
   if (ytest == "ZNS" .and. ith == 1) then
      open (unit=6868, file='E_temp_t.dat', iostat=ios)
   end if
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!      LECTURE FICHIER ZNS 1D ou warrick          C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
   if (ytest == "ZNS" .or. ytest == "WAR") then
      if (iclchgt == 1) then
         select case (icl_bas)
         case (-1)
            open (unit=45, file='E_debit_bas_t.dat', iostat=iop)
         case (-2)
            open (unit=45, file='E_charge_bas_t.dat', iostat=iop)
         end select
         select case (icl_haut)
         case (-1)
            open (unit=68, file='E_debit_haut_t.dat', iostat=ios)
         case (-2)
            open (unit=68, file='E_charge_haut_t.dat', iostat=ios)
         end select
      end if
   end if
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!      FICHIERS ERREUR VERIF                  C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!       open(46,file='S_erreur_t.dat')
!       open(20,file='S_pb_resolution.dat')
!       open(55,file='S_saturation.dat')
!       open(6001,file='S_permeabilite_longitudinale_verticale.dat')
!       open(75,file='S_charge_initiale.dat')
!       open(1000,file='S_pb_cfl.dat')
!       open(97,FILE='S_verif_drain_t.dat')

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!    INITIALISATION DU PAS DE TEMPS                   C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

   inano = 0
   inan = 0
   dta = dt

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!   LECTURE FILE CDT VS. TEMPS                C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC....INITIALISATION COMPTEUR LIGNES
   select case (ytest)
   case("ZHR", "ZHZ")
      ligne4 = 0

!CC....lecture des données
      call count_file(68, ios, ligne4)


      if (allocated(chgbot)) deallocate(chgbot)

      allocate(chgbot(ligne4))
      allocate(chgsurf(ligne4))
      allocate(tempbot(ligne4))
      allocate(tempsurf(ligne4))
      rewind (68)
      do j = 1, ligne4
         read (45, *, iostat=iop) chgsurf(j), chgbot(j)
         read (68, *, iostat=ios) tempsurf(j), tempbot(j)
      end do

   CASE('1DS')

      ligne4 = 0
!CC....lecture des données temperature time
      call count_file(68, ios, ligne4)


      if (allocated(chgbot)) deallocate(chgbot)

      allocate(chgbot(ligne4))
      allocate(chgsurf(ligne4))
      allocate(tempsurf(ligne4))
      allocate(tempbot(ligne4))
      rewind (68)
      do j = 1, ligne4
         read (45, *, iostat=iop) chgsurf(j), chgbot(j)
         read (68, *, iostat=ios) tempsurf(j), tempbot(j)
      end do

   case("1DJ", "ZND")
      ligne4 = 0
!CC....lecture des données temperature time
      do while (ios == 0)
         read (45, *, iostat=ios)
         if (ios == 0) then
            ligne4 = ligne4 + 1
         end if
      end do
      if (allocated(chgbot)) deallocate(chgbot)

      allocate(chgbot(ligne4))
      allocate(chgsurf(ligne4))

      rewind (45)
      do j = 1, ligne4
         read (45, *, iostat=iop) chgsurf(j), chgbot(j)
      end do


!CC....ZNS 1D ou Warrick
   case ("ZNS", "WAR")

!CCCC... ajout nicolas Radic
      ligne_temp_t = 0
      if (ith ==1) then
         call count_file(6868, ios, ligne_temp_t)

         allocate(tempbot(ligne_temp_t))
         allocate(tempsurf(ligne_temp_t))
         rewind (6868)
         do j = 1, ligne_temp_t
            read (6868, *, iostat=ios) tempsurf(j), tempbot(j)
         end do
      end if
!CCCC... fin ajout nicolas Radic
      if (iclchgt == 1) then
         ligne4 = 0
!CC....lecture des données
         do while (ios == 0)
            read (68, *, iostat=ios)
            if (ios == 0) then
               ligne4 = ligne4 + 1
            end if
         end do
         rewind (68)

         select case (icl_bas)
         case (-1)
            allocate(qbot(ligne4))
            do j = 1, ligne4
               read (45, *, iostat=iop) qbot(j)
            end do
         case (-2)
            allocate(chgbot(ligne4))
            do j = 1, ligne4
               read (45, *, iostat=iop) chgbot(j)
            end do

         end select
         select case (icl_haut)
         case (-1)
            allocate(qsurf(ligne4))
            do j = 1, ligne4
               read (68, *, iostat=iop) qsurf(j)
            end do
         case (-2)
            allocate(chgsurf(ligne4))
            do j = 1, ligne4
               read (68, *, iostat=iop) chgsurf(j)

            end do
         end select
      end if


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! MARINE
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!   LECTURE FILE CDT VS. TEMPS        C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC....INITIALISATION COMPTEUR LIGNES

   case ("AVA", "TEX")


      ligne = 0
      ligne1 = 0
      ligne2 = 0

!CC....lecture des données
!ccc...tempsol
      call count_file(225, iTsol, ligneT)  !E_TempSol.dat
      allocate(tempsol(ligneT))


      rewind (225)
      do j = 1, ligneT
         read (225, *, iostat=iTsol) tempsol(j)
      end do
      ligne=ligneT
!ccc...qpluie
      call count_file(226, i226, ligne_pluie)  !E_PluieR.dat

      ligne1=ligne_pluie
      allocate(qpluie(ligne_pluie))
      rewind (226)
      do j = 1, ligne_pluie
         read (226, *, iostat=i226) qpluie(j)
      end do
!ccc...boundaries wall
      call count_file(227, iCRD, ligneCL)  !E_chargeT_RD.dat'
      ligne2=ligneCL
      
      allocate(chgbot(ligneCL))
      allocate(chgsurf(ligneCL))
      allocate(tempbot(ligneCL))
      allocate(tempsurf(ligneCL))
      allocate(chgRD(ligneCL))
      allocate(chgRG(ligneCL))
      allocate(tempRD(ligneCL))
      allocate(tempRG(ligneCL))
      allocate(tempriver(ligneCL))
      allocate(chgriver(ligneCL))
      allocate(qsurf(ligneCL))
      allocate(qbot(ligneCL))
      rewind (227)
      do j = 1, ligneCL
         read (228, *) chgRG(j)
         read (227, *, iostat=iCRD) chgRD(j)
!      read(229,*)tempRD(j)
!      read(230,*)tempRG(j)
!      read(231,*)tempriver(j)
         read (232, *) chgriver(j)
      end do

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!   LECTURE FILE IBUILT   IDENTIFIANT      C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!CC....INITIALISATION COMPTEUR LIGNES
      ligne3 = 0
      ligne4 = 0
      ligne5 = 0
      ligne6 = 0
      ligne7 = 0

!CC....lecture des données
!ccc...bord droit
      call count_file(222, iidRD, ligne3)  ! 222 = E_BordRD.dat

      allocate(id_RD(ligne3))
      rewind (222)
      do j = 1, ligne3
         read (222, *, iostat=iidRD) id_RD(j)
      end do

!ccc...bord gauche
      call count_file(223, iidRG, ligne4) ! 223 = E_BordRG.dat
      allocate(id_RG(ligne4))
      rewind (223)
      do j = 1, ligne4
         read (223, *, iostat=iidRG) id_RG(j)
      end do

!ccc...riviere
      call count_file(224, iidrv, ligne5) ! 224 = E_Id_river.dat
      rewind (224)

      allocate(id_river(ligne5))

      do j = 1, ligne5
         read (224, *, iostat=iidrv) id_river(j)
      end do

!ccc...riviere a tester
      call count_file(2244, iidrvtest, ligne6) ! 2244 = E_Id_river_max.dat
      rewind (2244)

      allocate(id_rivert(ligne6))

      do j = 1, ligne6
         read (2244, *, iostat=iidrvtest) id_rivert(j)
      end do

!CC if ZH pour flux
      call count_file(2242, iidZH, ligne7) ! 2242 = E_bottomZH.dat
      rewind (2242)
      allocate(id_ZH(ligne7))

      do j = 1, ligne7
         read (2242, *, iostat=iidZH) id_ZH(j)
      end do


     case ("R2D")


      ligne = 0
      ligne1 = 0


!CC....lecture des données
!ccc...boundaries wall
      call count_file(227, iCRD, ligneCL)  !E_chargeT_RD.dat'
      ligne2=ligneCL
      
      allocate(chgRD(ligneCL))
      allocate(chgRG(ligneCL))
      allocate(tempRD(ligneCL))
      allocate(tempRG(ligneCL))
      allocate(tempriver(ligneCL))
      allocate(chgriver(ligneCL))
      rewind (227)
      do j = 1, ligneCL
         read (228, *) chgRG(j)
         read (227, *, iostat=iCRD) chgRD(j)
         read(229,*)tempRD(j)
         read(230,*)tempRG(j)
         read(231,*)tempriver(j)
         read (232, *) chgriver(j)
      end do

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!   LECTURE FILE IBUILT   IDENTIFIANT      C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!CC....INITIALISATION COMPTEUR LIGNES
      ligne3 = 0
      ligne4 = 0
      ligne5 = 0
      ligne6 = 0
      ligne7 = 0

!CC....lecture des données
!ccc...bord droit
      call count_file(222, iidRD, ligne3)  ! 222 = E_BordRD.dat

      allocate(id_RD(ligne3))
      rewind (222)
      do j = 1, ligne3
         read (222, *, iostat=iidRD) id_RD(j)
      end do

!ccc...bord gauche
      call count_file(223, iidRG, ligne4) ! 223 = E_BordRG.dat
      allocate(id_RG(ligne4))
      rewind (223)
      do j = 1, ligne4
         read (223, *, iostat=iidRG) id_RG(j)
      end do

!ccc...riviere
      call count_file(224, iidrv, ligne5) ! 224 = E_Id_river.dat
      rewind (224)

      allocate(id_river(ligne5))

      do j = 1, ligne5
         read (224, *, iostat=iidrv) id_river(j)
      end do

!ccc...riviere a tester
      call count_file(2244, iidrvtest, ligne6) ! 2244 = E_Id_river_max.dat
      rewind (2244)

      allocate(id_rivert(ligne6))

      do j = 1, ligne6
         read (2244, *, iostat=iidrvtest) id_rivert(j)
      end do

!CC if ZH pour flux
      call count_file(2242, iidZH, ligne7) ! 2242 = E_bottomZH.dat
      rewind (2242)
      allocate(id_ZH(ligne7))

      do j = 1, ligne7
         read (2242, *, iostat=iidZH) id_ZH(j)
      end do
   end select



!       ligne5=0
!CC....TEST POUR MAQUETTE DECONNECTION
!       if (icldrain == 1.or.iclriviere == 1) then
!       ligne=0
!       ligne1=0
!       ligne2=0
!       ligne3=0
!       ligne4=0
!       ligne5=0
!       eo1=0D+00
!       eo2=0D+00
!       eo3=0D+00
!       eo4=0D+00
!       eo5=0D+00
!       eo6=0D+00
!       eoinf=0D+00
!       seo1=0D+00
!       seo2=0D+00
!       seo3=0D+00
!       seo4=0D+00
!       seo5=0D+00
!       seo6=0D+00
!       seoinf=0D+00
!       endif

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
! Initialisation des variables d'enregistrementC
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

   dtrecord = 0
   dtrecord_mem = 0
   irecord = 0
   compteur_div = 0

! Regime Permanent
!       if(irp == 0) then
!       nitt=1.
!       endif
! Regime Permanent
!       if(irpth == 0) then
!       nitt=1.
!       endif
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!        Indice CONVERCENCE           C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC....ECOULEMENT
   amaxp = 0.D+00
   if (iec == 1) amaxp = 1D+05
   if (iec == 0) amaxp = crconvp
!CC....TRANSPORT
   amaxc = 0.D+00
   if (itr == 1) amaxc = 1D+05
   if (itr == 0) amaxc = crconvc
!CC....THERMIQUE
   amaxt = 0.D+00
   if (ith == 1) amaxt = 1D+05
   if (ith == 0) amaxt = crconvt

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!              MAILLAGE           C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC....CALCUL LIGNES COLONNES

   select case (imaille)
   CASE (0)

      call count_file(11, iocell, linecell)
      if (linecell .ne. nmi) then
         print *, 'le fichier de coordonnees et nm ne correspondent pas'
         print *, 'nb ligne E_coordonnee.dat =', linecell
         print *, 'nmi E_parametre.dat =', nmi
         stop
      end if
      ncnr = linecell
   CASE DEFAULT
      nc = nint(abs(al)/dx)
      nr = nint(abs(az)/dz)
      ncnr = nc*nr
   END SELECT

   allocate(x(ncnr))
   allocate(z(ncnr))
   allocate(inum(ncnr))
   allocate(inum2(ncnr))
   allocate(icol(ncnr))
   allocate(am(ncnr))
   allocate(bm(ncnr))
   allocate(ivois(ncnr, 4))

   select case (imaille)
   CASE (1)
!CC....MAILLAGE AUTO
!CCC....DEFINITION DU MAILLAGE BRUT
!cccc....Premiere possibilite maillage a pas d espace dx variable
!cccc....Progression logarithmique --> plus de precision au contact
!cccc....reseau de surface
!CC....INTERFROST TEST TH2
      if (ytest == "TH2") then
         dx = dble(0.1D+00/12D+00)
         dz = dble(0.1D+00/12D+00)
      end if

!cccc....Progression log du pas d espace en X ILOG+1
      select case (ilog)
      CASE (1)
         sum = 0
         nc = nclog
         do ii = 1, nc
            iat = ii
            sum = sum + log(2.*iat)
         end do

         pas = abs(al)/sum

         do i = 1, nc
            iat = nc - i + 1
            am(i) = dble(pas*log(2.*iat))
            bm(i) = dble(dz)
         end do
      CASE (0)
         do i = 1, nc*nr
            am(i) = dble(dx)
            bm(i) = dble(dz)
         end do
      END SELECT

      k = 0
      do i = 1, nr
      do j = 1, nc
         k = k + 1
         am(k) = dble(am(j))
         inum(k) = 0
         inum2(k) = 0
!cccc....initialisation tableau voisinage
         do ij = 1, 4
            ivois(k, ij) = -99
         end do
!cccc....definition des coordonnees
!cccc....modif  23/03/2011
!cccc....log en z
         select case (ilog)
         CASE (1)
            if (j == 1) then
               x(k) = (am(k)/2)
            else
               x(k) = (x(k - 1) + (am(k - 1) + am(k))/2)
            end if
         CASE (0)
            x(k) = (dx/2 + (j - 1)*dx)
         END SELECT

         if (i == 1) then
            z(k) = dble(az - bm(k)/2)
!cccc....Modif 19-09-2014- point repere haut de colonne pour modele 1D
            if (itopo == 0) then
            if (abs(reptop - repbot) + 1D-08 > az .and. &
                abs(reptop - repbot) - 1D-08 < az) then
               z(k) = dble(az - bm(k)/2.+reptop - az)
            else
               print *, 'probleme dans votre constuction de modele'
               print *, 'reptop-repbot=', reptop - repbot + 1D-08, 'et az=', az
               print *, reptop, repbot, reptop - repbot
               stop
            end if

            end if
         else

            z(k) = dble(z(k - nc) - (bm(k - nc) + bm(k))/2)
         end if
!CC....CREATION TABLEAU DE VOISINAGE
!cccc....ivois(ik,1)= voisin droite
!cccc....ivois(ik,2)= voisin gauche
!cccc....ivois(ik,3)= voisin haut
!cccc....ivois(ik,4)= voisin ibas
         if (j .ne. nc) ivois(k, 1) = k + 1
         if (j .ne. 1) ivois(k, 2) = k - 1
         if (i .ne. 1) ivois(k, 3) = k - nc
         if (i .ne. nr) ivois(k, 4) = k + nc
      end do
      end do

!CC....RENUMEROTATION en FONCTION DE TOPO ET BOTTOM

      kr = 0
      do i = 1, nr
      do j = 1, nc
         ii = (i - 1)*nc + j
         if (z(ii) <= topo(j) .and. z(ii) >= bot(j)) then
            kr = kr + 1
!CC....NOUVEAU NUMERO de MAILLE
            inum(kr) = ii
            inum2(ii) = kr
!CC...ajout nicolas radic
            select case (ivois(ii,3))
            case(-99)
            case default
            if (i > 1 .and. z(ivois(ii, 3)) > topo(j)) ivois(ii, 3) = -99
            end select
            select case (ivois(ii,4))
            case(-99)
            case default               
            if (i < nr .and. z(ivois(ii, 4)) < bot(j)) ivois(ii, 4) = -99
            end select
!CC...fin ajout nicolas radic
         end if
      end do
      end do

!CC....NM nombre de mailles reelles!!
      nm = kr

!CC....RECALCUL DE LA GEOMETRIE ET DU TABLEAU DE VOISINAGE
      do i = 1, nm
         x(i) = x(inum(i))
         z(i) = z(inum(i))
         am(i) = am(inum(i))
         bm(i) = bm(inum(i))
         do j = 1, 4
         if (ivois(inum(i), j) .ne. -99) then
            ivois(i, j) = inum2(ivois(inum(i), j))
!ccc....si inum2(ivois(inum(i),j))=0 c est que le voisin ds lancien maillage n est pas actif!!!
            if (ivois(i, j) == 0) ivois(i, j) = -99.
         else
            ivois(i, j) = -99
         end if
         end do
         write (11, *) i, x(i), z(i)
         write (12, *) am(i), bm(i)
         write (13, *) ivois(i, 1), ivois(i, 2), ivois(i, 3), ivois(i, 4)
      end do
      nm = nc*nr

   case (0)
!C imaille =0
      nc = nci
      nr = nri
      nm = linecell

      rewind (11)
      do i = 1, nm
         read (11, *) x(i), z(i)
         read (12, *) am(i), bm(i)
         read (13, *) ivois(i, 1), ivois(i, 2), ivois(i, 3), ivois(i, 4)
         if (icolone == 1 .or. ibuilt == 1) read (40, *) icol(i)
      end do

      close (11)
      close (12)
      close (13)
      close (40)

      do i = 1, nm
      if (ivois(i, 1) .ne. -99) then
      if (abs((x(ivois(i, 1)) - x(i))/((am(ivois(i, 1)) + am(i))/2)) &
          > 1.1 .or. &
          abs((x(ivois(i, 1)) - x(i))/((am(ivois(i, 1)) + am(i))/2)) < 0.99) then

         print *, 'pb coordonnee def maille x'
         print *, "maille=", i
         print *, "difference largeur voisine droite=", &
            (am(ivois(i, 1)) + am(i))/2
         print *, "difference coordonnee voisine droite=", &
            (x(ivois(i, 1)) - x(i))
         print *, "x coordonnee maille=", x(i)
         print *, "x coordonnee maille voisine=", x(ivois(i, 1))
         print *, "dx width maille=", bm(i)
         print *, "dx width voisine=", bm(ivois(i, 3))
         print *, "Ratio voisine droite=", &
            (x(ivois(i, 1)) - x(i))/((am(ivois(i, 1)) + am(i))/2)
         stop
      end if
      end if
      if (ivois(i, 3) .ne. -99) then

         if (abs((z(ivois(i, 3)) - z(i))/((bm(ivois(i, 3)) + bm(i))/2)) &
             > 1.1 .or. abs((z(ivois(i, 3)) - z(i))/ &
                            ((bm(ivois(i, 3)) + bm(i))/2)) < 0.99) then
            print *, 'pb coordonnee def maille z'
            print *, "maille:", i, "voisin haut:", ivois(i, 3)
            print *, "difference largeur voisine haute=", &
               (bm(ivois(i, 3)) + bm(i))/2
            print *, "difference coordonnee voisine  haute=", &
               (z(ivois(i, 3)) - z(i))
            print *, "z coordonnee maille=", z(i)
            print *, "z coordonnee maille voisine=", z(ivois(i, 3))
            print *, "dz thick maille=", bm(i)
            print *, "dz thick voisine=", bm(ivois(i, 3))
            print *, "Ratio voisine  haute=", &
               (z(ivois(i, 3)) - z(i))/((bm(ivois(i, 3)) + bm(i))/2)
            stop
         end if
      end if
      end do

      if (ibuilt == 1) then
         kr = 0
         do ii = 1, nm
         if (z(ii) - topo(icol(ii)) < 0.0001 .and. &
             z(ii) - bot(icol(ii)) > 0.0001) then
            kr = kr + 1
            inum(kr) = ii
            inum2(ii) = kr
            if (ivois(ii, 3) .ne. -99) then
               if (i > 1 .and. z(ivois(ii, 3)) > topo(icol(ii))) ivois(ii, 3) = -99
            end if
!c endif ivoi haut
            if (ivois(ii, 4) .ne. -99) then
               if (i < nr .and. z(ivois(ii, 4)) < bot(icol(ii))) ivois(ii, 4) = -99
            end if
!c endif ivoi bas
         end if
!c endif dans domain
         end do

!CC....NM nombre de mailles reelles!!
         nm = kr

!CC....RECALCUL DE LA GEOMETRIE ET DU TABLEAU DE VOISINAGE
         do i = 1, nm
            x(i) = x(inum(i))
            z(i) = z(inum(i))
            am(i) = am(inum(i))
            bm(i) = bm(inum(i))
            do j = 1, 4
            if (ivois(inum(i), j) .ne. -99) then
               ivois(i, j) = inum2(ivois(inum(i), j))
!ccc....si inum2(ivois(inum(i),j))=0 c est que le voisin ds lancien maillage n est pas actif!!!
               if (ivois(i, j) == 0) ivois(i, j) = -99
            else
               ivois(i, j) = -99
            end if
            end do
! c      write(11,*) x(i),z(i)
!     write(12,*)am(i),bm(i)
!     write(13,*)ivois(i,1),ivois(i,2),ivois(i,3),ivois(i,4)
         end do
!     call flush(11)
!ccc      call flush(12)
!ccc      call flush(13)

      end if

   END SELECT

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!   NUMEROS BAS DE COLONNE            C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

   kcol = 0
   do i = 1, nm
   if (ivois(i, 4) == -99) then
      kcol = kcol + 1
      ibas(kcol) = i
   end if
   end do

   if (icolone .ne. 1 .or. ibuilt .ne. 1) then
   do i = 1, nm

      do kkcol = 1, nc
      if (x(i) == x(ibas(kkcol))) then
         icol(i) = kkcol
      end if
      end do
   end do
   end if

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!              PARAMETRES_ALL     C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

   nmax = 5*nm + 1

   nmax1 = nm + 2

   allocate(icol_ind(nmax))
   allocate(irow_ptr(nmax1))
   allocate(val(nmax))
!         allocate(zaqui(nc))
!         allocate(zbot(nm))
   allocate(b(nm))
   allocate(pr(nm))
   allocate(pro(nm))
   allocate(conco(nm))

   allocate(izone(nm))

   allocate(vxm(nm))
   allocate(vxp(nm))
   allocate(vzp(nm))
   allocate(vzm(nm))

   allocate(ss(nm))
   allocate(dswdp(nm))
   allocate(om(nm))
   allocate(rhold(nm))
   allocate(akr(nm))
   allocate(ak(nm))
   allocate(rho(nm))
   allocate(sw(nm))
   allocate(asun(nm))
   allocate(ansun(nm))
   allocate(swresz(nm))
   allocate(akrv(nm))
   allocate(akv(nm))

   allocate(swp(nm))
   allocate(dswpdp(nm))
   allocate(dswdt(nm))
   allocate(iclt(nm, 4))
   allocate(valclt(nm, 4))
   allocate(cps(nm))
   if (iaquitard == 1) allocate(irow(nm))

   allocate(prk(nm))
   allocate(zs(nc, 2))
   allocate(zso(nc, 2))
   allocate(zsoo(nc, 2))
   allocate(rhos(nm))
   allocate(alandas(nm))
   allocate(temp(nm))
   allocate(tempo(nm))
   allocate(tempk(nm))
   allocate(tempoo(nm))
   allocate(valcl(nm, 4))
   allocate(icl(nm, 4))
   if (itr == 1) then
      allocate(iclc(nm, 4))
      allocate(valclc(nm, 4))
   end if
   if (ith == 1) then
      allocate(valclto(nm, 4))
      allocate(zl(nc, 2))
      allocate(zlo(nc, 2))
      allocate(zloo(nc, 2))
      allocate(alanda(nm))
      allocate(qtherm(nm))
      allocate(qad(nm))
      allocate(qcondu(nm))
   end if

   if (itr == 1) then
      allocate(conck(nm))
      allocate(conc(nm))
   end if

   if (ichi == 1) then
      select case (ytest)
      case ("1DS")
         allocate(chg(nm))
      case ("1DJ")
         allocate(chg(nm))
      case ("ZNS")
         allocate(chg(nm))
      case ("WAR")
         allocate(chg(nm))
      case default
         allocate(chg(nc))
      end select

   end if

   if (icycle == 1) then
      allocate(alph(nm))
      allocate(dsidtempoo(nm))
      allocate(dsidtempo(nm))
      allocate(dl(nc))
      allocate(def(nc))
      allocate(defo(nc))
      allocate(dsidp(nm))
      allocate(dsipdp(nm))
      allocate(dsipdtemp(nm))
      allocate(siceoo(nm))
   end if
   allocate(sicep(nm))
   allocate(siceo(nm))
   allocate(sice(nm))
   allocate(dsidtemp(nm))
   allocate(rhoi(nm))

   do ik = 1, nm
!ccc....MODIF AVRIL 2011
!ccc....Variation spatiale 05/04/2011
!ccc....pression initiale, alph, moy landa 05/04/2011
      ss(ik) = (sss)
      akr(ik) = akrx
      ak(ik) = akx
      akrv(ik) = akrz
      akv(ik) = akz
      asun(ik) = asp
      ansun(ik) = ans
      swresz(ik) = swres
      akc = akx
      akcv = akz
      rho(ik) = rho1
      om(ik) = omp
      sw(ik) = 1d+00
      dswdp(ik) = 0d+00
      swp(ik) = 1d+00
      dswpdp(ik) = 0d+00
      if (iss == 1) read (16, *) ss(ik)
      if (iom == 1) read (15, *) om(ik)
      if (irho == 1) read (25, *) rho(ik)
      if (ipermh == 1) read (14, *) ak(ik)
      if (ipermv == 1) read (140, *) akv(ik)
      sice(ik) = 0d+00
      dsidtemp(ik) = 0d+00
      rhoi(ik) = rhoii
      sice(ik) = 0d+00
   end do

   if (ith == 1) then
   do ik = 1, nm
      rhos(ik) = rhosi
      if (irhomi == 1) read (26, *) rhos(ik)
      alandas(ik) = (alandami)
      if (ilanda == 1) read (17, *) alandas(ik)
      cps(ik) = (cpm)
      if (icpm == 1) read (27, *) cps(ik)
   end do
   end if

   if (icycle == 1) then
   do ik = 1, nm

      siceo(ik) = 0D+00

      sicep(ik) = 0d+00
      dsipdtemp(ik) = 0d+00
      dsidp(ik) = 0d+00
      dsipdp(ik) = 0D+00
      dsidtempo(ik) = 0D+00
      dsidtempoo(ik) = 0D+00
      alph(ik) = (ss(ik)/(rho1*g))
   end do
   end if

   ito = 0
   ita = 0
   paso = 0.D+00

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!    READ DTS                             C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

   if (ytest == "DTS") then

      if (allocated(timeDTS)) deallocate(timeDTS)
      if (allocated(xDTS)) deallocate(xDTS)
      if (allocated(tempDTS)) deallocate(tempDTS)
      allocate(timeDTS(ligne))
      allocate(xDTS(nc))
      allocate(tempDTS(ligne, nc))
      if (allocated(chgRG)) deallocate(chgRG)
      if (allocated(cRivG)) deallocate(cRivG)
      if (allocated(timeG)) deallocate(timeG)
      if (allocated(tempRG)) deallocate(tempRG)
      allocate(chgRG(ligne1))
      allocate(cRivG(ligne1))
      allocate(timeG(ligne1))
      allocate(tempRG(ligne1))
      if (allocated(chgRD)) deallocate(chgRD)
      if (allocated(cRivD)) deallocate(cRivD)
      if (allocated(timeD)) deallocate(timeD)
      if (allocated(tempRD)) deallocate(tempRD)
      allocate(chgRD(ligne2))
      allocate(cRivD(ligne2))
      allocate(timeD(ligne2))
      allocate(tempRD(ligne2))

      allocate(slopeRH(2, ligne3))

      if (allocated(xriffle)) deallocate(xriffle)
      if (allocated(qout_wR)) deallocate(qout_wR)
      if (allocated(qin_wR)) deallocate(qin_wR)
      if (allocated(qout_hR)) deallocate(qout_hR)
      if (allocated(qin_hR)) deallocate(qin_hR)
      if (allocated(qcondout_hR)) deallocate(qcondout_hR)
      if (allocated(qadvout_hR)) deallocate(qadvout_hR)
      if (allocated(qcondin_hR)) deallocate(qcondin_hR)
      if (allocated(qadvin_hR)) deallocate(qadvin_hR)

      allocate(xriffle(2, ligne6))
      allocate(qout_wR(ligne6))
      allocate(qin_wR(ligne6))
      allocate(qout_hR(ligne6))
      allocate(qin_hR(ligne6))
      allocate(qcondout_hR(ligne6))
      allocate(qadvout_hR(ligne6))
      allocate(qcondin_hR(ligne6))
      allocate(qadvin_hR(ligne6))
      if (allocated(xpool)) deallocate(xpool)
      if (allocated(qout_w)) deallocate(qout_w)
      if (allocated(qin_w)) deallocate(qin_w)
      if (allocated(qout_h)) deallocate(qout_h)
      if (allocated(qin_h)) deallocate(qin_h)
      if (allocated(qcondout_h)) deallocate(qcondout_h)
      if (allocated(qadvout_h)) deallocate(qadvout_h)
      if (allocated(qcondin_h)) deallocate(qcondin_h)
      if (allocated(qadvin_h)) deallocate(qadvin_h)
      allocate(xpool(2, ligne4))
      allocate(qout_w(ligne4))
      allocate(qin_w(ligne4))
      allocate(qout_h(ligne4))
      allocate(qin_h(ligne4))
      allocate(qcondout_h(ligne4))
      allocate(qadvout_h(ligne4))
      allocate(qcondin_h(ligne4))
      allocate(qadvin_h(ligne4))

      call DTS_read(nc, xpool, ligne, ligne1, &
                    ligne2, ligne3, ligne4, &
                    ligne6, timeDTS, tempDTS, &
                    chgRG, cRivG, timeG, tempRG, &
                    chgRD, cRivD, timeD, tempRD, &
                    xriffle, slopeRH)

   end if
      
   !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
   !                                                 C
   !         Zonage des parametres         C
   !                                                 C
   !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


   SELECT CASE (ytest)

   CASE ('ZNS')

      nzone = 0
      do i = 1, nm
         read (32, *) izone(i)

         nzone = max(nzone, izone(i))

   !CC...Calcul le nombre de zone
      end do

      allocate(jzone(nzone))
      allocate(akzone(nzone))
      allocate(omzone(nzone))
      allocate(anszone(nzone))
      allocate(aspzone(nzone))
      allocate(swreszone(nzone))
      allocate(rhomzone(nzone))
   

      do j = 1, nzone
         read (321, *) jzone(j), akzone(j), omzone(j), aspzone(j), anszone(j), &
            swreszone(j)

      end do
      close(321)  ! ← garantit qu'on ne lira plus rien
      do i = 1, nm
      do j = 1, nzone
      if (izone(i) == jzone(j)) then
         ak(i) = akzone(j)
         akv(i) = akzone(j)
         om(i) = omzone(j)
         ss(i) = omzone(j)
         ansun(i) = anszone(j)
         asun(i) = aspzone(j)
         swresz(i) = swreszone(j)
      end if
!cc test zone
      end do
!cc loop zone
      end do
   CASE ('ZND')

      nzone = 0
      do i = 1, nm
         read (32, *) izone(i)

         nzone = max(nzone, izone(i))

   !CC...Calcul le nombre de zone
      end do

      allocate(jzone(nzone))
      allocate(akzone(nzone))
      allocate(omzone(nzone))
      allocate(anszone(nzone))
      allocate(aspzone(nzone))
      allocate(swreszone(nzone))
      allocate(rhomzone(nzone))

      do j = 1, nzone
         read (321, *) jzone(j), akzone(j), omzone(j), aspzone(j), anszone(j), &
            swreszone(j)

      end do
      close(321)  ! ← garantit qu'on ne lira plus rien
      do i = 1, nm
      do j = 1, nzone
      if (izone(i) == jzone(j)) then
         ak(i) = akzone(j)
         akv(i) = akzone(j)
         om(i) = omzone(j)
         ss(i) = omzone(j)
         ansun(i) = anszone(j)
         asun(i) = aspzone(j)
         swresz(i) = swreszone(j)
      end if
!cc test zone
      end do

!cc loop zone
      end do

   CASE ('1DS')
      nzone = 0
      do i = 1, nm
         read (32, *) izone(i)
         nzone = max(nzone, izone(i))
   !CC...Calcul le nombre de zone
      end do
      
      allocate(alandazone(nzone))
      allocate(rhomzone(nzone))
      allocate(akzone(nzone))
      allocate(omzone(nzone))
      allocate(aspzone(nzone))
      allocate(jzone(nzone))
      do i = 1, nm
         read (32, *) izone(i)
         nzone = max(nzone, izone(i))
!CC...Calcul le nombre de zone
      end do

      allocate(alandazone(nzone))
      allocate(cpmzone(nzone))
      allocate(rhomzone(nzone))
      allocate(jzone(nzone))
      allocate(akzone(nzone))
      allocate(omzone(nzone))
      allocate(anszone(nzone))
      allocate(aspzone(nzone))
      allocate(swreszone(nzone))
      do j = 1, nzone
         read (321, *) jzone(j), akzone(j), omzone(j), anszone(j), aspzone(j), &
            swreszone(j), alandazone(j), cpmzone(j), rhomzone(j)
      end do
      close(321)  ! ← garantit qu'on ne lira plus rien
      do i = 1, nm
      do j = 1, nzone
      if (izone(i) == jzone(j)) then
         alandas(i) = alandazone(j)
         cps(i) = cpmzone(j)
         rhos(i) = rhomzone(j)
         ak(i) = akzone(j)
         akv(i) = akzone(j)
         om(i) = omzone(j)
         ss(i) = omzone(j)
         ansun(i) = anszone(j)
         asun(i) = aspzone(j)
         swresz(i) = swreszone(j)
         rhos(i) = rhomzone(j)
      end if


!cc test zone
      end do
!cc loop zone
      end do

   CASE ('1DJ')
      nzone = 0
      do i = 1, nm
         read (32, *) izone(i)
         nzone = max(nzone, izone(i))
!CC...Calcul le nombre de zone
      end do

      allocate(jzone(nzone))
      allocate(akzone(nzone))
      allocate(omzone(nzone))
      allocate(anszone(nzone))
      allocate(aspzone(nzone))
      allocate(swreszone(nzone))      
      do j = 1, nzone
         read (321, *) jzone(j), akzone(j), omzone(j), anszone(j), aspzone(j), &
            swreszone(j)

      end do
      close(321)  ! ← garantit qu'on ne lira plus rien
      do i = 1, nm

         do j = 1, nzone
         if (izone(i) == jzone(j)) then
            ak(i) = akzone(j)
            akv(i) = akzone(j)
            om(i) = omzone(j)
            ss(i) = omzone(j)
            ansun(i) = anszone(j)
            asun(i) = aspzone(j)
            swresz(i) = swreszone(j)
         end if
!cc test zone
         end do
!cc loop zone
      end do
!cc loop element
 
   CASE ('AVA', 'TEX')
      do i = 1, nm
         read (32, *) izone(i)
         nzone = max(nzone, izone(i))
!CC...Calcul le nombre de zone
      end do
      allocate(cpmzone(nzone))
      allocate(anszone(nzone))
      allocate(swreszone(nzone))

      do j = 1, nzone
         read (321, *) jzone(j), akzone(j), omzone(j), anszone(j), aspzone(j), &
            swreszone(j), alandazone(j), cpmzone(j), rhomzone(j)
      end do
      close(321)  ! ← garantit qu'on ne lira plus rien
      do j = 1, nzone
            do i = 1, nm
      if (izone(i) == jzone(j)) then
         ak(i) = akzone(j)
         akv(i) = akzone(j)
         om(i) = omzone(j)
         ss(i) = omzone(j)
         ansun(i) = anszone(j)
         asun(i) = aspzone(j)
         swresz(i) = swreszone(j)
        alandas(i) = alandazone(j)
         cps(i) = cpmzone(j)
         rhos(i) = rhomzone(j)
      end if
!cc test zone
      end do
!cc loop zone
      end do


   CASE ('R2D')
      do i = 1, nm
         read (32, *) izone(i)
         nzone = max(nzone, izone(i))
!CC...Calcul le nombre de zone
      end do

      allocate(cpmzone(nzone))
      allocate(jzone(nzone))
      allocate(akzone(nzone))
      allocate(omzone(nzone))
      allocate(rhomzone(nzone))
      allocate(alandazone(nzone))
      print *, "nzone=", nzone
      do j = 1, nzone
         read(321, *, iostat=ierr) jzone(j), akzone(j), omzone(j), &
                                    alandazone(j), cpmzone(j), rhomzone(j)
         if (ierr /= 0) then
            print *, " Fichier trop court ou mal formé à la ligne ", j
            exit
         end if
      end do

      close(321)  ! ← garantit qu'on ne lira plus rien
      

      do j = 1, nzone
            do i = 1, nm
      if (izone(i) == jzone(j)) then
         ak(i) = akzone(j)
         akv(i) = akzone(j)
         om(i) = omzone(j)
         ss(i) = omzone(j)
        alandas(i) = alandazone(j)
         cps(i) = cpmzone(j)
         rhos(i) = rhomzone(j)
      end if
!cc test zone
      end do
!cc loop zone
      end do      
!cc loop element
   CASE ('DTS')
         nzone = 0
      do i = 1, nm
         read (32, *) izone(i)
         nzone = max(nzone, izone(i))
   !CC...Calcul le nombre de zone
      end do
      
      allocate(alandazone(nzone))
      allocate(rhomzone(nzone))
      allocate(akzone(nzone))
      allocate(omzone(nzone))
      allocate(aspzone(nzone))
      allocate(jzone(nzone))


      do j = 1, nzone
         read (321, *) jzone(j), akzone(j), omzone(j), &
            alandazone(j), rhomzone(j)

      end do
      close(321)  ! ← garantit qu'on ne lira plus rien
      do i = 1, nm
      do j = 1, nzone
      if (izone(i) == jzone(j)) then
         ak(i) = akzone(j)
         akv(i) = akzone(j)
         om(i) = omzone(j)
         ss(i) = omzone(j)
         alandas(i) = alandazone(j)
         rhos(i) = rhomzone(j)
      end if
!cc test zone
      end do

!cc loop zone

!     if(i == 10001) print*,ak(i)
!     if(i == 12101) print*,ak(i)
!     if(i == 47405) print*,ak(i)
      end do
!cc loop element

   CASE ('ZHZ')
      nzone = 0
      do i = 1, nm
         read (32, *) izone(i)
         nzone = max(nzone, izone(i))
   !CC...Calcul le nombre de zone
      end do
      if (allocated(alandazone)) deallocate(alandazone)
      if (allocated(rhomzone)) deallocate(rhomzone)
      if (allocated(akzone)) deallocate(akzone) 
      if (allocated(omzone)) deallocate(omzone)
      if (allocated(aspzone)) deallocate(aspzone)
      if (allocated(jzone)) deallocate(jzone)

      allocate(alandazone(nzone))
      allocate(rhomzone(nzone))
      allocate(akzone(nzone))
      allocate(omzone(nzone))
      allocate(aspzone(nzone))
      allocate(jzone(nzone))


      do j = 1, nzone
         read (321, *) jzone(j), akzone(j), omzone(j), &
            alandazone(j), rhomzone(j)
      end do
      close(321)  ! ← garantit qu'on ne lira plus rien
      do i = 1, nm
      do j = 1, nzone
      if (izone(i) == jzone(j)) then
         ak(i) = akzone(j)
         akv(i) = akzone(j)
         om(i) = omzone(j)
         ss(i) = omzone(j)
         alandas(i) = alandazone(j)
         rhos(i) = rhomzone(j)
      end if
!cc test zone
      end do
!cc loop zone
      end do
!cc loop element
   END SELECT

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!        INITIALISATION DES PARAMETRES        C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC....ISOTHERMES
   do kkcol = 1, nc
   do kg = 1, 2
      zs(kkcol, kg) = -99
      zsoo(kkcol, kg) = -99
      zso(kkcol, kg) = -99
   end do
   end do

!CC....INITIALISATION PRESSURE
   do ik = 1, nm
   if (ichi2 .ne. 1 .and. ichi .ne. 1) then
      pr(ik) = (rho1*g*(chgi - z(ik)))
!       write(75,*) x(ik),z(ik),pr(ik)
   end if
   end do

!CC....INITIALISATION CONCENTRATION
   if (itr == 1) then
   do ik = 1, nm
      conc(ik) = (conci)
      if (iconci == 1) read (22, *) conc(ik)
   end do
   end if

!CC....INITIALISATION TEMPERATURE
   do ik = 1, nm
      tempoo(ik) = (tempini)
      tempo(ik) = (tempini)
      temp(ik) = (tempini)
   end do

   if (ith == 1) then
   if (itempi == 1) then
   do ik = 1, nm
      read (23, *) temp(ik)
      tempo(ik) = temp(ik)
      tempoo(ik) = temp(ik)
   end do
   end if
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!         Thermal conductivity      C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!CC....THERMIQUE CONDUCTIVITE MOYENNE SANS GEL SANS ZONE NON SATUREE
   SELECT CASE (ymoycondtherm)
   CASE ("WOODS")
      do ik = 1, nm
         alanda(ik) = ((sqrt(alandae)*om(ik) + &
                        sqrt(alandas(ik))*(1D+00 - om(ik)))**2)
      end do
   CASE ("GEOME")
      do ik = 1, nm
         alanda(ik) = (alandae**(om(ik))*alandas(ik)**(1 - om(ik)))
      end do
   CASE ("ARITH")
      do ik = 1, nm
         alanda(ik) = (alandae*(om(ik)) + alandas(ik)*(1 - om(ik)))
      end do
   CASE ("LUNAR")
      do ik = 1, nm
         alanda(ik) = (2.417196D+00)
      end do
   CASE ("NEUMA")
      do ik = 1, nm
         alanda(ik) = 2.619D+00
      end do
   END SELECT

!CC....ISOTHERMES
   do kkcol = 1, nc
   do kg = 1, 2
      zl(kkcol, kg) = -99
      zloo(kkcol, kg) = -99
      zlo(kkcol, kg) = -99
   end do
   end do

   if (icycke == 1) then
   do kkcol = 1, nc
      dl(kkcol) = 0.D00
      def(kkcol) = 0.D00
!     zaqui(kkcol)=-99
   end do
   end if

   end if

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!    VARIATION DE LA CHARGE INITIALE        C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC...HYDROSTATIQUE
   if (ichi == 1) then
      SELECT CASE (ytest)
      CASE ("ZNS")
         do i = 1, nm
            read (24, *) chg(i)
            if (abs(chg(j)) < 10D-9) chg(j) = 0D+00
            pr(i) = dble((rho(i)*g*(chg(i) - z(i))))
         end do
      CASE ("ZND")
         do i = 1, nm
            read (24, *) chg(i)
            if (abs(chg(j)) < 10D-9) chg(j) = 0D+00
            pr(i) = dble((rho(i)*g*(chg(i) - z(i))))
         end do
      CASE ("WAR")
         do i = 1, nm
            read (24, *) chg(i)
            if (abs(chg(j)) < 10D-9) chg(j) = 0D+00
            pr(i) = dble((rho(i)*g*(chg(i))))
         end do
      CASE ("1DS")
         do i = 1, nm
            read (24, *) chg(i)
            if (abs(chg(j)) < 10D-9) chg(j) = 0D+00
            pr(i) = dble((rho(i)*g*(chg(i) - z(i))))
         end do

      CASE ("1DJ")
         do i = 1, nm
            read (24, *) chg(i)
            if (abs(chg(j)) < 10D-9) chg(j) = 0D+00
            pr(i) = dble((rho(i)*g*(chg(i) - z(i))))
         end do

      case default
         do j = 1, nc
            read (24, *) chg(j)
         end do
         do j = 1, nc
         do ii = 1, nm
         if (icol(ii) == j) then
            if (abs(chg(j)) < 10D-9) chg(j) = 0D+00
            pr(ii) = ((rho(ii)*g*(chg(j) - z(ii))))
            if (abs(pr(ii)) < 10D-9) pr(ii) = 0D+00
         end if
         end do
         end do
      END SELECT
   end if

!CC...VARIATION DE LA CHARGE INITIALE sur tout le model
   if (ichi2 == 1) then
   do i = 1, nm
      read (242, *) pr(i)
   end do
   end if

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!              BASSINS                C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC...AQUITARD TEST PETIT BASSIN 1 COUCHE AQUITARD
!       if(iaquitard == 1) then
!       do i=1,nm
!       read(21,*) irow(i)
!       if(irow(i) == nrowtard) then
!       ak(i)=(aktardx)
!       akv(i)=(aktardz)
!       om(i)=(omptard)
!       alandas(i)=(alandatard)
!       ss(i)=(sstard)
!       do j=1,nc
!       if(icol(i) == j) zaqui(j)=z(i)
!       enddo
!       endif
!       do jk=nrowtard-nmailleaqui+1,nrowtard-1
!          if(irow(i) == jk) then
!          ak(i)=(aktardx)
!          akv(i)=(aktardz)
!          om(i)=(omptard)
!          alandas(i)=(alandatard)
!          ss(i)=(sstard)
!          endif
!          enddo
!       enddo
!       endif
!CC...Terme de surpression
!     if (ysupdp == "SPWAL") then
!         if(iaquitard == 0) then
!         do i=1,nm
! c            do j=1,nc
!         if(icol(i) == j) zbot(i)=bot(j)
!         enddo
!         enddo
!         endif
!       if(iaquitard == 1) then
!          do j=1,nc
!          do i=1,nm
!          if(icol(i) == j) then
!          if (zaqui(j).ne.-99.and.z(i) >= zaqui(j)) zbot(i)=zaqui(j)
!          if (zaqui(j) == -99) zbot(i)=bot(j)
!          if (z(i) < zaqui(j)) zbot(i)=bot(j)
!              endif
!          enddo
!          enddo
!          endif
!     endif
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!       temperature liquidus solidus        C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
   ts = tsg
   tl = tlg

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
! CONDITIONS LIMITES SANS VARIATION SPATIALE  C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC...ECOULEMENT
!CC...Tableau ICL et VALCL
!ccc....-1 Flux impose sur une face
!ccc....-2 potentiel impose sur une face
!ccc.... 1 Mailles 'normale'
!ccc.... ICL(ik,1)=-1, -2 ou 1 ik=1,...4: face droite,gauche,haute et BASSE
!ccc....Conditions aux limites appliquees sur les faces!!!!
!ccc....VALCL(ik,1), valeur de la condition limite sur ces meme faces..
!CC....INITIALISATION DES CONDITIONS AUX LIMITES=FLUX NUL
   do i = 1, nm

      do j = 1, 4
         icl(i, j) = 1
!CC....FLUX NUL PAR DEFAUT
         if (ivois(i, j) == -99) icl(i, j) = -1
         valcl(i, j) = 0.D+00
      end do
   end do

!CC....ECOULEMENT
   do i = 1, nm
   if (ivois(i, 1) == -99) then
      icl(i, 1) = icl_droite
      valcl(i, 1) = dble(valcl_droite)
   end if

   if (ivois(i, 2) == -99) then
      icl(i, 2) = icl_gauche
      valcl(i, 2) = dble(valcl_gauche)
   end if
   if (ivois(i, 3) == -99) then
      icl(i, 3) = icl_haut
      valcl(i, 3) = dble(valcl_haut)
   end if
   if (ivois(i, 4) == -99) then
      icl(i, 4) = icl_bas
      valcl(i, 4) = dble(valcl_bas)
   end if

!ccc....VARIATION SPATIALE : une condition par maille
   if (iclect == 1) then
      read (37, *) icl(i, 1), icl(i, 2), icl(i, 3), icl(i, 4), valcl(i, 1), &
         valcl(i, 2), valcl(i, 3), valcl(i, 4)
!       else
!       if(iclect.ne.1) then
!       write(37,*) icl(i,1),icl(i,2),icl(i,3),icl(i,4),valcl(i,1),
!     &valcl(i,2),valcl(i,3),valcl(i,4)
!       endif
   end if

!ccc....message erreur
   do j = 1, 4
   if (ivois(i, j) == -99 .and. icl(i, j) == 1) then
      print *, i, &
         "erreur pas de voisin et pas de conditions aux limites"
   end if
   end do
   end do
   do i = 1, nm
!ccc....CHARGE IMPOSEE aux faces (attention pas au centre des mailles !!!)
      if (icl(i, 1) == -2) then
         valcl(i, 1) = (rho(i)*g*(valcl(i, 1) - z(i)))
!ccc....TEST INTERFROST TH2 TH3
         if (ytest == "TH2" .or. ytest == "TH3") then
            valcl(i, 1) = (rho(i)*g*valcl_droite)
         end if
!ccc....VALEURS NULLES
         if (abs(valcl(i, 1)) < 10D-10) valcl(i, 1) = 0D+00
      end if

      if (icl(i, 2) == -2) then
         valcl(i, 2) = (rho(i)*g*(valcl(i, 2) - z(i)))
!ccc....TEST INTERFROST TH2 TH
         if (ytest == "TH2" .or. ytest == "TH3") then
            valcl(i, 2) = (rho(i)*g*valcl_gauche)
         end if
!ccc....VALEURS NULLES
         if (abs(valcl(i, 2)) < 10D-10) valcl(i, 2) = 0
      end if

      if (icl(i, 3) == -2) then
         valcl(i, 3) = dble(rho(i)*g*(valcl(i, 3) - z(i) - bm(i)/2))
!ccc....VALEURS NULLES
         if (ytest == "WAR") valcl(i, 3) = valcl_haut*rho1*g
         if (abs(valcl(i, 3)) < 10D-10) valcl(i, 3) = 0D+00
      end if

      if (icl(i, 4) == -2) then
         valcl(i, 4) = (rho(i)*g*(valcl(i, 4) - z(i) + bm(i)/2))

         if (ytest == "WAR") valcl(i, 4) = valcl_bas*rho1*g

!ccc....VALEURS NULLES
         if (abs(valcl(i, 4)) < 10D-10) valcl(i, 4) = 0
      end if

      if (ytest == "ZNS" .or. ytest == "ZND") then
         if (icl(i, 4) == -2) then
            zbas = z(nm) - bm(nm)/2
            if (abs(zbas) < 1e-6) zbas = 0D+00
            valcl(nm, 4) = (rho(nm)*g*(valcl_bas - zbas))
         end if

         if (icl(i, 3) == -2) then
            zhaut = z(1) + bm(1)/2
            if (abs(zhaut) < 1e-6) zhaut = 0D+00
            valcl(1, 3) = rho(1)*g*(valcl_haut - zhaut)
         end if
      end if
   end do

!ccc....CDT MAILLES RIVIERE
   if (iclriviere == 1) then
   do i = 1, nm
   if (ivois(i, 3) == -99) then
      icl(i, 3) = -1
      valcl(i, 3) = qre
   end if
   end do
   end if

!CC....TRANSPORT
!ccc....Tableau ICL et VALCL
!ccc....TABLEAU ICLC 1 NORMAL, -1 FLUX NUL, -2 CONCENTRATION IMPOSEE SUR LES FACES CORRESPONDANTES
!ccc....TABLEAU VALCLC VALEURS CORRESPONDANTES
!ccc....N.B.: LES FLUX IMPOSES (ICLC(i,j)=-1) SONT TRAITES VIA LE TERME ADVECTIF
   if (itr == 1) then
      do ik = 1, nm
      do j = 1, 4
         iclc(ik, j) = 1
!ccc....FLUX NUL PAR DEFAUT
         if (ivois(ik, j) == -99) iclc(ik, j) = -1
         valclc(ik, j) = 0.D+00
      end do
      end do
   end if

!CC....THERMIQUE
!ccc....Conditions limites
!ccc....TABLEAU ICLT 1 NORMAL, -1 FLUX NUL, -2 CONCENTRATION IMPOSEE SUR LES FACES CORRESPONDANTES
!ccc....TABLEAU VALCLT VALEURS CORRESPONDANTES
!ccc....N.B.: LES FLUX IMPOSES (ICLC(i,j)=-1) SONT TRAITES VIA LE TERME ADVECTIF
   if (ith == 1) then
      do ik = 1, nm
      do j = 1, 4
         iclt(ik, j) = 1
      end do
      if (ivois(ik, 2) == -99) then
         iclt(ik, 2) = iclt_gauche
         valclt(ik, 2) = dble(valclt_gauche)
      end if
      if (ivois(ik, 1) == -99) then
         iclt(ik, 1) = iclt_droite
         valclt(ik, 1) = dble(valclt_droite)
      end if
      if (ivois(ik, 3) == -99) then
         iclt(ik, 3) = iclt_haut
         valclt(ik, 3) = dble(valclt_haut)
      end if
      if (ivois(ik, 4) == -99) then
         iclt(ik, 4) = iclt_bas
         valclt(ik, 4) = dble(valclt_bas)
      end if

!ccc....VARIATION SPATIALE : une condition par maille
      if (icltherm == 1) then
         read (38, *) iclt(ik, 1), iclt(ik, 2), iclt(ik, 3) &
            , iclt(ik, 4), valclt(ik, 1), valclt(ik, 2), valclt(ik, 3), valclt(ik, 4)
!     else
!     write(38,*) iclt(ik,1),iclt(ik,2),iclt(ik,3)
!     &,iclt(ik,4),valclt(ik,1),valclt(ik,2),valclt(ik,3),valclt(ik,4)
      end if
      do j = 1, 4
      if (ivois(ik, j) == -99 .and. iclt(ik, j) == 1) then
         print *, ik, &
            "erreur pas de voisin et pas de conditions aux limites"
      end if
      end do
      end do
   end if

!CC....EXUTOIRE
!       if(iexutoire == 1) then
!       do i=1,nm
!       if(x(i) < xexutoire.and.ivois(i,3) == -99) then
!       valclt(i,3)=tempexutoite
!       iclt(i,3)=-2
!       endif
!       enddo
!       endif
!       qruis=0.D+00
!      if (iinfil == 1) then
!       call infiltration(pr,nm,n,icol,nc,om,ss,bm,rho,pore,
!     &am,valcl,icl,qtpore,ruis,ivois,z,g,dt,xberg,x,qre,qruis)
!       endif
!       if (iriv == 1) then
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!         RIVIERE               C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC....Conditions limites
!ccc....CONDITION LIMITES BORDS BERGE RIVIERE TEST CELLULES VOISINES DE LA RIVIERE
!ccc....CONDITION TEMPERATURE IMPOSE SUR BERGE ET LIT DE RIVIERE
!ccc....MAILLE VOISINE RIVIERE PERMEABILITE DU LIT DE LA RIVIERE
!       if(iqriv == 1) then
!       hriv=((qre*xberg*1+qriva+qruis)**(3./5.))*(rug**(-3./5.))*
!     &((al-xberg)**(-3./5.))*(pent**(-3./10.))
!       endif
!ccc....CDT LIMITE RIVIERE
!       call cdt_riviere(hriv,g,nm,z,hbot,
!     &bm,ivois,am,rho,xberg,x,al,n,icl,valcl,iclt,valclt,
!     &tempriv,aklit,aklitv,ak,akv,akr,akrv,dsw,sw,elit,
!     &akc,akcv,it,ita,tempo,ts)
!       endif
!ccc....Variation conditions limites lecture fichier
!       if(iclriviere == 1) then
!       do while (iot == 0)
!       read(42,*,iostat=iot)
!       if (iot == 0) then
!       ligne1=ligne1+1
!       endif
!       enddo
!       rewind(42)
!       do j=1,ligne1-1
!       read(42,*,iostat=iot),qrivent1(j),qpluie1(j)
!       enddo
!       do while (io2 == 0)
!       read(43,*,iostat=io2)
!       if (io2 == 0) then
!       ligne2=ligne2+1
!       endif
!       enddo
!       rewind(43)
!       do j=1,ligne2
!       read(43,*,iostat=io2),qrivent2(j),qpluie2(j)
!       enddo
!       do while (io3 == 0)
!       read(44,*,iostat=io3)
!       if (io3 == 0) then
!       ligne3=ligne3+1
!       ENDIF
!       enddo
!       rewind(44)
!       do j=1,ligne3
!       read(44,*,iostat=io3),qrivent3(j),qpluie3(j)
!       enddo
!       endif
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!      Topo en marche d escalier          C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!       if (itopomch == 1) then
!     do i=1,nm
!     if(x(i) < xberg-am(i).and.ivois(i,1) == -99
!     &.and.ivois(i,3) == -99) then
!       valcl(i,1)=qre/2
!       valcl(i,3)=qre/2
!     icl(i,1)=-1
!     icl(i,3)=-1
!       valclt(i,1)=10D+00
!     valclt(i,3)=10D+00
!     iclt(i,1)=-2
!     iclt(i,3)=-2
!     endif
!     write(37,*) icl(i,1),icl(i,2),icl(i,3),icl(i,4),valcl(i,1),
!     &valcl(i,2),valcl(i,3),valcl(i,4)
!     enddo
!     endif

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
! CDT INTERFROST INITIALE ET GEOMETRIE    C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!CC....CDT TEST TH2
   if (ytest == "TH2") then
   do i = 1, nm
   if (x(i) >= 1 - 0.333/2 .and. x(i) <= 1 + 0.333/2.) then
   if (z(i) >= 0.5 - 0.333/2 .and. z(i) <= 0.5 + 0.333/2.) then
      temp(i) = -5D+00
      sice(i) = 0.95D+00
      tempo(i) = -5D+00
      siceo(i) = 0.95D+00
      sw(i) = 0.05D+00
      akr(i) = 1D-06
      akrv(i) = 1D-06

   end if
   end if
   tempo(i) = (temp(i))
   end do
   end if
!CC....CDT TEST TH3 demi cercle
   if (ytest == "TH3") then
   do i = 1, nm
      ab = (x(i) - 0.5D+00)**2 + (z(i) + 0.1D+00)**2
      if (ab <= (0.5099)**2) then
         temp(i) = -5D+00
         sice(i) = 0.95D+00
         siceo(i) = 0.95D+00
         sw(i) = 0.05D+00
         akr(i) = 1D-06
         akrv(i) = 1D-06
      else
         temp(i) = 5D+00
         sice(i) = 0D+00
         siceo(i) = 0D+00
         sw(i) = 1D+00
      end if
      tempo(i) = (temp(i))
   end do
   end if

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!     MESSAGE UTILISATEUR              C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      if (iec == 0) then
!       print*,'attention vous avez annule l ecoulement'
!       print*,'vous devez indiquer des flux imposes en ecoulement'
!       print*,'dans le cas contraire calcul de conduction'
!      endif
!      if  (irp == 1) then
!      Print*,'nb de lignes dans les fichiers temporaires ='
!      print*,nitt*unitsim/(itsortie*us)
!      Print*,'attention cela doit etre un entier'
!      endif
!     if (iec == 1) then
!      do i=1,nm
!       if (icl(i,1) == -2) then
!       if (ivois(i,2) == -99) then
!       print*,'pb de configuration 2 cellule par ligne'
!       print*,'avec une charge imposee'
!       endif
!       endif
!       if (icl(i,2) == -2) then
!       if (ivois(i,1) == -99) then
!       print*,'pb de configuration 2 cellule par ligne'
!       print*,'avec une charge imposee'
!       endif
!       endif
!       if (icl(i,3) == -2) then
!       if (ivois(i,4) == -99) then
!       print*,'pb de configuration 2 cellule par ligne'
!       print*,'avec une charge imposee'
!       endif
!       endif
!       if (icl(i,4) == -2) then
!       if (ivois(i,3) == -99) then
!       print*,'pb de configuration 2 cellule par ligne'
!       print*,'avec une charge imposee'
!       endif
!       endif
!       enddo
!       endif
!       if(ithec == 1) then
!       do j=1,4
!       do i=1,nm
!       if (icl(i,j) == -1.and.iclt(i,j) == -1) then
!       if(abs(valcl(i,j)) > 0) then
!c      print*,'Attention si vous imposez un flux nul en temperature'
!       print*,'vous allez convecter une temperature nulle'
!c      endif
!       endif
!       enddo
!       enddo
!       endif

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!      TERME dswdp  CDT INITIALE          C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC....DEF PARAMETRES VS CDT INITIALES
!CC....ZNS
   if (ivg == 1 .or. yunconfined == "UNS") then
      Call unsaturated(pr, swp, dswpdp, swresz, ans, akr, &
                       nm, akrv, rho1, g, ansun, asun)

!CC....THERMIQUE CONDUCTIVITE MOYENNE SANS GEL AVEC ZONE NON SATUREE
      if (igelzns == 0 .and. ith == 1) then
      if (ymoycondtherm == "WOODS") then
         alanda(i) = DBLE(sqrt(alandae)*om(i)*sw(i) + &
                          sqrt(alandas(i))*(1D+00 - om(i)) + &
                          sqrt(alandag)*om(i)*(1D+00 - sw(i)))**2
      else if (ymoycondtherm == "GEOME") then
         alanda(i) = DBLE(alandae**(om(i)*sw(i))*alandas(i)**(1 - om(i))* &
                          alandag**(om(i)*(1D+00 - sw(i))))
      else if (ymoycondtherm == "ARITH") then
         alanda(i) = DBLE(alandae*(om(i)*sw(i)) + alandas(i)*(1 - om(i)) &
                          + alandag*(om(i)*(1D+00 - sw(i))))
      end if

      end if

   else if (yunconfined == "CAP") then
      do i = 1, nm
!CC....nappe captive
!ccc....COEFFICIENT D EMMAGASINEMENT
         swp(i) = 1.D00
         dswpdp(i) = ss(i)/(rho1*g*om(i))
!ccc....Test Interfrost
         if (ytest == "TH2" .or. ytest == "TH3") then
            dswpdp(i) = sw(i)*ss(i)
         end if
      end do
!CC....nappe libre
!ccc....simplification solution analytique
   else if (yunconfined == "UNC") then
      do i = 1, nm
         akr(i) = 1D+00
         akrv(i) = 1D+00
         swp(i) = 1D+00
         if (irp == 1) dswpdp(i) = 1D+00/(rho1*g*bm(i))
!       if (irp == 1) dswpdp(i)=1D+00/(rho1*g)
         if (pr(i) < 0d+00) then
!     akr(i)=0D+00
            akrv(i) = 1D+00
            swp(i) = 0D+00
         end if
         if (pr(i) >= 0d+00) then
            akr(i) = 1D+00
            akrv(i) = 1D+00
            swp(i) = 1D+00
         end if
      end do
   end if

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!   PARAMETRE VS. TEMPERATURE INITIALE    C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC....Indice gel/degel
!ccc....test cycle gel valeur de igel!!!!! Attention test seulement sur la maille 1 à améliorer

   if (icycle .ne. 0) then
      if (temp(1) < 0 .and. tempo(1) - temp(1) > 0) igel = 1
      if (temp(1) < 0 .and. tempo(1) - temp(1) < 0) igel = 2
      if (ytest == "TH1" .or. ytest == "TH2" .or. ytest == "TH3") igel = 2
      if (ytest == "MAQ") then
         if (paso <= 49.1*86400) igel = 1
         if (paso > 49.1*86400) igel = 2
      end if
!       if(ytest == "MAQ") then
!     if (paso <= 1180*3600) igel=1
!     if (paso > 1180*3600) then
!     igel=2
!     endif
   else
      igel = 0
   end if

   if (icycle == 1) then
   if (ytest == "TH2" .or. ytest == "TH3") then
   do i = 1, nm
      dswpdp(i) = ss(i)*sw(i)
   end do
   if (ytest == "MAQ") then
      do i = 1, nm
!CC....Changement coef emmagasinement
         if (temp(i) < tl) dswpdp(i) = ss(i)/g/rho(i)/om(i)
         if (temp(i) > tl) dswpdp(i) = 1
         do kcol = 1, nc
!CC...GEL....dswpd subpermafrost
            if (igel == 1 .and. zs(kcol, 1) .ne. -99 .and. icol(i) == kcol .and. &
                z(i) < zs(kcol, 1) .and. pro(i)/(rho(i)*g) + z(i) > zs(kcol, 1)) then
               dswpdp(i) = ss(i)/g/rho(i)/om(i)
            end if
!CC...GEL....dswpd subpermafrost
            if (igel == 1 .and. zs(kcol, 1) .ne. -99 .and. icol(i) == kcol .and. &
                z(i) < zs(kcol, 1) .and. pro(i)/(rho(i)*g) + z(i) > topo(kcol)) &
               then
               dswpdp(i) = ss(i)/g/rho(i)/om(i)
            end if
!CC....DEGEL....dswpd suprapermafrost
            if (igel == 2 .and. zs(kcol, 1) .ne. -99 .and. icol(i) == kcol .and. &
                z(i) > zs(kcol, 1) .and. temp(i) > tl) then
               dswpdp(i) = 1
            end if
!CC....DEGEL....dswpd subpermafrost
            if (igel == 2 .and. zs(kcol, 2) .ne. -99 .and. icol(i) == kcol .and. &
                z(i) < zs(kcol, 2) .and. temp(i) > tl) then
               dswpdp(i) = ss(i)/(rho1*g*om(i))
            end if
         end do
      end do
   end if
   end if
!CC....Calcul saturation permeabilite
!ccc....Cas completement sature en eau
!ccc....derivation de la saturation en glace en fonction de la temperature dsice
!ccc....sw (liquide) = (1D0 - sice)
!CC....Implicite
   if (icycle .ne. 0 .and. ibigridice .ne. 1 .and. iparo == 1) then
      CALL icesatperm(nm, tl, ts, akr, akrv, dk, tempo, &
                      dsipdtemp, igelzns, ytypakrice, &
                      sicep, om, sw, swressi, ytypsice, &
                      cimp, tempoo, omega, siceo)
   else
      do i = 1, nm
         sicep(i) = 0.D0
         dsipdtemp(i) = 0.D0
      end do
   end if
!CC....Explicite
   if (icycle .ne. 0 .and. ibigridice .ne. 1 .and. iparo == 0) then
      CALL icesatperm(nm, tl, ts, akr, akrv, dk, temp, &
                      dsipdtemp, igelzns, ytypakrice, &
                      sicep, om, sw, swressi, ytypsice, &
                      cimp, tempo, omega, siceo)
   else
      do i = 1, nm
         sicep(i) = 0.D0
         dsipdtemp(i) = 0.D0
      end do
   end if
!CC....grosse maille saturation en glace et permeabilite relative calculees par rapport a la position des isothermes dans les mailles.
   if (icycle .ne. 0 .and. ibigridice == 1) then
      CALL biggridice(nm, tl, ts, akr, akrv, dk, temp, igel, &
                      nc, z, zl, zs, bm, dsipdtemp, valclt, ivois, &
                      sicep, swressi, siceo, tempo)
   end if

!ccc.....SWTOT est la quantite total d'eau (liquide+glace)
   do i = 1, nm
      swtot = swp(i)
      sice(i) = MAX(0.D0, swtot - (1.D0 - sicep(i)))
      sw(i) = swtot - sice(i)
      if (dsipdtemp(i) == 0) then
         dswdp(i) = dswpdp(i)
         dswdt(i) = 0.D0
         dsidp(i) = 0.D0
         dsidtemp(i) = 0.D0
      else
         dswdp(i) = ss(i)/(rho1*g*om(i))
!ccc....nappe captive
         dswdt(i) = -dsipdtemp(i)
         dsidp(i) = dswpdp(i)
         dsidtemp(i) = dsipdtemp(i)

      end if
      if (irp == 0) dswdp(i) = 0D+00

   end do
   else
   do i = 1, nm
      dswdp(i) = dswpdp(i)
   end do
   end if


 !!!!!!!!!!!!!!!!!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                               C
!     OUVERTURE FICHIERS SORTIES                                C
!                                                               C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!CC....FICHIER FIN DE SIMULATION
      open (74, file='S_pression_charge_temperature.dat')
!CC....FICHIERS A TRIER
!       open(49,file='S_profil_temp_colonne_choisie_t.dat')
!       open(91,file='S_pression_charge_temperature_t.dat')
!       open(62,file='S_bound_permafrost_t.dat')
!       open(63,file='S_soulevemnt_t.dat')
!       open(67,file='S_debit_riviere_t.dat')
!       open(782,file='S_moy_vitesses.dat')
!       open(786,file='S_flux.dat')
!       open(52,file='S_pression_charge_temperature_maillechoisie.dat')
!       open(53,file='S_charge_capteurs.dat')
!       OPEN(94,FILE='S_permeabilite_t.dat')
!       open(96,FILE='S_surpression_t.dat')
      select case(ytest)
!CC....FICHIERS these texier
      case("TEX") 
      if (irpth == 1) then
         open (7782, file='S_fichier_tot.dat')
      end if

!CC....FICHIERS ZNS 1D
      case ("ZNS") 
         open (1818, file='S_hydraulic_conductivities_profil_t.dat')
         open (18181, file='S_saturation_profil_t.dat')
         open (18182, file='S_pressure_profil_t.dat')
!CC...ajout nicolas radic
         open (942, FILE='S_temperature_t.dat')
         open (9455, file='S_pressure_flux__last_maille.dat')
         open (944944, file='S_bilan.dat')
         open (9456, file='S_vitesse_profil.dat')
!CC...fin ajout nicolas radic

      case ("ZND") 
         open (1818, file='S_hydraulic_conductivities_profil_t.dat')
         open (18181, file='S_saturation_profil_t.dat')
         open (18182, file='S_pressure_profil_t.dat')

!CC....FICHIERS ZNS 1D
      case ("WAR") 
      if (irp == 1) then
         open (1818, file='S_hydraulic_conductivities_profil_t.dat')
         open (18181, file='S_saturation_profil_t.dat')
         open (18182, file='S_pressure_profil_t.dat')
      end if


!CC....FICHIERS KARINA SCRIPT INVERSION COUPLAGE MAD
      case ("ZHR", "ZHZ")

         if (irpth == 1 .and. irp .ne. 0) then
            open (unit=59, file='Sim_temperature_maille1_t.dat')
            open (unit=60, file='Sim_temperature_maille2_t.dat')
            open (unit=61, file='Sim_temperature_maille3_t.dat')
            open (18181, file='Sim_velocity_profil_t.dat')
            open (18182, file='Sim_heat_flux_profil_t.dat')
            open (18183, file='Sim_temperature_profil_t.dat')
!       open(unit=64,file='Sim_temperature_maille4_t.dat')
!       open(unit=65,file='Sim_temperature_maille5_t.dat')
         end if

!CC....FICHIERS KARINA SCRIPT EMMANUEL LEGER
      case ("1DS") 

         if (irpth == 1 .and. irp .ne. 0) then
            open (7782, file='S_vitesse_nmaille2_hb.dat')
            open (unit=59, file='Sim_temperature_maille1_t.dat')
            open (unit=60, file='Sim_temperature_maille2_t.dat')
            open (unit=61, file='Sim_temperature_maille3_t.dat')
            open (unit=65, file='Sim_temperature_maille5_t.dat')
            open (unit=66, file='Sim_temperature_maille6_t.dat')
            open (unit=67, file='Sim_temperature_maille7_t.dat')
            open (unit=591, file='Sim_temperature_maille8_t.dat')
            open (unit=69, file='Sim_temperature_maille9_t.dat')
            open (unit=591, file='Sim_temperature_maille10_t.dat')
            open (7782, file='S_vitesse_nmaille2_hb.dat')
            open (181822, file='Sim_pressure_profil_t.dat')
            open (1818, file='S_flux_therm_velocity_1_t.dat')
            open (18181, file='Sim_velocity_profil_t.dat')
            open (18182, file='Sim_heat_flux_profil_t.dat')
            open (18183, file='Sim_temperature_profil_t.dat')
            open (unit=64, file='Sim_temperature_maille4_t.dat')
            OPEN (94, FILE='S_permeabilite_t.dat')
            OPEN (941, FILE='S_saturation_profil_t.dat')
            OPEN (942, FILE='S_conduct_therm_t.dat')
            OPEN (943, FILE='S_cap_term_t.dat')

         end if


!CC....FICHIERS KARINA SCRIPT EMMANUEL LEGER
      case("1DJ") 

         if (irpth == 1 .and. irp .ne. 0) then
            open (7782, file='S_vitesse_nmaille2_hb.dat')
            open (181822, file='Sim_pressure_profil_t.dat')
            open (18181, file='Sim_velocity_profil_t.dat')
            OPEN (94, FILE='S_permeabilite_t.dat')
            OPEN (941, FILE='S_saturation_profil_t.dat')

         end if


!CC....SORTIE MAQ GEL/DEGEL
      case ("MAQ") 
         open (19, file='S_soulevement_1_t.dat')
         open (53, file='S_Charges_cap_t.dat')
         open (54, file='S_temp_cap_t.dat')
         open (56, file='S_temp_PT100_t.dat')
         open (18, file='S_bound_permaf_1_t.dat')

!CC...INTERFROST LUNARDINI
      case ("THL") 
         open (743, file='S_pression_charge_temperature_day_3.dat')
         open (742, file='S_pression_charge_temperature_day_2.dat')
         open (741, file='S_pression_charge_temperature_day_1.dat')
!CC...INTERFROST TEST NEUMAN
      case ("TH1") 
         open (18, file='S_bound_permaf_1_t.dat')
         open (53, file='S_Temp_id400_t.dat')
         open (1818, file='S_bilan_therm_t.dat')

!CC....INTERFROST TEST TH3
      case ("TH3") 
         OPEN (181818, FILE='S_TH3.dat', form='unformatted')
         OPEN (181819, FILE='S_TH31.dat')
         OPEN (181820, FILE='S_TH3nm.dat')


!CC....INTERFROST TEST TH2
      case ("TH2") 
         OPEN (181819, FILE='S_TH21.dat')
         OPEN (181820, FILE='S_sortie_TH2nm.dat')
         OPEN (181818, FILE='S_TH2.dat', form='unformatted')


      case ("VIT") 
         open (unit=59, file='Sim_temperature_maille1_t.dat')
         open (unit=60, file='Sim_temperature_maille2_t.dat')
         open (unit=61, file='Sim_temperature_maille3_t.dat')


!CC...AVAV Marine Dangead
      case ("AVA") 
         if (ith == 1) then
            OPEN (181819, FILE='S_temp_zns1.dat')
            OPEN (181820, FILE='S_temp_zns2.dat')
            OPEN (181818, FILE='S_temp_zns3.dat')
            OPEN (181821, FILE='S_temp_TRes1.dat')
            OPEN (181822, FILE='S_temp_TRes2.dat')
            OPEN (181823, FILE='S_temp_HoboRD.dat')
            OPEN (181824, FILE='S_temp_HoboRG.dat')
         end if
         OPEN (94, FILE='S_permeabilite_t.dat')
         OPEN (941, FILE='S_saturation_t.dat')
         OPEN (942, FILE='S_temperature_t.dat')
         OPEN (943, FILE='S_pression_t.dat')                     
         OPEN (181825, FILE='S_piezoB_RD.dat')
         OPEN (181826, FILE='S_piezoB_RG.dat')
         OPEN (181827, FILE='S_surf_piezo_avril.dat')
         OPEN (181831, FILE='S_surf_piezo_juin_bf_rain.dat')
         OPEN (181828, FILE='S_surf_piezo_juin.dat')
         OPEN (181829, FILE='S_surf_piezo_aout.dat')
         OPEN (62, file='S_surf_piezo_finale.dat')
         open (751, file='S_charge_initiale_10days.dat')
         open (181830, FILE='S_flux_ZH_aq.dat')

      case ("R2D") 
         if (ith == 1) then
         open (56, file='S_temp_PT100_t.dat')
         end if
!         OPEN (942, FILE='S_temperature_t.dat')
!         OPEN (943, FILE='S_pression_t.dat')                     
         open (181830, FILE='S_flux_ZH_aq.dat')

!CC...DTS
      case ("DTS") 
!       OPEN(181819,FILE='S_velocity.dat')
!       OPEN(181820,FILE='S_temperature.dat')
!       OPEN(181821,FILE='S_pressure.dat')
         open (751, file='S_charge_initiale_10days.dat')
         OPEN (181822, FILE='S_MolTempP40.dat')
         OPEN (181823, FILE='S_MolTempP42.dat')
         OPEN (181840, FILE='S_MolTempP41.dat')
         OPEN (181841, FILE='S_MolTempP43.dat')
         OPEN (181842, FILE='S_MolTempP44.dat')
         OPEN (181843, FILE='S_MolTempP45.dat')
!     OPEN(181845,FILE='S_MolFluxpP41.dat')
!     OPEN(181846,FILE='S_MolFluxP43.dat')
!     OPEN(181847,FILE='S_MolFluxP44.dat')
!     OPEN(181848,FILE='S_MolFluxP45.dat')
!     OPEN(329,FILE='S_Hd_P_month1.dat')
!     OPEN(330,FILE='S_Hd_P_month2.dat')
!     OPEN(331,FILE='S_Hd_P_month3.dat')
      end select

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                          C
!          BOUCLE TEMPS               C
!                          C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC....COMPTEUR

   icpt = 0
   it = 0
!CC....indice regime permanant thermique
   irptha = irpth
!CC....Boucle temps jusqu a fin de simulation
   do while (nitt*unitsim - paso > 0)
!         print*,dt,paso
      it = it + 1

!CC....Compteur iteration calcul PICARD
     nk = 0
!CC...Retour pas de temps initial impose par l utilisateur
      dt = dble(dta)
      dtreco = dble(dtrecord)
      !print*, "on recommence dans la boucle while avec dt =", int(dt),"paso=",int(paso), "et dtreco =", dtreco!, "et irecord =", irecord
!CC....ENREGISTREMENT A DES PAS DE TEMPS CONSTANT
!ccc....Changement du pas de temps pour enregistrer au pas de temps constant itsortie*unitsortie
!ccc....dt record = temps ecoulé depuis le dernier enregistrement
!ccc....irecord = booléen si vrai ecriture sinon rien

      if (dtreco + dt < itsortie*unitsortie) then
         dtrecord = dble(dtreco) + dble(dt)
         irecord = 0

      elseif (dtreco + dt == itsortie*unitsortie) then
         irecord = 1
         dtrecord = 0
         dtreco = 0
      elseif (dtreco + dt > itsortie*unitsortie) then
         dt = dble(itsortie)*dble(unitsortie) - dble(dtreco)
         irecord = 1
         dtrecord = 0
         dtreco = 0
         !print *, '1 On est passé par la donc dt = ',dt
      end if
      !print*, "ici dtreco = ", dtreco, "et dtrecord =", dtrecord

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!        Indice CONVERCENCE           C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC....ECOULEMENT
      amaxp = 0.D+00
      if (iec == 1) amaxp = 1D+05
      if (iec == 0) amaxp = crconvp
!CC....TRANSPORT
      amaxc = 0.D+00
      if (itr == 1) amaxc = 1D+05
      if (itr == 0) amaxc = crconvc
!CC....THERMIQUE
      amaxt = 0.D+00
      if (ith == 1) amaxt = 1D+05
      if (ith == 0) amaxt = crconvt

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!      Compteur temps simulation          C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      if (irp == 0 .and. irptha == 0 .and. it > 3) paso = nitt*unitsim
      if (irptha == 1 .or. irp == 1) paso = dble(dt) + dble(paso)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!       Variables et cdt lim old          C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      do i = 1, nm
         if (pr(i) + 1 .ne. pr(i)) pro(i) = pr(i) ! Dans le cas où la pression n'a pas convergé on reprend la pression du temps précédent
         rhold(i) = rho(i)
      end do

      if (itr == 1) then
         do i = 1, nm
            conco(i) = conc(i)
         end do
      end if

      if (ith == 1) then
         do i = 1, nm

            if (temp(i) + 1 .ne. temp(i)) tempoo(i) = tempo(i)
            if (temp(i) + 1 .ne. temp(i)) tempo(i) = temp(i)
            do k = 1, 4
               valclto(i, K) = valclt(i, k)
            end do
         end do
         valclto1 = valclt(1, 3)
      end if

      if (icycle == 1) then
      do i = 1, nm
         if (temp(i) + 1 .ne. temp(i)) siceo(i) = sice(i)
      end do
      do kcol = 1, nc
      do iiso = 1, 2
         zloo(kcol, iiso) = zlo(kcol, iiso)
         zlo(kcol, iiso) = zl(kcol, iiso)
      end do
      defo(kcol) = def(kcol)
      end do
      end if

      do kcol = 1, nc
      do iiso = 1, 2
         zsoo(kcol, iiso) = zso(kcol, iiso)
         zso(kcol, iiso) = zs(kcol, iiso)
      end do
      end do

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                          C
!  Variation des conditions aux limites vs. temps C
!                          C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      if (iclchgt == 1) then
         select case (ytest)
         case ("WAR")
            ntsortie = ligne4
         case ("ZNS")
            ntsortie = ligne4
         case ("ZND")
            ntsortie = ligne4
         case ("ZHR")
            ntsortie = ligne4
         case ("1DS")
            ntsortie = ligne4
            if (ligne == 0) ligne = 1
            if (ligne1 == 0) ligne1 = 1
            if (ligne2 == 0) ligne2 = 1
            if (ligne3 == 0) ligne3 = 1
            if (ligne4 == 0) ligne4 = 1
            if (ligne5 == 0) ligne5 = 1
            if (ligne6 == 0) ligne6 = 1
         case ("1DJ")
            ntsortie = ligne4
            if (ligne == 0) ligne = 1
            if (ligne1 == 0) ligne1 = 1
            if (ligne2 == 0) ligne2 = 1
            if (ligne3 == 0) ligne3 = 1
            if (ligne4 == 0) ligne4 = 1
            if (ligne5 == 0) ligne5 = 1
            if (ligne6 == 0) ligne6 = 1
         case ("ZHZ")
            ntsortie = ligne4
         case default
            ntsortie = ligne2
         end select
!CC...ajout nicolas radic
         iecriture_pluie = 0
         ! print*, "Avant calcule vitesse maille n 203 :", vzm(203)
         ! print*, "icl(nm,4)=", icl(nm,4)
         if (it .ne. 1) then
            swo = sw(1)
         else 
            swo = 0
         end if
!CC...fin ajout nicolas radic
         if (.not.allocated(tempsurf)) then
            allocate(tempsurf(1))
         end if
         if (.not.allocated(tempbot)) then
            allocate(tempbot(1))
         end if
         if (.not.allocated(chgsurf)) then
            allocate(chgsurf(1))
         end if
         if (.not.allocated(chgbot)) then
            allocate(chgbot(1))
         end if
         if (.not.allocated(chgriver)) then
            allocate(chgriver(1))
         end if
         if (.not.allocated(chgRD)) then
            allocate(chgRD(1))
         end if
         if (.not.allocated(tempRD)) then
            allocate(tempRD(1))
         end if
         if (.not.allocated(tempRG)) then
            allocate(tempRG(1))
         end if
         if (.not.allocated(tempriver)) then
            allocate(tempriver(1))
         end if
         if (.not.allocated(chgRG)) then
            allocate(chgRG(ntsortie))
         end if

         if (.not.allocated(qbot)) then
            allocate(qbot(1))
         end if
         if (.not.allocated(qsurf)) then
            allocate(qsurf(1))
         end if


         call variation_cdt_limites(nm, paso, itlecture, ytest, &
                                    ligne, ligne1, ligne2, ligne3, ligne4, ligne5, ligne6, &
                                    icl, valcl, iclt, valclt, ivois, &
                                    z, g, ntsortie, bm, irptha, &
                                    rho, qpluie, chgriver, &
                                    chgRD, chgRG, tempRD, tempRG, &
                                    tempriver, &
                                    id_RD, id_RG, id_river, id_rivert, tempsurf, &
                                    tempbot, chgsurf, chgbot, &
                                    x, qsurf, qbot,iecriture_pluie,ak(nm),akr(nm),akr(1),ak(1),vzm(nm-1),swo,pro)


      end if
      !print*, "Rentrer dans la boucle de picard"
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                          C
!       BOUCLE DE PICARD               C
!                          C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      do while (amaxp > crconvp .or. pr(1) + 1 == pr(1) &
                .or. amaxt > crconvt .or. temp(1) + 1 == temp(1) &
                .or. temp(1) - temp(1) .ne. 0 &
                .or. pr(1) - pr(1) .ne. 0 .or. nk == iteration - 1)
         nk = nk + 1

!CC...Sauvegarde de literation picard precedante
         do i = 1, nm
            prk(i) = pr(i)
         end do
         if (itr == 1) then
         do i = 1, nm
            conck(i) = conc(i)
         end do
         end if
         if (ith == 1) then
         do i = 1, nm
            tempk(i) = temp(i)
         end do
         end if

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                          C
!         non convergence ou NaN           C
!         pas de temps adaptatif           C
!                          C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            !print*,sw(1),akv(1),akrv(1),pr(1)/rho1/g+z(1),pr(200)/rho1/g+z(200)
         if (nk == iteration - 1 .or. pr(1) + 1 == pr(1) &
             .or. dt <= 1D-5) then
            !EXPLICATION DE POURQUOI IL N Y A PAS DE CONVERGENCE
            ! if ( nk == iteration -1 ) then
            !    print*, 'NON CONVERGENCE CAR NK EGALE ITERATION -1'
            ! else if ( pr(1) + 1 == pr(1) ) then
            !    print*, 'NON CONVERGENCE CAR LA PRESSION A EXPLOSE'
            ! else if (abs(dt) <= 1D-5) then
            !    print*, 'NON CONVERGENCE CAR DT TROP PETIT'
            !    print*, 'paso = ',paso ,"dt =", dt, "valcl(1,3) =", valcl(1,3)
            !    stop
            ! end if
            if (dtreco >= 0) then
               dtrecord = dble(dtreco) - dble(dt)
               irecord = 0
            end if

            if (dtreco <= 0) then
               dtrecord = 0
               irecord = 0
            end if

            dto = dble(dt)
            dt = dble(dto)/2
            compteur_div=compteur_div+1
            if (dt > dta) dt = dta

            
            !print *, 'ici on a dt ==========',dt,'et dta =============',dta
            !print*,"dans la boucle picard ===================================", modulo(dta, dt) .ne. 0
            ! print*, 'et compteur div =', compteur_div
            if (dt>=1) then
               if (modulo(dta, dt) .ne. 0) dt = dble(dta)/10
            else 
               if (modulo(dta/10**CEILING(log10(dta)),dt) > 1.0d-15) then
                  
                  dt = dble(dta)/10 ! Ajout de ce if car sinon ginette tourne dans le vide quand dt<1
               end if
                  ! print *, "On est dans le cas ou dt<1 avec dt =", dt
            end if
            it = int(it - 1)
!CC....ENREGISTREMENT A DES PAS DE TEMPS CONSTANT
!ccc....Changement du pas de temps pour enregistrer au pas de temps constant itsortie*unitsortie
!ccc....dt record = temps ecoulé depuis le dernier enregistrement
!ccc....irecord = booléen si vrai ecriture sinon rien
            if (dtreco + dt < itsortie*unitsortie) then
               dtrecord = dble(dtreco) + dble(dt)
               irecord = 0
            elseif (dtreco + dt == itsortie*unitsortie) then
               irecord = 1
               dtrecord = 0
               dtreco = 0
            elseif (dtreco + dt > itsortie*unitsortie) then
               dt = dble(itsortie)*dble(unitsortie) - dble(dtreco)
               irecord = 1
               dtrecord = 0
               dtreco = 0
            end if
            !print *, "On change de dt passant de", dto, dt
            !print*, "ici dtreco = ", dtreco, "et dtrecord =", dtrecord
            it = it + 1
            nk = 1

            paso = dble(dt) + dble(paso) - dble(dto)

            if (irp == 0 .and. irptha == 0 .and. it > 3) paso = nitt*unitsim

            if (iclchgt == 1) then
               select case (ytest)
               case ("WAR")
                  ntsortie = ligne4
               case ("ZNS")
                  ntsortie = ligne4
               case ("ZND")
                  ntsortie = ligne4
               case ("ZHR")
                  ntsortie = ligne4
               case ("1DS")
                  ntsortie = ligne4
               case ("1DJ")
                  ntsortie = ligne4
               case ("ZHZ")
                  ntsortie = ligne4
               case default
                  ntsortie = ligne2
               end select

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                          C
!  Variation des conditions aux limites vs. temps C
!                          C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!c aaaaaaaaaaaaaaaaaaaaaaaaaaa
               if (ligne == 0) ligne = 1
               if (ligne1 == 0) ligne1 = 1
               if (ligne2 == 0) ligne2 = 1
               if (ligne3 == 0) ligne3 = 1
               if (ligne4 == 0) ligne4 = 1
               if (ligne5 == 0) ligne5 = 1
               if (ligne6 == 0) ligne6 = 1

!ccc....recalcul car chgt pas temps
               call variation_cdt_limites(nm, paso, itlecture, ytest, &
                                    ligne, ligne1, ligne2, ligne3, ligne4, ligne5, ligne6, &
                                    icl, valcl, iclt, valclt, ivois, &
                                    z, g, ntsortie, bm, irptha, &
                                    rho, qpluie, chgriver, &
                                    chgRD, chgRG, tempRD, tempRG, &
                                    tempriver, &
                                    id_RD, id_RG, id_river, id_rivert, tempsurf, &
                                    tempbot, chgsurf, chgbot, &
                                    x,  &
                                    qsurf, qbot,iecriture_pluie,ak(nm),akr(nm),akr(1),ak(1),vzm(nm-1),swo,pro)

            end if

            do i = 1, nm
               pr(i) = pro(i)
               if (ith == 1) tempo(i) = tempoo(i)
               if (icycle == 1) siceo(i) = siceoo(i)
               if (ith == 1) temp(i) = tempo(i)
               if (icycle == 1) sice(i) = siceo(i)
            end do
         end if

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                          C
!         PAS DE CONVERGENCE           C
!                          C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC....arret du calcul
!CC....impression du dernier pas de temps

         if (dt < 1D-5 .and. nk == iteration) then
         if (iec == 1) then
         if (amaxp > crconvp .or. pr(1) + 1 == pr(1) .or. &
             pr(1) .ne. pr(1) .or. abs(pr(1)) > 1D+20) then
            print *, 'non convergence Picard apres', iteration, 'iterations'
            print *, 'pas de temps =', dt
            print *, 'ecoulement non convergence'
            print *, 'amaxp =', amaxp, 'crit=', crconvp
            print *, 'maille a probeme = ', ipb
            print *, 'pr(ipb)+1 =', pr(ipb) + 1
            print *, 'pr(ipb) = ', pr(ipb)
            print *, 'pr(ipb) nk-1 = ', prk(ipb)
            print *, 'charge(ipb) = ', pr(ipb)/(rho(ipb)*g) + z(ipb)
            do i = 1, nm
               write (74, *) i, x(i), z(i), pro(i), pro(i)/(rho1*g) + z(i), tempo(i)
               write (778, *) i, x(i), z(i), vxp(i), vxm(i), vzp(i), vzm(i)
            end do
            stop
         end if
         end if

         if (ith == 1) then
         if (amaxt > crconvt .or. temp(1) + 1 == temp(1) .or. &
             temp(1) .ne. temp(1) .or. abs(temp(1)) > 1D+16 .or. amazcfl > 1 &
             .or. amaxcfl > 1) then
            print *, 'non convergence Picard apres', iteration, 'iterations'
            print *, 'pas de temps =', dt
            print *, 'thermique non convergence'
            Print *, 'amaxt =', amaxt, 'crit=', crconvt
            if (amaxcfl > 1) print *, 'nombre courant horizontal', amaxcfl
            if (amazcfl > 1) print *, 'nombre courant vertical', amazcfl
            print *, 'temp(1)+1 =', temp(ipb) + 1
            print *, 'temp(1) = ', temp(ipb)
            print *, 'amaxcfl = ', amaxcfl
            print *, 'amazcfl = ', amazcfl
            do i = 1, nm
               write (74, *) i, x(i), z(i), pro(i), pro(i)/(rho1*g) + z(i), tempo(i)
               write (778, *) i, x(i), z(i), vxp(i), vxm(i), vzp(i), vzm(i)
            end do
            stop
         end if
         end if

!      if(rtrho == 1) then
!       print*,'non convergence Picard apres', iteration, 'iterations'
!       print*,'pas de temps =',dt
!       print*, 'rtrho',rtrho
!       do i=1,nm
!       write(74,*) i,x(i),z(i),pro(i),pro(i)/(rho1*g)+z(i),tempo(i)
!       enddo
!       stop
! c     endif
!CCCCCCCCC fin  non convergence
         goto 15
         end if

!CCCCCCCCCCCCCCC

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!   subroutine a netoyer          C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC Topo en marche d escalier
!     if (itopomch == 1) then
!     do i=1,nm
!     if(x(i) < xberg-am(i).and.ivois(i,1) == -99
!     &.and.ivois(i,3) == -99) then
!     valcl(i,1)=qre/2
!       valcl(i,3)=qre/2
!     icl(i,1)=-1
!     icl(i,3)=-1
!     valclt(i,1)=10D+00
!     valclt(i,3)=10+00
!     iclt(i,1)=-2
!     iclt(i,3)=-2
!     endif
!     enddo
!     endif
!        qruis=0.
!      if (iinfil == 1) then
!       call infiltration(pr,nm,n,icol,nc,om,ss,bm,rho,pore,
!     &am,valcl,icl,qtpore,ruis,ivois,z,g,dt,xberg,x,qre,qruis)
!       endif
!       if(iriv == 1) then
!CCriviere
!       n2=n
!       qinf=0
!       qriv=0
!       call debit_riviere(hriv,qriva,qriv,nm,z,hbot,
!     &bm,ivois,qinf,am,xberg,x,al,n,
!     &rug,pent,iqriv,vxp,vzp,qruis)
!       endif
!       if(iriv == 0) goto 79
! RIVIERE
!       call cdt_riviere(hriv,g,nm,z,hbot,
!     &bm,ivois,am,rho,xberg,x,al,n,icl,valcl,iclt,valclt,
!     &tempriv,aklit,aklitv,ak,akv,akr,akrv,dswdp,sw,elit,
!     &akc,akcv,it,ita,tempo,ts)
!79     continue
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!         TERME dswdp  temporel         C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!         ZNS                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

         if (ivg == 1 .or. yunconfined == "UNS") then
            Call unsaturated(pr, swp, dswpdp, swresz, ans, akr, &
                             nm, akrv, rho1, g, ansun, asun)
            if (igelzns == 0 .and. ith == 1) then
            if (ymoycondtherm == "WOODS") then
               alanda(i) = DBLE(sqrt(alandae)*om(i)*sw(i) + &
                                sqrt(alandas(i))*(1D+00 - om(i)) + &
                                sqrt(alandag)*om(i)*(1D+00 - sw(i)))**2
            else if (ymoycondtherm == "GEOME") then
               alanda(i) = DBLE(alandae**(om(i)*sw(i))*alandas(i)**(1 - om(i))* &
                                alandag**(om(i)*(1D+00 - sw(i))))
            else if (ymoycondtherm == "ARITH") then
               alanda(i) = DBLE(alandae*(om(i)*sw(i)) + alandas(i)*(1 - om(i)) &
                                + alandag*(om(i)*(1D+00 - sw(i))))
            end if

            end if
         else if (yunconfined == "CAP" .and. ivg .ne. 1) then
            do i = 1, nm
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          nappe captive             C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!ccc....COEFFICIENT D EMMAGASINEMENT
               swp(i) = 1.D00
               dswpdp(i) = ss(i)/(rho1*g*om(i))
               if (ytest == "TH2" .or. ytest == "TH3") then
                  dswpdp(i) = sw(i)*ss(i)
               end if
            end do

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!         nappe libre           C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         else if (yunconfined == "UNC" .and. ivg .ne. 1) then
            do i = 1, nm
               akr(i) = 1D+00
               akrv(i) = 1D+00
               swp(i) = 1D+00
               dswpdp(i) = 1D+00/(rho1*g*bm(i))
!     dswpdp(i)=1D+00/(rho1*g)
               if (pr(i) < 0d+00) then
!     akr(i)=0D+00
                  akrv(i) = 1D+00
                  swp(i) = 0D+00
               end if
               if (pr(i) >= 0d+00) then
                  akr(i) = 1D+00
                  akrv(i) = 1D+00
                  swp(i) = 1D+00
               end if
            end do
         end if

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!         PARAMETRE VS. TEMPERATURE       C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC....Indice gel/degel
!ccc....test cycle gel valeur de igel!!!!! Attention test seulement sur la maille 1 à améliorer
         if (icycle == 0) igel = 0

         if (icycle == 1) then
         if (ytest .ne. "TH2" .or. ytest .ne. "TH3" .and. &
             ytest .ne. "MAQ" .and. ytest .ne. "TH1") then
            if (temp(1) < 0 .and. tempo(1) - temp(1) > 0) igel = 1
            if (temp(1) < 0 .and. tempo(1) - temp(1) < 0) igel = 2
         end if
         if (ytest == "TH1" .or. ytest == "TH2" .or. ytest == "TH3") igel = 2
!CC... SPECIAL CASES

!CC.... TESTS Interfrost
         if (ytest == "TH2" .or. ytest == "TH3") then
         do i = 1, nm
            dswpdp(i) = ss(i)*sw(i)
         end do
         end if

!       if(ytest == "MAQ") then
!     if (paso <= 1180*3600) igel=1
!     if (paso > 1180*3600) then
!       igel=2
!     cimp=0006d-01
!     endif
         if (ytest == "MAQ") then
            if (paso <= 49.1*86400) igel = 1
            if (paso > 49.1*86400) then
               cimp = 3D-2
               igel = 2
            end if

!CC....Changement coef emmagasinement
            do i = 1, nm
               if (temp(i) < tl) dswpdp(i) = ss(i)/g/rho(i)/om(i)
               if (temp(i) > tl) dswpdp(i) = 1
               do kcol = 1, nc
!CC...GEL....dswpd subpermafrost
                  if (igel == 1 .and. zs(kcol, 1) .ne. -99 .and. icol(i) == kcol .and. &
                      z(i) < zs(kcol, 1) .and. pro(i)/(rho(i)*g) + z(i) > zs(kcol, 1)) then
                     dswpdp(i) = ss(i)/g/rho(i)/om(i)
                  end if
!CC...GEL....dswpd subpermafrost
                  if (igel == 1 .and. zs(kcol, 1) .ne. -99 .and. icol(i) == kcol .and. &
                      z(i) < zs(kcol, 1) .and. pro(i)/(rho(i)*g) + z(i) > topo(kcol)) then
                     dswdp(i) = ss(i)/g/rho(i)/om(i)
                  end if
!CC....DEGEL....dswpd suprapermafrost
                  if (igel == 2 .and. zs(kcol, 1) .ne. -99 .and. icol(i) == kcol .and. &
                      z(i) > zs(kcol, 1) .and. temp(i) > tl) then
                     dswpdp(i) = 1
                  end if
!CC....DEGEL....dswpd subpermafrost
                  if (igel == 2 .and. zs(kcol, 2) .ne. -99 .and. icol(i) == kcol .and. &
                      z(i) < zs(kcol, 2)) then
                     dswpdp(i) = ss(i)/(rho1*g*om(i))
                  end if
                  if (igel == 2 .and. zs(kcol, 1) .ne. -99 .and. icol(i) == kcol .and. &
                      z(i) < zs(kcol, 1) .and. pro(i)/(rho(i)*g) + z(i) < topo(kcol)) then
                     dswdp(i) = 1
                  end if
                  if (igel == 2 .and. zs(kcol, 1) .ne. -99 .and. icol(i) == kcol .and. &
                      z(i) < zs(kcol, 1) .and. pro(i)/(rho(i)*g) + z(i) < zs(kcol, 1)) then
                     dswpdp(i) = 1
                  end if
                  if (igel == 2 .and. zs(kcol, 1) == -99 .and. icol(i) == kcol) then
                     dswpdp(i) = 1
                  end if
               end do
            end do

         end if

!CC....Calcul saturation permeabilite
!CC....Cas completement sature en eau
!CC....Derivation de la saturation en glace en fonction de la temperature dsice
!CC....sw (liquide) = (1D0 - sice)
!CC....Implicite
         if (ibigridice .ne. 1 .and. iparo == 1) then
            CALL icesatperm(nm, tl, ts, akr, akrv, dk, tempo, &
                            dsipdtemp, igelzns, ytypakrice, &
                            sicep, om, sw, swressi, ytypsice, &
                            cimp, tempoo, omega, siceo)
         else
            do i = 1, nm
               sicep(i) = 0.D0
               dsipdtemp(i) = 0.D0
            end do
         end if

!CC....Explicite
         if (ibigridice .ne. 1 .and. iparo == 0) then
            CALL icesatperm(nm, tl, ts, akr, akrv, dk, temp, &
                            dsipdtemp, igelzns, ytypakrice, &
                            sicep, om, sw, swressi, ytypsice, &
                            cimp, tempo, omega, siceo)
         else
            do i = 1, nm
               sicep(i) = 0.D0
               dsipdtemp(i) = 0.D0
            end do
         end if
!CC....grosse maille saturation en glace et permeabilite relative calculees par rapport a la position des isothermes dans les mailles.
         if (ibigridice == 1) then
            CALL biggridice(nm, tl, ts, akr, akrv, dk, temp, igel, &
                            nc, z, zl, zs, bm, dsipdtemp, valclt, ivois, &
                            sicep, swressi, siceo, tempo)
         end if

!CC.....SWTOT est la quantite total d'eau (liquide+glace)
         do i = 1, nm

            swtot = swp(i)
            sice(i) = MAX(0.D0, swtot - (1.D0 - sicep(i)))
            sw(i) = swtot - sice(i)
            if (dsipdtemp(i) == 0) then
               dswdp(i) = dswpdp(i)
               dswdt(i) = 0.D0
               dsidp(i) = 0.D0
               dsidtemp(i) = 0.D0
            else
               dswdp(i) = ss(i)/(rho1*g*om(i))
!CC....nappe captive
               dswdt(i) = -dsipdtemp(i)
               dsidp(i) = dswpdp(i)
               dsidtemp(i) = dsipdtemp(i)
            end if
!CC....Regime permanant
            if (irp == 0) dswdp(i) = 0D+00
         end do
         else
!cc fin icycle=1
!cc debut icycle=0

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                         c
!         PAS DE GEL                                  c
!                                         c
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         do i = 1, nm
            sw(i) = swp(i)
!      sice(i)=0
            dswdp(i) = dswpdp(i)
            dswdt(i) = 0.D0
!      dsidp(i)=0.D0
!      dsidtemp(i)=0.D0

!CC....Regime permanant
            if (irp == 0) dswdp(i) = 0D+00

         end do

         end if

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!             ECOULEMENT           C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC....CONSTRUCTION DE LA MATRICE PRESSION
         if (iec == 1) then
            nz = nm
            nmaxz = nmax
            nmaxzz = nmax1


            call matp(val, icol_ind, irow_ptr, x, z, b, am, ivois, &
                      rho, ak, akr, amu, dt, ia2, g, icl, valcl, rhold, om, pro, dswdp, &
                      sw, nz, irp, nmax1, nmax, bm, akv, akrv, &
                      igel, ysupdp, rhoi, ixy, dsidtemp, temp, tempo)
!CC....Resolution
            n1 = nm
            nmaxz = nmax
            nmaxzz = nmax
         select case(ysolv)
            case("BIC")
               call bicg(pr, b, n1, k, val, icol_ind, irow_ptr, nmaxz, nmaxzz)
            case("CGS")
               call cgs(pr, b, n1, k, val, icol_ind, irow_ptr, nmaxz, nmaxzz)
         end select

!     if (ysolv == "LIB") then
!       do i=1,nm+1
!       irow_ptr(i)=irow_ptr(i)-1
!       enddo
!     call GC_init_sys(irow_ptr,icol_ind,n1)
!           call GC_solve (val,irow_ptr(nmax1),b,n1,pr,pr,sw_int,sw_reel)
!       k=0
!     endif

!CC....calcul de critere amaxp
!CC....Test du Picard
!     if(ivg == 1.or.iriv == 1.or.iqriv == 0.or.
!     &icycle == 0.or.igelzns == 0.or.iomdegel == 0.or.irp == 0) then
            amaxp = 0D+00
            do ii = 1, nm
               if (abs(prk(ii) - pr(ii)) >= amaxp) amaxp = abs(prk(ii) - pr(ii))
               if (abs(prk(ii) - pr(ii)) >= amaxp) ipb = ii


      
      !   endif
            end do
            if (k > 999) amaxp = 1000
            if (k > 999) ipb = -9999

!     else
!     amaxp=0D00
!     endif

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!         Calcul de vitesse           C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            do i = 1, nm
!CC....VITESSE EN X A TRAVERS LES FACES
!ccc....vxm(i) face gauche
!ccc....vxp(i) face droite
!ccc....vzm(i) face du dessous
!ccc....vzp(i) face du dessus
               if (iec == 1 .and. idecouplage == 0) then
!ccc....Horizontal INITIALISATION
                  vxm(i) = 0.D+00
                  vxp(i) = 0.D+00
!ccc....cacul vitesse face DROITE
                  if (icl(i, 1) == 1) then
                     rhom = dble((rho(i) + rho(ivois(i, 1)))/2D+00)
                     iv = ivois(i, 1)
                     vxp(i) = dble((am(i) + am(iv))/((am(i)/akr(i)/ak(i)) + &
                                                     (am(iv)/akr(iv)/ak(iv))))
                     vxp(i) = dble(-vxp(i)/amu*((pr(i) - pr(iv))/(x(i) - x(iv)) &
                                                + ixy*g*rhom*(z(iv) - z(i))/(x(iv) - x(i))))
                  end if
!ccc....cacul vitesse face GAUCHE
                  if (icl(i, 2) == 1) then
                     rhom = dble((rho(i) + rho(ivois(i, 2)))/2D+00)
                     iv = ivois(i, 2)
                     vxm(i) = dble((am(i) + am(iv))/((am(i)/akr(i)/ak(i)) + &
                                                     (am(iv)/akr(iv)/ak(iv))))
                     vxm(i) = dble(-vxm(i)/amu*((pr(iv) - pr(i))/(x(iv) - x(i)) &
                                                + ixy*g*rhom*(z(iv) - z(i))/(x(iv) - x(i))))
                  end if
!ccc....vertical INITIALISATION
                  vzm(i) = 0.D+00
                  vzp(i) = 0.D+00
!ccc....cacul vitesse face BAS
                  if (icl(i, 4) == 1) then
                     rhom = dble((rho(i) + rho(ivois(i, 4)))/2D+00)
                     iv = ivois(i, 4)
                     vzm(i) = dble((bm(i) + bm(iv))/((bm(i)/akrv(i)/akv(i)) + &
                                                     (bm(iv)/akrv(iv)/akv(iv))))
                     vzm(i) = dble((-vzm(i)/amu*((pr(i) - pr(iv))/(z(i) - z(iv)) &
                                                 + ixy*g*rhom)))
                  end if
!ccc....cacul vitesse face HAUT
                  if (icl(i, 3) == 1) then
                     rhom = dble((rho(i) + rho(ivois(i, 3)))/2D+00)
                     iv = ivois(i, 3)
                     vzp(i) = dble(((bm(i) + bm(iv))/((bm(i)/akrv(i)/akv(i)) + &
                                                      (bm(iv)/akrv(iv)/akv(iv)))))
                     vzp(i) = dble((-vzp(i)/amu*((pr(iv) - pr(i))/(z(iv) - z(i)) &
                                                 + ixy*g*rhom)))
                  end if
!ccc....cacul vitesse avec CHARGE IMPOSEE
                  if (icl(i, 1) == -2) then
                     vxp(i) = dble(vxm(i))
                  end if
                  if (icl(i, 2) == -2) then
                     vxm(i) = dble(vxp(i))
                  end if
                  if (icl(i, 3) == -2) then
                     vzp(i) = dble(vzm(i))
                  end if
                  if (icl(i, 4) == -2) then
                     vzm(i) = dble(vzp(i))
                  end if
!ccc....cacul vitesse avec FLUX IMPOSE
                  if (icl(i, 1) == -1) then
                     vxp(i) = dble(-valcl(i, 1))
                  end if
                  if (icl(i, 2) == -1) then
                     vxm(i) = dble(valcl(i, 2))
                  end if
                  if (icl(i, 3) == -1) then
                     vzp(i) = dble(-valcl(i, 3))
                     vzm(i) = vzp(i)
                  end if
                  if (icl(i, 4) == -1) then
                     vzm(i) = dble(valcl(i, 4))
                     vzp(i) = vzm(i)
                  end if
               end if

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!         DECOUPLAGE               C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC....cacul vitesse face  decouplage code  =1
!CC....utilisation des pression n-1 pour le terme adveectif
               if (iec == 1 .and. idecouplage == 1) then
!CC....INITIALISATION
                  vxm(i) = 0.D+00
                  vxp(i) = 0.D+00
!ccc....cacul vitesse face DROITE
                  if (icl(i, 1) == 1) then
                     rhom = (rho(i) + rho(ivois(i, 1)))/2
                     iv = ivois(i, 1)
                     vxp(i) = (am(i) + am(iv))/((am(i)/akr(i)/ak(i)) + &
                                                (am(iv)/akr(iv)/ak(iv)))
                     vxp(i) = -vxp(i)/amu*((pro(i) - pro(iv))/(x(i) - x(iv)) &
                                           + ixy*g*rhom*(z(iv) - z(i))/(x(iv) - x(i)))
                  end if
                  if (icl(i, 1) == -2) then
                     vxp(i) = -ak(i)*akr(i)/amu*(valcl(i, 1) - pro(i))/am(i)*2
                  end if
                  if (icl(i, 1) == -1) vxp(i) = -valcl(i, 1)
!ccc....cacul vitesse face GAUCHE
                  if (icl(i, 2) == 1) then
                     rhom = (rho(i) + rho(ivois(i, 2)))/2
                     iv = ivois(i, 2)
                     vxm(i) = (am(i) + am(iv))/((am(i)/akr(i)/ak(i)) + &
                                                (am(iv)/akr(iv)/ak(iv)))
                     vxm(i) = -vxm(i)/amu*((pro(iv) - pro(i))/(x(iv) - x(i)) &
                                           + ixy*g*rhom*(z(iv) - z(i))/(x(iv) - x(i)))
                  end if
                  if (icl(i, 2) == -2) then
                     vxm(i) = -ak(i)*akr(i)/amu*(pro(i) - valcl(i, 2))/am(i)*2
                  end if
                  if (icl(i, 2) == -1) vxm(i) = valcl(i, 2)
                  vzm(i) = 0.D+00
                  vzp(i) = 0.D+00
!ccc....cacul vitesse face BAS
                  if (icl(i, 4) == 1) then
                     rhom = (rho(i) + rho(ivois(i, 4)))/2
                     iv = ivois(i, 4)
                     vzm(i) = (bm(i) + bm(iv))/((bm(i)/akrv(i)/akv(i)) + &
                                                (bm(iv)/akrv(iv)/akv(iv)))
                     vzm(i) = -vzm(i)/amu*((pro(i) - pro(iv))/(z(i) - z(iv)) &
                                           + ixy*g*rhom)
                  end if
                  if (icl(i, 4) == -2) then
                     vzm(i) = -akv(i)*akrv(i)/amu*((pro(i) - valcl(i, 4))/bm(i)*2 &
                                                   + ixy*g*rho(i))
                  end if
                  if (icl(i, 4) == -1) vzm(i) = valcl(i, 4)
!ccc....cacul vitesse face HAUT
                  if (icl(i, 3) == 1) then
                     rhom = (rho(i) + rho(ivois(i, 3)))/2
                     iv = ivois(i, 3)
                     vzp(i) = (bm(i) + bm(iv))/((bm(i)/akrv(i)/akv(i)) + &
                                                (bm(iv)/akrv(iv)/akv(iv)))
                     vzp(i) = -vzp(i)/amu*((pro(iv) - pro(i))/(z(iv) - z(i)) &
                                           + ixy*g*rhom)
                  end if
                  if (icl(i, 3) == -2) then
                     vzp(i) = -akv(i)*akrv(i)/amu*((valcl(i, 3) - pro(i))/bm(i)*2 &
                                                   + ixy*g*rho(i))
                  end if
                  if (icl(i, 3) == -1) vzp(i) = -valcl(i, 3)
!ccc....cacul vitesse avec CHARGE IMPOSEE
                  if (icl(i, 1) == -2) then
                     vxp(i) = vxm(i)
                  end if
                  if (icl(i, 2) == -2) then
                     vxm(i) = vxp(i)
                  end if
                  if (icl(i, 3) == -2) then
                     vzp(i) = vzm(i)
                  end if
                  if (icl(i, 4) == -2) then
                     vzm(i) = vzp(i)
                  end if
               end if
!CC....VITESSE NULLE
               if (abs(vxp(i)) < 1D-17) vxp(i) = 0D+00
               if (abs(vxm(i)) < 1D-17) vxm(i) = 0D+00
               if (abs(vzm(i)) < 1D-17) vzm(i) = 0D+00
               if (abs(vzp(i)) < 1D-17) vzp(i) = 0D+00
            end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          Fin alcul de vitesse           C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!         fin ecoulement calcule         C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         end if

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                         C
!    Advection sans calcul de l'ecoulement     C
!                         C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC....avec des vitesses donnees par l'utilisateur
!CC.... ou annulation de l'advection C
         if (iec == 0) then
         do i = 1, nm
         if (ithec == 1) then
!ccc....Vitesse face DROITE
            vxp(i) = dble(-valcl_droite)
!ccc....Vitesse face GAUCHE
            vxm(i) = dble(valcl_gauche)
!ccc....Vitesse face BAS
            vzm(i) = dble(valcl_bas)
!ccc....Vitesse face HAUT
            vzp(i) = dble(-valcl_haut)
!CC....VITESSE NULLE
            if (abs(vxp(i)) < 1D-17) vxp(i) = 0D+00
            if (abs(vxm(i)) < 1D-17) vxm(i) = 0D+00
            if (abs(vzm(i)) < 1D-17) vzm(i) = 0D+00
            if (abs(vzp(i)) < 1D-17) vzp(i) = 0D+00
         else
!CC....PAS D'ADVECTION
            vxp(i) = 0.D+00
            vxm(i) = 0.D+00
            vzm(i) = 0.D+00
            vzp(i) = 0.D+00
         end if
         end do
         end if
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!        TRANSPORT             C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC....SI ITR=0, ON ZAPPE LE TRANSPORT!!!!

         if (itr == 1) then
            nz = nm
            dtc = dt
            call matc(val, icol_ind, irow_ptr, x, b, am, ivois, conco, &
                      dtc, iclc, valclc, om, nm, allg, alt, &
                      nmax, nmax1, bm, z, vxp, vxm, vzp, vzm)
            n1 = nm
            nmaxz = nmax
            nmaxzz = nmax1
!CC....resolution
            call cgs(conc, b, n1, k, val, icol_ind, irow_ptr, nmaxz, nmaxzz)
!CC....Test du Picard transport
            amaxc = 0.D+00
            do ii = 1, nm
               rho(ii) = rho1 + 1/38.D+00*25.D+00*conc(ii)
               if (abs(conck(ii) - conc(ii)) >= amaxc) then
                  amaxc = abs(conck(ii) - conc(ii))
               end if
            end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!         fin transport calcule           C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         end if

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!     THERMIQUE                   C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC....SI ITh=0, ON ZAPPE LE Thermique!!!!
         if (ith == 1) then
            nth = nm
            nmaxth = nmax
            nmaxthh = nmax1

            call matt(val, icol_ind, irow_ptr, x, b, am, ivois, tempo, &
                      rho, dt, iclt, valclt, om, nm, bll, blt, &
                      cpe, cps, rhos, alanda, chlat, nmaxth, igel, nmaxthh, &
                      z, bm, temp, ithec, irpth, ts, tl, &
                      alandae, alandas, alandai, rhoi, cpice, sice, sw, &
                      icycle, rhog, alandag, cpg, igelzns, &
                      vxp, vxm, vzp, vzm, &
                      ymoycondtherm, dsidtemp, ytest, ysolv)
!CC....resolution
            n1 = nm
            nmaxz = nmax
            nmaxzz = nmax1
            if (ysolv == "BIC") then
               call bicg(temp, b, n1, k, val, icol_ind, irow_ptr, nmaxz, nmaxzz)
            end if
            if (ysolv == "CGS") then
               call cgs(temp, b, n1, k, val, icol_ind, irow_ptr, nmaxz, nmaxzz)
            end if

!     if (ysolv == "LIB") then
!       do i=1,nm+1
!       irow_ptr(i)=irow_ptr(i)-1
!       enddo
!     do ii=1,10
!     print*,irow_ptr(ii),val(ii),icol_ind(ii)
!     enddo
!     sw_int(12)=1
!     call GC_init_sys(irow_ptr,icol_ind,n1)
! C         call GC_solve (val,irow_ptr(nmax1),b,n1,temp,temp,sw_int,sw_reel)
!       k=0
!     endif

!CC....Test du Picard thermique
            if (iriv == 1 .or. iqriv == 0 .or. &
                icycle == 0 .or. igelzns == 0 .or. iomdegel == 0) then
               amaxt = 0D+00
               do ii = 1, nm
               if (abs(tempk(ii) - temp(ii)) >= amaxt) then
                  amaxt = abs(tempk(ii) - temp(ii))
                  ipb = ii
               end if
               end do
               if (k > 999) amaxt = 1000
               if (k > 999) ipb = -9999
            else
               amaxt = 0D00
            end if

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!     calcul du soulevement           C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            if (icycle == 1) then
               nz = n
               nmz = nm
               ncz = nc
               call upheaval(ncz, nmz, igel, zs, zsoo, ivois, pr, pro, &
                             alph, dl, def, defo, icol)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!     ISOTHERMES                 C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
               call interpol(nc, zs, nm, icol, temp, ivois, igel, &
                             bm, z, ncmax, valclt, ts, iclt, topo)
               call interpol(nc, zl, nm, icol, temp, ivois, igel, &
                             bm, z, ncmax, valclt, tl, iclt, topo)
            end if
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!         fin thermique calcule           C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         end if

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!     nombre de courant               C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         amaxcfl = -1
         amazcfl = -1
         if (ithec == 1) then
         do i = 1, nm
!C      CCCCCC
            if ((abs(vxm(i) + vxp(i))/2*dt/am(i)) > amaxcfl) then
               amaxcfl = (abs(vxm(i) + vxp(i))/2*dt/am(i))
            end if
            if ((abs(vzm(i) + vzp(i))/2*dt/bm(i)) > amazcfl) then
               amazcfl = (abs(vzm(i) + vzp(i))/2*dt/bm(i))

            end if
         end do
         end if

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  fin boucle while (tests convergence    C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!CC....CATTENTION ENDDO A NE SURTOUT PAS VIRER BOUCLE DE PICARD TEST DE CONVERGENCE
15       continue
!cccc debug print
 !        if ( amaxp <= crconvp .and. pr(1) + 1 .ne. pr(1) &
 !        .and. amaxt <= crconvt .and. temp(1) + 1 .ne. temp(1) &
  !       .and. temp(1) - temp(1) == 0 &
  !       .and. pr(1) - pr(1) == 0 .and. nk .ne. iteration - 1 ) then
         ! print*,"IL Y A EU CONVERGENCE et dt =", dt
         !print*, 'et dtrecord = ', dtrecord, "et paso = ",paso
   !   end if
      end do
      if (irecord == 1) then
         compteur_div = 0
         dtrecord_mem = 0
      else 
         dtrecord_mem = dtrecord_mem + dtrecord
         dtrecord = dtrecord_mem
         !dtrecord = dtrecord * compteur_div
      end if
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!              VERIF           C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!   bilan flux ecoulement              C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC....VITESSE EN X A TRAVERS LES FACES
!ccc....vxm(i) face gauche
!ccc....vxp(i) face droite
!ccc....vzm(i) face du dessous
!ccc....vzp(i) face du dessus
!CC....ivois(ik,1)= voisin droite
!CC....ivois(ik,2)= voisin gauche
!CC....ivois(ik,3)= voisin haut
!CC....ivois(ik,4)= voisin ibas
      sumflux = 0D+00
      do i = 1, nm
         sumflux = sumflux - vxp(i)*bm(i) + vxm(i)*bm(i) + vzm(i)*am(i) &
                   - vzp(i)*am(i)
      end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  CDT DE DRAIN                    C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!       if (icldrain == 1.or.iriviere == 1) then
!       call mse(n,nmaille1,nmaille2,nmaille3,nmaille4,nmaille5,
!     &nmaille6,pr,qinf,qinfobs,paso,itlecture,chgobs1,chgobs2,
!     &chgobs3,chgobs4,rho,g,chgobs5,chgobs6,eo1,eo2,eo3,eo4,eo5,
!     &eo6,eoinf,z,seo1,seo2,seo3,seo4,seo5,
!     &seo6,seoinf)
!       if(modulo(paso,itsortie*unitsortie) == 0.and.paso >= itsortie*unitsortie) then
!       write(46,*) paso/unitsortie,ans,asp,seo6,qinf,seoinf
!       endif
!       endif
!CC....TEST COND PARTICULIERES IMPOSEES SUR CELLULE VOISINE DRAIN
!       if (icldrain == 1) then
!          drain=0
!          do i=1,nm
!       if (z(i) > 0.05.and.z(i) <= 0.05+bm(i)
!     &.and.ivois(i,4) == -99) then
!       drain=drain+vzm(i)*am(i)*1
!       endif
!       if (ivois(i,2) == -99
!     &.and.x(i) <= 0.05+am(i).and.
!     &x(i) > 0.05) then
!       drain=drain+vxm(i)*bm(i)*1
!       endif
!       enddo
!CC....somme flux en fonction du temps
!       sumflux=0
!       do i=1,nm
!       sumflux=sumflux-vxp(i)*bm(i)+vxm(i)*bm(i)+vzm(i)*am(i)
!     &-vzp(i)*am(i)
!       enddo
!       endif
!       if (pr(1)+1 == pr(1).or.pr(1).ne.pr(1).or.
!     &temp(1) == temp(1)+1.or.temp(1).ne.temp(1)) then
!       print*,'ca ne marche pas pr ou temp=nan ou infiny '
!       stop
!       endif

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!          TEST PERFORMANCE           C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      if (ytest == "ZHR" .or. ytest == "ZHZ"&
    &.or. ytest == "1DS") then
         do i = 1, nm
            qad(i) = 0D+00
            qcondu(i) = 0D+00
            qtherm(i) = 0D+00
         end do
!CC....flux THERMIQUE
         do i = 1, nm
            if (ymoycondtherm == "WOODS") then
               alanda(i) = DBLE(sqrt(alandae)*om(i)*sw(i) &
                                + sqrt(alandai)*(om(i)*sice(i)) + &
                                sqrt(alandas(i))*(1D+00 - om(i)) + &
                                sqrt(alandag)*om(i)*(1D+00 - sw(i) - sice(i)))**2
            else if (ymoycondtherm == "GEOME") then
               alanda(i) = DBLE(alandae**(om(i)*sw(i))*alandas(i)**(1 - om(i))* &
                                alandai**(om(i)*sice(i))*alandag**(om(i)*(1D+00 - sw(i) - sice(i))))
            else if (ymoycondtherm == "ARITH") then
               alanda(i) = DBLE(alandae*(om(i)*sw(i)) + alandas(i)*(1 - om(i)) &
                                + alandai*(om(i)*sice(i)) + alandag*(om(i)*(1D+00 - sw(i) - sice(i))))
            end if

!CC....FACE INF
!ccc....vzp(i) face du dessous entrante kg/m3 J⋅kg−1⋅K−1 m/s K  > J /s  > W/m2
            if (vzm(i) > 0 .and. iclt(i, 4) == 1) then
               qad(i) = dble(rho(ivois(i, 4))*cpe*vzm(i))*temp(ivois(i, 4))
            end if
!   if(vzm(i) > 0.and.iclt(i,4) == -2) then
!   qad(i)=dble(rho(i)*cpe*vzm(i)*valclt(i,4))
!   endif
            if (vzm(i) < 0) then
               qad(i) = dble(rho(i)*cpe*vzm(i)*temp(i))
            end if
!CC...conduction dispersion
            dtjm = dble(blt*rho(i)*cpe*((vxm(i) + vxp(i))**2./4. &
                                        +(vzm(i) + vzp(i))**2./4.)**0.5 + alanda(i))

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!         FACE BASSE              C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            if (iclt(i, 4) == 1) then
               qcondu(i) = -dble(dtjm/(z(i) - z(ivois(i, 4)))* &
                                 (temp(i) - temp(ivois(i, 4))))

            end if
!CC....temperature imposee W·m−1·K−1 K /m  > W/m2
            if (iclt(i, 4) == -2) then
               qcondu(i) = -dble(dtjm*(temp(i) - valclt(i, 4))/bm(i)*2.D00)
            end if
            qtherm(i) = qcondu(i) + qad(i)
         end do
      end if

      if (ytest == "TH2") then
         PF2 = 0
         qtin = 0D+00
         qtout = 0D+00
         qthermtot = 0D+00
         swtotal = 0D+00
         sicetotal = 0D+00
         tempmin = temp(1924)
!CC....BILAN THERMIQUE
         do i = 1, nm
            if (ivois(i, 2) == -99) then
               qtin = qtin + (rho(i)*cpe*(temp(i) - 5)*vxp(i) - dble(alanda(i)* &
                                                                     (temp(ivois(i, 1)) - temp(i))/(x(ivois(i, 1)) - x(i))))*bm(i)
            end if

            if (ivois(i, 1) == -99) then
               qtout = qtout + (rho(i)*cpe*(temp(i) - 5)*vxm(i) - &
                                dble(alanda(i)*(temp(i) - temp(ivois(i, 2))) &
                                     /(x(i) - x(ivois(i, 2)))))*bm(i)
            end if

            swtotal = dble(sw(i)*om(i)*am(i)*bm(i) + swtotal)
            sicetotal = dble(sice(i)*om(i)*am(i)*bm(i) + sicetotal)
            tempmin = Min(tempmin, temp(i))
         end do
         PF2 = (qtout - qtin)/az

      end if
!CC....Critere performance Interfrost
      if (ytest == "TH3") then
         qtcol = 0D+00
         qeout = 0D+00
         akeq = 0D+00
         qthermtot = 0D+00
         pt1 = 0D+00
         pt2 = 0D+00

!ccc....BILAN THERMIQUE
         do i = 1, nm
            qtherm(i) = (temp(i) + 273.15)*(om(i)*sw(i)*rho(i)*cpe + &
                                            om(i)*sice(i)*rhoi(i)*cpice + &
                                            (1.-om(i))*rhos(i)*cps(i))*bm(i)*am(i)
            qthermtot = qthermtot + qtherm(i)
!ccc....Terme advectif
!ccc....VITESSE EN X A TRAVERS LES FACES
!ccc....vxm(i) face gauche
!ccc....vxp(i) face droite
!ccc....vzm(i) face du dessous
!ccc....vzp(i) face du dessus
            if (ivois(i, 4) == -99) then
               qtcol = qtcol + alanda(i)*(temp(i) - valclt_bas) &
                       /z(i)*am(i)
            end if
            if (ivois(i, 1) == -99) then
               qeout = qeout + vxp(i)*bm(i)
            end if
         end do
         qthermtot = qthermtot*2
!     grad=0.06D+00
         grad = valcl_gauche
         akeq = (qeout/grad)*2
         qtcol = qtcol/al*2
         pt1 = (temp(nmaille1) + temp(nmaille2))/2.
         pt2 = (temp(nmaille5) + temp(nmaille6) + &
                temp(nmaille7) + temp(nmaille8))/4.
      end if
!CC....TEST NEUMAN
      if (ytest == "TH1") then
         gxl = 2D0 - zl(1, 1)
         gxs = 2D0 - zs(1, 1)
      end if
      
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!                                                 C
!                     OUTPUTS                          C
!                                                 C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!       ECRITURE  FICHIERS SORTIES        C
!                                                 C

!CC....SORTIE these texier
         select case(ytest)
         
         case ("TEX") 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!     SURFACE PIEZO  INTERPOLATION LINEAR  C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         if (irp == 1) then
         if (irecord == 1) then
            ts = 0
            call interpolsurf(nc, zs, nm, icol, pr, ivois, &
                              bm, z, valcl, icl, id_river, id_rivert, chgriver, &
                              ligne5, ligne6, ligne2, paso, itlecture, rho1, g)
            do i = 1, nm
               write (7782, *) paso/unitsortie, paso/86400, x(i), &
                  z(i), pr(i)/rho1/g + z(i), sw(i), vxp(i), &
                  vxm(i), vzp(i), vzm(i)
            end do
         end if
         end if


!CC....SORTIE ZNS 1D
      case ("ZNS") 
         print *, "out", paso/86400,"time", paso, "dt", dt,'prnm', pr(nm)/rho1/g + z(nm), &
            valcl(1, 3), sw(1), sw(50), sw(2), akr(1)
         if (modulo(paso, itsortie*1.0) == 0) then
            irecord = 1
         else
            print *, "############################### PAS D ECRITURE DE FICHIER##############################################"
         end if
!.... option debug
         if (irp == 1 .and. irecord == 1) then
                     ! if (modulo(int(paso), itsortie) .ne. 0) then
                     !       print*,"PAS DE PAS REGULIER"
                     ! end if
         ! print *, "writing Output TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
         teneur_bilan_sable = 0 ! A Enlever pour le commit
         do i = 1, nm
            write (1818, *) paso/unitsortie, z(i), akr(i)*ak(i)
            write (18181, *) int(paso/unitsortie), z(i), sw(i)
            write (18182, *) int(paso/unitsortie), z(i), pr(i), pr(i)/rho1/g + z(i)
            write(942, *) int(paso/unitsortie), z(i), temp(i)
            write (9456, *) int(paso/unitsortie), z(i), vzm(i)
            teneur_bilan_sable = teneur_bilan_sable + sw(i)*0.43
         end do
         write (9455, *) int(paso/unitsortie), z(nm), pr(nm), vzm(nm)
         write (944944, *) int(paso/unitsortie), vzp(1), vzm(nm)*900, teneur_bilan_sable ! A Enlever pour le commit
         end if


      case ("ZND") 
         print *, "time", paso, "dt", dt, pr(1)/rho1/g + z(1), &
            valcl(1, 3)/rho1/g + z(1), sw(100), ak(100)
         if (modulo(int(paso), itsortie) == 0) then
            irecord = 1
         end if

         if (irp == 1 .and. irecord == 1) then
         do i = 1, nm
            write (1818, *) paso/unitsortie, z(i), akr(i)*ak(i)
            write (18181, *) paso/unitsortie, z(i), sw(i)
            write (18182, *) paso/unitsortie, z(i), pr(i), pr(i)/rho1/g + z(i)
         end do
         end if


!CC....SORTIE ZNS 1D
      case ("WAR") 
         print *, "time", paso, "dt", dt, vzm(1), pr(1), "K", ak(1)*akr(1)

         if (modulo(int(paso), itsortie) == 0) then
            irecord = 1
         end if

         if (irp == 1 .and. irecord == 1) then
         do i = 1, nm
            write (1818, *) paso/unitsortie, z(i), akr(i)*ak(i)
            write (18181, *) paso/unitsortie, z(i), sw(i)
            write (18182, *) paso/unitsortie, z(i), pr(i), pr(i)/rho1/g + z(i)
         end do
         end if


!CC...SORTIES DTS
      case ("DTS") 
      if (irecord == 1) then

         if (paso < 86400 + 10 .and. paso > 86400 - 10) then
         do i = 1, nm
            write (751, *) pr(i)
         end do
         end if

!     if(paso == 1559500) then
!     do i=1,nm
!        write(329,*)paso/unitsortie,paso/86400,x(i),
!     &z(i),pr(i)/rho1/g+z(i),temp(i),vxp(i),
!     &vxm(i),vzp(i),vzm(i)
!      enddo
!     endif

!     if(paso == 4238000) then
!     do i=1,nm
!        write(330,*)paso/unitsortie,paso/86400,x(i),
!     &z(i),pr(i)/rho1/g+z(i),temp(i),vxp(i),
!     &vxm(i),vzp(i),vzm(i)
!      enddo
!     endif

!     if(paso == 6830000) then
!     do i=1,nm
!        write(331,*)paso/unitsortie,paso/86400,x(i),
!     &z(i),pr(i)/rho1/g+z(i),temp(i),vxp(i),
!     &vxm(i),vzp(i),vzm(i)
!      enddo
!     endif

!CCC MOLONARI PRINT start time record
!ccc.... to add an end time paso < XXXX)

!ccc....molonari40
!   print*, izone(), x(55), z(55)
!        stop

!ccc....molonari40
!CCC            print*,z(60469),izone(60469)
!CCC            stop

         if (paso > 86400*2) then
!cc...Temperature for the surface and each depth of the molonari
!ccc....molonari40
!       write(181822,*)paso/unitsortie,(temp(60445)+temp(60446))/2,
!     &(temp(23960)+temp(23961))/2,(temp(84873)+temp(84874))/2,
!     &(temp(91895)+temp(91896))/2
!ccc....molonari42
            write (181823, *) paso/unitsortie, &
               (temp(98524) + temp(98525))/2, (temp(98505) + temp(98506))/2, &
               (temp(98483) + temp(98484))/2, (temp(98445) + temp(98446))/2
!ccc....molonari41
            write (181840, *) paso/unitsortie, &
               !     &(valclt(91896,3)+valclt(92463,3))/2,&
               (temp(91297) + temp(91298))/2, (temp(90670) + temp(90671))/2, &
               (temp(90031) + temp(90032))/2, (temp(89379) + temp(89380))/2
!ccc....molonari43
            write (181841, *) paso/unitsortie, &
               !     &(valclt(70796,3)+valclt(69588,3))/2,&
               (temp(69588) + temp(68589))/2, (temp(68372) + temp(67128))/2, &
               (temp(67128) + temp(67129))/2, (temp(65866) + temp(65867))/2
!ccc....molonari44
            write (181842, *) paso/unitsortie, &
               !     &(valclt(65723,3)+valclt(65722,3))/2,&
               (temp(64445) + temp(64446))/2, (temp(63145) + temp(63146))/2, &
               (temp(61830) + temp(61831))/2, (temp(60469) + temp(60470))/2
!ccc....molonari45
            write (181843, *) paso/unitsortie, &
               !     &(valclt(33221,3)+valclt(33220,3))/2,&
               (temp(23960) + temp(23961))/2, (temp(22055) + temp(22056))/2, &
               (temp(20118) + temp(20119))/2, (temp(18148) + temp(18149))/2
         end if
      end if

!ccccccc Find Location of Molonari
!         do i=1,nm
!         if(x(i) >= 773.5.and.x(i) <= 774.5) then
!          print*,i,x(i),z(i),ivois(i,3)
!         endif
!         if(x(i) >= 436.5.and.x(i) <= 437.5) then
!c        print*,i,x(i),z(i),ivois(i,3)
!         endif
!         if(x(i) >= 364.5.and.x(i) <= 365.5) then
!         print*,i,x(i),z(i),ivois(i,3)
!         endif
!         if(x(i) >= 74.5.and.x(i) <= 75.5) then
!         print*,i,x(i),z(i),ivois(i,3)
!         endif
!         enddo
      print *, "cell ", "time ", "hh", "temp "
      print *, "cell 1", paso/86400, pr(98445)/rho1/g + z(98445), temp(98445)
      print *, "cell 1", paso/86400, pr(98525)/rho1/g + z(98525), temp(98525)
      print *, "cell 1", izone(98524)

!CC...SORTIES VAUCLIN parametre clement
!       print*,paso/86400,dtCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      case( "VAU") 
         ts = 0
         call interpolsurf(nc, zs, nm, icol, pr, ivois, &
                           bm, z, valcl, icl, id_river, id_rivert, chgriver, &
                           ligne5, ligne6, ligne2, paso, itlecture, rho1, g)
         print *, "out", paso, zs(1, 1), zs(nc, 1)
         if (irecord == 1) then

            do i = 1, nm
               write (91, *) paso/us, i, x(i), z(i), pr(i), pr(i)/rho(i)/g + z(i), &
                  sw(i)
            end do
            OPEN (181819, FILE='S_2H.dat')
            OPEN (181821, FILE='S_3H.dat')
            OPEN (181820, FILE='S_4H.dat')
            OPEN (181818, FILE='S_8H.dat')
            print *, zs(1, 1)
            write (18, *) paso/unitsortie, zs(1, 1), zs(1, 2)
            do kkcol = 1, nc
               write (622, *) paso/unitsortie, kkcol, x(ibas(kkcol)), &
                  zs(kkcol, 1), zs(kkcol, 2)
            end do

            if (paso < 2*3600 + 10 .and. paso > 2*3600 - 10) then
               do kkcol = 1, nc
                  write (181819, *) x(ibas(kkcol)), zs(kkcol, 1)
               end do
            end if

            if (paso < 3*3600 + 10 .and. paso > 3*3600 - 10) then
               do kkcol = 1, nc
                  write (181821, *) x(ibas(kkcol)), zs(kkcol, 1)
               end do
            end if

            if (paso < 4*3600 + 10 .and. paso > 4*3600 - 10) then
               do kkcol = 1, nc
                  write (181820, *) x(ibas(kkcol)), zs(kkcol, 1)
               end do
            end if

            if (paso < 8*3600 + 10 .and. paso > 8*3600 - 10) then
               do kkcol = 1, nc
                  write (181818, *) x(ibas(kkcol)), zs(kkcol, 1)
               end do
            end if

         end if


!CC...SORTIES AVAV
!CCC     id_ZH(j)

      case ("AVA")
      
      if (irecord == 1) then
         ts = 0
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!     SURFACE PIEZO               C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         call interpolsurf(nc, zs, nm, icol, pr, ivois, &
                           bm, z, valcl, icl, id_river, id_rivert, chgriver, &
                           ligne5, ligne6, ligne2, paso, itlecture, rho1, g)
         qflux = 0D+00
         do j = 1, ligne7
            qflux = vzm(id_ZH(j))*am(id_ZH(j)) + qflux
         end do
         write (181830, *) paso/unitsortie, qflux
         print *, "out", paso/86400, paso,dt,am(id_ZH(1)), qflux,irecord
         print *, modulo(int(paso), itsortie) == 0
         if (modulo(int(paso), itsortie) .ne. 0) then
            print*,"PAS DE PAS REGULIER"
            stop
         end if
!     &pr(2706)/rho(2706)/g+z(2706),pr(2858)/rho1/g+z(2858),
!     &zs(65,1),zs(217,1)
         if (ith == 1) then
         do i=1,nm
         write (942, *) paso/unitsortie,x(i), z(i), temp(i)
         enddo
!ccc....S_temp_zns1.dat') au milieu RD
            write (181819, *) paso/unitsortie, temp(257), temp(578), &
               temp(866), temp(1335)
!ccc....S_temp_zns2.dat') proche riviere RD
            write (181820, *) paso/unitsortie, temp(1301), temp(2146), &
               temp(3036), temp(3958)
!ccc....S_temp_zns3.dat') loin riv RD
            write (181818, *) paso/unitsortie, temp(137), temp(511), &
               temp(1035), temp(1856)
!ccc....S_temp_TRes1.dat') proche riv
            write (181821, *) paso/unitsortie, temp(187), temp(495), temp(763), &
               temp(1175), temp(1840), temp(2358)
!ccc....S_temp_TRes2.dat') loin riv
            write (181822, *) paso/unitsortie, temp(125), temp(420), temp(674), &
               temp(1023), temp(1675), temp(2188)
!ccc....S_temp_HoboRD.dat')
            write (181823, *) paso/unitsortie, temp(6683), temp(7863), temp(8748), &
               temp(9633)
!ccc....S_temp_HoboRG.dat')
            write (181824, *) paso/unitsortie, temp(6927), temp(8107), temp(8992), &
               temp(9877)
         end if
            do i=1,nm
            write (94, *) paso/unitsortie,x(i), z(i), akv(i)*akrv(i)*rho(i)*g/amu
            write (941, *) paso/unitsortie,x(i), z(i), sw(i)

            write (943, *) paso/unitsortie,x(i), z(i), pr(i)
            end do              
!ccc....S_piezoB_RG.dat') zone 4
         write (181826, *) paso/unitsortie, pr(2706)/rho(2706)/g + z(2706), &
            zs(65, 1)
!ccc....piezoB_RD.dat' zone 7
         write (181825, *) paso/unitsortie, pr(2858)/rho(2858)/g + z(2858), &
            zs(217, 1)

         if (paso < 57*86400 + 1000 .and. paso > 57*86400 - 1000) then
         do kkcol = 1, nc
            write (181827, *) x(ibas(kkcol)), &
               zs(kkcol, 1)
         end do
         end if

         if (paso < 133*86400 + 1000 .and. paso > 133*86400 - 1000) then
         do kkcol = 1, nc
            write (181831, *) x(ibas(kkcol)), &
               zs(kkcol, 1)
         end do
         end if

         if (paso < 134*86400 + 1000 .and. paso > 134*86400 - 1000) then
         do kkcol = 1, nc
            write (181828, *) x(ibas(kkcol)), &
               zs(kkcol, 1)
         end do
         end if

         if (paso < 195*86400 + 1000 .and. paso > 195*86400 - 1000) then
         do kkcol = 1, nc
            write (181829, *) x(ibas(kkcol)), &
               zs(kkcol, 1)
         end do
         end if

         if (paso < 10*86400 + 1000 .and. paso > 10*86400 - 1000) then
         do i = 1, nm
            write (751, *) pr(i)
         end do
         end if
      end if

      case ("R2D")
         
      if (irecord == 1) then
         qflux = 0D+00
         do j = 1, ligne6
            qflux = vzm(id_river(j))*am(id_river(j)) + qflux
         end do
        write (181830, *) paso/unitsortie, qflux
        print *, "out", paso/86400, temp(nmaille1),temp(1),vzm(100),vzp(100)
!        write (18181, *) paso/unitsortie, x(i),z(i), vzm(i),vzp(i)
!     &pr(2706)/rho(2706)/g+z(2706),pr(2858)/rho1/g+z(2858),
!     &zs(65,1),zs(217,1)
         if (ith == 1) then
!         do i=1,nm
!         write (942, *) paso/unitsortie,x(i), z(i), temp(i)
!         enddo
            write (56, *) paso/unitsortie, temp(nmaille1), temp(nmaille2), temp(nmaille3), &
            temp(nmaille4), temp(nmaille5), temp(nmaille6), temp(nmaille7), &
            temp(nmaille8)
            write (57, *) paso/unitsortie, pr(nmaille4)/rho1/g+z(nmaille4), &
            pr(nmaille8)/rho1/g+z(nmaille8)
         end if
 !           do i=1,nm
 !           write (943, *) paso/unitsortie,x(i), z(i), pr(i)
 !           end do              

      end if

!CC....SORTIE ZH Karina
       case( "ZHR" ,"ZHZ")
    !           print*,paso,vzm(1),valclt(1,3),temp(nmaille1)
      !         print*,temp(nmaille2),temp(nmaille3)
         if (irecord == 1 .and. irpth == 1) then
            write (59, *) paso/unitsortie, temp(nmaille1)
            write (60, *) paso/unitsortie, temp(nmaille2)
            write (61, *) paso/unitsortie, temp(nmaille3)
            do i = 1, nm
               write (18181, *) paso/unitsortie, z(i), vzm(i)
               write (18182, *) paso/unitsortie, z(i), qcondu(i), qad(i), qtherm(i)
               write (18183, *) paso/unitsortie, z(i), temp(i)
            end do
         end if


!CC....SORTIE sol Emmanuel leger
      case ( "1DS") 

         print *, paso/unitsortie, dt, pr(1)/9800 + z(1), pr(nm)/9800 + z(nm)
         print *, 'bas=', valcl(nm, 4)/9800 + z(nm)
         print *, 'pluie=', valcl(1, 3)*1000*86400, 'mm/j'
         print *, 'sat=', sw(1)

         if (irecord == 1 .and. irpth == 1) then
            write (59, *) paso/unitsortie, temp(nmaille1)
            write (60, *) paso/unitsortie, temp(nmaille2)
            write (61, *) paso/unitsortie, temp(nmaille3)
            write (64, *) paso/unitsortie, temp(nmaille4)
            write (65, *) paso/unitsortie, temp(nmaille5)
            write (66, *) paso/unitsortie, temp(nmaille6)
            write (67, *) paso/unitsortie, temp(nmaille7)
            write (592, *) paso/unitsortie, temp(nmaille8)
            write (69, *) paso/unitsortie, temp(nmaille9)
            write (591, *) paso/unitsortie, temp(nmaille10)
            write (7782, *) paso/unitsortie, vzp(nmaille2), vzm(nmaille2)
            write (1818, *) paso/unitsortie, qcondu(1), qad(1), qtherm(1), vzm(1), &
               temp(1), temp(2)
            do i = 1, nm
               write (18181, *) paso/unitsortie, z(i), vzm(i)
               write (18182, *) paso/unitsortie, z(i), qcondu(i), qad(i), qtherm(i)
               write (18183, *) paso/unitsortie, z(i), temp(i)
               write (94, *) paso/unitsortie, z(i), akv(i)*akrv(i)*rho(i)*g/amu
               write (941, *) paso/unitsortie, z(i), sw(i)
               write (942, *) paso/unitsortie, z(i), alanda(i)
               write (943, *) paso/unitsortie, z(i), om(i)*sw(i)*rho(i)*cpe + &
                  om(i)*sice(i)*rhoi(i)*cpice + &
                  (1.-om(i))*rhos(i)*cps(i)
               write (181822, *) paso/unitsortie, z(i), pr(i), pr(i)/rho1/g + z(i)
            end do
         end if


!CC....SORTIE sol Jose
      case ( "1DJ") 

         print *, paso/unitsortie, dt, pr(1)/9800 + z(1), pr(nm)/9800 + z(nm)
         print *, 'bas=', valcl(nm, 4)/9800 + z(nm)
         print *, 'pluie=', valcl(1, 3)*1000*86400, 'mm/j'
         print *, 'sat=', sw(1)

         if (irecord == 1 .and. irpth == 1) then
            write (7782, *) paso/unitsortie, vzp(nmaille2), vzm(nmaille2)
            do i = 1, nm
               write (18181, *) paso/unitsortie, z(i), vzm(i)
               write (94, *) paso/unitsortie, z(i), akv(i)*akrv(i)*rho(i)*g/amu
               write (941, *) paso/unitsortie, z(i), sw(i)
               write (181822, *) paso/unitsortie, z(i), pr(i), pr(i)/rho1/g + z(i)
            end do
         end if
      

!CC....SORTIES TEST LUNARDINI
      case( "THL") 
         print *, "OUT", dt, paso
         if (irecord == 1) then
            print *, "OUT", paso/unitsortie, pr(nm)/rho(1)/g + z(nm), &
               temp(nmaille1), valclt(nm, 4), temp(1), pr(1), akr(1)
            if (paso < 86400 + 1000 .and. paso > 86400 - 1000) then
               do i = 1, nm
                  write (741, *) i, x(i), z(i), pr(i), pr(i)/rho1/g + z(i), temp(i), sw(i)
               end do
            end if
            if (paso < 2*86400 + 1000 .and. paso > 2*86400 - 1000) then
               do i = 1, nm
                  write (742, *) i, x(i), z(i), pr(i), pr(i)/rho1/g + z(i), temp(i), sw(i)
               end do
            end if
            if (paso < 3*86400 + 1000 .and. paso > 3*86400 - 1000) then
               do i = 1, nm
                  write (743, *) i, x(i), z(i), pr(i), pr(i)/rho1/g + z(i), temp(i), sw(i)
               end do
            end if
         end if

!CC...SORTIES MAQUETTE EXP CYCLE 2
!       print*,paso/86400,dt
!       print*,pr(nm)/rho(nm)/g+z(nm)
!     &,akr(7666),cimp,alanda(10321),temp(1)

      case ("MAQ") 
         if (irecord == 1) then
         if (paso < 10*86400 + 1000 .and. paso > 10*86400 - 1000) then
            do i = 1, nm
               write (742, *) i, x(i), z(i), pr(i), pr(i)/rho1/g + z(i), temp(i), sw(i)
            end do
         end if
         print *, "OUT", paso/86400
         write (53, *) paso/unitsortie, &
            pr(nmaille1)/rho(nmaille1)/g + z(nmaille1), &
            pr(nmaille2)/rho(nmaille2)/g + z(nmaille2), &
            pr(nmaille3)/rho(nmaille3)/g + z(nmaille3), &
            pr(nmaille4)/rho(nmaille4)/g + z(nmaille4), &
            pr(nmaille5)/rho(nmaille5)/g + z(nmaille5), &
            pr(nmaille6)/rho(nmaille6)/g + z(nmaille6), &
            pr(nmaille7)/rho(nmaille7)/g + z(nmaille7), &
            pr(nmaille8)/rho(nmaille8)/g + z(nmaille8), &
            pr(nmaille9)/rho(nmaille9)/g + z(nmaille9)
         write (54, *) paso/unitsortie, temp(nmaille1), &
            temp(nmaille2), temp(nmaille3), temp(nmaille4), temp(nmaille5), &
            temp(nmaille6), temp(nmaille7), temp(nmaille8), temp(nmaille9)
!CC...POSITION PT100 MILIEU MAQUETTE ID X Z
!ccc....81       1.0062500000000001      1.0062499999999999
!ccc....1361       1.0062500000000001       0.90625000000000022
!ccc....2641       1.0062500000000001       0.80625000000000058
!ccc....3921       1.0062500000000001       0.70625000000000093
!ccc....5201       1.0062500000000001       0.60625000000000129
!ccc....6481       1.0062500000000001       0.50625000000000164
!ccc....7761       1.0062500000000001       0.40625000000000155
!ccc....9041       1.0062500000000001       0.30625000000000147
!ccc....10321       1.0062500000000001       0.20625000000000138
!ccc....11601       1.0062500000000001       0.10625000000000132
!ccc....12881       1.0062500000000001      6.2500000000013309E-003
         write (56, *) paso/unitsortie, temp(81), temp(1361), temp(2641), &
            temp(3921), temp(5201), temp(6481), temp(7761), &
            temp(9041), temp(10321), temp(11601), temp(12881)
!CC....Position isotherm 0 et liquidus
         write (18, *) paso/unitsortie, zs(20, 1), zs(20, 2), zl(20, 1), zl(20, 2)
!CC....Soulevement
         write (19, *) paso/unitsortie, dl(20), def(20)
         call flush (19)
         call flush (56)
         call flush (53)
!      do kkcol=1,nc
!      write(62,*) paso/unitsortie,kkcol,x(ibas(kkcol)),zs(kkcol,1),zs(kkcol,2),
!     &zl(kkcol,1),zl(kkcol,2)
!      write(63,*) paso/unitsortie,dl(kkcol),def(kkcol)
!      enddo
         end if

!CC....sortie tests performances
!CC...TEST NEUMAN
      case ( "TH1") 
         print *, "OUT", paso/unitsortie, irecord
         if (irecord == 1) then

            write (18, *) paso/unitsortie, zs(1, 1), zs(1, 2), zl(1, 1), zl(1, 2)
            write (53, *) paso/unitsortie, temp(400)
            write (1818, *) paso/unitsortie, qtherm(1)
            write (181818) real(paso/unitsortie), real(gxl), real(gxs)
            call flush (181818)
            open (74, file='S_pts', position='rewind', &
                  form='unformatted')
            write (74) real(paso/unitsortie)
            do i = 1, nm
               write (74) real(pr(i)), real(temp(i)), real(sw(i))
               call flush (74)
            end do
            rewind (unit=74)
         end if

!CC....TEST TH2
      case ("TH2") 
      if (irecord == 1) then
         write (181818) real(paso/unitsortie), real(tempmin), &
            real(swtotal), real(PF2)
         call flush (181818)
         open (74, file='S_pts', position='rewind', &
               form='unformatted')
         do i = 1, nm
            write (74) real(pr(i)), real(temp(i)), real(sw(i))
         end do
         rewind (unit=74)
      end if

!CC....TEST TH3
      case ("TH3") 
         print *, real(paso/unitsortie), real(akeq), real(qtcol), &
            real(qthermtot), real(pt1), real(pt2)
         if (irecord == 1) then
            write (181818) real(paso/unitsortie), real(akeq), &
               real(qtcol), real(qthermtot), real(pt1), real(pt2)
!              t,PM1,PM2,PM3,PM4_Pt1,PM4_Pt2
            call flush (181818)
            open (74, file='S_pts', position='rewind', &
                  form='unformatted')
            write (74) real(paso/unitsortie)
            do i = 1, nm
               write (74) real(pr(i)), real(temp(i)), real(sw(i))
               call flush (74)
            end do
            rewind (unit=74)
         end if
      end select

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                      C
!       Compteur temps simul        C
!                      C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!       CALL CPU_Time(Seconds)
!       trest=(nitt-it)*(Seconds-seicold)
!       if(it > 1) then
!       seicold=Seconds
!       print*,it,' temps restant: ',trest/60.,dt,pr(1),paso
!       endif
   end do
      print *, "simu done"
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                      C
!       FIN BOUCLE TEMPS          C
!                      C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!       CALL CPU_Time(Seconds)
!...FICHIERS FIN DE SIMULATION
!ccc....Dernier pas de temps ou regime permanant
   do i = 1, nm
      write (74, FMT='(i6,BN,F14.2,F14.2,F14.3,F15.5,F14.2,E14.4,E14.4,E14.4,E14.4)') &
         i, x(i), z(i), pr(i), &
         pr(i)/(rho1*g) + z(i), temp(i),vxm(i), vzm(i),vxp(i),vzp(i)
   end do

   if (ytest == "AVA") then
      rewind (62)
      do kkcol = 1, nc
         write (62, *) x(ibas(kkcol)), zs(kkcol, 1)
      end do
   end if

   if (ytest == "VIT") then
      if (it == 1) then
         print *, 'thermal conductivity=', alanda(1)
         print *, 'apparent heat capacity=', om(1)*sw(1)*rho(1)*cpe + &
            om(1)*(1 - sw(1))*rhog*cpg + (1.-om(1))*rhos(1)*cps(1)
         print *, 'hydraulic conductivity=', akr(1)*ak(1)*rho(1)/(amu)
         do i = 1, nm
         if (icl(i, 4) == -2) then
            print *, 'h(nm) [m]', pr(i)/(rho(i)*g) + z(i), 'pr(nm) [Pa]', pr(i)
            prbas = pr(i)
            zbas = z(i)
         end if
         if (icl(i, 3) == -2) then
            print *, 'h(1)[m]', pr(i)/(rho(i)*g) + z(i), 'pr(1) [Pa]', pr(i)
            prhaut = pr(i)
            zhaut = z(i)
         end if
         end do
         print *, 'Z_Obs1', z(nmaille1)
         print *, 'Z_Obs2', z(nmaille2)
         print *, 'Z_Obs3', z(nmaille3)
         print *, 'imposed hydraulic gradient [m]', (prhaut/(rho1*g) + zhaut &
                                                     - (prbas/(rho1*g) + zbas))/(zhaut - zbas)
!CC vitessse face haute si sortant négatif, entrant positif
         print *, 'velocity=', -vzp(1)
      end if

      if (irecord == 1) then
         print *, 't(s)', paso, 'T_obs [°C]', temp(nmaille1), temp(nmaille2), &
            temp(nmaille3), 'velocity=', -vzp(1)
         do i = 1, nm
         if (icl(i, 4) == -2) then
            print *, 'h(nm) [m]', pr(i)/(rho(i)*g) + z(i), 'pr(nm) [Pa]', pr(i)
            prbas = pr(i)
            zbas = z(i)
         end if
         if (icl(i, 3) == -2) then
            print *, 'h(1)[m]', pr(i)/(rho(i)*g) + z(i), 'pr(1) [Pa]', pr(i)
            prhaut = pr(i)
            zhaut = z(i)
         end if
         end do
         print *, 'imposed hydraulic gradient [m]', (prhaut/(rho1*g) + zhaut &
                                                     - (prbas/(rho1*g) + zbas))/(zhaut - zbas)
!CC vitessse face haute si sortant négatif, entrant positif
         print *, 'velocity=', -vzp(1)
         write (59, *) paso/us, temp(nmaille1)
         write (60, *) paso/us, temp(nmaille2)
         write (61, *) paso/us, temp(nmaille3)
      end if
   end if

!CC....Fermeture des fichiers!!!
!     call purge_noms_fichiers

   close (96)
   close (91)
   close (94)
   close (42)
   close (41)
   close (49)
   close (11)
   close (12)
   close (13)
   close (14)
   close (140)
   close (15)
   close (16)
   close (17)
   close (20)
   close (22)
   close (23)
   close (24)
   close (25)
   close (26)
   close (27)
   close (28)
   close (29)
   close (31)
   close (32)
   close (33)
   close (34)
   close (35)
   close (40)
   close (55)
   close (54)
   close (56)
   close (57)
   close (58)
   close (59)
   close (60)
   close (61)
   close (62)
   close (18)
   close (19)
   close (63)
   close (67)
   close (53)
   close(321)
   close (6001)
   ! close (74)
   close (778)
   close (782)
   close (1000)
   close (37)
   close (38)
   close (784)
   close (786)
   close (45)
   close (68)
   close (74)
   close (49)
   close (91)
   close (62)
   close (63)
   close (67)
   close (782)
   close (786)
   close (52)
   close (53)
   close (94)
   close (96)
   close (59)
   close (60)
   close (61)
   close (7782)
   close (64)
   close (65)
   close (19)
   close (53)
   close (54)
   close (56)
   close (18)
   close (743)
   close (742)
   close (741)
   close (18)
   close (53)
   close (1818)
   close (181818)
   close (181819)
   close (181820)
   close (181821)
   close (181822)
   close (181823)
   close (181824)
   close (181825)
   close (181826)
   close (62)
   close (751)
   close (181827)
   close (181840)
   close (181841)
   close (181842)
   close (181843)
   close (181845)
   close (181846)
   close (181847)
   close (181848)
   close (329)
   close (330)
   close (32)
   close (321)
   close (331)

   if (allocated(irow)) deallocate(irow)
   if (allocated(icol_ind)) deallocate(icol_ind)
   if (allocated(irow_ptr)) deallocate(irow_ptr)
   if (allocated(b)) deallocate(b)

   if (allocated(val)) deallocate(val)
   if (allocated(ss)) deallocate(ss)
   if (allocated(dswdp)) deallocate(dswdp)
   if (allocated(om)) deallocate(om)
   if (allocated(rhold)) deallocate(rhold)
   if (allocated(akr)) deallocate(akr)
   if (allocated(ak)) deallocate(ak)
   if (allocated(akv)) deallocate(akv)
   if (allocated(akrv)) deallocate(akrv)

   if (allocated(rho)) deallocate(rho)
   if (allocated(sw)) deallocate(sw)

   if (allocated(x)) deallocate(x)
   if (allocated(z)) deallocate(z)
   if (allocated(inum)) deallocate(inum)
   if (allocated(inum2)) deallocate(inum2)
   if (allocated(icol)) deallocate(icol)
   if (allocated(am)) deallocate(am)
   if (allocated(bm)) deallocate(bm)

   if (allocated(qbot)) deallocate(qbot)
   if (allocated(qsurf)) deallocate(qsurf)
   if (allocated(icl)) deallocate(icl)
   if (allocated(valcl)) deallocate(valcl)

   if (allocated(ivois)) deallocate(ivois)

   if (allocated(topo)) deallocate(topo)
   if (allocated(bot)) deallocate(bot)

   if (allocated(ibas)) deallocate(ibas)

!     if    (allocated(zbot))   DEALLOCATE(zbot)

   if (allocated(vxm)) deallocate(vxm)
   if (allocated(vxp)) deallocate(vxp)
   if (allocated(vzp)) deallocate(vzp)
   if (allocated(vzm)) deallocate(vzm)
   if (allocated(pr)) deallocate(pr)
   if (allocated(prk)) deallocate(prk)
   if (allocated(pro)) deallocate(pro)

   if (icycle == 1) then
      if (allocated(alph)) deallocate(alph)
      if (allocated(dl)) deallocate(dl)
      if (allocated(def)) deallocate(def)
      if (allocated(defo)) deallocate(defo)
      if (allocated(dsipdtemp)) deallocate(dsipdtemp)
      if (allocated(dsidp)) deallocate(dsidp)
      if (allocated(dsipdp)) deallocate(dsipdp)
      if (allocated(dsidtempoo)) deallocate(dsidtempoo)
      if (allocated(dsidtempo)) deallocate(dsidtempo)
      if (allocated(siceoo)) deallocate(siceoo)
      if (allocated(dsidtempoo)) deallocate(dsidtempoo)
      if (allocated(dsidtempo)) deallocate(dsidtempo)
      if (allocated(zlo)) deallocate(zlo)
      if (allocated(zloo)) deallocate(zloo)
      if (allocated(zl)) deallocate(zl)

   end if
   if (itr == 1) then
      if (allocated(iclc)) deallocate(iclc)
      if (allocated(valclc)) deallocate(valclc)
      if (allocated(conco)) deallocate(conco)
      if (allocated(conck)) deallocate(conck)
      if (allocated(conc)) deallocate(conc)
   end if
   if (allocated(chgbot)) deallocate(chgbot)
   if (allocated(chgsurf)) deallocate(chgsurf)
   if (allocated(tempbot)) deallocate(tempbot)
   if (allocated(tempsurf)) deallocate(tempsurf)
   if (allocated(asun)) deallocate(asun)
   if (allocated(ansun)) deallocate(ansun)
   if (allocated(dswpdp)) deallocate(dswpdp)
   if (allocated(dswdt)) deallocate(dswdt)

   if (allocated(zs)) deallocate(zs)
   if (allocated(zso)) deallocate(zso)
   if (allocated(zsoo)) deallocate(zsoo)
   if (allocated(dsidtemp)) deallocate(dsidtemp)
   if (allocated(sice)) deallocate(sice)
   if (allocated(rhoi)) deallocate(rhoi)
   if (allocated(siceo)) deallocate(siceo)
   if (allocated(sicep)) deallocate(sicep)

   if (ichi == 1) then
      if (allocated(chg)) deallocate(chg)
   end if

   if (allocated(tempoo)) deallocate(tempoo)
   if (allocated(temp)) deallocate(temp)
   if (allocated(tempo)) deallocate(tempo)
   if (allocated(tempk)) deallocate(tempk)
   if (ith == 1) then
      if (allocated(valclt)) deallocate(valclt)
      if (allocated(rhos)) deallocate(rhos)
      if (allocated(cps)) deallocate(cps)
      if (allocated(alanda)) deallocate(alanda)
      if (allocated(qad)) deallocate(qad)
      if (allocated(alandas)) deallocate(alandas)
      if (allocated(iclt)) deallocate(iclt)

      if (allocated(qtherm)) deallocate(qtherm)
      if (allocated(valclto)) deallocate(valclto)

      if (allocated(qcondu)) deallocate(qcondu)

   end if

   if (ytest == "AVA" .or. ytest == "ZHZ"&
  &.or. ytest == "DTS" .or. ytest == "TEX"&
  &.or. ytest == "R2D") then
      if (allocated(alandazone)) deallocate(alandazone)
      if (allocated(cpmzone)) deallocate(cpmzone)
      if (allocated(rhomzone)) deallocate(rhomzone)
      if (allocated(izone)) deallocate(izone)
      if (allocated(jzone)) deallocate(jzone)
      if (allocated(akzone)) deallocate(akzone)
      if (allocated(omzone)) deallocate(omzone)
   end if

   if (ytest == "AVA" .or. ytest == "TEX"&
  &.or. ytest == "R2D") then
      if (allocated(tempRD)) deallocate(tempRD)
      if (allocated(tempRG)) deallocate(tempRG)
      if (allocated(id_RD)) deallocate(id_RD)
      if (allocated(id_RG)) deallocate(id_RG)
      if (allocated(id_ZH)) deallocate(id_ZH)
      if (allocated(chgRD)) deallocate(chgRD)
      if (allocated(chgRG)) deallocate(chgRG)
      if (allocated(id_river)) deallocate(id_river)
      if (allocated(id_rivert)) deallocate(id_rivert)
      if (allocated(chgriver)) deallocate(chgriver)
      if (allocated(tempriver)) deallocate(tempriver)
      if (allocated(tempsol)) deallocate(tempsol)
      if (allocated(qpluie)) deallocate(qpluie)
      if (allocated(anszone)) deallocate(anszone)
      if (allocated(swreszone)) deallocate(swreszone)
      if (allocated(aspzone)) deallocate(aspzone)
      if (allocated(swresz)) deallocate(swresz)
   end if
   if (ytest == "DTS") then
      if (allocated(tempRD)) deallocate(tempRD)
      if (allocated(tempRG)) deallocate(tempRG)
      if (allocated(timeG)) deallocate(timeG)
      if (allocated(timeD)) deallocate(timeD)
      if (allocated(timeDTS)) deallocate(timeDTS)
      if (allocated(cRivG)) deallocate(cRivG)
      if (allocated(cRivD)) deallocate(cRivD)
      if (allocated(xDTS)) deallocate(xDTS)
      if (allocated(tempDTS)) deallocate(tempDTS)
      if (allocated(xDTS)) deallocate(xDTS)
      if (allocated(slopeRH)) deallocate(slopeRH)
      if (allocated(xpool)) deallocate(xpool)
      if (allocated(xriffle)) deallocate(xriffle)
      if (allocated(qout_w)) deallocate(qout_w)
      if (allocated(qout_wR)) deallocate(qout_wR)
      if (allocated(qin_w)) deallocate(qin_w)
      if (allocated(qin_wR)) deallocate(qin_wR)
      if (allocated(qout_h)) deallocate(qout_h)
      if (allocated(qout_hR)) deallocate(qout_hR)
      if (allocated(qin_h)) deallocate(qin_h)
      if (allocated(qin_hR)) deallocate(qin_hR)
      if (allocated(qadvout_h)) deallocate(qadvout_h)
      if (allocated(qadvout_hR)) deallocate(qadvout_hR)
      if (allocated(qadvin_h)) deallocate(qadvin_h)
      if (allocated(qadvin_hR)) deallocate(qadvin_hR)
      if (allocated(qcondout_h)) deallocate(qcondout_h)
      if (allocated(qcondout_hR)) deallocate(qcondout_hR)
      if (allocated(qcondin_h)) deallocate(qcondin_h)
      if (allocated(qcondin_hR)) deallocate(qcondin_hR)
   end if

!     if    (allocated(zaqui))  DEALLOCATE(zaqui)

end program

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                      C
!                                      C
!         SUBROUTINES APPELEES DS Le PG PCPAL             C
!                                      C
!                                      C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                         C
!    open file            C
!                         C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC...PURPOSE :
!ccc...
!CC...Modification
!ccc...
!CC....TODOLIST
!ccc....
SUBROUTINE open_file(FileName, An_Error, My_LUN)
   character(*), INTENT(IN) :: FileName
   LOGICAL, INTENT(OUT) :: An_Error

   Integer, INTENT(IN)       :: My_LUN
   LOGICAL                         :: It_Exists

   INQUIRE (File=FileName, EXIST=It_Exists)
   IF (It_Exists) THEN
      An_Error = .FALSE.
   ELSE
      An_Error = .TRUE.
      RETURN
   END IF

   OPEN (UNIT=My_LUN, FILE=FileName)

   return
end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                         C
!    compter ligne            C
!                         C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC...PURPOSE :
!ccc...
!CC...Modification
!ccc...
!CC....TODOLIST
!ccc....
SUBROUTINE count_file(My_LUN, iosnew, Num)
   integer, INTENT(OUT) :: num
   integer :: iosnew
   Integer                         :: My_LUN
        num=0
   do while (iosnew == 0)
      read (My_LUN, *, iostat=iosnew)
      if (iosnew == 0) then
         num = num + 1
      end if
   end do
   rewind (My_LUN)
       

   return
end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                         C
!    Bi-gradient conjugue             C
!                         C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC...PURPOSE :
!ccc...
!CC...Modification
!ccc...
!CC....TODOLIST
!ccc....
subroutine bicg(x, b, n, k, val, icol_ind, irow_ptr, nmax, nmax1)
   implicit double precision(a - h, o - x, z), integer(I - N)
   implicit CHARACTER*5(y)
   integer k
   dimension val(nmax), icol_ind(nmax), irow_ptr(nmax1)
   dimension x(n), r(n), z(n), p(n)
   dimension b(n), q(n), rt(n), qt(n), zt(n), pt(n)
!CC....Residu INITIAL
!ccc....produit matrice vecteur
   do i = 1, n
      sum = 0.D+00
      do j = irow_ptr(i), irow_ptr(i + 1) - 1
         sum = sum + val(j)*x(icol_ind(j))
      end do
      r(i) = b(i) - sum
      rt(i) = r(i)
   end do
!ccc....meth gradient conjugue
   amaxr = 1.
   k = 0
   do while (amaxr >= 1.e-10 .and. k <= 1000)
      k = k + 1
      if (k > 1) rhoo = rho
      rho = 0.D+00
      tt = 1
      do i = 1, n
      do j = irow_ptr(i), irow_ptr(i + 1) - 1
         if (i == j) tt = val(j)
      end do
      z(i) = r(i)/tt
      zt(i) = rt(i)/tt
      rho = rho + rt(i)*z(i)
      end do
      if (rho == 0.) then
!       print*,'attention rho=0'
!       write (20,*) 'rho=0'
         go to 11
      end if
      if (k == 1) then
      do i = 1, n
         p(i) = z(i)
         pt(i) = zt(i)
      end do
      else
      beta = rho/rhoo
      if (beta == 0.) then
!       print*,'attention beta=0'
!       write (20,*) 'sum=0'
      end if
      do i = 1, n
         p(i) = z(i) + beta*p(i)
         pt(i) = zt(i) + beta*pt(i)

      end do
      end if

      sum = 0.D+00
      do i = 1, n
         q(i) = 0.D+00
!       sum=0.
!ccc....produit matrice vecteur
         do j = irow_ptr(i), irow_ptr(i + 1) - 1
            q(i) = q(i) + val(j)*p(icol_ind(j))
         end do
         sum = sum + q(i)*pt(i)
      end do
      alfa = rho/sum
      do i = 1, n
         qt(i) = 0.D+00
      end do
      do j = 1, n
      do i = irow_ptr(j), (irow_ptr(j + 1) - 1)
         qt(icol_ind(i)) = qt(icol_ind(i)) + val(i)*pt(j)
      end do
      end do

      do i = 1, n
         x(i) = x(i) + alfa*p(i)
         r(i) = r(i) - alfa*q(i)
         rt(i) = rt(i) - alfa*qt(i)
         if (r(i) == 0.) rt(i) = 0.D+00

      end do

!ccc....limitation du nombre d iteration
      amaxr = 0.D+00
      do i = 1, n
      if (abs(r(i)) > amaxr) then
         amaxr = abs(r(i))
         iii = i
      end if
      end do
   end do
11 continue
end
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                         C
!   CGS                       C
!                         C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC...PURPOSE :
!ccc...
!CC...Modification
!ccc...
!CC....TODOLIST
!ccc....
subroutine cgs(x, b, n, k, val, icol_ind, irow_ptr, nmax, nmax1)
   implicit double precision(A - H, O - Z), integer*4(I - N)
   integer k
   dimension val(nmax), icol_ind(nmax), irow_ptr(nmax1)
   dimension x(n), r(n), u(n), p(n), vc(n)
   dimension b(n), q(n), rt(n), qc(n), uc(n), pc(n)
   double precision :: sum, rho, rhoo, beta, alfa, amaxr
   integer :: i, j, iii

!ccccccResidu INITIAL
! produit matrice vecteur
   !$omp parallel do private(i, j, sum) shared(val, icol_ind, irow_ptr, x, b, r, rt)
   do i = 1, n
      sum = 0.D0
      do j = irow_ptr(i), irow_ptr(i + 1) - 1
         sum = sum + val(j) * x(icol_ind(j))
      end do
      r(i) = b(i) - sum
      rt(i) = r(i)
   end do
   !$omp end parallel do

!ccccccmeth gradient conjugue
   amaxr = 1.D0
   k = 0

   do while (amaxr >= 1.e-10 .and. k <= 1000)
      k = k + 1
      if (k > 1) rhoo = rho
      rho = 0.D0

      !$omp parallel for reduction(+:rho) private(i) shared(rt, r)
      do i = 1, n
         rho = rho + rt(i) * r(i)
      end do
      !$omp end parallel for

      if (rho == 0.D0) then
         exit
      end if

      if (k == 1) then
         !$omp parallel do private(i) shared(u, p, r)
         do i = 1, n
            u(i) = r(i)
            p(i) = u(i)
         end do
         !$omp end parallel do
      else
         beta = rho / rhoo
         if (beta == 0.D0) then
            ! print*,'attention beta=0'
         end if
         !$omp parallel do private(i) shared(u, p, r, q, beta)
         do i = 1, n
            u(i) = r(i) + beta * q(i)
            p(i) = u(i) + beta * (q(i) + beta * p(i))
         end do
         !$omp end parallel do
      end if

      !$omp parallel do private(i) shared(pc, p)
      do i = 1, n
         pc(i) = p(i)
      end do
      !$omp end parallel do

      sum = 0.D0
      !$omp parallel do private(i, j) reduction(+:sum) shared(vc, val, icol_ind, irow_ptr, pc, rt)
      do i = 1, n
         vc(i) = 0.D0
         do j = irow_ptr(i), irow_ptr(i + 1) - 1
            vc(i) = vc(i) + val(j) * pc(icol_ind(j))
         end do
         sum = sum + vc(i) * rt(i)
      end do
      !$omp end parallel do

      alfa = rho / sum

      !$omp parallel do private(i) shared(q, u, vc, alfa, uc, x)
      do i = 1, n
         q(i) = u(i) - alfa * vc(i)
         uc(i) = u(i) + q(i)
         x(i) = x(i) + alfa * uc(i)
      end do
      !$omp end parallel do

      !$omp parallel do private(i, j) shared(qc, val, icol_ind, irow_ptr, uc, r, alfa)
      do i = 1, n
         qc(i) = 0.D0
         do j = irow_ptr(i), irow_ptr(i + 1) - 1
            qc(i) = qc(i) + val(j) * uc(icol_ind(j))
         end do
         r(i) = r(i) - alfa * qc(i)
      end do
      !$omp end parallel do

      amaxr = 0.D0
      !$omp parallel do private(i) shared(r) reduction(max:amaxr)
      do i = 1, n
         if (abs(r(i)) > amaxr) then
            amaxr = abs(r(i))
            iii = i
         end if
      end do
      !$omp end parallel do

   end do

end


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!                       ECOULEMENT                      C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC...PURPOSE :
!ccc...
!CC...Modification
!ccc...
!CC....TODOLIST
!ccc....
subroutine matp(val, icol_ind, irow_ptr, x, z, b, am, ivois, &
                rho, ak, akr, amu, dt, ia2, g, icl, valcl, rhold, om, pro, dswdp, sw, nm, &
                irp, nmax, nmax1, bm, akv, akrv, igel, &
                ysupdp, rhoi, ixy, dsidtemp, temp, tempo)
   implicit double precision(a - h, o - x, z), integer(I - N)
   implicit CHARACTER*5(y)
   dimension b(nm), x(nm), z(nm)
   dimension ivois(nm, 4), dsidtemp(nm), temp(nm)
   dimension val(nmax1), icol_ind(nmax1), irow_ptr(nmax1)
   dimension dswdp(nm), om(nm), rhold(nm)
   dimension akr(nm), ak(nm), rho(nm), sw(nm), pro(nm)
   dimension icl(nm, 4), am(nm), valcl(nm, 4)
   dimension bm(nm), tempo(nm)
   dimension akrv(nm), akv(nm)
   dimension rhoi(nm)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                         C
!   Definition de la matrice pr lecoulement    C
!                         C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

   irow_ptr(1) = 1
   il = 0
   do ik = 1, nm

      b(ik) = 0D+00
      il = il + 1
      val(il) = 0D+00

!CC....Stockage du numero de la diagonale
      ilik = il
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                         C
!   transient state: storage yield         C
!                         C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC....irp=0 steady state annule term
      icol_ind(il) = ik
      val(il) = dble(-rho(ik)*om(ik)*dswdp(ik)/dt*irp)
      irow_ptr(ik + 1) = irow_ptr(ik) + 1
!CC....terme b
      b(ik) = dble(b(ik) - sw(ik)*ia2 &
                   *om(ik)*(rho(ik) - rhold(ik))/dt*irp &
                   - rho(ik)*om(ik)*dswdp(ik)*pro(ik)/dt*irp)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                         C
!         FACE HAUTE              C
!                         C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      if (icl(ik, 3) == 1) then
!ccc....Maille du haut existe terme gravite J+1/2 negatif annuler avec ixy
         il = il + 1
         jp = ivois(ik, 3)
         icol_ind(il) = jp
         dzh = dble(bm(ik)*abs(z(jp) - z(ik)))
         val(il) = dble(1./dzh*(bm(jp) + bm(ik)) &
                        /(bm(jp)/(akrv(jp)*akv(jp)*rho(jp)/amu) + &
                          bm(ik)/(akrv(ik)*akv(ik)*rho(ik)/amu)))
         val(ilik) = dble(val(ilik) - val(il))
         b(ik) = dble(b(ik) - ixy*(bm(jp) + bm(ik)) &
                      /(bm(jp)/(akrv(jp)*akv(jp)*rho(jp)/amu) + &
                        bm(ik)/(akrv(ik)*akv(ik)*rho(ik)/amu)) &
                      *1*g/bm(ik)*(rho(jp) + rho(ik))/2)
         irow_ptr(ik + 1) = irow_ptr(ik + 1) + 1
      end if

      if (icl(ik, 3) == -1) then
!ccc....Face haute a flux impose....
         b(ik) = dble(b(ik) - valcl(ik, 3)/bm(ik)*rho(ik))
      end if

      if (icl(ik, 3) == -2) then
!ccc....Face haute a potentiel impose....
         dzhi = dble((bm(ik)*bm(ik)/2))
         val(ilik) = dble(val(ilik) - 1/dzhi*rho(ik)*akv(ik)*akrv(ik)/amu)
         b(ik) = dble(b(ik) - ixy*1/dzhi*rho(ik)*akv(ik)*akrv(ik)/amu &
                      *valcl(ik, 3)) - 1/bm(ik)*rho(ik)*akv(ik)*akrv(ik)/amu*rho(ik)*g
      end if

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                         C
!         FACE BASSE              C
!                         C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!ccc....terme gravite j-1/2 positif
      if (icl(ik, 4) == 1) then
!ccc....Maille du ibas existe
         il = il + 1
         jm = ivois(ik, 4)
         icol_ind(il) = jm
         dzb = dble(bm(ik)*abs(z(ik) - z(jm)))
         val(il) = dble(1/dzb*(bm(jm) + bm(ik)) &
                        /(bm(jm)/(akrv(jm)*akv(jm)*rho(jm)/amu) + &
                          bm(ik)/(akrv(ik)*akv(ik)*rho(ik)/amu)))
         val(ilik) = dble(val(ilik) - val(il))
         b(ik) = dble(b(ik) + ixy*(bm(jm) + bm(ik)) &
                      /(bm(jm)/(akrv(jm)*akv(jm)*rho(jm)/amu) + &
                        bm(ik)/(akrv(ik)*akv(ik)*rho(ik)/amu)) &
                      *1*g/bm(ik)*(rho(jm) + rho(ik))/2)
         irow_ptr(ik + 1) = irow_ptr(ik + 1) + 1
      end if
      if (icl(ik, 4) == -1) then
!ccc....Face BASSE a flux impose....
         b(ik) = dble(b(ik) - valcl(ik, 4)/bm(ik)*rho(ik))
      end if
      if (icl(ik, 4) == -2) then
!ccc....Face basse a potentiel impose....
         dzbi = dble((bm(ik)*bm(ik)/2))
         val(ilik) = dble(val(ilik) - 1/dzbi*rho(ik)*akv(ik)*akrv(ik)/amu)
         b(ik) = dble(b(ik) - 1/dzbi*rho(ik)*akv(ik)*akrv(ik)/amu &
                      *valcl(ik, 4)) + ixy*1/bm(ik)*rho(ik)*akv(ik)*akrv(ik)/amu*rho(ik)*g
      end if

      if (icl(ik, 4) == -3) then
!ccc....Face basse a free drainage..
         b(ik) = dble(b(ik) - rho(ik)*akv(ik)*akrv(ik)/amu/bm(ik)*rho(ik))
      end if

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                         C
!         FACE GAUCHE              C
!                         C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      if (icl(ik, 2) == 1) then

!ccc....Maille de gauche existe CORRECTION TERME DZ/DX*RHO*G/AM terme i-1/2 positif

         il = il + 1
         im = ivois(ik, 2)
         icol_ind(il) = im

         dxg = dble(am(ik)*abs(x(ik) - x(im)))
         val(il) = dble(1/dxg*(am(im) + am(ik)) &
                        /(am(im)/(akr(im)*ak(im)*rho(im)/amu) + &
                          am(ik)/(akr(ik)*ak(ik)*rho(ik)/amu)))

         b(ik) = dble(b(ik) + (am(im) + am(ik))/(am(im)/ &
                                                 (akr(im)*ak(im)*rho(im)/amu) + &
                                                 am(ik)/(akr(ik)*ak(ik)*rho(ik)/amu)) &
                      *(rho(im) + rho(ik))/2*g*ixy*(z(im) - z(ik))/(x(im) - x(ik))/am(ik))

         val(ilik) = dble(val(ilik) - val(il))
         irow_ptr(ik + 1) = irow_ptr(ik + 1) + 1
      end if

      if (icl(ik, 2) == -1) then
!ccc....Face gauche a flux impose....
         b(ik) = dble(b(ik) - valcl(ik, 2)/am(ik)*rho(ik))
      end if
      if (icl(ik, 2) == -2) then
!cccccc Face gauchee a potentiel impose....
         b(ik) = dble(b(ik) - 2/am(ik)**2*rho(ik)*ak(ik)*akr(ik)/amu* &
                      valcl(ik, 2))
         val(ilik) = dble(val(ilik) - 2/am(ik)**2*rho(ik)*ak(ik)*akr(ik)/amu)
      end if

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                         C
!         FACE DROITE              C
!                         C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      if (icl(ik, 1) == 1) then
!ccc....Maille de droite existe CORRECTION TERME DZ/DX*RHO*G/AM terme i+1/2 negatif
         il = il + 1
         ip = ivois(ik, 1)
         icol_ind(il) = ip
         dxd = dble(am(ik)*abs(x(ip) - x(ik)))
         val(il) = dble(1/dxd*(am(ip) + am(ik))/(am(ip)/(akr(ip)*ak(ip) &
                                                         *rho(ip)/amu) + am(ik)/(akr(ik)*ak(ik)*rho(ik)/amu)))
         b(ik) = dble(b(ik) - (am(ip) + am(ik))/(am(ip) &
                                                 /(akr(ip)*ak(ip)*rho(ip)/amu) + am(ik)/ &
                                                 (akr(ik)*ak(ik)*rho(ik)/amu))*(rho(ip) + rho(ik))/2*g*ixy* &
                      (z(ip) - z(ik))/(x(ip) - x(ik))/am(ik))
         val(ilik) = dble(val(ilik) - val(il))
         irow_ptr(ik + 1) = irow_ptr(ik + 1) + 1
      end if

      if (icl(ik, 1) == -1) then
!ccc....Face droite a flux impose....
         b(ik) = dble(b(ik) - valcl(ik, 1)/am(ik)*rho(ik))
      end if

      if (icl(ik, 1) == -2) then
!ccc....Face droite a potentiel impose....
         b(ik) = dble(b(ik) - 2/am(ik)**2*rho(ik)*ak(ik)*akr(ik)/amu* &
                      valcl(ik, 1))
         val(ilik) = dble(val(ilik) - 2/am(ik)**2*rho(ik)*ak(ik)*akr(ik)/amu)
      end if

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                         C
!  Terme de supression expansion glace     C
!                         C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC....MCKENZIE CALCUL TAUX DE PRODUCTION DE GLACE
      if (ysupdp == "SPMCK") then
!ccc....GEL TERME SURPRESSION
         if (igel == 1) then
            b(ik) = b(ik) - &
                    abs((rho(ik) - rhoi(ik))*om(ik)*(dsidtemp(ik)) &
                        *(tempo(ik) - temp(ik))/dt)
         end if
!ccc....fin gel

!ccc....DEGEL  TERME SURPRESSION
         if (igel == 2) then
            b(ik) = b(ik) + &
                    abs((rho(ik) - rhoi(ik))*om(ik)*(dsidtemp(ik)) &
                        *(tempo(ik) - temp(ik))/dt)
         end if
      end if
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                         C
!       Ordonner Vecteurs CRS                 C
!                         C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC....ORDONNER LES TABLEAUX icol_IND ET VAL
      ideb = irow_ptr(ik)
      ifin = irow_ptr(ik + 1) - 1
      if (ifin == ideb) goto 4
!     if (ysolv == "BIC".or.ysolv == "CGS") then
2     continue
      do ii = ideb, ifin - 1
      if (icol_ind(ii) > icol_ind(ii + 1)) then
         idum = icol_ind(ii + 1)
         dum = val(ii + 1)
         icol_ind(ii + 1) = icol_ind(ii)
         val(ii + 1) = val(ii)
         icol_ind(ii) = idum
         val(ii) = dum
         goto 2
      end if
      end do
!     endif
!     if (ysolv == "LIB") then
!     do ii=ideb,ifin-1
!     if(icol_ind(ii) > icol_ind(ii+1)) then
!     idum=icol_ind(ii+1)
!     dum=val(ii+1)
!     icol_ind(ii+1)=icol_ind(ii)
!     val(ii+1)=val(ii)
!     icol_ind(ii)=idum
!     val(ii)=dum
!     endif
!     enddo
!     do ii=ideb,ifin-1
!     if(icol_ind(ii) == ik) then
!     idum=icol_ind(ifin)
!     dum=val(ifin)
!     icol_ind(ifin)=icol_ind(ii)
!     val(ifin)=val(ii)
!     icol_ind(ii)=idum
!     val(ii)=dum
!     endif
!     enddo
!     endif
4     continue
   end do
   return
end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!                 TRANSPORT                                C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC...PURPOSE :
!ccc...
!CC...Modification
!ccc...
!CC....TODOLIST
!ccc....
subroutine matc(val, icol_ind, irow_ptr, x, b, am, ivois, conco, &
                dt, iclc, valclc, om, nm, allg, alt, &
                nmax, nmax1, bm, z, vxp, vxm, vzp, vzm)
   implicit double precision(a - h, o - x, z), integer(I - N)
   implicit CHARACTER*5(y)
   dimension b(nm), x(nm), z(nm), am(nm)
   dimension ivois(nm, 4)
   dimension val(nmax), icol_ind(nmax), irow_ptr(nmax1)
   dimension om(nm)
   dimension bm(nm)
   dimension iclc(nm, 4), valclc(nm, 4), conco(nm)
   dimension vxm(nm), vxp(nm), vzp(nm), vzm(nm)

!CC....Definition de la matrice pr le transport

   irow_ptr(1) = 1
   il = 0
   do i = 1, nm
!CC....VITESSE EN X A TRAVERS LES FACES thermique
!ccc....vxm(i) face gauche
!ccc....vxp(i) face droite
!ccc....vzm(i) face du dessous
!ccc....vzp(i) face du dessus

!CC....TERME DISPERSIF
      dlip = allg*(vxp(i)**2 + (vzm(i) + vzp(i))**2/4)**0.5
      dlim = allg*(vxm(i)**2 + (vzm(i) + vzp(i))**2/4)**0.5
      dtjp = alt*(vzp(i)**2 + (vxm(i) + vxp(i))**2/4)**0.5
      dtjm = alt*(vzm(i)**2 + (vxm(i) + vxp(i))**2/4)**0.5
      dlip = allg*((vxm(i) + vxp(i))**2/4.+(vzm(i) + vzp(i))**2/4.)**0.5
      dlim = dlip
      dtjp = alt*((vxm(i) + vxp(i))**2/4.+(vzm(i) + vzp(i))**2/4.)**0.5
      dtjm = dtjp

      il = il + 1
      ili = il
      icol_ind(il) = i
      irow_ptr(i + 1) = irow_ptr(i) + 1
      val(il) = -1
      b(i) = 0.D+00

!CC....FACE DROITE
      if (iclc(i, 1) == 1) then
         val(ili) = val(ili) - dlip/am(i)/(x(ivois(i, 1)) - x(i))*dt/om(i)
         il = il + 1
         icol_ind(il) = ivois(i, 1)
         val(il) = dlip/am(i)/(x(ivois(i, 1)) - x(i))*dt/om(i)
         irow_ptr(i + 1) = irow_ptr(i + 1) + 1
      end if
!ccc....concentration imposee
      if (iclc(i, 1) == -2) then
         val(ili) = val(ili) - dlip/am(i)**2*2*dt/om(i)
         b(i) = b(i) - valclc(i, 1)*dlip/am(i)**2*2*dt/om(i)
      end if
!CC....FACE GAUCHE
      if (iclc(i, 2) == 1) then
         val(ili) = val(ili) - dlim/am(i)/(x(i) - x(ivois(i, 2)))*dt/om(i)
         il = il + 1
         icol_ind(il) = ivois(i, 2)
         val(il) = dlim/am(i)/(x(i) - x(ivois(i, 2)))*dt/om(i)
         irow_ptr(i + 1) = irow_ptr(i + 1) + 1
      end if
!ccc....concentration imposee
      if (iclc(i, 2) == -2) then
         val(ili) = val(ili) - dlim/am(i)**2*2*dt/om(i)
         b(i) = b(i) - valclc(i, 2)*dlim/am(i)**2*2*dt/om(i)
      end if
!CC....FACE HAUTE
      if (iclc(i, 3) == 1) then
         val(ili) = val(ili) - dtjp/bm(i)/(z(ivois(i, 3)) - z(i))*dt/om(i)
         il = il + 1
         icol_ind(il) = ivois(i, 3)

         val(il) = dtjp/bm(i)/(z(ivois(i, 3)) - z(i))*dt/om(i)
         irow_ptr(i + 1) = irow_ptr(i + 1) + 1
      end if
!ccc....concentration imposee
      if (iclc(i, 3) == -2) then
         val(ili) = val(ili) - dtjp/bm(i)**2*2*dt/om(i)
         b(i) = b(i) - valclc(i, 3)*dtjp/bm(i)**2*2*dt/om(i)
      end if
!CC....FACE BASSE
      if (iclc(i, 4) == 1) then
         val(ili) = val(ili) - dtjm/bm(i)/(z(i) - z(ivois(i, 4)))*dt/om(i)
         il = il + 1
         icol_ind(il) = ivois(i, 4)
         val(il) = dtjm/bm(i)/(z(i) - z(ivois(i, 4)))*dt/om(i)
         irow_ptr(i + 1) = irow_ptr(i + 1) + 1
      end if
!ccc....concentration imposee
      if (iclc(i, 4) == -2) then
         val(ili) = val(ili) - dtjm/bm(i)**2*2*dt/om(i)
         b(i) = b(i) - valclc(i, 4)*dtjm/bm(i)**2*2*dt/om(i)
      end if

!CC....ORDONNER LES TABLEAUX icol_IND ET VAL
      ideb = irow_ptr(i)
      ifin = irow_ptr(i + 1) - 1
      if (ifin == ideb) goto 11
10    continue
      do ii = ideb, ifin - 1
      if (icol_ind(ii) > icol_ind(ii + 1)) then
         idum = icol_ind(ii + 1)
         dum = val(ii + 1)
         icol_ind(ii + 1) = icol_ind(ii)
         val(ii + 1) = val(ii)
         icol_ind(ii) = idum
         val(ii) = dum
         goto 10
      end if
      end do
11    continue

!CC....REPERAGE DES VOISINS DS irow_ptr
      ii = 0
      id = 0
      ig = 0
      ih = 0
      ib = 0
      do ij = irow_ptr(i), irow_ptr(i + 1) - 1
         if (icol_ind(ij) == i) ii = ij
         if (icol_ind(ij) == ivois(i, 1)) id = ij
         if (icol_ind(ij) == ivois(i, 2)) ig = ij
         if (icol_ind(ij) == ivois(i, 3)) ih = ij
         if (icol_ind(ij) == ivois(i, 4)) ib = ij
      end do

!CC....TERMES ADVECTIFS

!ccc....FACE DROITE
      if (vxp(i) > 0) then
         val(ii) = val(ii) - vxp(i)*dt/am(i)/om(i)
      end if
      if (iclc(i, 1) == 1 .and. vxp(i) < 0) then
         val(id) = val(id) - vxp(i)*dt/am(i)/om(i)
      end if
!ccc....si face droite imposee et vxp(i)<0!!!!
      if (iclc(i, 1) == -2 .and. vxp(i) < 0) then
         b(i) = b(i) + vxp(i)*dt/am(i)/om(i)*valclc(i, 1)
      end if
!ccc....face droite flux impose (>0 U*C mol/s ou g/s)
      if (iclc(i, 1) == -1) then
         b(i) = b(i) - valclc(i, 1)*dt/am(i)/om(i)
      end if

!CC....FACE GAUCHE
      if (vxm(i) > 0 .and. iclc(i, 2) == 1) then
         val(ig) = val(ig) + vxm(i)*dt/am(i)/om(i)
      end if
!ccc....si face gauche imposee et vxm(i)>0!!!!
      if (vxm(i) > 0 .and. iclc(i, 2) == -2) then
         b(i) = b(i) - vxm(i)*dt/am(i)/om(i)*valclc(i, 2)
      end if
      if (vxm(i) < 0.) then
         val(ii) = val(ii) + vxm(i)*dt/am(i)/om(i)
      end if
!ccc....face gauche flux impose (U*C mol/s ou g/s)
      if (iclc(i, 2) == -1) then
         b(i) = b(i) - valclc(i, 2)*dt/am(i)/om(i)
      end if

!CC....FACE SUP
      if (vzp(i) > 0) then
         val(ii) = val(ii) - vzp(i)*dt/bm(i)/om(i)
      end if
      if (vzp(i) < 0 .and. iclc(i, 3) == 1) then
         val(ih) = val(ih) - vzp(i)*dt/bm(i)/om(i)
      end if
      if (vzp(i) < 0 .and. iclc(i, 3) == -2) then
         b(i) = b(i) + vzp(i)*dt/bm(i)/om(i)*valclc(i, 3)
      end if
!ccc....face sup flux impose (U*C mol/s ou g/s)
      if (iclc(i, 3) == -1) then
         b(i) = b(i) - valclc(i, 3)*dt/bm(i)/om(i)
      end if

!CC....FACE INF
      if (vzm(i) > 0 .and. iclc(i, 4) == 1) then
         val(ib) = val(ib) + vzm(i)*dt/bm(i)/om(i)
      end if
      if (vzm(i) > 0 .and. iclc(i, 4) == -2) then
         b(i) = b(i) - vzm(i)*dt/bm(i)/om(i)*valclc(i, 4)
      end if
!ccc....face droite flux impose (U*C mol/s ou g/s)
      if (iclc(i, 4) == -1) then
         b(i) = b(i) - valclc(i, 4)*dt/bm(i)/om(i)
      end if
      if (vzm(i) < 0.) then
         val(ii) = val(ii) + vzm(i)*dt/bm(i)/om(i)
      end if

!CC....Terme B
      b(i) = b(i) - conco(i)

   end do

   return
end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                         C
!          SOULEVEMENT          C
!                         C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC...PURPOSE :
!ccc...
!CC...Modification
!ccc...
!CC....TODOLIST
!ccc....
subroutine upheaval(nc, nm, igel, zs, zsoo, ivois, pr, pro, &
                    alph, dl, def, defo, icol)
   implicit double precision(a - h, o - x, z), integer(I - N)
   implicit CHARACTER*5(y)
   dimension zs(nc, 2), zsoo(nc, 2), ivois(nm, 4), icol(nc)
   dimension pr(nm), pro(nm), dl(nc), def(nc), defo(nc), alph(nm)

   if (igel == 1) then
   do kkcol = 1, nc
   do i = 1, nm
   if (zs(kkcol, 1) .ne. -99 .and. zsoo(kkcol, 1) .ne. -99) then
   if (icol(i) == kkcol .and. ivois(i, 4) == -99 &
       .and. pr(i) > pro(i)) then
      dl(kkcol) = (pr(i) - pro(i))*alph(i)*zs(kkcol, 1)
   end if
   if (icol(i) == kkcol .and. ivois(i, 4) == -99 &
       .and. pr(i) < pro(i)) then
      dl(kkcol) = 00D+00
   end if
   def(kkcol) = dl(kkcol) + defo(kkcol)
   end if
   end do
   end do
   end if
   if (igel == 2) then
   do kkcol = 1, nc
   do i = 1, nm
   if (zs(kkcol, 1) .ne. -99 .and. zsoo(kkcol, 1) .ne. -99) then
   if (icol(i) == kkcol .and. ivois(i, 4) == -99 &
       .and. pr(i) > pro(i)) then
      dl(kkcol) = (pr(i) - pro(i))*alph(i)*zs(kkcol, 1)
   end if
   if (icol(i) == kkcol .and. ivois(i, 4) == -99 &
       .and. pr(i) < pro(i)) then
      dl(kkcol) = -(pr(i) - pro(i))*alph(i)*zs(kkcol, 1)
   end if
   def(kkcol) = dl(kkcol) + defo(kkcol)
   end if
   end do
   end do
   end if
   return
end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                         C
!          ZNS                C
!                         C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

subroutine unsaturated(pr, sw, dswdp, swresz, ans, akr, &
                       nm, akrv, rho1, g, ansun, asun)
   implicit double precision(a - h, o - x, z), integer(I - N)
   implicit CHARACTER*5(y)
   dimension pr(nm), sw(nm)
   dimension akr(nm), dswdp(nm), akrv(nm)
   dimension ansun(nm), asun(nm), swresz(nm)

   double precision :: prc, as, ans, swres, swt

!CCCC ZNS SSSSSSSSSSSSSSSS


!$omp parallel do private(prc, as, ans, swres, swt)
   do i = 1, nm
      ! Recalcul des parametres dependant de P
      prc = 0.D+00
      if (pr(i) <= 0d+00) prc = -pr(i)

      ! Calcul des paramètres de Van Genuchten
      as = asun(i) / (rho1 * g)
      ans = ansun(i)
      swres = swresz(i)

      ! Calcul de la saturation en eau
      sw(i) = swres + (1 - swres) * (1.D0 / (1.D0 + (as * prc)**ans))**((ans - 1.D0) / ans)
      if (sw(i) - swres <= 1d-05) sw(i) = swres + 1d-05

      ! Calcul de la dérivée de la saturation par rapport à la pression
      dswdp(i) = as * (ans - 1.D0) * (1.D0 - swres) * (as * prc)**(ans - 1.D0) / &
                 (1.D0 + (as * prc)**ans)**((2.D0 * ans - 1.D0) / ans)

      ! Calcul de la perméabilité relative
      swt = (sw(i) - swres) / (1.D0 - swres)
      akr(i) = (swt**0.5D0) * ((1.D0 - (1.D0 - swt**(ans / (ans - 1.D0)))**((ans - 1.D0) / ans))**2)
      akrv(i) = akr(i)
   end do
!$omp end parallel do

   return
end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                         C
!          RIVIERE            C
!                         C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC...PURPOSE :
!ccc...
!CC...Modification
!ccc...
!CC....TODOLIST
!ccc....
subroutine cdt_riviere(hriv, g, nm, z, hbot, &
                       bm, ivois, am, rho, xberg, x, icl, valcl, iclt, valclt, &
                       tempriv, aklit, aklitv, ak, akv, elit, &
                       akc, akcv, it, ita, tempo, ts)
   implicit double precision(a - h, o - x, z), integer(I - N)
   implicit CHARACTER*5(y)
   dimension z(nm), bm(nm), ivois(nm, 4), akv(nm), ak(nm), tempo(nm)
   dimension am(nm), rho(nm), x(nm), icl(nm, 4), valcl(nm, 4)
   dimension iclt(nm, 4), valclt(nm, 4)
   if (it == ita) goto 22
!       TEST COND PARTICULIERES IMPOSEES SUR CELLULE VOISINE RIVIERE
!       Bords en contact avec riviere!!!
   do i = 1, nm
   if (z(i) > hbot - bm(i) .and. z(i) <= hbot &
       .and. ivois(i, 3) == -99 .and. tempo(i) > ts) then
      icl(i, 3) = -2
      valcl(i, 3) = hriv*rho(i)*g
      iclt(i, 3) = -2
      valclt(i, 3) = tempriv
      akv(i) = bm(i)/((elit)/aklitv + &
                      (bm(i) - elit)/akcv)
      ak(i) = ((elit)*aklit + &
               (am(i) - elit)*akc)/am(i)
      if (ivois(i, 1) == -99) then
         icl(i, 1) = -1
         valcl(i, 1) = 0.D+00
      end if
   end if

   if (ivois(i, 1) == -99 &
       .and. x(i) <= xberg .and. &
       x(i) > (xberg - am(i)) .and. tempo(i) > ts) then
      icl(i, 1) = -2
      valcl(i, 1) = rho(i)*g*(hbot + hriv - z(i))
      iclt(i, 1) = -2
      valclt(i, 1) = tempriv
      ak(i) = am(i)/((elit)/aklitv + &
                     (am(i) - elit)/akcv)
      akv(i) = ((elit)*aklit + &
                (bm(i) - elit)*akc)/bm(i)
   end if

   if (z(i) > hbot - bm(i) .and. z(i) <= hbot &
       .and. ivois(i, 3) == -99 .and. tempo(i) <= ts) then
      icl(i, 3) = -1
      valcl(i, 3) = 0D+00
      iclt(i, 3) = -2
      valclt(i, 3) = tempriv
      akv(i) = bm(i)/((elit)/aklitv + &
                      (bm(i) - elit)/akcv)
      ak(i) = ((elit)*aklit + &
               (am(i) - elit)*akc)/am(i)
      if (ivois(i, 1) == -99) then
         icl(i, 1) = -1
         valcl(i, 1) = 0.D+00
      end if
   end if

   if (ivois(i, 1) == -99 &
       .and. x(i) <= xberg .and. &
       x(i) > (xberg - am(i)) .and. tempo(i) <= ts) then
      icl(i, 1) = -1
      valcl(i, 1) = 0D+00
      iclt(i, 1) = -2
      valclt(i, 1) = tempriv
      ak(i) = am(i)/((elit)/aklitv + &
                     (am(i) - elit)/akcv)
      akv(i) = ((elit)*aklit + &
                (bm(i) - elit)*akc)/bm(i)
   end if

   end do
   ita = it

22 continue
   return
end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                         C
!         DEBIT RIVIERE              C
!                         C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC...PURPOSE :
!ccc...
!CC...Modification
!ccc...
!CC....TODOLIST
!ccc....

subroutine debit_riviere(hriv, qriva, qriv, nm, z, hbot, &
                         bm, ivois, qinf, am, xberg, x, al, &
                         rug, pent, iqriv, vxp, vzp, qruis)
   implicit double precision(a - h, o - x, z), integer(I - N)
   implicit CHARACTER*5(y)
   dimension z(nm), bm(nm), ivois(nm, 4)
   dimension am(nm), x(nm), vzp(nm), vxp(nm)
!c RIVIERE
! Qriv=debit maille riviere
! qrivent= debit entrant.
! qinf= debit d infiltration riviere vers nappe
! hypo pression berge=pression dans riviere
! calcul intermediare

!CC....calcul debit de la riviere
   qinf = 0D+00
   qriv = 0D+00
   if (iqriv == 1) qriv = qriva + qruis

!CC....cellule sous rivière
   do ik = 1, nm
   if (z(ik) > (hbot - bm(ik)) &
       .and. z(ik) < hbot .and. &
       ivois(ik, 3) == -99) then
      qriv = qriv + am(ik)*1*vzp(ik)

      qinf = qinf - am(ik)*1*vzp(ik)
   end if

   if (ivois(ik, 1) == -99 &
       .and. x(ik) < xberg .and. &
       x(ik) > (xberg - am(ik))) then
!       if (hriv+hbot >= z(ik)-bm(ik)/2.) then
      qriv = qriv + bm(ik)*1*vxp(ik)

      qinf = qinf - bm(ik)*1*vxp(ik)
!       endif
   end if
   end do

!CC....DETERMINATION hauteur DANS RIVIERE hriv
   if (iqriv == 1) then
      if (qriv < qriva) qriv = qriva
      hriv = ((qriv)**(3./5.))*(rug**(-3./5.))*((al - xberg)**(-3./5.))* &
             (pent**(-3./10.))
   end if

   return
end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!                     Thermique                            C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC...PURPOSE :
!ccc...
!CC...Modification
!ccc...
!CC....TODOLIST
!ccc....
subroutine matt(val, icol_ind, irow_ptr, x, b, am, ivois, tempo, &
                rho, dt, iclt, valclt, om, nm, bll, blt, &
                cpe, cps, rhos, alanda, chlat, nmax, igel, nmax1, &
                z, bm, temp, ithec, irpth, ts, tl, &
                alandae, alandas, alandai, rhoi, cpice, sice, sw, &
                icycle, rhog, alandag, cpg, igelzns, &
                vxp, vxm, vzp, vzm, &
                ymoycondtherm, dsidtemp, ytest, ysolv)
   implicit double precision(a - h, o - x, z), integer(I - N)
   implicit CHARACTER*5(y)
   CHARACTER(3) :: ysolv
   dimension b(nm), x(nm)
   dimension ivois(nm, 4), sice(nm)
   dimension val(nmax), icol_ind(nmax), irow_ptr(nmax1)
   dimension om(nm), rho(nm), alandas(nm), rhoi(nm)
   dimension am(nm), dsidtemp(nm)
   dimension iclt(nm, 4), valclt(nm, 4), tempo(nm)
   dimension z(nm), bm(nm), rhos(nm), alanda(nm), temp(nm)
   dimension sw(nm), cps(nm)
   dimension vxm(nm), vxp(nm), vzp(nm), vzm(nm)
   CHARACTER(3) :: ytest

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                         C
!CC....Definition de la matrice pr le transport
   irow_ptr(1) = 1
   il = 0

   do i = 1, nm
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                         C
!     INITIALISATION VARIABLES GEL/DEGEL     C
!                         C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      if (icycle .ne. 1) then
         dsidtemp(i) = 0D0
         sice(i) = 0D0
         ap = 0D0
         chlat = 0D+00

      else
         ap = -abs(dsidtemp(i))
      end if
      if (sw(i) + sice(i) > 1 .and. icycle == 1) print *, "probleme sat"
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                         C
!      Termal conductivity average      C
!                         C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC... arythmetic geometric wooside
!CC
      if (igelzns == 1) then
      if (ymoycondtherm == "WOODS") then
         alanda(i) = DBLE(sqrt(alandae)*om(i)*sw(i) &
                          + sqrt(alandai)*(om(i)*sice(i)) + &
                          sqrt(alandas(i))*(1D+00 - om(i)) + &
                          sqrt(alandag)*om(i)*(1D+00 - sw(i) - sice(i)))**2
      else if (ymoycondtherm == "GEOME") then
         alanda(i) = DBLE(alandae**(om(i)*sw(i))*alandas(i)**(1 - om(i))* &
                          alandai**(om(i)*sice(i))*alandag**(om(i)*(1D+00 - sw(i) - sice(i))))
      else if (ymoycondtherm == "ARITH") then
         alanda(i) = DBLE(alandae*(om(i)*sw(i)) + alandas(i)*(1 - om(i)) &
                          + alandai*(om(i)*sice(i)) + alandag*(om(i)*(1D+00 - sw(i) - sice(i))))
      end if
      end if
      if (ymoycondtherm == "NEUMA") then
         if (temp(i) < ts) alanda(i) = 2.619D+00
         if (temp(i) >= ts) alanda(i) = 1.839D+00
      end if
      if (ymoycondtherm == "LUNAR") then
!cccc....solution avec conductivite constante sur 3 zones
!      if(temp(i) <= ts) alanda(i)=3.462696D+00
!      if(temp(i) >= ts.and.temp(i) < tl) then
!      alanda(i)=293994600.D-08
         if (temp(i) <= ts) alanda(i) = 3.4644D+00
         if (temp(i) >= ts .and. temp(i) < tl) then
            alanda(i) = 294110000.D-08
         end if
         if (temp(i) >= tl) alanda(i) = 2.4184D+00
!    if(i.ne.1.and.ap.ne.0) print*,dsidtemp(i),ap,temp(i),i
!    if(i.ne.1.and.ap.ne.0) print*,"pb",alanda(i)
      else
         if (sw(i) < 1) then
            if (ymoycondtherm == "WOODS") then
               alanda(i) = DBLE(sqrt(alandae)*om(i)*sw(i) + &
                                sqrt(alandas(i))*(1D+00 - om(i)) + &
                                sqrt(alandag)*om(i)*(1D+00 - sw(i)))**2
            else if (ymoycondtherm == "GEOME") then
               alanda(i) = DBLE(alandae**(om(i)*sw(i))*alandas(i)**(1 - om(i))* &
                                alandag**(om(i)*(1D+00 - sw(i))))
            else if (ymoycondtherm == "ARITH") then
               alanda(i) = DBLE(alandae*(om(i)*sw(i)) + alandas(i)*(1 - om(i)) &
                                + alandag*(om(i)*(1D+00 - sw(i) - sice(i))))
            end if

         end if
      end if

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                         C
!         TERME DISPERSIF          C
!                         C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      dlip = dble(bll*rho(i)*cpe*(vxp(i)**2 &
                                  + (vzm(i) + vzp(i))**2/4)**0.5 + alanda(i))
      dlim = dble(bll*rho(i)*cpe*(vxm(i)**2 &
                                  + (vzm(i) + vzp(i))**2/4)**0.5 + alanda(i))
      dtjp = dble(blt*rho(i)*cpe*(vzp(i)**2 &
                                  + (vxm(i) + vxp(i))**2/4)**0.5 + alanda(i))
      dtjm = dble(blt*rho(i)*cpe*(vzm(i)**2 &
                                  + (vxm(i) + vxp(i))**2/4)**0.5 + alanda(i))
      dlip = dble(bll*rho(i)*cpe*((vxm(i) + vxp(i))**2/4 &
                                  + (vzm(i) + vzp(i))**2/4.)**0.5 + alanda(i))
      dlim = dlip
      dtjp = dble(blt*rho(i)*cpe*((vxm(i) + vxp(i))**2/4 &
                                  + (vzm(i) + vzp(i))**2/4)**0.5 + alanda(i))
      dtjm = dtjp

      il = il + 1
      ili = il
      icol_ind(il) = i
      irow_ptr(i + 1) = irow_ptr(i) + 1

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!        INITIALISATION TERME B            C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC....chaleur latente
!ccc....modification de rho dpar rhoi chaleur latente 19-11-2014
!ccc.... la densite de la glace en degel et la densite de l eau en gel
!ccc.... TH1 utilisation rho eau degel
      if (ytest .ne. "THL" .and. igel == 0) then
         val(il) = dble(-(om(i)*sw(i)*rho(i)*cpe + &
                          om(i)*(1 - sw(i))*rhog*cpg + &
                          (1.-om(i))*rhos(i)*cps(i))/dt*irpth)
      else if (ytest .ne. "THL" .and. igel == 1) then
         val(il) = dble(-((om(i)*sw(i)*rho(i)*cpe + &
                           om(i)*sice(i)*rhoi(i)*cpice + om(i)*(1 - sw(i) - sice(i))*rhog*cpg + &
                           (1.-om(i))*rhos(i)*cps(i)) - rho(i)*om(i)*chlat*ap)/dt*irpth)
      else if (ytest .ne. "THL" .and. igel == 2) then
         val(il) = dble(-((om(i)*sw(i)*rho(i)*cpe + &
                           om(i)*sice(i)*rhoi(i)*cpice + om(i)*(1 - sw(i) - sice(i))*rhog*cpg + &
                           (1.-om(i))*rhos(i)*cps(i)) - rhoi(i)*om(i)*chlat*ap)/dt*irpth)
      end if
      if (ytest == "THL") then
         val(il) = dble(-(690360D+00 - rhoi(i)*om(i)*chlat*ap)/dt)
      end if

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!        INITIALISATION TERME B            C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      b(i) = 0.D00

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          FACE DROITE          C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      if (iclt(i, 1) == 1) then
         val(ili) = dble(val(ili) - dlip/am(i)/(x(ivois(i, 1)) - x(i)))
         il = il + 1
         icol_ind(il) = ivois(i, 1)
         val(il) = dble(dlip/am(i)/(x(ivois(i, 1)) - x(i)))
         irow_ptr(i + 1) = irow_ptr(i + 1) + 1
      end if
!ccc  temperature imposee
      if (iclt(i, 1) == -2) then
         val(ili) = dble(val(ili) - dlip/am(i)**2.D0*2.D00)
         b(i) = dble(b(i) - valclt(i, 1)*dlip/am(i)**2.D0*2.D0)
      end if
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!            FACE GAUCHE          C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      if (iclt(i, 2) == 1) then
         val(ili) = dble(val(ili) - dlim/am(i)/(x(i) - x(ivois(i, 2))))
         il = il + 1
         icol_ind(il) = ivois(i, 2)
         val(il) = dble(dlim/am(i)/(x(i) - x(ivois(i, 2))))
         irow_ptr(i + 1) = irow_ptr(i + 1) + 1
      end if
!CC....temperature imposee
      if (iclt(i, 2) == -2) then
         val(ili) = dble(val(ili) - dlim/am(i)**2.D00*2.D00)
         b(i) = dble(b(i) - valclt(i, 2)*dlim/am(i)**2.D00*2.D00)
      end if
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          FACE HAUTE            C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      if (iclt(i, 3) == 1) then
         val(ili) = dble(val(ili) - dtjp/bm(i)/(z(ivois(i, 3)) - z(i)))
         il = il + 1
         icol_ind(il) = ivois(i, 3)
         val(il) = dble(dtjp/bm(i)/(z(ivois(i, 3)) - z(i)))
         irow_ptr(i + 1) = irow_ptr(i + 1) + 1
      end if
!CC....temperature imposee
      if (iclt(i, 3) == -2) then
         val(ili) = dble(val(ili) - dtjp/bm(i)**2.D00*2.D00)
         b(i) = dble(b(i) - valclt(i, 3)*dtjp/bm(i)**2.D00*2.D00)
      end if
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!         FACE BASSE              C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      if (iclt(i, 4) == 1) then
         val(ili) = dble(val(ili) - dtjm/bm(i)/(z(i) - z(ivois(i, 4))))
         il = il + 1
         icol_ind(il) = ivois(i, 4)
         val(il) = dble(dtjm/bm(i)/(z(i) - z(ivois(i, 4))))
         irow_ptr(i + 1) = irow_ptr(i + 1) + 1
      end if
!CC....temperature imposee
      if (iclt(i, 4) == -2) then
         val(ili) = dble(val(ili) - dtjm/bm(i)**2.D00*2.D00)
         b(i) = dble(b(i) - valclt(i, 4)*dtjm/bm(i)**2.D00*2.D00)
      end if
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                         C
!     ORDONNER LES TABLEAUX icol_IND ET VAL    C
!                         C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      ideb = irow_ptr(i)
      ifin = irow_ptr(i + 1) - 1

      if (ifin == ideb) goto 11
      if (ysolv == "BIC" .or. ysolv == "CGS") then
10       continue
         do ii = ideb, ifin - 1
         if (icol_ind(ii) > icol_ind(ii + 1)) then
            idum = icol_ind(ii + 1)
            dum = val(ii + 1)
            icol_ind(ii + 1) = icol_ind(ii)
            val(ii + 1) = dble(val(ii))
            icol_ind(ii) = idum
            val(ii) = dum
            goto 10
         end if
         end do
      end if
!     if (ysolv == "LIB") then
!     do ii=ideb,ifin-1
!     if(icol_ind(ii) > icol_ind(ii+1)) then
!     idum=icol_ind(ii+1)
!     dum=val(ii+1)
!     icol_ind(ii+1)=icol_ind(ii)
!     val(ii+1)=val(ii)
!     icol_ind(ii)=idum
!     val(ii)=dum
!     endif
!     enddo
!     do ii=ideb,ifin-1
!     if(icol_ind(ii) == ik) then
!     idum=icol_ind(ifin)
!     dum=val(ifin)
!     icol_ind(ifin)=icol_ind(ii)
!     val(ifin)=val(ii)
!     icol_ind(ii)=idum
!     val(ii)=dum
!     endif
!     enddo
!     endif
11    continue
!CC....REPERAGE DES VOISINS DS irow_ptr
      ii = 0
      id = 0
      ig = 0
      ih = 0
      ib = 0
      do ij = irow_ptr(i), irow_ptr(i + 1) - 1
         if (icol_ind(ij) == i) ii = ij
         if (icol_ind(ij) == ivois(i, 1)) id = ij
         if (icol_ind(ij) == ivois(i, 2)) ig = ij
         if (icol_ind(ij) == ivois(i, 3)) ih = ij
         if (icol_ind(ij) == ivois(i, 4)) ib = ij
      end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                         C
!      TERMES ADVECTIFS                C
!                         C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC....VITESSE EN X A TRAVERS LES FACES
!ccc....vxm(i) face gauche
!ccc....vxp(i) face droite
!ccc....vzm(i) face du dessous
!ccc....vzp(i) face du dessus
      if (ithec == 0) then
         vxp(i) = 00D+00
         vxm(i) = 00D+00
         vzm(i) = 00D+00
         vzp(i) = 00D+00
      end if
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          FACE DROITE          C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      if (vxp(i) > 0) then
         val(ii) = dble(val(ii) - rho(i)*cpe*vxp(i)/am(i))
      end if
      if (iclt(i, 1) == 1 .and. vxp(i) < 0) then
         val(id) = dble(val(id) - rho(ivois(i, 1))*cpe*vxp(i)/am(i))
      end if
!CC....si face droite imposee et vxp(i)<0!!!!
      if (iclt(i, 1) == -2 .and. vxp(i) < 0) then
         b(i) = dble(b(i) + rho(i)*cpe*vxp(i)/am(i)*valclt(i, 1))
         if (abs(vxp(i) - vxm(i)) > abs(vxm(i))/1D9) then
            print *, 'pb de vitesse x maille ='
            print *, i, vxp(i), vxm(i)
            print *, 'verifiez votre etat initial'
            stop
         end if
      end if
!cccc....face droite flux impose (>0 U*C mol/s ou g/s)
      if (iclt(i, 1) == -1) then
         b(i) = dble(b(i) - valclt(i, 1)/am(i))
      end if

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!            FACE GAUCHE          C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      if (vxm(i) > 0 .and. iclt(i, 2) == 1) then
         val(ig) = dble(val(ig) + rho(ivois(i, 2))*cpe*vxm(i)/am(i))
      end if
!CC....si face gauche imposee et vxm(i)>0!!!!
      if (vxm(i) > 0 .and. iclt(i, 2) == -2) then
         b(i) = dble(b(i) - rho(i)*cpe*vxm(i)/am(i)*valclt(i, 2))
         if (abs(vxp(i) - vxm(i)) > abs(vxm(i))/1D9) then
            print *, 'pb de vitesse x'
            print *, i, vxp(i), vxm(i)
            print *, 'verifiez votre etat initial'
            stop
         end if
      end if
      if (vxm(i) < 0) then
         val(ii) = dble(val(ii) + rho(i)*cpe*vxm(i)/am(i))
      end if
!CC....face gauche flux impose (U*C mol/s ou g/s)
      if (iclt(i, 2) == -1) then
         b(i) = dble(b(i) - valclt(i, 2)/am(i))
      end if

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          FACE SUP                C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      if (vzp(i) > 0) then
         val(ii) = dble(val(ii) - rho(i)*cpe*vzp(i)/bm(i))
      end if
      if (vzp(i) < 0 .and. iclt(i, 3) == 1) then
         val(ih) = dble(val(ih) - rho(ivois(i, 3))*cpe*vzp(i)/bm(i))
      end if
      if (vzp(i) < 0 .and. iclt(i, 3) == -2) then
      if (abs(vzp(i) - vzm(i)) >= abs(vzm(i))/1D9) then
         print *, 'pb de vitesse z', vzp(i), vzm(i), i, abs(vzp(i) - vzm(i)), dt
         print *, 'verifiez votre etat initial sup'
         stop
      end if
      b(i) = dble(b(i) + rho(i)*cpe*vzp(i)/bm(i)*valclt(i, 3))
      end if
!CC....face sup flux impose (U*C mol/s ou g/s)
      if (iclt(i, 3) == -1) then
         b(i) = dble(b(i) - valclt(i, 3)/bm(i))

      end if

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          FACE INF                C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      if (vzm(i) > 0 .and. iclt(i, 4) == 1) then
         val(ib) = dble(val(ib) + rho(ivois(i, 4))*cpe*vzm(i)/bm(i))
      end if
      if (vzm(i) > 0 .and. iclt(i, 4) == -2) then
         b(i) = dble(b(i) - rho(i)*cpe*vzm(i)/bm(i)*valclt(i, 4))
         if (abs(vzp(i) - vzm(i)) >= abs(vzm(i))/1D9) then
            print *, 'pb de vitesse z', abs(vzp(i) - vzm(i)), i
            print *, 'verifiez votre etat initial'
            stop
         end if
      end if
!CC....face basse flux impose (U*C mol/s ou g/s)
      if (iclt(i, 4) == -1) then
         b(i) = dble(b(i) - valclt(i, 4)/bm(i))
      end if
      if (vzm(i) < 0.) then
         val(ii) = dble(val(ii) + rho(i)*cpe*vzm(i)/bm(i))
      end if
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                         C
!      Terme B                C
!                         C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC....Chaleur latente+capacite calorifique
!CC....Calcul differents en fonction des tests INTERFOST
!cccc....THL Lunardi
!cccc....TH1 Neuman
!cccc....TH2/TH3
!CC....chaleur latente
!ccc....modification de rho dpar rhoi chaleur latente 19-11-2014
!ccc.... la densite de la glace en degel et la densite de l eau en gel.
      if (ytest .ne. "THL" .and. igel == 0) then
         b(i) = dble(b(i) - tempo(i)*(om(i)*sw(i)*rho(i)*cpe + &
                                      om(i)*(1 - sw(i))*rhog*cpg + &
                                      (1.-om(i))*rhos(i)*cps(i))/dt*irpth)
      else if (ytest .ne. "THL" .and. igel == 2) then
         b(i) = dble(b(i) - tempo(i)*((om(i)*sw(i)*rho(i)*cpe + &
                                       om(i)*sice(i)*rhoi(i)*cpice + om(i)*(1 - sw(i) - sice(i))*rhog*cpg + &
                                       (1.-om(i))*rhos(i)*cps(i)) - rhoi(i)*om(i)*chlat*ap)/dt*irpth)
      else if (ytest .ne. "THL" .and. igel == 1) then
         b(i) = dble(b(i) - tempo(i)*((om(i)*sw(i)*rho(i)*cpe + &
                                       om(i)*sice(i)*rhoi(i)*cpice + om(i)*(1 - sw(i) - sice(i))*rhog*cpg + &
                                       (1.-om(i))*rhos(i)*cps(i)) - rho(i)*om(i)*chlat*ap)/dt*irpth)
      end if
      if (ytest == "THL") then
         b(i) = dble(b(i) - tempo(i)*(690360D+00 - rhoi(i)*om(i)*chlat*ap)/dt)
      end if

   end do

   return
end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!                INFILTRATION                              C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC...PURPOSE :
!ccc...
!CC...Modification
!ccc...
!CC....TODOLIST
!ccc....

!       subroutine infiltration(pr,nm,n,icol,nc,om,ss,bm,rho,pore,
!     &am,valcl,icl,qtpore,ruis,ivois,z,g,dt,xberg,x,qre,qruis)
!CC...PURPOSE :
!ccc...
!CC...Modification
!ccc...
!CC....TODOLIST
!ccc....
!       implicit double precision(a-h,o-x,z)
!       implicit integer (i-n)
!       dimension pr(n),icol(n),om(n),ss(n),bm(n),x(n)
!       dimension rho(n),pore(n),qtpore(nc),am(n)
!       dimension valcl(n,4),icl(n,4),ruis(nc),ivois(n,4),z(n)
!       do j=1,nc
!       qtpore(j)=0.
!       do i=1,nm
!       pore(i)=0.
!C test colonne
!       if(icol(i) == j) then
!test charge < top maille
!       if(pr(i)/(rho(i)*g)+z(i) < z(i)+bm(i)/2) then
!       pore(i)=(bm(i)/2.D00-pr(i)/(rho(i)*g))*om(i)*am(i)
!       qtpore(j)=qtpore(j)+pore(i)
!       endif
!       endif
!       enddo
!      enddo
!       do j=1,nc
!       ruis(j)=0
!       do i=1,nm
!          if(icol(i) == j) then
!      if (ivois(i,3) == -99.and.icl(i,3) == -1.and.
! c    &valcl(i,3) > 0.and.x(i) < xberg) then
!     if (qtpore(j) < a*am(i)*dt) then
!       a=qre
!       ruis(j)=qre*am(i)*dt-qtpore(j)
!       valcl(i,3)=qtpore(j)/am(i)/dt
!       qruis=qruis+ruis(j)/dt
!       endif
!       endif
!       if (ivois(i,3) == -99.and.x(i) < xberg) then
!       if(pr(i)/(rho(i)*g)+z(i) > z(i)+bm(i)/2) then
!       ruis(j)=ruis(j)+pr(i)/(rho(i)*g)-bm(i)/2
!       endif
!       endif
!C fin test colonne
!       endif
!       enddo
!       enddo
!       return
!       end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!                     INTERPOLATION                            C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC...PURPOSE :
!ccc...Calcul de la position d'une isovaleur d'une variable d'état
!CC...Modification
!ccc...
!CC....TODOLIST
!ccc....
subroutine interpol(nc, zo, nm, icol, temp, ivois, igel, &
                    bm, z, ncmax, valclt, to, iclt, topo)
   implicit double precision(a - h, o - x, z), integer(I - N)
   implicit CHARACTER*5(y)
   dimension z(nm)
   dimension ivois(nm, 4), valclt(nm, 4)
   dimension zo(ncmax, 2)
   dimension bm(nm), temp(nm), topo(ncmax)
   dimension icol(nm), iclt(nm, 4)
!CCC position permafrost,ice
   do kkcol = 1, nc
      zo(kkcol, 1) = -99D+00
      zo(kkcol, 2) = -99D+00
!CC PASSAGE SUR TOUT LES MAILLES
      do i = 1, nm
! TEST COLONNE
         if (icol(i) == kkcol) then
!DEGEL
!CC PAS DE VOISIN EN HAUT CC CAS DEGEL SUP>0 I<0
            if (ivois(i, 3) == -99 .and. igel == 2 &
                .and. temp(i) < to .and. valclt(i, 3) > to .and. &
                iclt(i, 3) == -2) then
               zo(kkcol, 1) = z(i) + &
                              (bm(i)/2)/(valclt(i, 3) - temp(i)) &
                              *(to - temp(i))
               if (zo(kkcol, 1) > topo(kkcol)) then
                  zo(kkcol, 1) = -99
               end if
            end if

!CC PAS DE VOISIN EN HAUT CC CAS GEL SUP>0 I<0
            if (ivois(i, 3) == -99 .and. igel == 1 .and. &
                temp(i) > to .and. &
                valclt(i, 3) < to .and. &
                iclt(i, 3) == -2) then
               zo(kkcol, 1) = z(i) + &
                              (bm(i)/2)/(valclt(i, 3) - temp(i)) &
                              *(to - temp(i))
               if (zo(kkcol, 1) > topo(kkcol)) then
                  zo(kkcol, 1) = -99
               end if
            end if

!CPAS DE VOISIN EN BAS Degel
            if (ivois(i, 4) == -99 .and. zo(kkcol, 2) == -99 .and. &
                igel == 2 .and. temp(i) > to .and. &
                iclt(i, 4) == -2) then
               zo(kkcol, 2) = z(i) + &
                              (bm(i)/2)/(valclt(i, 4) - temp(i)) &
                              *(to - temp(i))
               if (zo(kkcol, 2) > z(i) + bm(i)/2 .or. &
                   zo(kkcol, 2) < z(i) - bm(i)/2) then
                  zo(kkcol, 2) = -99
               end if

            end if
!CCC cellule interpolation par le bas DEGEL
            if (igel == 2 .and. zo(kkcol, 2) == -99 .and. &
                temp(i) <= to .and. temp(ivois(i, 4)) > to .and. &
                ivois(i, 4) .ne. -99) then
               zo(kkcol, 2) = z(i) + &
                              ((z(i) - z(ivois(i, 4)))/(temp(i) - temp(ivois(i, 4))) &
                               *(to - temp(i)))
               if (zo(kkcol, 2) > z(i) + bm(i)/2 .or. &
                   zo(kkcol, 2) < z(i) - bm(i)/2) then
                  zo(kkcol, 2) = -99
               end if
            end if
!CC ISOTHERME 1
!ccc cellule milieu modele VOISIN DU HAUT EXISTE  MILIEU MODELE CCCCCCCC
!CCCC TEMPERATURE VOISIN SUP >0 CAS DEGEL
            if (ivois(i, 3) .ne. -99 .and. temp(i) <= to .and. &
                temp(ivois(i, 3)) > to .and. igel == 2) then
               zo(kkcol, 1) = z(ivois(i, 3)) + &
                              ((z(ivois(i, 3)) - z(i))/(temp(ivois(i, 3)) - temp(i)) &
                               *(to - temp(ivois(i, 3))))
!CC FIN DE TEST VOISIN SUP >0
            end if
!C ISOTHERME  1 Gel
!CCCC TEMP VOISIN SUP<0 CAS DEGEL SOUS PERMAFROST / CAS GEL
            if (igel == 1 .and. ivois(i, 3) .ne. -99 .and. temp(i) >= to .and. &
                temp(ivois(i, 3)) < to) then
               zo(kkcol, 1) = z(ivois(i, 3)) + &
                              ((z(ivois(i, 3)) - z(i))/(temp(ivois(i, 3)) - temp(i)) &
                               *(to - temp(ivois(i, 3))))
            end if
!C ISOTHERME  2 degel
            if (igel == 2 .and. ivois(i, 3) .ne. -99 .and. temp(i) >= to .and. &
                temp(ivois(i, 3)) < to) then
               zo(kkcol, 2) = z(ivois(i, 3)) + &
                              ((z(ivois(i, 3)) - z(i))/(temp(ivois(i, 3)) - temp(i)) &
                               *(to - temp(ivois(i, 3))))
            end if
!CCC FIN DE test colonne
         end if
      end do
   end do

   return
end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!                     icesatperm                          C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC...PURPOSE :
!CC...saturation en glace en fonction de la temperature
!CC...permeabilite  relative en fonction de la temperature
!CC...les saturations sont caluclees en fonction de la saturation en eau initiale dans le cas d'un sol non sature
!CC...actuellement une seule fonction pourra etre utilisee a l'ensemeble du domaines
!CC...Modification
!ccc...
!CC....TODOLIST
!ccc....
subroutine icesatperm(nm, tl, ts, akr, akrv, dk, temp, &
                      dsidtemp, igelzns, ytypakrice, &
                      sice, om, sw, swressi, ytypsice, &
                      cimp, tempo, omega, siceo)
   implicit double precision(a - h, o - x, z), integer(I - N)
   implicit CHARACTER*5(y)
   dimension akr(nm), akrv(nm), temp(nm), tempo(nm)
   dimension dsidtemp(nm)
   dimension om(nm)
   dimension sice(nm), sw(nm), siceo(nm)
   CHARACTER(5) ::  ytypakrice, ytypsice
   do i = 1, nm
!     ---------------------------------------------
!       pas de gel
!     ---------------------------------------------
      sice(i) = 0D+00
      akr(i) = 1D+00
      akrv(i) = 1D+00
      dsidtemp(i) = 0D+00
!     ---------------------------------------------
!     EXPONENTIAL MODEL OF LUNARDINI
!........sice VS. TEMP
!     ---------------------------------------------
      if (ytypsice == "EXPON") then
!       if(tempo(i) < tl.and.temp(i) > tl) then
!     dsidtemp(i)=((1.D0-swressi)*(2.D0*(temp(i)
!     &/(cimp*cimp)))*siexp)
!     endif
!     if(tempo(i) > tl.and.temp(i) < tl) then
!     dsidtemp(i)=((1.D0-swressi)*(2.D0*(temp(i)
!     &/(cimp*cimp)))*siexp)
!     endif
         if (temp(i) > 0) then
            sice(i) = 0.D0
         else
            tld = -0.00
            siexp = DEXP(-((temp(i) - tld)/cimp)**2)
            sice(i) = (1.D0 - ((1.D0 - swressi) &
                               *siexp + swressi))

!..........derivee de la saturation en glace en fonction de la temperature
            dsidtemp(i) = ((1.D0 - swressi)*(2.D0*(temp(i) &
                                                   /(cimp*cimp)))*siexp)

         end if

!     ---------------------------------------------
!     POWER MODEL Anderson and Tice
!........sice VS. TEMP
!     ---------------------------------------------
      else if (ytypsice == "POWER") then
         if (temp(i) > tl) then
            sice(i) = 0.D0
         else if (temp(i) > ts) then
            sice(i) = 1.D0 - swressi
         else
            tld = 0
            xi = 0.3
            beta = 3
            siexp = DEXP(xi*((tl - temp(i))/(temp(i) - ts))**beta)
            sice(i) = (1.D0 - ((1.D0 - swressi) &
                               *siexp + swressi))

!..........derivee de la saturation en glace en fonction de la temperature
            dsidtemp(i) = (siceo(i) - sice(i))/(tempo(i) - temp(i))

         end if

!     ---------------------------------------------
!     LINEAR FUNCTION
!........sice VS. TEMP
!     ---------------------------------------------              !
      else if (ytypsice == "LINEA") then

!     if(tempo(i) < tl.and.temp(i) > tl) then
!     dsidtemp(i)=(1.D0-swressi)/(ts-tl)
!     endif
!     if(tempo(i) > tl.and.temp(i) < tl) then
!     dsidtemp(i)=(1.D0-swressi)/(ts-tl)
!     endif
         if (temp(i) < ts) then
            sice(i) = 1.D0 - swressi
            dsidtemp(i) = 0.D0
         else if (temp(i) > tl) then
            sice(i) = 0D+00
            dsidtemp(i) = 0.D0
         else if (temp(i) >= ts .and. temp(i) <= tl) then
            sice(i) = (1.D0 - swressi)/(ts - tl)*(temp(i) - tl)
            dsidtemp(i) = (1.D0 - swressi)/(ts - tl)
         end if

!--------------------------
!     USER-DEFINED MODEL
!     -----------------
      else if (ytypsice == "UDEF") then
!     *  ASSIGN  VALUES TO sice(i) and  dsidtemp(i)
!     if(tempo(i) < ts.and.temp(i) > tl) then
!     dsidtemp(i)=(siceo(i)-sice(i))/dt
!     endif
!     if(tempo(i) > tl.and.temp(i) < ts) then
!     dsidtemp(i)=(siceo(i)-sice(i))/dt
!     endif

         if (temp(i) <= ts) then
            sice(i) = 1D00
         end if
         if (temp(i) >= tl) then
            sice(i) = 0D+00
         end if

         if (temp(i) >= -0.2 .and. temp(i) <= -0.1) then
            sice(i) = -1*temp(i) + 0.8
            dsidtemp(i) = -1
         end if
         if (temp(i) >= -0.1 .and. temp(i) <= 0) then
            sice(i) = -9*temp(i)
            dsidtemp(i) = -9

         end if

      end if
!ccc   COMPLEMENT SATURE
      if (igelzns .ne. 1) then
!ccc         ---------------------------------------------
!ccc          Fonction exponentielle brutale en fonction de la temperature (Lundin (19??)
!ccc         ---------------------------------------------
         if (ytypakrice == "IMPED") then
!ccc...........permeabilite relative calcule en fonction de la teneur en glace
!CC! dans le cas non sature on aura akr=(10**(-1*om*(SI/(SI+SW)))) car SI+SW<1
            akr(i) = 1.D1**(-1*omega*om(i)*sice(i))
!ccc         limitation de la permeabilite relative pour eviter les valeurs nulles
!c        if(akr(i) <= dk) akr(i)=dk
!ccc      ---------------------------------------------
!ccc         Fonction lineaire brutale en fonction de la temperature
!ccc      ---------------------------------------------
         elseif (ytypakrice == "LINET") then
!ccc..........permeabilite relative calcule en fonction de la temperature (akr VS. TEMP)
            if (temp(i) < ts) akr(i) = dk
            if (temp(i) <= tl .and. temp(i) >= ts) then
!          akr(i)=(((dk-1.D0)*(temp(i)-tl)/(tl-ts))+1.D0)
               akr(i) = 10**((log(dk))/(tl - ts)*(tl - temp(i)))
            end if
            if (temp(i) > tl) akr(i) = 1
!ccc      ---------------------------------------------
!ccc        fonction lineaire en fonction de la saturation
!ccc      ---------------------------------------------
         elseif (ytypakrice == "LINES") then
            slopek = -(1.D0 - dk)/(1.D0 - swressi)
            if (sw(i) <= swressi) then
               akr(i) = dk
            else
               akr(i) = dk + slopek*(sw(i) - swressi)
            end if
!ccc      ---------------------------------------------
!ccc        UTILISATEUR DEFINIES LA FONCTION
!ccc      ---------------------------------------------
         elseif (ytypakrice == "UDEF") then
!C    ecrire fonction akr(i) et akrv(i)

            if (temp(i) <= ts) then
               akr(i) = dk
            end if
            if (temp(i) >= tl) then
               akr(i) = 1
            end if

            if (temp(i) >= -0.2 .and. temp(i) <= -0.1) then
               akr(i) = (1D-5 - dk)/(-0.1 - ts)*temp(i) + (-0.1*dk - ts*1D-5)/ &
                        (-0.1 - ts)
            end if
            if (temp(i) >= -0.1 .and. temp(i) <= 0) then
               akr(i) = (1 - 1D-5)/(0.1)*temp(i) + 1D0

            end if
         end if
!c       if(akr(i) <= dk) akr(i)=dk
      end if
!C permeabilite verticale = permeabilite horizontale
      akrv(i) = akr(i)
   end do
!cc
   return
end

!
!     SUBROUTINE       PERMEABILITE NON SATUREE GEL
! --- PURPOSE :
!CCC SUBROUTINE DATANT ET NON VERIFIEE A NE PAS UTILISEE
!
SUBROUTINE unsatpermice(nm, akr, akrv, dk, &
                        ytypakrice, ans, &
                        sice, om, sw, swressi, rlamb)
   implicit double precision(a - h, o - x, z), integer(I - N)
   implicit CHARACTER*5(y)
   dimension akr(nm), akrv(nm)
   dimension om(nm)
   dimension sice(nm), sw(nm)
   Do i = 1, nm
   if (ytypakrice == "VGEN") then
!     --------------------------------
!     MODEL OF VANCENUCHTEN
!     --------------------------------
      swt = (sw(i) - swressi)/(1.-swressi)
      akr(i) = ((swt**0.5)*((1 - (1 - swt**(ans/(ans - 1))) &
                             **((ans - 1)/ans))**2))
   end if
   if (ytypakrice == "BCOR") then
!     --------------------------------
!     MODEL OF BROOKS AND COREY (1964)
!     --------------------------------
      swrel = (sw(i) - swressi)/(1D+00 - swressi)
!........RELATIVE PERMEABILITY AS A FUNCTION OF SATURATION (akr VS. SW)
!        (CALCULATED ONLY WHEN IUNSAT=2)
      akr(i) = (swrel**(3D0 + 2D0/RLAMB))
   end if
   if (ytypakrice == "PLIN") then
!     ----------------------
!      LINEAR MODEL
!     ----------------------
!........akr VS. SW
      if (sw(i) <= swressi) then
         akr(i) = dk
      else
         SLOPEK = -(1.D0 - dk)/(1.D0 - swressi)
         akr(i) = (1.D0 + SLOPEK*(1.D0 - sw(i)))
         if (akr(i) <= dk) akr(i) = dk
      end if

      if (ytypakrice == "IMPE") then
!     ----------------------------------------------
!     AD-HOC EXPONENTIAL EXPRESSION OF LUNDIN (19??)
!     ----------------------------------------------
         akr(i) = (10**(-1*om(i)*omega*(sice(i)/(sice(i) + sw(i)))))
         if (akr(i) <= dk) akr(i) = dk

      end if
      if (ytypakrice == "UDEF") then
!     ------------------
         dkswr = dk/swressi
         SLKRSW = (1.D0 - dkswr)/(1.D0 - swressi)
         akr(i) = sw(i)*(dkswr + SLKRSW)*(sw(i) - swressi)
         if (akr(i) <= dk) akr(i) = dk
      end if
      akrv(i) = akr(i)
   end if
   END DO
   RETURN
END

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!                Subroutine Biggridice                        C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC...PURPOSE :
!ccc...Calcul la teneur en glace sice et permeabilite relative akr akrv en fonction de la position des isothermes dans les mailles
!ccc...fonction linear
!ccc...1) Calcul de la portion de glace (<ts), front de gel (ts<t<tl), eau liquide (>tl)
!ccc...2) Calcul de la temperature du milieu du front de gel tifrt et de la permeabilite relative akrifrt
!ccc...3) Moyenne de saturation en glace sur la maille
!ccc...4) Moyenne de la permeabilite en glace sur la maille
!ccc...ATTENTION NE FONCTIONNE PAS EN NON SATURE
!CC...Modification
!ccc...
!CC....TODOLIST
!ccc....
subroutine biggridice(nm, tl, ts, akr, akrv, dk, temp, igel, &
                      nc, z, zl, zs, bm, dsidtemp, valclt, ivois, &
                      sice, swressi, siceo, tempo)
   implicit double precision(a - h, o - x, z), integer(I - N)
   implicit CHARACTER*5(y)
   dimension akr(nm), akrv(nm), temp(nm)
   dimension z(nm)
   dimension tempo(nm)
   dimension zs(nc, 2), zl(nc, 2), icol(nm), bm(nm), siceo(nm)
   dimension ivois(nm, 4), valclt(nm, 4), sice(nm), dsidtemp(nm)
   do i = 1, nm
      sice(i) = 0D+00
      akr(i) = 1D+00
      akrv(i) = 1D+00
      dsidtemp(i) = 0D+00
!C POUR ETRE SURE QUI NE SE BALLADE PAS AVEC DES VALEURS NIMP
      if (temp(i) < tl .and. temp(i) > ts) then
         akr(i) = (((dk - 1.D0)*temp(i)/ts) + 1.D0)
         sice(i) = ((1.D0 - swressi)/(ts - tl)*(temp(i) - tl))
         dsidtemp(i) = ((1.D0 - swressi)/(ts - tl))
      end if
      if (temp(i) <= ts) then
         akr(i) = dk
         sice(i) = (1D+00 - swressi)
         dsidtemp(i) = 0D+00
      end if
      if (temp(i) >= tl) then
         akr(i) = 1D+00
         sice(i) = 0D+00
         dsidtemp(i) = 0D+00
      end if

      do j = 1, nc
      if (icol(i) == j) then
!CC Gel maille audessus liquidus
!CCCCCCCCCCCCCCC GEL attention codé pour une seule aire interface frozen unfrozen
         if (zs(j, 1) == -99 .and. zl(j, 1) .ne. -99 .and. igel == 1) then
         if (zl(j, 1) < z(i) + bm(i)/2 &
             .and. zl(j, 1) > z(i) - bm(i)/2) then
         if (ivois(i, 3) == -99) then
            tifrt = ((temp(i) - valclt(i, 3)) &
                     /(-bm(i)/2D00) &
                     *(bm(i)/2D00 + zl(j, 1) - z(i))/2.D00 + temp(i))
         else
            tifrt = ((temp(i) - temp(ivois(i, 3)))/ &
                     (z(i) - z(ivois(i, 3))) &
                     *(bm(i)/2.D00 + zl(j, 1) - z(i))/2.D00 + temp(i))
         end if
         akrifrt = ((1.D00 - dk)/(tl - ts)*tifrt + (tl*dk - ts*1)/(tl - ts))
         sifrt = ((1 - swressi)/(ts - tl)*(tifrt - tl))
         akr(i) = (((z(i) + bm(i)/2.D00 - zl(j, 1))*akrifrt + &
                    (zl(j, 1) - z(i) + bm(i)/2.D00)*1)/bm(i))
         sice(i) = ((sifrt*(z(i) + bm(i)/2.D00 - zl(j, 1)))/bm(i))
         if (tempo(i) - temp(i) .ne. 0) then
            dsidtemp(i) = (abs((siceo(i) - sice(i))/(tempo(i) - temp(i))))
         end if
         end if
         end if
!CCCCCCCCCCCCCCfin pas d'isotherme solidus'

         if (zs(j, 1) .ne. -99 .and. zl(j, 1) .ne. -99 .and. igel == 1) then
! maille avec un liquidus mais pas de solidus
            if (zs(j, 1) > z(i) + bm(i)/2 .and. zl(j, 1) < z(i) + bm(i)/2 &
                .and. zl(j, 1) > z(i) - bm(i)/2) then
            if (ivois(i, 3) == -99) then
               tifrt = ((temp(i) - valclt(i, 3))/ &
                        (-bm(i)/2D00) &
                        *(bm(i)/2.D00 + zl(j, 1) - z(i))/2.D00 + temp(i))
            else
               tifrt = ((temp(i) - temp(ivois(i, 3)))/(z(i) - z(ivois(i, 3))) &
                        *(bm(i)/2.D00 + zl(j, 1) - z(i))/2.D00 + temp(i))
            end if
            akrifrt = ((1.D00 - dk)/(tl - ts)*tifrt + (tl*dk - ts*1)/(tl - ts))
            sifrt = ((1 - swressi)/(ts - tl)*(tifrt - tl))
            akr(i) = (((z(i) + bm(i)/2.D00 - zl(j, 1))*akrifrt + (zl(j, 1) - &
                                                                  z(i) + bm(i)/2)*1)/bm(i))
            sice(i) = ((sifrt*(z(i) + bm(i)/2.D00 - zl(j, 1)))/bm(i))
            if (tempo(i) - temp(i) .ne. 0) then
               dsidtemp(i) = (abs((siceo(i) - sice(i))/(tempo(i) - temp(i))))
            end if

            end if

!CCCfin 1ier cas
! maille avec solidus et liquidus
            if (zs(j, 1) < z(i) + bm(i)/2 .and. zs(j, 1) > z(i) - bm(i)/2 &
                .and. zl(j, 1) < z(i) + bm(i)/2 .and. zl(j, 1) > z(i) - bm(i)/2) &
               then
            if (ivois(i, 3) == -99) then
               tifrt = ((temp(i) - valclt(i, 3)) &
                        /(-bm(i)/2.D00) &
                        *((zs(j, 1) + zl(j, 1))/2.D00 - z(i)) + temp(i))
            else
               tifrt = ((temp(i) - temp(ivois(i, 3)))/(z(i) - z(ivois(i, 3))) &
                        *((zs(j, 1) + zl(j, 1))/2.D00 - z(i)) + temp(i))
            end if
            akrifrt = ((1.D00 - dk)/(tl - ts)*tifrt + (tl*dk - ts*1)/(tl - ts))
            sifrt = ((1 - swressi)/(ts - tl)*(tifrt - tl))
            akr(i) = (((z(i) + bm(i)/2.D00 - zs(j, 1))*dk + (zs(j, 1) - zl(j, 1)) &
                       *akrifrt + (zl(j, 1) - z(i) + bm(i)/2)*1)/bm(i))
            sice(i) = ((sifrt*(zs(j, 1) - zl(j, 1)) + &
                        (z(i) + bm(i)/2.D00 - zs(j, 1))*1)/bm(i))
            if (tempo(i) - temp(i) .ne. 0) then
               dsidtemp(i) = (abs((siceo(i) - sice(i))/(tempo(i) - temp(i))))
            end if
            end if
!cccc fin 2 eme cas

!C maille avec solidus
            if (zs(j, 1) < z(i) + bm(i)/2 .and. &
                zs(j, 1) > z(i) - bm(i)/2 .and. &
                zl(j, 1) < z(i) - bm(i)/2) then
            if (ivois(i, 3) == -99) then
               tifrt = ((temp(i) - valclt(i, 3)) &
                        /(-bm(i)/2.D00) &
                        *((zs(j, 1) + z(i) - bm(i)/2)/2.D00 - z(i)) + temp(i))
            else
               tifrt = ((temp(i) - temp(ivois(i, 3)))/(z(i) - z(ivois(i, 3))) &
                        *((zs(j, 1) + z(i) - bm(i)/2)/2.D00 - z(i)) + temp(i))
            end if
            akrifrt = ((1.D00 - dk)/(tl - ts)*tifrt + (tl*dk - ts*1)/(tl - ts))
            sifrt = ((1 - swressi)/(ts - tl)*(tifrt - tl))
            akr(i) = (((z(i) + bm(i)/2.D00 - zs(j, 1))*dk + (zs(j, 1) - &
                                                             z(i) + bm(i)/2)*akrifrt)/bm(i))
            sice(i) = ((sifrt*(zs(j, 1) - z(i) + bm(i)/2) + &
                        (z(i) + bm(i)/2.D00 - zs(j, 1))*1)/bm(i))
            dsidtemp(i) = (abs((siceo(i) - sice(i))/(tempo(i) - temp(i))))
            end if
!cc fin cas 3
         end if
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC fin gel
!CCCCCCCCCCCCCCC DEGEL attention codé pour deux aires interfaces frozen/unfrozen
!cc interface SUBPERMAFROST
         if (zl(j, 2) == -99 .and. zs(j, 2) .ne. -99 .and. igel == 2) then
! pas de liquidus
            if (zs(j, 2) < z(i) + bm(i)/2 .and. &
                zs(j, 2) > z(i) - bm(i)/2) then
            if (ivois(i, 3) == -99) then
               tifrt = ((temp(i) - valclt(i, 3))/(-bm(i)/2.D00) &
                        *((zs(j, 2) - z(i) - bm(i)/2)/2) + temp(i))
            else
               tifrt = ((temp(i) - temp(ivois(i, 3)))/(z(i) - z(ivois(i, 3))) &
                        *((zs(j, 2) - z(i) - bm(i)/2)/2) + temp(i))
            end if
            akrifrt = ((1.D00 - dk)/(tl - ts)*tifrt + (tl*dk - ts*1)/(tl - ts))
            sifrt = ((1 - swressi)/(ts - tl)*(tifrt - tl))
            akr(i) = (((z(i) + bm(i)/2.D00 - zs(j, 2))*dk + (zs(j, 2) - z(i) + &
                                                             bm(i)/2)*akrifrt)/bm(i))
            sice(i) = ((sifrt*(zs(j, 2) - (z(i) - bm(i)/2)) + &
                        (z(i) + bm(i)/2.D00 - zs(j, 2))*1)/bm(i))
            if (tempo(i) - temp(i) .ne. 0) then
               dsidtemp(i) = (abs((siceo(i) - sice(i))/(tempo(i) - temp(i))))
            end if
            end if
         end if
!cc fin pas d'isotherme liquidus'
         if (zs(j, 2) .ne. -99 .and. zl(j, 2) .ne. -99 .and. igel == 2) then
!C maille avec solidus sans liquidus
            if (zs(j, 2) < z(i) + bm(i)/2 .and. zs(j, 2) > z(i) - bm(i)/2 &
                .and. zl(j, 2) < z(i) - bm(i)/2) then
            if (ivois(i, 3) == -99) then
               tifrt = ((temp(i) - valclt(i, 3))/(-bm(i)/2.D00) &
                        *((zs(j, 2) + z(i) - bm(i)/2)/2.D00 - z(i)) + temp(i))
            else
               tifrt = ((temp(i) - temp(ivois(i, 3)))/(z(i) - z(ivois(i, 3))) &
                        *((zs(j, 2) + z(i) - bm(i)/2)/2.D00 - z(i)) + temp(i))
            end if
            akrifrt = ((1.D00 - dk)/(tl - ts)*tifrt + (tl*dk - ts*1)/(tl - ts))
            sifrt = ((1 - swressi)/(ts - tl)*(tifrt - tl))
            akr(i) = (((z(i) + bm(i)/2.D00 - zs(j, 2))*dk + (zs(j, 2) - z(i) + &
                                                             bm(i)/2)*akrifrt)/bm(i))
            sice(i) = ((sifrt*(zs(j, 2) - (z(i) - bm(i)/2)) + &
                        (z(i) + bm(i)/2.D00 - zs(j, 2))*1)/bm(i))
            if (tempo(i) - temp(i) .ne. 0) then
               dsidtemp(i) = (abs((siceo(i) - sice(i))/(tempo(i) - temp(i))))
            end if
            end if

!C fin cas 1
! Maille avec liquidus et avec solidus
            if (zs(j, 2) < z(i) + bm(i)/2 .and. &
                zs(j, 2) > z(i) - bm(i)/2 .and. zl(j, 2) < z(i) + bm(i)/2 .and. &
                zl(j, 2) > z(i) - bm(i)/2) then
            if (ivois(i, 3) == -99) then
               tifrt = ((temp(i) - valclt(i, 3))/(-bm(i)/2.D00) &
                        *((zs(j, 2) + zl(j, 2))/2.D00 - z(i)) + temp(i))
            else
               tifrt = ((temp(i) - temp(ivois(i, 3)))/(z(i) - z(ivois(i, 3))) &
                        *((zs(j, 2) + zl(j, 2))/2.D00 - z(i)) + temp(i))
            end if
            akrifrt = ((1.D00 - dk)/(tl - ts)*tifrt + (tl*dk - ts*1)/(tl - ts))
            sifrt = ((1 - swressi)/(ts - tl)*(tifrt - tl))
            akr(i) = (((z(i) + bm(i)/2.D00 - zs(j, 2))*dk + &
                       (zs(j, 2) - zl(j, 2))*akrifrt + (zl(j, 2) - z(i) + bm(i)/2)*1)/bm(i))
            sice(i) = ((sifrt*(zs(j, 2) - zl(j, 2)) + &
                        (z(i) + bm(i)/2.D00 - zs(j, 2))*1)/bm(i))
            if (tempo(i) - temp(i) .ne. 0) then
               dsidtemp(i) = (abs((siceo(i) - sice(i))/(tempo(i) - temp(i))))
            end if
            end if
!cc fin cas 2
! Maille avec liquidus et sans solidus
            if (zl(j, 2) < z(i) + bm(i)/2 .and. zl(j, 2) > z(i) - bm(i)/2 &
                .and. zs(j, 2) > z(i) + bm(i)/2) then
            if (ivois(i, 3) == -99) then
               tifrt = ((temp(i) - valclt(i, 3))/(-bm(i)/2.D00) &
                        *((zl(j, 2) + z(i) - bm(i)/2)/2.D00 - z(i)) + temp(i))
            else
               tifrt = ((temp(i) - temp(ivois(i, 3)))/(z(i) - z(ivois(i, 3))) &
                        *((zl(j, 2) + z(i) - bm(i)/2)/2.D00 - z(i)) + temp(i))
            end if
            akrifrt = ((1.D00 - dk)/(tl - ts)*tifrt + (tl*dk - ts*1)/(tl - ts))
            sifrt = ((1 - swressi)/(ts - tl)*(tifrt - tl))
            akr(i) = (((z(i) + bm(i)/2.D00 - zl(j, 2))*akrifrt + &
                       (zl(j, 2) - z(i) + bm(i)/2)*1)/bm(i))
            sice(i) = ((sifrt*(z(i) + bm(i)/2.D00 - zl(j, 2))) &
                       /bm(i))
            if (tempo(i) - temp(i) .ne. 0) then
               dsidtemp(i) = (abs((siceo(i) - sice(i))/(tempo(i) - temp(i))))
            end if
            end if
! fin cas 3
         end if
!cc interface sur permafrost SUPRAPERMAFROST SENS INVERSE
         if (zl(j, 1) == -99 .and. zs(j, 1) .ne. -99 .and. igel == 2) then
!CC pas de liquidus
            if (zs(j, 1) < z(i) + bm(i)/2 .and. &
                zs(j, 1) > z(i) - bm(i)/2) then
            if (ivois(i, 3) == -99) then
               tifrt = ((temp(i) - valclt(i, 3))/(-bm(i)/2.D00) &
                        *(bm(i)/2.D00 + zs(j, 1) - z(i))/2.D00 + temp(i))
            else
               tifrt = ((temp(i) - temp(ivois(i, 3)))/(z(i) - z(ivois(i, 3))) &
                        *(bm(i)/2.D00 + zs(j, 1) - z(i))/2.D00 + temp(i))
            end if
            akrifrt = ((1.D00 - dk)/(tl - ts)*tifrt + (tl*dk - ts*1)/(tl - ts))
            sifrt = ((1 - swressi)/(ts - tl)*(tifrt - tl))
            akr(i) = (((z(i) + bm(i)/2.D00 - zs(j, 1))*akrifrt + &
                       (zs(j, 1) - z(i) + bm(i)/2)*dk)/bm(i))
            sice(i) = ((sifrt*(z(i) + bm(i)/2.D00 - zs(j, 1)) + &
                        (-z(i) + bm(i)/2.D00 + zs(j, 1))*1)/bm(i))
            if (tempo(i) - temp(i) .ne. 0) then
               dsidtemp(i) = (abs((siceo(i) - sice(i))/(tempo(i) - temp(i))))
            end if
            end if
         end if
!CC fin pas de liquidus

! liquidus au dessus haut cellule
         if (zl(j, 1) .ne. -99 .and. zs(j, 1) .ne. -99 .and. igel == 2) then
         if (zl(j, 1) > z(i) + bm(i)/2 .and. zs(j, 1) < z(i) + bm(i)/2 &
             .and. zs(j, 1) > z(i) - bm(i)/2) then
         if (ivois(i, 3) == -99) then
            tifrt = ((temp(i) - valclt(i, 3))/(-bm(i)/2.D00) &
                     *(bm(i)/2.D00 + zs(j, 1) - z(i))/2.D00 + temp(i))
         else
            tifrt = ((temp(i) - temp(ivois(i, 3)))/(z(i) - z(ivois(i, 3))) &
                     *(bm(i)/2.D00 + zs(j, 1) - z(i))/2.D00 + temp(i))
         end if
         akrifrt = ((1.D00 - dk)/(tl - ts)*tifrt + (tl*dk - ts*1)/(tl - ts))
         sifrt = ((1 - swressi)/(ts - tl)*(tifrt - tl))
         akr(i) = (((z(i) + bm(i)/2.D00 - zs(j, 1))*akrifrt + &
                    (zs(j, 1) - z(i) + bm(i)/2)*dk)/bm(i))
         sice(i) = ((sifrt*(z(i) + bm(i)/2.D00 - zs(j, 1)) + &
                     (-z(i) + bm(i)/2.D00 + zs(j, 1))*1)/bm(i))
         if (tempo(i) - temp(i) .ne. 0) then
            dsidtemp(i) = (abs((siceo(i) - sice(i))/(tempo(i) - temp(i))))
         end if
!       print*,'degel4',akr(i),sice(i),temp(i),zl(j,1),zs(j,1)
         end if
! fin cas 1
         if (zl(j, 1) < z(i) + bm(i)/2 .and. zl(j, 1) > &
             z(i) - bm(i)/2 .and. &
             zs(j, 1) < z(i) + bm(i)/2 .and. zs(j, 1) > z(i) - bm(i)/2) then
         if (ivois(i, 3) == -99) then
            tifrt = ((temp(i) - valclt(i, 3))/(-bm(i)/2.D00) &
                     *((zl(j, 1) + zs(j, 1))/2.D00 - z(i)) + temp(i))
         else
            tifrt = ((temp(i) - temp(ivois(i, 3)))/(z(i) - z(ivois(i, 3))) &
                     *((zl(j, 1) + zs(j, 1))/2.D00 - z(i)) + temp(i))
         end if
         akrifrt = ((1.D00 - dk)/(tl - ts)*tifrt + (tl*dk - ts*1)/(tl - ts))
         sifrt = ((1 - swressi)/(ts - tl)*(tifrt - tl))
         akr(i) = (((z(i) + bm(i)/2.D00 - zl(j, 1))*1.D00 + &
                    (zl(j, 1) - zs(j, 1))*akrifrt + (zs(j, 1) - z(i) + bm(i)/2)*dk)/bm(i))
         sice(i) = ((sifrt*(zl(j, 1) - zs(j, 1)) + &
                     (-z(i) + bm(i)/2.D00 + zs(j, 1))*1)/bm(i))
         if (tempo(i) - temp(i) .ne. 0) then
            dsidtemp(i) = (abs((siceo(i) - sice(i))/(tempo(i) - temp(i))))
         end if
         end if
!cccc fin cas 2
         if (zl(j, 1) < z(i) + bm(i)/2 .and. &
             zl(j, 1) > z(i) - bm(i)/2 .and. &
             zs(j, 1) < z(i) - bm(i)/2) then
         if (ivois(i, 3) == -99) then
            tifrt = ((temp(i) - valclt(i, 3))/(-bm(i)/2.D00) &
                     *((zl(j, 1) + z(i) - bm(i)/2)/2.D00 - z(i)) + temp(i))
         else
            tifrt = ((temp(i) - temp(ivois(i, 3)))/(z(i) - z(ivois(i, 3))) &
                     *((zl(j, 1) + z(i) - bm(i)/2)/2.D00 - z(i)) + temp(i))
         end if
         akrifrt = ((1.D00 - dk)/(tl - ts)*tifrt + (tl*dk - ts*1)/(tl - ts))
         sifrt = ((1 - swressi)/(ts - tl)*(tifrt - tl))
         akr(i) = (((z(i) + bm(i)/2.D00 - zl(j, 1))*1 + &
                    (zl(j, 1) - z(i) + bm(i)/2)*akrifrt)/bm(i))
         sice(i) = ((sifrt*(-z(i) + bm(i)/2.D00 + zl(j, 1))) &
                    /bm(i))
         if (tempo(i) - temp(i) .ne. 0) then
            dsidtemp(i) = (abs((siceo(i) - sice(i))/(tempo(i) - temp(i))))
         end if

         end if
         end if
      end if
      end do
      if (akr(i) < dk) akr(i) = dk
      akrv(i) = akr(i)
   end do
   return
end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!                Subroutine lecture parametres                    C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC...PURPOSE :
!ccc...
!CC...Modification
!ccc...
!CC....TODOLIST
!ccc....
subroutine lecture_parametre( &
   ilog, ipermh, ipermv, &
   iec, irp, ith, dt, nitt, unitsim, al, &
   az, ixy, imaille, itopo, reptop, repbot, icolone, nmi, nci, &
   nri, dx, dz, nclog, rho1, irho, g, amu, &
   akrx, akx, akrz, akz, iom, omp, iss, sss, ia2, yunconfined, ivg, &
   ans, asp, swres, itr, allg, alt, iriv, iqriv, qre, hbot, xberg, &
   rug, pent, qriva, hriv, aklit, aklitv, tempriv, elit, akdrain, &
   edrain, crconvp, crconvc, &
   iteration, itsortie, unitsortie, icalvit, mcol, nmaille1, nmaille2, &
   nmaille3, nmaille4, nmaille5, &
   nmaille6, nmaille7, nmaille8, nmaille9, nmaille10, iparo, ibuilt, &
   ysolv)
   implicit double precision(a - h, o - x, z), integer(I - N)
   implicit CHARACTER*5(y)
   CHARACTER(3) :: ysolv, yunconfined
   open (3, file='E_parametre.dat', form='formatted', status='old')
   READ (3, '(4x,i1)') iec    !!  l1
   READ (3, '(3x,i1)') irp    !!  l2
   READ (3, '(4x,i1)') ith    !!  l3
   READ (3, '(3x,d9.0)') dt   !!  l4
   READ (3, '(5x,i10)') nitt  !!  l5
   READ (3, '(6x,d9.0)') unitsim  !!  l6
   READ (3, '(3x,d8.0)') al   !!  l7
   READ (3, '(3x,d9.0)') az   !!  l8
   READ (3, '(4x,i1)') ixy    !!  l9
   READ (3, '(8x,i1)') imaille    !!  l10
   READ (3, '(6x,i1)') itopo  !!  l11
   READ (3, '(7x,i1)') ibuilt     !!  l12
   READ (3, '(7x,d9.0)') reptop   !!  l13
   READ (3, '(7x,d9.0)') repbot   !!  l14
   READ (3, '(9x,i1)') icolone    !!  l15
   READ (3, '(4x,i6)') nmi    !!  l16
   READ (3, '(4x,i5)') nci    !!  l17
   READ (3, '(4x,i5)') nri    !!  l18
   READ (3, '(3x,d8.0)') dx   !!  l19
   READ (3, '(3x,d8.0)') dz   !!  l20
   READ (3, '(6x,i5)') nclog  !!  l21
   READ (3, '(5x,i1)') ilog   !!  l22
   READ (3, '(5x,d8.0)') rho1     !!  l23
   READ (3, '(5x,i1)') irho   !!  l24
   READ (3, '(2x,d8.0)') g    !!  l25
   READ (3, '(4x,d8.0)') amu  !!  l26
   READ (3, '(7x,i1)') ipermh     !!  l27
   READ (3, '(7x,i1)') ipermv     !!  l28
   READ (3, '(5x,d8.0)') akrx     !!  l29
   READ (3, '(4x,d8.0)') akx  !!  l30
   READ (3, '(5x,d8.0)') akrz     !!  l31
   READ (3, '(4x,d8.0)') akz  !!  l32
   READ (3, '(4x,i1)') iom    !!  l33
   READ (3, '(4x,f5.3)') omp  !!  l34
   READ (3, '(4x,i1)') iss    !!  l35
   READ (3, '(3x,d8.0)') sss  !!  l36
   READ (3, '(4x,i1)') ia2    !!  l37
   READ (3, '(12x,A3)') yunconfined   !!  l38
   READ (3, '(4x,i1)') ivg    !!  l39
   READ (3, '(4x,d8.0)') ans  !!  l40
   READ (3, '(4x,d8.0)') asp  !!  l41
   READ (3, '(6x,f6.4)') swres    !!  l42
   READ (3, '(4x,i1)') itr    !!  l43
   READ (3, '(4x,d8.0)') allg     !!  l44
   READ (3, '(4x,d8.0)') alt  !!  l45
   READ (3, '(5x,i1)') iriv   !!  l46
   READ (3, '(6x,i1)') iqriv  !!  l47
   READ (3, '(4x,d8.0)') qre  !!  l48
   READ (3, '(5x,d8.0)') hbot     !!  l49
   READ (3, '(6x,d8.0)') xberg    !!  l50
   READ (3, '(4x,d8.0)') rug  !!  l51
   READ (3, '(5x,d8.0)') pent     !!  l52
   READ (3, '(6x,d8.0)') qriva    !!  l53
   READ (3, '(5x,d8.0)') hriv     !!  l54
   READ (3, '(6x,d8.0)') aklit    !!  l55
   READ (3, '(7x,d8.0)') aklitv   !!  l56
   READ (3, '(8x,d8.0)') tempriv  !!  l57
   READ (3, '(5x,d8.0)') elit     !!  l58
   READ (3, '(8x,d8.0)') akdrain  !!  l59
   READ (3, '(7x,d8.0)') edrain   !!  l60
   READ (3, '(8x,d8.0)') crconvp  !!  l61
   READ (3, '(8x,d8.0)') crconvc  !!  l62
   READ (3, '(10x,i8)') iteration     !!  l63
   READ (3, '(9x,i8)') itsortie   !!  l64
   READ (3, '(11x,d9.0)') unitsortie  !!  l65
   READ (3, '(8x,i1)') icalvit    !!  l66
   READ (3, '(5x,i5)') mcol   !!  l67
   READ (3, '(9x,i5)') nmaille1   !!  l68
   READ (3, '(9x,i5)') nmaille2   !!  l69
   READ (3, '(9x,i5)') nmaille3   !!  l70
   READ (3, '(9x,i5)') nmaille4   !!  l71
   READ (3, '(9x,i5)') nmaille5   !!  l72
   READ (3, '(9x,i5)') nmaille6   !!  l73
   READ (3, '(9x,i5)') nmaille7   !!  l74
   READ (3, '(9x,i5)') nmaille8   !!  l75
   READ (3, '(9x,i5)') nmaille9   !!  l76
   READ (3, '(10x,i5)') nmaille10     !!  l77
   READ (3, '(6x,i1)') iparo  !!  l78
   READ (3, '(6x,A3)') ysolv  !!  l79
   close (3)
   return
end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!            lecture parametres thermique                          C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC...PURPOSE :
!ccc...
!CC...Modification
!ccc...
!CC....TODOLIST
!ccc....
subroutine lecture_parametre_thermique(irpth, crconvt, &
                                       ithec, idecouplage, ysupdp, &
                                       bll, blt, alandae, cpe, ilanda, &
                                       alandami, irhomi, rhosi, icpm, cpm, ymoycondtherm, &
                                       icycle, igel, igelzns, iomdegel, &
                                       ytypakrice, ytypsice, &
                                       alandai, rhoii, cpice, &
                                       chlat, rhog, alandag, cpg, dk, tsg, tlg, tsd, tld, cimp, swressi, &
                                       iaquitard, aktardx, aktardz, omptard, alandatard, sstard, &
                                       nrowtard, nmailleaqui, iexutoire, xexutoire, tempexutoire, &
                                       iinfil, itopomch, ibigridice, ytest, omega, cimpt)
   implicit double precision(a - h, o - x, z), integer(I - N)
   implicit CHARACTER*5(y)
   CHARACTER(5) :: ytypakrice, ytypsice
   CHARACTER(3) :: ytest
   open (unit=7, file='E_p_therm.dat', &
         form='formatted', status='old')
   READ (7, '(5x,i1)') irpth  !!  l1
   READ (7, '(8x,d8.0)') crconvt  !!  l2
   READ (7, '(6x,i1)') ithec  !!  l3
   READ (7, '(12x,i1)') idecouplage   !!  l4
   READ (7, '(7x,A5)') ysupdp     !!  l5
   READ (7, '(4x,d8.0)') bll  !!  l6
   READ (7, '(4x,d8.0)') blt  !!  l7
   READ (7, '(8x,d8.0)') alandae  !!  l8
   READ (7, '(4x,d8.0)') cpe  !!  l9
   READ (7, '(7x,i1)') ilanda     !!  l10
   READ (7, '(8x,d8.0)') alandami     !!  l11
   READ (7, '(7x,i1)') irhomi     !!  l12
   READ (7, '(6x,d8.0)') rhosi    !!  l13
   READ (7, '(5x,i1)') icpm   !!  l14
   READ (7, '(4x,d8.0)') cpm  !!  l15
   READ (7, '(14x,A5)') ymoycondtherm     !!  l16
   READ (7, '(7x,i1)') icycle     !!  l17
   READ (7, '(5x,i1)') igel   !!  l18
   READ (7, '(8x,i1)') igelzns    !!  l19
   READ (7, '(9x,i1)') iomdegel   !!  l20
   READ (7, '(11x,A5)') ytypakrice    !!  l21
   READ (7, '(9x,A5)') ytypsice   !!  l22
   READ (7, '(11x,i1)') ibigridice    !!  l23
   READ (7, '(8x,d10.0)') alandai     !!  l24
   READ (7, '(5x,d8.0)') rhoii    !!  l25
   READ (7, '(6x,d10.0)') cpice   !!  l26
   READ (7, '(4x,d10.0)') chlat   !!  l27
   READ (7, '(5x,d8.0)') rhog     !!  l28
   READ (7, '(8x,d8.0)') alandag  !!  l29
   READ (7, '(4x,d8.0)') cpg  !!  l30
   READ (7, '(3x,d8.0)') dk   !!  l31
   READ (7, '(4x,d8.0)') tsg  !!  l32
   READ (7, '(4x,d8.0)') tlg  !!  l33
   READ (7, '(4x,d8.0)') tsd  !!  l34
   READ (7, '(4x,d8.0)') tld  !!  l35
   READ (7, '(5x,d8.0)') cimp     !!  l36
   READ (7, '(5x,d8.0)') cimpt    !!  l37
   READ (7, '(6x,d8.0)') omega    !!  l38
   READ (7, '(6x,d8.0)') rlamb    !!  l39
   READ (7, '(8x,f6.4)') swressi  !!  l40
   READ (7, '(10x,i1)') iaquitard     !!  l41
   READ (7, '(8x,d8.0)') aktardx  !!  l42
   READ (7, '(8x,d8.0)') aktardz  !!  l43
   READ (7, '(8x,f5.3)') omptard  !!  l44
   READ (7, '(11x,d8.0)') alandatard  !!  l45
   READ (7, '(7x,d8.0)') sstard   !!  l46
   READ (7, '(9x,i5)') nrowtard   !!  l47
   READ (7, '(12x,i5)') nmailleaqui   !!  l48
   READ (7, '(10x,i1)') iexutoire     !!  l49
   READ (7, '(10x,d8.0)') xexutoire   !!  l50
   READ (7, '(13x,d8.0)') tempexutoire    !!  l51
   READ (7, '(7x,i1)') iinfil     !!  l52
   READ (7, '(9x,i1)') itopomch   !!  l53
   READ (7, '(6x,A3)') ytest  !!  l54
   close (7)
   return
end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!                lecture CDT INITIALE                          C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC...PURPOSE :
!ccc...
!CC...Modification
!ccc...
!CC....TODOLIST
!ccc....

subroutine lecture_cdt_ini(chgi, conci, tempini, ichi, ichi2, iconci, &
                           itempi)
   implicit double precision(a - h, o - x, z), integer(I - N)
   implicit CHARACTER*5(y)
   open (unit=5, file='E_cdt_initiale.dat', form='formatted', &
         status='old')
   READ (5, '(5x,d9.0)') chgi
   READ (5, '(6x,d8.0)') conci
   READ (5, '(6x,d8.0)') tempini
   READ (5, '(5x,i1)') ichi
   READ (5, '(6x,i1)') ichi2
   READ (5, '(7x,i1)') iconci
   READ (5, '(7x,i1)') itempi
   close (5)
   return
end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!                lecture CDT LIMITES                          C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC...PURPOSE :
!ccc...
!CC...Modification
!ccc...
!CC....TODOLIST
!ccc....
subroutine lecture_cdt_lt(icl_gauche, valcl_gauche, icl_droite, &
                          valcl_droite, icl_haut, valcl_haut, icl_bas, valcl_bas, iclc_gauche, &
                          valclc_gauche, iclc_droite, valclc_droite, iclc_haut, valclc_haut, &
                          iclc_bas, valclc_bas, &
                          iclt_gauche, valclt_gauche, iclt_droite, valclt_droite, iclt_haut, &
                          valclt_haut, iclt_bas, valclt_bas, iclect, icltherm, iclectchgt &
                          , iclchgt, icldrain, iclriviere, itlecture)
   implicit double precision(a - h, o - x, z), integer(I - N)
   implicit CHARACTER*5(y)
   open (unit=4, file='E_cdt_aux_limites.dat', form='formatted', &
         status='old')
   READ (4, '(11x,i2)') icl_gauche
   READ (4, '(13x,d12.0)') valcl_gauche
   READ (4, '(11x,i2)') icl_droite
   READ (4, '(13x,d12.0)') valcl_droite
   READ (4, '(9x,i2)') icl_haut
   READ (4, '(11x,d12.0)') valcl_haut
   READ (4, '(8x,i2)') icl_bas
   READ (4, '(10x,d12.0)') valcl_bas
   READ (4, '(12x,i2)') iclc_gauche
   READ (4, '(14x,d12.0)') valclc_gauche
   READ (4, '(12x,i2)') iclc_droite
   READ (4, '(14x,d12.0)') valclc_droite
   READ (4, '(10x,i2)') iclc_haut
   READ (4, '(12x,d12.0)') valclc_haut
   READ (4, '(9x,i2)') iclc_bas
   READ (4, '(11x,d12.0)') valclc_bas
   READ (4, '(12x,i2)') iclt_gauche
   READ (4, '(14x,d12.0)') valclt_gauche
   READ (4, '(12x,i2)') iclt_droite
   READ (4, '(14x,d12.0)') valclt_droite
   READ (4, '(10x,i2)') iclt_haut
   READ (4, '(12x,d12.0)') valclt_haut
   READ (4, '(9x,i2)') iclt_bas
   READ (4, '(11x,d12.0)') valclt_bas
   READ (4, '(7x,i1)') iclect
   READ (4, '(9x,i1)') icltherm
   READ (4, '(11x,i1)') iclectchgt
   READ (4, '(8x,i1)') iclchgt
   READ (4, '(9x,i1)') icldrain
   READ (4, '(11x,i1)') iclriviere
   READ (4, '(10x,i8)') itlecture
   close (4)
   return
end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!                     INTERPOLATION  surf piezo                    C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC...PURPOSE :
!ccc...Calcul de la position d'une isovaleur d'une variable d'état
!CC...Modification
!ccc...
!CC....TODOLIST
!ccc....
subroutine interpolsurf(nc, zo, nm, icol, pr, ivois, &
                        bm, z, valcl, icl, id_river, id_rivert, chgriver, &
                        ligne5, ligne6, ligne2, paso, itlecture, rho1, g)
   implicit double precision(a - h, o - x, z), integer(I - N)
   implicit CHARACTER*5(y)
   dimension z(nm)
   dimension ivois(nm, 4), valcl(nm, 4)
   dimension zo(nc, 2)
   dimension bm(nm), pr(nm)
   dimension icol(nm), icl(nm, 4), id_rivert(ligne6)
   dimension id_river(ligne5), chgriver(ligne2)
   to=0
   kimp = int(paso/itlecture) + 1
   if (kimp > ligne2) kimp = ligne2
!CCC position surf piezo
   do kkcol = 1, nc
      zo(kkcol, 1) = -99D+00
!CC PASSAGE SUR TOUT LES MAILLES
      do i = 1, nm
! TEST COLONNE
         if (icol(i) == kkcol) then
            ili = ivois(i, 3)
            select case (ivois(i, 3))
            case (-99)
         !CC PAS DE VOISIN EN HAUT
            if (zo(kkcol, 1) == -99 .and.  &
                pr(i) > to .and. &
                valcl(i, 3) < to .and. &
                icl(i, 3) == -2) then
               zo(kkcol, 1) = z(i) + &
                              (bm(i)/2)/(valcl(i, 3) - pr(i)) &
                              *(to - valcl(i, 3))
            end if
            case default
            if (ili < nm + 1 .and. &
                pr(i) >= to .and. pr(ili) < to) then
               zo(kkcol, 1) = z(ivois(i, 3)) + &
                              ((z(ili) - z(i))/(pr(ili) - pr(i)) &
                               *(to - pr(ili)))
            end if
            end select



!CPAS DE VOISIN EN BAS
            if (ivois(i, 4) == -99 .and. zo(kkcol, 1) == -99 .and. &
                pr(i) < to .and. icl(i, 4) == -2 .and. valcl(i, 4) > to) then
               zo(kkcol, 1) = z(i) - &
                              (bm(i)/2)/(pr(i) - valcl(i, 4)) &
                              *(to - pr(i))
            end if

            if (zo(kkcol, 1) == -99 .and. ivois(i, 3) == -99) then
               zo(kkcol, 1) = pr(i)/(rho1*g) + z(i)
            end if

            do j = 1, ligne5
            if (i == id_river(j)) then
               zo(kkcol, 1) = chgriver(kimp)
            end if
            end do

            do j = 1, ligne6
            if (i == id_rivert(j)) then
               zo(kkcol, 1) = chgriver(kimp)
            end if
            end do

!CCC FIN DE test colonne
         end if
      end do
   end do

!CCC position surf piezo
   do kkcol = 1, nc

!CC PASSAGE SUR TOUT LES MAILLES
      do i = 1, nm
! TEST COLONNE
         if (icol(i) == kkcol) then
         do j = 1, ligne5
         if (i == id_river(j)) then
            zo(kkcol, 1) = chgriver(kimp)
         end if
         end do
         do j = 1, ligne6
         if (i == id_rivert(j)) then
            zo(kkcol, 1) = chgriver(kimp)
         end if
         end do

!CCC FIN DE test colonne
         end if
      end do
   end do

   return
end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!                CALCUL MSE                                C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC...PURPOSE :
!ccc...
!CC...Modification
!ccc...
!CC....TODOLIST
!ccc....
subroutine mse(n, nmaille1, nmaille2, nmaille3, nmaille4, nmaille5, &
               nmaille6, pr, qinf, qinfobs, paso, itlecture, chgobs1, chgobs2, &
               chgobs3, chgobs4, rho, g, chgobs5, chgobs6, eo1, eo2, eo3, eo4, eo5, &
               eo6, eoinf, z, seo1, seo2, seo3, seo4, seo5, &
               seo6, seoinf)
   implicit double precision(a - h, o - x, z), integer(I - N)
   implicit CHARACTER*5(y)
   dimension pr(n), rho(n), z(n)
   dimension chgobs1(1), chgobs2(1), chgobs3(1), chgobs4(1)
   dimension chgobs5(1), chgobs6(1), qinfobs(1)

   if (modulo(int(paso), itlecture) == 0 .and. paso >= itlecture) then
      k = int((paso)/itlecture) + 1
      eo1 = 0D+00
      eo2 = 0D+00
      eo3 = 0D+00
      eo4 = 0D+00
      eo5 = 0D+00
      eo6 = 0D+00
      eoinf = 0D+00
      eo1 = (pr(nmaille1)/(rho(nmaille1)*g) + z(nmaille1) - chgobs1(k))**2
      eo2 = (pr(nmaille2)/(rho(nmaille2)*g) + z(nmaille2) - chgobs2(k))**2
      eo3 = (pr(nmaille3)/(rho(nmaille3)*g) + z(nmaille3) - chgobs3(k))**2
      eo4 = (pr(nmaille4)/(rho(nmaille4)*g) + z(nmaille4) - chgobs4(k))**2
      eo5 = (pr(nmaille5)/(rho(nmaille5)*g) + z(nmaille5) - chgobs5(k))**2
      eo6 = (pr(nmaille6)/(rho(nmaille6)*g) + z(nmaille6) - chgobs6(k))**2
      eoinf = (qinf - qinfobs(k))**2
      seo1 = seo1 + eo1
      seo2 = seo2 + eo2
      seo3 = seo3 + eo3
      seo4 = seo4 + eo4
      seo5 = seo5 + eo5
      seo6 = seo6 + eo6
      seoinf = seoinf + eoinf

   end if
   return
end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!              Funcions to open DTS files                     C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine DTS_open(icolone, imaille, itopo, ligne, ligne1, &
                    ligne2, ligne3, ligne4, ligne6)
   implicit integer(i - n)
   LOGICAL :: An_Error
   icolone = 1
   imaille = 0
   itopo = 1

   call open_file('E_zone.dat', An_Error, 32)
   call open_file('E_zone_parameter.dat', An_Error, 321)

   open (unit=222, file='E_p_boundaryB.dat', iostat=ios, action='read')
   open (unit=2233, file='E_t_boundaryB.dat')
   call count_file(222, ios, ligne1)

   open (unit=2244, file='E_p_boundaryPP.dat', iostat=ioPP)
   open (unit=225, file='E_t_boundaryPP.dat', iostat=ios3)
   call count_file(2244, ioPP, ligne2)

   open (unit=226, file='FODTSginette.dat', iostat=i226)
   call count_file(226, i226, ligne)
   open (unit=227, file='E_heads.dat', iostat=ioslope)
   call count_file(227, ioslope, ligne3)

   open (unit=300, file='E_pool_coord.dat', iostat=ipool)
   call count_file(300, ipool, ligne4)

   open (unit=300400, file='E_riffle_coord.dat', iostat=iriff)
   call count_file(300400, iriff, ligne6)

   return
end

subroutine DTS_read(nc, xpool, ligne, ligne1, &
                    ligne2, ligne3, ligne4, &
                    ligne6, timeDTS, tempDTS, &
                    chgRG, cRivG, timeG, tempRG, &
                    chgRD, cRivD, timeD, tempRD, &
                    xriffle, slopeRH)

   implicit double precision(a - h, o - z)
   implicit integer(i - n)
   dimension timeDTS(ligne), tempDTS(ligne, nc)
   dimension chgRG(ligne1), cRivG(ligne1), timeG(ligne1)
   dimension tempRG(ligne1), chgRD(ligne2), cRivD(ligne2)
   dimension timeD(ligne2), tempRD(ligne2), xriffle(2, ligne6)
   dimension xpool(2, ligne4), slopeRH(2, ligne3)

   rewind (226)
   do j = 1, ligne
      read (226, *, iostat=ios4) timeDTS(j) &
         , (tempDTS(j, i), i=1, nc)
   end do

   do j = 1, ligne1
      read (222, *, iostat=ios) timeG(j), chgRG(j), cRivG(j)
      read (2233, *, iostat=ios1) a, tempRG(j), a
   end do
   do j = 1, ligne2
      read (2244, *, iostat=ios2) timeD(j), chgRD(j), cRivD(j)
      read (225, *, iostat=ios3) a, tempRD(j), a
   end do

   do j = 1, ligne3
      read (227, *, iostat=ioslope) slopeRH(1, j), slopeRH(2, j)
   end do

   do j = 1, ligne6
      read (300400, *, iostat=iriffle) xriffle(1, j), xriffle(2, j)
   end do
   do j = 1, ligne4
      read (300, *, iostat=ipool) xpool(1, j), xpool(2, j)
   end do

   return
end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!              Funcions for output                    C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC...PURPOSE :
!ccc...
!CC...Modification
!ccc...
!CC....TODOLIST
!ccc....

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                 C
!              VARIATION CDT LIMITES VS. TEMPS                    C
!                                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC...PURPOSE :
!ccc...
!CC...Modification
!ccc...
!CC....TODOLIST
!ccc....

subroutine variation_cdt_limites(nm, paso, itlecture, ytest, &
                                 ligne, ligne1, ligne2, ligne3, ligne4, ligne5, ligne6, &
                                 icl, valcl, iclt, valclt, ivois, &
                                 z, g, ntsortie, bm, irptha, &
                                 rho, qpluie, chgriver, &
                                 chgRD, chgRG, tempRD, tempRG, &
                                 tempriver, &
                                 id_RD, id_RG, id_river, id_rivert, tempsurf, &
                                 tempbottom, chgsurf, chgbottom, &
                                 x,&
                                 qsurf, qbottom,iecriture_pluie,akr_bottom,ak_bottom,akr_surf,ak_surf,vzm_nm1,swo,pro)
!     nc nb de colonnes
!     nr nb de ligne
!       qre debit pluie
!       ivois(ik,1)= voisin droite
!       ivois(ik,2)= voisin gauche
!       ivois(ik,3)= voisin haut
!       ivois(ik,4)= voisin ibas
!       nm nombre de mailles reelles
! ECOULEMENT icl condition valcl valeur
!       ICL=-1 Flux impose sur une face
!       ICL=-2 potentiel impose sur une face
!       ICL=1 Mailles 'normale'
!       ICL(ik,1)=-1, -2 ou 1 ik=1,...4: face droite,gauche,haute et BASSE
! THERMIQUE iclt condition valclt valeur
!       ICLT=-1 Flux impose sur une face
!       ICLT=-2 potentiel impose sur une face
!       ICLT=1 Mailles 'normale'
!       ICLT(ik,1)=-1, -2 ou 1 ik=1,...4: face droite,gauche,haute et BASSE
!       akr(i)=permeabilite relative
!       ak(i)=permeabilite intrinseque
!MODIF 18/03/2025 nradic Pour gerer la pluie lorsque que le sol est saturé.

   implicit double precision(a - h, o - z)
   implicit integer(i - n)
   dimension ivois(nm, 4), valcl(nm, 4), icl(nm, 4)
   dimension valclt(nm, 4), iclt(nm, 4), z(nm)
   dimension pro(nm)
   dimension rho(nm), bm(nm), x(nm)
   dimension chgbottom(ntsortie), chgsurf(ntsortie)
   dimension qbottom(ntsortie), qsurf(ntsortie)
   dimension qpluie(ntsortie), chgRD(ntsortie)
   dimension chgRG(ntsortie), tempRD(ntsortie)
   dimension chgriver(ntsortie), tempriver(ntsortie)
   dimension id_RD(ligne3), id_RG(ligne4)
   dimension id_river(ligne5)
   dimension tempRG(ntsortie), id_rivert(ligne6)
   dimension tempsurf(ntsortie), tempbottom(ntsortie)
   CHARACTER(3) :: ytest

   kimp = int(paso/itlecture)
   if (kimp >= ntsortie) kimp = ntsortie
   
   select case (ytest)
   case ("R2D") 



      do i = 1, nm
!ccc....CDT LIMITE SOL
         if (ivois(i, 3) == -99) then

            icl(i, 3) = -1
            valcl(i, 3) = 0
            valclt(i, 3) = 0
           iclt(i, 3) = -1

         end if

!ccc....CDT LIMITE BOTTOM
         if (ivois(i, 4) == -99) then
            iclt(i, 4) = -1
            icl(i, 4) = -1
            valclt(i,4)=0
            valcl(i,4)=0
         end if

!ccc....CDT LIMITE RD
         do j = 1, ligne3  ! # E_BordRD.dat
         if (i == id_RD(j)) then
            icl(i, 1) = -2
            valcl(i, 1) = (chgRD(kimp) - z(i))*rho(i)*g
     iclt(i,1)=-2
     valclt(i,1)=tempRD(kimp)
         end if
         end do

!ccc....CDT LIMITE RG
         do j = 1, ligne4 ! 223 = E_BordRG.dat
            
         if (i == id_RG(j)) then
            icl(i, 2) = -2
            valcl(i, 2) = (chgRG(kimp) - z(i))*rho(i)*g
        iclt(i,2)=-2
        valclt(i,2)=tempRG(kimp)
!        print*,tempRG(kimp),kimp
!     &-(tempRG(kimp)-tempsol(kheure))
!     &/(76.83-79.81)*(76.83-z(i))
         end if
         end do
!ccc....CDT LIMITE RIVER TOUJOURS SOUS L'eau
         do j = 1, ligne5 ! 224 = E_Id_river.dat
         if (i == id_river(j)) then
            icl(i, 3) = -2
            zhaut = z(i) + bm(i)/2
            valcl(i, 3) = (chgriver(kimp) - zhaut)*rho(i)*g
           iclt(i,3)=-2
           valclt(i,3)=tempriver(kimp)

         end if
         end do

!ccc....CDT LIMITE RIVER a tester
         do j = 1, ligne6 ! 2244 = E_Id_river_max.dat

         if (i == id_rivert(j)) then
            if (z(i) + bm(i)/2 < chgriver(kimp)) then
               icl(i, 3) = -2
               zhaut = z(i) + bm(i)/2
               valcl(i, 3) = (chgriver(kimp) - zhaut)*rho(i)*g
              iclt(i,3)=-2
               valclt(i,3)=tempriver(kimp)
            end if
         end if
         end do

      end do
   case ("WAR")
      if (kimp > ligne4) kimp = ligne4
      select case (icl(nm, 4))
      case (-2)
         valcl(nm, 4) = chgbottom(kimp)*rho(nm)*g
      case (-1)
         valcl(nm, 4) = qbottom(kimp)
      end select
      select case (icl(1, 3))
      case (-2)
         valcl(1, 3) = chgsurf(kimp)*rho(1)*g
      case (-1)
         valcl(1, 3) = qsurf(kimp)
      end select
   case ("ZNS")

      kheure = int(paso/3600) + 1 ! Modif test
      if (kheure > ligne4) kheure = ligne4
      if (kimp > ligne4) kimp = ligne4
      select case (icl(nm, 4))
      case (-2)
         zbas = z(nm) - bm(nm)/2
         if (abs(zbas) < 1e-6) zbas = 0D+00
         valcl(nm, 4) = (rho(nm)*g*(chgbottom(ligne4) - zbas))
      case (-1)
         valcl(nm, 4) = qbottom(kimp)
!....option debug
!         print *, vzm_nm1,valcl(nm, 4), qbottom(kimp), akr_bottom, ak_bottom,'bottom'


      case (-4) !cas où l'on prends la vitesse de la maille précédente pour calculer le débit
         valcl(nm, 4) = pro(nm-1)-rho(nm)*g*(z(nm-1)-z(nm)) ! Le flux sortant correspond au à la vitesse  de la face du bas de l'avant dernière maille
         print *, valcl(nm, 4), akr_bottom, ak_bottom,'bottom'

      end select

      if (ivois(1, 3) == -99) then
         select case (icl(1, 3))
         case (-2)
            zhaut = z(1) + bm(1)/2
            if (abs(zhaut) < 1e-6) zhaut = 0D+00
            valcl(1, 3) = rho(1)*g*(chgsurf(kimp) - zhaut)
         case (-1)
            !Ajouter le ruisselement
            !ccc.....MODIF AVEC PLUIE nicolas radic
            qpossible = akr_surf*1000*9.81/0.001            !akr_surf*ak_surf*1000*9.81/0.001 ! Calcul du flux max possible à l'aide de la permeabilité relative.
            if (qpossible <= qsurf(kimp)) then ! Le débit de pluie est plus grand que le debit max possible dans le sol
                 valcl(1,3) = qpossible ! L'infiltration est plafonné par la permeabilté relative
                 print *, "ON A DIMINUE LE FLUX EN SURFACE en passant de ", qsurf(kimp), qpossible, akr_surf
                 print *, 'HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH'
                 print *, paso 
               !   stop
            else ! Dans le cas ou la pluie est plus faible que la permabilité relative 
            valcl(1, 3) = qsurf(kimp)
            end if 
            if (swo .gt. 0.985) then
               valcl(1,3) = 0
               print *, "Profil saturée le flux d'infiltration en surface est nul. 100% ruisselement valcl(1, 3) = 0"
               ! stop
            end if
            !ccc....END MODIF
            !valcl(1, 3) = qsurf(kimp)
         end select
      end if
      if (irptha == 1 ) then
      valclt(1, 3) = tempsurf(kimp)
      valclt(nm, 4) = tempbottom(kimp)
      end if

   case ("ZND")

      if (kimp > ligne4) kimp = ligne4
      select case (icl(nm, 4))
      case (-2)
         zbas = z(nm) - bm(nm)/2
         if (abs(zbas) < 1e-6) zbas = 0D+00
         valcl(nm, 4) = (rho(nm)*g*(chgbottom(ligne4) - zbas))
      case (-1)
         valcl(nm, 4) = 0
      end select

      if (ivois(1, 3) == -99) then
         zhaut = z(1) + bm(1)/2
         if (chgsurf(kimp) - zhaut >= 0.001) then
            icl(1, 3) = -2
         else
            icl(1, 3) = -1
         end if
         select case (icl(1, 3))
         case (-2)
            print *, 'coucou', chgsurf(kimp)
            zhaut = z(1) + bm(1)/2
            if (abs(zhaut) < 1e-6) zhaut = 0D+00
            valcl(1, 3) = rho(1)*g*(chgsurf(kimp) - zhaut)
         case (-1)
            valcl(1, 3) = 0
         end select
      end if



   case ("ZHR","ZHZ") 
      if (kimp > ligne4) kimp = ligne4
      iclt(1, 3) = -2
      valclt(1, 3) = tempsurf(kimp)
      icl(1, 3) = -2
      zhaut = z(1) + bm(1)/2
      if (abs(zhaut) < 1e-6) zhaut = 0D+00
      valcl(1, 3) = rho(1)*g*(chgsurf(kimp) - zhaut)
      icl(nm, 4) = -2
      zbas = z(nm) - bm(nm)/2
      if (abs(zbas) < 1e-6) zbas = 0D+00
      valcl(nm, 4) = rho(nm)*g*(chgbottom(kimp) - zbas)
      iclt(nm, 4) = -2
      valclt(nm, 4) = tempbottom(kimp)
      if (kimp > ligne4) then
         iclt(1, 3) = -2
         valclt(1, 3) = tempsurf(ligne4)
         icl(1, 3) = -2
         zhaut = z(1) + bm(1)/2
         if (abs(zhaut) < 1e-6) zhaut = 0D+00
         valcl(1, 3) = (rho(1)*g*(chgsurf(ligne4) - zhaut))
         icl(nm, 4) = -2
         zbas = z(nm) - bm(nm)/2
         if (abs(zbas) < 1e-6) zbas = 0D+00
         valcl(nm, 4) = (rho(nm)*g*(chgbottom(ligne4) - zbas))
         iclt(nm, 4) = -2
         valclt(nm, 4) = tempbottom(ligne4)
      end if
      if (irptha == 0) then
         iclt(1, 3) = -2
         valclt(1, 3) = tempsurf(1)
         icl(1, 3) = -2
         zhaut = z(1) + bm(1)/2
         if (abs(zhaut) < 1e-6) zhaut = 0D+00
         valcl(1, 3) = rho(1)*g*(chgsurf(1) - zhaut)
         icl(nm, 4) = -2
         zbas = z(nm) - bm(nm)/2
         if (abs(zbas) < 1e-6) zbas = 0D+00
         valcl(nm, 4) = rho(nm)*g*(chgbottom(1) - zbas)
         iclt(nm, 4) = -2
         valclt(i, 4) = tempbottom(1)
      end if


   case("1DS") 
      if (kimp < 1) kimp = 1

      if (kimp > ligne4) kimp = ligne4
      valclt(1, 3) = tempsurf(kimp)
      zhaut = z(1) + bm(1)/2
      if (abs(zhaut) < 1e-6) zhaut = 0D+00
      if (icl(1, 3) == -2) valcl(1, 3) = rho(1)*g*(chgsurf(kimp) - zhaut)
      if (icl(1, 3) == -1) valcl(1, 3) = chgsurf(kimp)
      zbas = z(nm) - bm(nm)/2
      if (abs(zbas) < 1e-6) zbas = 0D+00
      if (icl(nm, 4) == -2) valcl(nm, 4) = rho(nm)*g*(chgbottom(kimp) - zbas)
      if (icl(nm, 4) == -1) valcl(nm, 4) = chgbottom(kimp)
      valclt(nm, 4) = tempbottom(kimp)

      if (kimp > ligne4) then
         valclt(1, 3) = tempsurf(ligne4)
         zhaut = z(1) + bm(1)/2
         if (abs(zhaut) < 1e-6) zhaut = 0D+00
         if (icl(1, 3) == -2) valcl(1, 3) = (rho(1)*g*(chgsurf(ligne4) - zhaut))
         if (icl(1, 3) == -1) valcl(1, 3) = chgsurf(ligne4)
         zbas = z(nm) - bm(nm)/2
         if (abs(zbas) < 1e-6) zbas = 0D+00
         if (icl(nm, 4) == -2) then
            valcl(nm, 4) = rho(nm)*g*(chgbottom(ligne4) - zbas)
         end if
         if (icl(nm, 4) == -1) then
            valcl(nm, 4) = chgbottom(ligne4)
         end if
         valclt(nm, 4) = tempbottom(ligne4)
      end if

      if (irptha == 0) then
         iclt(1, 3) = -2
         valclt(1, 3) = tempsurf(1)
         icl(1, 3) = -2
         zhaut = z(1) + bm(1)/2
         if (abs(zhaut) < 1e-6) zhaut = 0D+00
         if (icl(1, 3) == -2) valcl(1, 3) = rho(1)*g*(chgsurf(1) - zhaut)
         if (icl(1, 3) == -1) valcl(1, 3) = chgsurf(1)
         zbas = z(nm) - bm(nm)/2
         if (abs(zbas) < 1e-6) zbas = 0D+00
         if (icl(nm, 4) == -2) valcl(nm, 4) = rho(nm)*g*(chgbottom(1) - zbas)
         if (icl(nm, 4) == -1) valcl(nm, 4) = chgbottom(1)
      end if


   case ( "1DJ") 
      if (kimp < 1) kimp = 1

      if (kimp > ligne4) kimp = ligne4
      zhaut = z(1) + bm(1)/2
      if (abs(zhaut) < 1e-6) zhaut = 0D+00
      if (icl(1, 3) == -2) valcl(1, 3) = rho(1)*g*(chgsurf(kimp) - zhaut)
      if (icl(1, 3) == -1) valcl(1, 3) = chgsurf(kimp)
      zbas = z(nm) - bm(nm)/2
      if (abs(zbas) < 1e-6) zbas = 0D+00
      if (icl(nm, 4) == -2) valcl(nm, 4) = rho(nm)*g*(chgbottom(kimp) - zbas)
      if (icl(nm, 4) == -1) valcl(nm, 4) = chgbottom(kimp)

      if (kimp > ligne4) then
         zhaut = z(1) + bm(1)/2
         if (abs(zhaut) < 1e-6) zhaut = 0D+00
         if (icl(1, 3) == -2) valcl(1, 3) = (rho(1)*g*(chgsurf(ligne4) - zhaut))
         if (icl(1, 3) == -1) valcl(1, 3) = chgsurf(ligne4)
         zbas = z(nm) - bm(nm)/2
         if (abs(zbas) < 1e-6) zbas = 0D+00
         if (icl(nm, 4) == -2) then
            valcl(nm, 4) = rho(nm)*g*(chgbottom(ligne4) - zbas)
         end if
         if (icl(nm, 4) == -1) then
            valcl(nm, 4) = chgbottom(ligne4)
         end if
      end if

      if (irptha == 0) then
         iclt(1, 3) = -2
         icl(1, 3) = -2
         zhaut = z(1) + bm(1)/2
         if (abs(zhaut) < 1e-6) zhaut = 0D+00
         if (icl(1, 3) == -2) valcl(1, 3) = rho(1)*g*(chgsurf(1) - zhaut)
         if (icl(1, 3) == -1) valcl(1, 3) = chgsurf(1)
         zbas = z(nm) - bm(nm)/2
         if (abs(zbas) < 1e-6) zbas = 0D+00
         if (icl(nm, 4) == -2) valcl(nm, 4) = rho(nm)*g*(chgbottom(1) - zbas)
         if (icl(nm, 4) == -1) valcl(nm, 4) = chgbottom(1)
      end if


   case ( "MAQ") 
      if (kimp > ligne4) kimp = ligne4
      do i = 1, nm
!     if(paso >= 86400*1) then
!       if (ivois(i,4) == -99) then
!     icl(i,4)=-2
!     valcl(i,4)=(1-z(i))*rho(i)*g
!     endif
!     endif
!     if(paso > 86400*3) then
!     if (ivois(i,4) == -99) then
!     icl(i,4)=-1
!     valcl(i,4)=0
!     endif
!     endif

         if (ivois(i, 3) == -99) then
            iclt(i, 3) = -2
            valclt(i, 3) = tempsurf(kimp)
         end if
         if (ivois(i, 4) == -99) then
            iclt(i, 4) = -2
            valclt(i, 4) = tempbottom(kimp)
         end if
      end do


   case ("AVA") 

      if (kimp > ligne2) kimp = ligne2
!ccc.... 1 heure
      kheure = int(paso/3600) + 1
      if (kheure > ligne) kheure = ligne
!ccc.... 1 jour
      kjour = int(paso/86400) + 1
      if (kjour > ligne1) kjour = ligne1
      do i = 1, nm
!ccc....CDT LIMITE SOL
         if (ivois(i, 3) == -99) then
            iclt(i, 3) = -2
            icl(i, 3) = -1
            valcl(i, 3) = qpluie(kjour)
!     valclt(i,3)=tempsol(kheure)
         end if
!ccc....CDT LIMITE BOTTOM
         if (ivois(i, 4) == -99) then
            iclt(i, 4) = -1
            icl(i, 4) = -1
!     valclt(i,4)=0
!     valcl(i,4)=0
         end if
!ccc....CDT LIMITE RD

         do j = 1, ligne3
         if (i == id_RD(j)) then
            if (kimp >= ligne2) kimp = ligne2
            icl(i, 1) = -2
            valcl(i, 1) = (chgRD(kimp) - z(i))*rho(i)*g
!     iclt(i,1)=-2
!     valclt(i,1)=tempRD(kimp)-
!     &(tempRD(kimp)-tempsol(kheure))
!     &/(76.83-79.81)*(76.83-z(i))
         end if
         end do

!ccc....CDT LIMITE RG
         do j = 1, ligne4
         if (i == id_RG(j)) then
            if (kimp >= ligne2) kimp = ligne2
            icl(i, 2) = -2
            valcl(i, 2) = (chgRG(kimp) - z(i))*rho(i)*g
!     print*,valcl(i,2),chgRG(kimp)
!     iclt(i,2)=-2
!     valclt(i,2)=tempRG(kimp)-
!     &(tempRG(kimp)-tempsol(kheure))
!     &/(76.83-79.81)*(76.83-z(i))
         end if
         end do
!ccc....CDT LIMITE RIVER TOUJOURS SOUS L'eau
         do j = 1, ligne5
         if (i == id_river(j)) then
            if (kimp >= ligne2) kimp = ligne2
            icl(i, 3) = -2
            zhaut = z(i) + bm(i)/2
            valcl(i, 3) = (chgriver(kimp) - zhaut)*rho(i)*g
!     iclt(i,3)=-2
!     valclt(i,3)=tempriver(kimp)
         end if
         end do

!ccc....CDT LIMITE RIVER a tester
         do j = 1, ligne6
         if (i == id_rivert(j)) then
            if (kimp >= ligne2) kimp = ligne2
            if (z(i) + bm(i)/2 < chgriver(kimp)) then
               icl(i, 3) = -2
               zhaut = z(i) + bm(i)/2
               valcl(i, 3) = (chgriver(kimp) - zhaut)*rho(i)*g
!     iclt(i,3)=-2
!     valclt(i,3)=tempriver(kimp)
            end if
         end if
         end do

      end do


                                              
   case ( "VAU") 
!CC....ECOULEMENT VAUCLIN
      do i = 1, nm
      if (ivois(i, 1) == -99) then
         icl(i, 1) = -2
         valcl(i, 1) = dble(0.65 - z(i))*rho(i)*g
      end if

      if (ivois(i, 2) == -99) then
         icl(i, 2) = -1
         valcl(i, 1) = dble(0)
      end if
      qre = 148D-03/3600
      if (ivois(i, 3) == -99) then
         icl(i, 3) = -1
         valcl(i, 3) = dble(0)
      end if
      if (ivois(i, 3) == -99 .and. x(i) <= 0.5) then
         icl(i, 3) = -1
         valcl(i, 3) = qre
      end if
      if (ivois(i, 4) == -99) then
         icl(i, 4) = -1
         valcl(i, 4) = dble(0)
      end if
      end do
   

   case ( "TEX") 

      if (kimp > ligne2) kimp = ligne2
!ccc.... 1 heure
      kheure = int(paso/3600) + 1
      if (kheure > ligne) kheure = ligne
!ccc.... 1 jour
      kjour = int(paso/86400) + 1
      if (kjour > ligne1) kjour = ligne1
      do i = 1, nm
!ccc....CDT LIMITE SOL
         if (ivois(i, 3) == -99) then
            iclt(i, 3) = -2
            icl(i, 3) = -1
            valcl(i, 3) = qpluie(kjour)
         end if
!ccc....CDT LIMITE BOTTOM
         if (ivois(i, 4) == -99) then
            iclt(i, 4) = -1
            icl(i, 4) = -1
         end if
!ccc....CDT LIMITE RD

         do j = 1, ligne3
         if (i == id_RD(j)) then
            if (kimp >= ligne2) kimp = ligne2
            icl(i, 1) = -2
            valcl(i, 1) = (chgRD(kimp) - z(i))*rho(i)*g
         end if
         end do

!ccc....CDT LIMITE RG
         do j = 1, ligne4
         if (i == id_RG(j)) then
            if (kimp >= ligne2) kimp = ligne2
            icl(i, 2) = -2
            valcl(i, 2) = (chgRG(kimp) - z(i))*rho(i)*g
         end if
         end do
!ccc....CDT LIMITE RIVER TOUJOURS SOUS L'eau
         do j = 1, ligne5
         if (i == id_river(j)) then
            if (kimp >= ligne2) kimp = ligne2
            icl(i, 3) = -2
            zhaut = z(i) + bm(i)/2
            valcl(i, 3) = (chgriver(kimp) - zhaut)*rho(i)*g
         end if
         end do

!ccc....CDT LIMITE RIVER variable dans le temps
         do j = 1, ligne6
         if (i == id_rivert(j)) then
            if (kimp >= ligne2) kimp = ligne2
            if (z(i) + bm(i)/2 < chgriver(kimp)) then
               icl(i, 3) = -2
               zhaut = z(i) + bm(i)/2
               valcl(i, 3) = (chgriver(kimp) - zhaut)*rho(i)*g
            end if
         end if
         end do

      end do
   end select

   return
end


subroutine variation_cdt_limitesDTS(nm, paso, itlecture, ytest, &
                                 ligne, ligne1, ligne2, ligne3, ligne4, ligne5, ligne6, &
                                 icl, valcl, iclt, valclt, ivois, &
                                 z, g, ntsortie, bm, irptha, &
                                 rho, qpluie, chgriver, &
                                 chgRD, chgRG, tempRD, tempRG, &
                                 id_RD, id_RG, id_river, id_rivert, tempsurf, &
                                 tempbottom, chgsurf, chgbottom, &
                                 cRivG, cRivD, &
                                 tempDTS, x, icol, nc, tempo, pro, slopeRH, &
                                 qsurf, qbottom)
!     nc nb de colonnes
!     nr nb de ligne
!       qre debit pluie
!       ivois(ik,1)= voisin droite
!       ivois(ik,2)= voisin gauche
!       ivois(ik,3)= voisin haut
!       ivois(ik,4)= voisin ibas
!       nm nombre de mailles reelles
! ECOULEMENT icl condition valcl valeur
!       ICL=-1 Flux impose sur une face
!       ICL=-2 potentiel impose sur une face
!       ICL=1 Mailles 'normale'
!       ICL(ik,1)=-1, -2 ou 1 ik=1,...4: face droite,gauche,haute et BASSE
! THERMIQUE iclt condition valclt valeur
!       ICLT=-1 Flux impose sur une face
!       ICLT=-2 potentiel impose sur une face
!       ICLT=1 Mailles 'normale'
!       ICLT(ik,1)=-1, -2 ou 1 ik=1,...4: face droite,gauche,haute et BASSE
!       akr(i)=permeabilite relative
!       ak(i)=permeabilite intrinseque

   implicit double precision(a - h, o - z)
   implicit integer(i - n)
   dimension ivois(nm, 4), valcl(nm, 4), icl(nm, 4)
   dimension valclt(nm, 4), iclt(nm, 4), z(nm)
   dimension rho(nm), bm(nm), x(nm)
   dimension icol(nm)
   dimension chgbottom(ntsortie), chgsurf(ntsortie)
   dimension qbottom(ntsortie), qsurf(ntsortie)
   dimension qpluie(ligne1), chgRD(ligne2)
   dimension chgRG(ligne2), tempRD(ligne2)
   dimension id_RD(ligne3), id_RG(ligne4)
   dimension id_river(ligne5), chgriver(ligne2)
   dimension tempRG(ligne1), id_rivert(ligne6)
   dimension tempsurf(ntsortie), tempbottom(ntsortie)
   dimension cRivG(ligne1)
   dimension cRivD(ligne2)
   dimension pro(nm)
   dimension tempo(nm)
   dimension tempDTS(ligne, nc)
   dimension slopeRH(2, ligne3)
   CHARACTER(3) :: ytest
   kimp = int(paso/itlecture)

   select case (ytest)

   case ("WAR")
      if (kimp > ligne4) kimp = ligne4
      select case (icl(nm, 4))
      case (-2)
         valcl(nm, 4) = chgbottom(kimp)*rho(nm)*g
      case (-1)
         valcl(nm, 4) = qbottom(kimp)
      end select
      select case (icl(1, 3))
      case (-2)
         valcl(1, 3) = chgsurf(kimp)*rho(1)*g
      case (-1)
         valcl(1, 3) = qsurf(kimp)
      end select
   case ("ZNS")

      if (kimp > ligne4) kimp = ligne4
      select case (icl(nm, 4))
      case (-2)
         zbas = z(nm) - bm(nm)/2
         if (abs(zbas) < 1e-6) zbas = 0D+00
         valcl(nm, 4) = (rho(nm)*g*(chgbottom(ligne4) - zbas))
      case (-1)
         valcl(nm, 4) = qbottom(kimp)
      end select

      if (ivois(1, 3) == -99) then
         select case (icl(1, 3))
         case (-2)
            zhaut = z(1) + bm(1)/2
            if (abs(zhaut) < 1e-6) zhaut = 0D+00
            valcl(1, 3) = rho(1)*g*(chgsurf(kimp) - zhaut)
         case (-1)
            valcl(1, 3) = qsurf(kimp)
         end select
      end if

   case ("ZND")

      if (kimp > ligne4) kimp = ligne4
      select case (icl(nm, 4))
      case (-2)
         zbas = z(nm) - bm(nm)/2
         if (abs(zbas) < 1e-6) zbas = 0D+00
         valcl(nm, 4) = (rho(nm)*g*(chgbottom(ligne4) - zbas))
      case (-1)
         valcl(nm, 4) = 0
      end select

      if (ivois(1, 3) == -99) then
         zhaut = z(1) + bm(1)/2
         if (chgsurf(kimp) - zhaut >= 0.001) then
            icl(1, 3) = -2
         else
            icl(1, 3) = -1
         end if
         select case (icl(1, 3))
         case (-2)
            print *, 'coucou', chgsurf(kimp)
            zhaut = z(1) + bm(1)/2
            if (abs(zhaut) < 1e-6) zhaut = 0D+00
            valcl(1, 3) = rho(1)*g*(chgsurf(kimp) - zhaut)
         case (-1)
            valcl(1, 3) = 0
         end select
      end if

   end select

   if (ytest == "ZHR" .or. ytest == "ZHZ") then
      if (kimp > ligne4) kimp = ligne4
      iclt(1, 3) = -2
      valclt(1, 3) = tempsurf(kimp)
      icl(1, 3) = -2
      zhaut = z(1) + bm(1)/2
      if (abs(zhaut) < 1e-6) zhaut = 0D+00
      valcl(1, 3) = rho(1)*g*(chgsurf(kimp) - zhaut)
      icl(nm, 4) = -2
      zbas = z(nm) - bm(nm)/2
      if (abs(zbas) < 1e-6) zbas = 0D+00
      valcl(nm, 4) = rho(nm)*g*(chgbottom(kimp) - zbas)
      iclt(nm, 4) = -2
      valclt(nm, 4) = tempbottom(kimp)
      if (kimp > ligne4) then
         iclt(1, 3) = -2
         valclt(1, 3) = tempsurf(ligne4)
         icl(1, 3) = -2
         zhaut = z(1) + bm(1)/2
         if (abs(zhaut) < 1e-6) zhaut = 0D+00
         valcl(1, 3) = (rho(1)*g*(chgsurf(ligne4) - zhaut))
         icl(nm, 4) = -2
         zbas = z(nm) - bm(nm)/2
         if (abs(zbas) < 1e-6) zbas = 0D+00
         valcl(nm, 4) = (rho(nm)*g*(chgbottom(ligne4) - zbas))
         iclt(nm, 4) = -2
         valclt(nm, 4) = tempbottom(ligne4)
      end if
      if (irptha == 0) then
         iclt(1, 3) = -2
         valclt(1, 3) = tempsurf(1)
         icl(1, 3) = -2
         zhaut = z(1) + bm(1)/2
         if (abs(zhaut) < 1e-6) zhaut = 0D+00
         valcl(1, 3) = rho(1)*g*(chgsurf(1) - zhaut)
         icl(nm, 4) = -2
         zbas = z(nm) - bm(nm)/2
         if (abs(zbas) < 1e-6) zbas = 0D+00
         valcl(nm, 4) = rho(nm)*g*(chgbottom(1) - zbas)
         iclt(nm, 4) = -2
         valclt(i, 4) = tempbottom(1)
      end if

   end if

   if (ytest == "1DS") then
      if (kimp < 1) kimp = 1

      if (kimp > ligne4) kimp = ligne4
      valclt(1, 3) = tempsurf(kimp)
      zhaut = z(1) + bm(1)/2
      if (abs(zhaut) < 1e-6) zhaut = 0D+00
      if (icl(1, 3) == -2) valcl(1, 3) = rho(1)*g*(chgsurf(kimp) - zhaut)
      if (icl(1, 3) == -1) valcl(1, 3) = chgsurf(kimp)
      zbas = z(nm) - bm(nm)/2
      if (abs(zbas) < 1e-6) zbas = 0D+00
      if (icl(nm, 4) == -2) valcl(nm, 4) = rho(nm)*g*(chgbottom(kimp) - zbas)
      if (icl(nm, 4) == -1) valcl(nm, 4) = chgbottom(kimp)
      valclt(nm, 4) = tempbottom(kimp)

      if (kimp > ligne4) then
         valclt(1, 3) = tempsurf(ligne4)
         zhaut = z(1) + bm(1)/2
         if (abs(zhaut) < 1e-6) zhaut = 0D+00
         if (icl(1, 3) == -2) valcl(1, 3) = (rho(1)*g*(chgsurf(ligne4) - zhaut))
         if (icl(1, 3) == -1) valcl(1, 3) = chgsurf(ligne4)
         zbas = z(nm) - bm(nm)/2
         if (abs(zbas) < 1e-6) zbas = 0D+00
         if (icl(nm, 4) == -2) then
            valcl(nm, 4) = rho(nm)*g*(chgbottom(ligne4) - zbas)
         end if
         if (icl(nm, 4) == -1) then
            valcl(nm, 4) = chgbottom(ligne4)
         end if
         valclt(nm, 4) = tempbottom(ligne4)
      end if

      if (irptha == 0) then
         iclt(1, 3) = -2
         valclt(1, 3) = tempsurf(1)
         icl(1, 3) = -2
         zhaut = z(1) + bm(1)/2
         if (abs(zhaut) < 1e-6) zhaut = 0D+00
         if (icl(1, 3) == -2) valcl(1, 3) = rho(1)*g*(chgsurf(1) - zhaut)
         if (icl(1, 3) == -1) valcl(1, 3) = chgsurf(1)
         zbas = z(nm) - bm(nm)/2
         if (abs(zbas) < 1e-6) zbas = 0D+00
         if (icl(nm, 4) == -2) valcl(nm, 4) = rho(nm)*g*(chgbottom(1) - zbas)
         if (icl(nm, 4) == -1) valcl(nm, 4) = chgbottom(1)
      end if

   end if

   if (ytest == "1DJ") then
      if (kimp < 1) kimp = 1

      if (kimp > ligne4) kimp = ligne4
      zhaut = z(1) + bm(1)/2
      if (abs(zhaut) < 1e-6) zhaut = 0D+00
      if (icl(1, 3) == -2) valcl(1, 3) = rho(1)*g*(chgsurf(kimp) - zhaut)
      if (icl(1, 3) == -1) valcl(1, 3) = chgsurf(kimp)
      zbas = z(nm) - bm(nm)/2
      if (abs(zbas) < 1e-6) zbas = 0D+00
      if (icl(nm, 4) == -2) valcl(nm, 4) = rho(nm)*g*(chgbottom(kimp) - zbas)
      if (icl(nm, 4) == -1) valcl(nm, 4) = chgbottom(kimp)

      if (kimp > ligne4) then
         zhaut = z(1) + bm(1)/2
         if (abs(zhaut) < 1e-6) zhaut = 0D+00
         if (icl(1, 3) == -2) valcl(1, 3) = (rho(1)*g*(chgsurf(ligne4) - zhaut))
         if (icl(1, 3) == -1) valcl(1, 3) = chgsurf(ligne4)
         zbas = z(nm) - bm(nm)/2
         if (abs(zbas) < 1e-6) zbas = 0D+00
         if (icl(nm, 4) == -2) then
            valcl(nm, 4) = rho(nm)*g*(chgbottom(ligne4) - zbas)
         end if
         if (icl(nm, 4) == -1) then
            valcl(nm, 4) = chgbottom(ligne4)
         end if
      end if

      if (irptha == 0) then
         iclt(1, 3) = -2
         icl(1, 3) = -2
         zhaut = z(1) + bm(1)/2
         if (abs(zhaut) < 1e-6) zhaut = 0D+00
         if (icl(1, 3) == -2) valcl(1, 3) = rho(1)*g*(chgsurf(1) - zhaut)
         if (icl(1, 3) == -1) valcl(1, 3) = chgsurf(1)
         zbas = z(nm) - bm(nm)/2
         if (abs(zbas) < 1e-6) zbas = 0D+00
         if (icl(nm, 4) == -2) valcl(nm, 4) = rho(nm)*g*(chgbottom(1) - zbas)
         if (icl(nm, 4) == -1) valcl(nm, 4) = chgbottom(1)
      end if

   end if

   if (ytest == "MAQ") then
      if (kimp > ligne4) kimp = ligne4
      do i = 1, nm
!     if(paso >= 86400*1) then
!       if (ivois(i,4) == -99) then
!     icl(i,4)=-2
!     valcl(i,4)=(1-z(i))*rho(i)*g
!     endif
!     endif
!     if(paso > 86400*3) then
!     if (ivois(i,4) == -99) then
!     icl(i,4)=-1
!     valcl(i,4)=0
!     endif
!     endif

         if (ivois(i, 3) == -99) then
            iclt(i, 3) = -2
            valclt(i, 3) = tempsurf(kimp)
         end if
         if (ivois(i, 4) == -99) then
            iclt(i, 4) = -2
            valclt(i, 4) = tempbottom(kimp)
         end if
      end do
   end if

   if (ytest == "AVA") then

      if (kimp > ligne2) kimp = ligne2
!ccc.... 1 heure
      kheure = int(paso/3600) + 1
      if (kheure > ligne) kheure = ligne
!ccc.... 1 jour
      kjour = int(paso/86400) + 1
      if (kjour > ligne1) kjour = ligne1
      do i = 1, nm
!ccc....CDT LIMITE SOL
         if (ivois(i, 3) == -99) then
            iclt(i, 3) = -2
            icl(i, 3) = -1
            valcl(i, 3) = qpluie(kjour)
!     valclt(i,3)=tempsol(kheure)
         end if
!ccc....CDT LIMITE BOTTOM
         if (ivois(i, 4) == -99) then
            iclt(i, 4) = -1
            icl(i, 4) = -1
!     valclt(i,4)=0
!     valcl(i,4)=0
         end if
!ccc....CDT LIMITE RD

         do j = 1, ligne3
         if (i == id_RD(j)) then
            if (kimp >= ligne2) kimp = ligne2
            icl(i, 1) = -2
            valcl(i, 1) = (chgRD(kimp) - z(i))*rho(i)*g
!     iclt(i,1)=-2
!     valclt(i,1)=tempRD(kimp)-
!     &(tempRD(kimp)-tempsol(kheure))
!     &/(76.83-79.81)*(76.83-z(i))
         end if
         end do

!ccc....CDT LIMITE RG
         do j = 1, ligne4
         if (i == id_RG(j)) then
            if (kimp >= ligne2) kimp = ligne2
            icl(i, 2) = -2
            valcl(i, 2) = (chgRG(kimp) - z(i))*rho(i)*g
!     print*,valcl(i,2),chgRG(kimp)
!     iclt(i,2)=-2
!     valclt(i,2)=tempRG(kimp)-
!     &(tempRG(kimp)-tempsol(kheure))
!     &/(76.83-79.81)*(76.83-z(i))
         end if
         end do
!ccc....CDT LIMITE RIVER TOUJOURS SOUS L'eau
         do j = 1, ligne5
         if (i == id_river(j)) then
            if (kimp >= ligne2) kimp = ligne2
            icl(i, 3) = -2
            zhaut = z(i) + bm(i)/2
            valcl(i, 3) = (chgriver(kimp) - zhaut)*rho(i)*g
!     iclt(i,3)=-2
!     valclt(i,3)=tempriver(kimp)
         end if
         end do

!ccc....CDT LIMITE RIVER a tester
         do j = 1, ligne6
         if (i == id_rivert(j)) then
            if (kimp >= ligne2) kimp = ligne2
            if (z(i) + bm(i)/2 < chgriver(kimp)) then
               icl(i, 3) = -2
               zhaut = z(i) + bm(i)/2
               valcl(i, 3) = (chgriver(kimp) - zhaut)*rho(i)*g
!     iclt(i,3)=-2
!     valclt(i,3)=tempriver(kimp)
            end if
         end if
         end do

      end do
   end if

   if (ytest == "VAU") then
!CC....ECOULEMENT VAUCLIN
      do i = 1, nm
      if (ivois(i, 1) == -99) then
         icl(i, 1) = -2
         valcl(i, 1) = dble(0.65 - z(i))*rho(i)*g
      end if

      if (ivois(i, 2) == -99) then
         icl(i, 2) = -1
         valcl(i, 1) = dble(0)
      end if
      qre = 148D-03/3600
      if (ivois(i, 3) == -99) then
         icl(i, 3) = -1
         valcl(i, 3) = dble(0)
      end if
      if (ivois(i, 3) == -99 .and. x(i) <= 0.5) then
         icl(i, 3) = -1
         valcl(i, 3) = qre
      end if
      if (ivois(i, 4) == -99) then
         icl(i, 4) = -1
         valcl(i, 4) = dble(0)
      end if
      end do
   end if

   if (ytest == "DTS") then

      if (kimp >= ligne2) kimp = ligne2
      if (kimp >= ligne1) kimp = ligne1
      if (kimp < 1) kimp = 1
!     headG=cRivG(kimp)
!     headPG=chgRG(kimp)
!     temperatureG=tempRG(kimp)

!cc side boundary petit paris and bertin
      do i = 1, nm

!ccPetit Paris
         if (ivois(i, 1) == -99 .and. x(i) > 999.5) then
         if (ivois(i, 4) .ne. -99 .or. ivois(i, 3) .ne. -99) then
            iclt(i, 1) = -2
            valclt(i, 1) = tempo(i)
         end if
         if (z(i) < 1) then
            icl(i, 1) = -2
            valcl(i, 1) = (chgRD(kimp) - z(i) - 4.188)*rho(i)*g
         else
            icl(i, 1) = -2
            valcl(i, 1) = pro(i)
         end if
         end if

!ccBertin
         if (ivois(i, 2) == -99 .and. x(i) < 0.5) then
         if (ivois(i, 4) .ne. -99 .or. ivois(i, 3) .ne. -99) then
            iclt(i, 2) = -2
            valclt(i, 2) = tempo(i)
         end if

         if (z(i) < 1) then
            icl(i, 2) = -2
            valcl(i, 2) = (chgRG(kimp) - z(i))*rho(i)*g
         else
            icl(i, 2) = -2
            valcl(i, 2) = pro(i)
         end if
         end if

!ccc....CDT LIMITE BOTTOM for the water zero flux
         if (ivois(i, 4) == -99) then
!cc heat transport
            iclt(i, 4) = -2
            valclt(i, 4) = tempRG(kimp)
!cc water flow
            icl(i, 4) = -1
            valcl(i, 4) = 0
            if (ivois(i, 1) == -99) then
!cc corner cell
               iclt(i, 1) = -2
               valclt(i, 1) = tempRG(kimp)
               icl(i, 1) = -2
               valcl(i, 1) = pro(i)
            end if
            if (ivois(i, 2) == -99) then
!cc corner cell
               iclt(i, 2) = -2
               valclt(i, 2) = tempRG(kimp)
               icl(i, 2) = -2
               valcl(i, 2) = pro(i)
            end if
         end if
      end do

!ccc....Boundary condition river/ZH for the heat transport
      if (kimp >= ligne) kimp = ligne
      do i = 1, nm
         if (icol(i) > nc) icol(i) = nc
         if (ivois(i, 3) == -99) then
            temperature = tempDTS(kimp, icol(i))
            icl(i, 3) = -2
            valclt(i, 3) = temperature
            if (ivois(i, 2) == -99) then
               iclt(i, 2) = -2
               valclt(i, 2) = temperature
            end if
            if (ivois(i, 1) == -99) then
               iclt(i, 1) = -2
               valclt(i, 1) = temperature
            end if
         end if
      end do

!CCC River/ZH hydraulic head
      do i = 1, nm
         headbk = 0
!ccc loop on the element
         if (ivois(i, 3) == -99) then
!ccc no neigh top
            icl(i, 3) = -2
!ccc boundary diricler for the water flow
            do j = 1, ligne3
!cc loop on the piece of the line of the water level of the river
!ccc slope for each element
               if (x(i) <= 392) slope = slopeRH(2, 1)
               if (x(i) <= 392) kj = 1
               if (x(i) >= 392 .and. x(i) < 397) slope = slopeRH(2, 2)
               if (x(i) >= 392 .and. x(i) < 397) kj = 2
               if (x(i) >= 397 .and. x(i) < 929) slope = slopeRH(2, 3)
               if (x(i) >= 397 .and. x(i) < 929) kj = 3
               if (x(i) >= 929 .and. x(i) < 931) slope = slopeRH(2, 4)
               if (x(i) >= 929 .and. x(i) < 931) kj = 4
!cc id = number of upstream part
               if (x(i) > 931) then
                  kj = 5
                  slope = slopeRH(2, 5)
               end if
            end do

!cc slope from bertin water level for the downstream part
            if (x(i) < 392) slope = slopeRH(2, 1)
            if (x(i) < 392) headbk = cRivG(kimp)

!cc between the downstream and the first cascade
            if (kj == 2) then
            do k = 1, kj
               if (k == 1) headbk = cRivG(kimp)
               if (k > 1) headbk = slopeRH(2, k - 1)* &
                                   (slopeRH(1, k) - slopeRH(1, k - 1)) + headbk
            end do
            end if
!ccc water level of the river = altitude of the first cascade
            if (kj == 3) headbk = 116.29999694824218D+00
!ccc water level of the river = altitude of the second cascade
            if (kj == 4) headbk = 120.5D+00
!ccc water level of the upstream part
            if (kj == ligne3) headbk = cRivD(kimp) - 2.16

!     if(kj == ligne3.or.kj == 3) then
            if (kj == ligne3) then
               slopeRH(2, kj) = (cRivD(kimp) - 2.16 - 120.5D+00)/ &
                                (1000 - slopeRH(1, ligne3))
               slope = slopeRH(2, kj)
            end if

!cc slope inside the first cascade
            if (kj == 2) slope = (headbk - 116.29999694824218D+00)/ &
                                 (slopeRH(1, 2) - slopeRH(1, 3))

!cc Inside the second cascade
!cc calculatation the slope
            if (kj == 4) then
               headbk = slopeRH(2, 3)* &
                        (slopeRH(1, kj) - slopeRH(1, kj - 1)) + 116.29999694824218D+00
               slope = (headbk - 120.5D+00)/(slopeRH(1, kj) - slopeRH(1, ligne3))
               slopeRH(2, kj) = slope
            end if

            if (kj <= 4) then
               head = slope*(x(i) - slopeRH(1, kj)) + headbk
            end if

            if (kj == ligne3) head = slope*(x(i) - 1000) + headbk
            valcl(i, 3) = (head - z(i) - bm(i)/2)*rho(i)*g

            if (ivois(i, 2) == -99) then
               icl(i, 2) = -2
               valcl(i, 2) = (head - z(i))*rho(i)*g
            end if
            if (ivois(i, 1) == -99 .and. x(i) < 999.5) then
               icl(i, 1) = -2
               valcl(i, 1) = (head - z(i))*rho(i)*g
            end if
            if (ivois(ivois(i, 4), 2) == -99 .and. x(i) > 0.5) then
               icl(i, 2) = -2
               valcl(i, 2) = (head - z(i))*rho(i)*g
            end if
         end if

!ccc....Boundary UPSTREAM
         if (ivois(i, 1) == -99 .and. x(i) > 999.5) then
            icl(i, 1) = -2
            iclt(i, 1) = -2
            valcl(i, 1) = (chgRD(kimp) - z(i))*rho(i)*g
            valclt(i, 1) = tempRD(kimp)
         end if

!ccc....Boundary DOWNSTREAM
         if (ivois(i, 2) == -99 .and. x(i) < 0.5) then
            icl(i, 2) = -2
            iclt(i, 2) = -2
            valcl(i, 2) = (chgRD(kimp) - z(i))*rho(i)*g
            valclt(i, 2) = tempRG(kimp)
         end if

      end do

   end if

   if (ytest == "TEX") then

      if (kimp > ligne2) kimp = ligne2
!ccc.... 1 heure
      kheure = int(paso/3600) + 1
      if (kheure > ligne) kheure = ligne
!ccc.... 1 jour
      kjour = int(paso/86400) + 1
      if (kjour > ligne1) kjour = ligne1
      do i = 1, nm
!ccc....CDT LIMITE SOL
         if (ivois(i, 3) == -99) then
            iclt(i, 3) = -2
            icl(i, 3) = -1
            valcl(i, 3) = qpluie(kjour)
         end if
!ccc....CDT LIMITE BOTTOM
         if (ivois(i, 4) == -99) then
            iclt(i, 4) = -1
            icl(i, 4) = -1
         end if
!ccc....CDT LIMITE RD

         do j = 1, ligne3
         if (i == id_RD(j)) then
            if (kimp >= ligne2) kimp = ligne2
            icl(i, 1) = -2
            valcl(i, 1) = (chgRD(kimp) - z(i))*rho(i)*g
         end if
         end do

!ccc....CDT LIMITE RG
         do j = 1, ligne4
         if (i == id_RG(j)) then
            if (kimp >= ligne2) kimp = ligne2
            icl(i, 2) = -2
            valcl(i, 2) = (chgRG(kimp) - z(i))*rho(i)*g
         end if
         end do
!ccc....CDT LIMITE RIVER TOUJOURS SOUS L'eau
         do j = 1, ligne5
         if (i == id_river(j)) then
            if (kimp >= ligne2) kimp = ligne2
            icl(i, 3) = -2
            zhaut = z(i) + bm(i)/2
            valcl(i, 3) = (chgriver(kimp) - zhaut)*rho(i)*g
         end if
         end do

!ccc....CDT LIMITE RIVER variable dans le temps
         do j = 1, ligne6
         if (i == id_rivert(j)) then
            if (kimp >= ligne2) kimp = ligne2
            if (z(i) + bm(i)/2 < chgriver(kimp)) then
               icl(i, 3) = -2
               zhaut = z(i) + bm(i)/2
               valcl(i, 3) = (chgriver(kimp) - zhaut)*rho(i)*g
            end if
         end if
         end do

      end do
   end if

   return


end
