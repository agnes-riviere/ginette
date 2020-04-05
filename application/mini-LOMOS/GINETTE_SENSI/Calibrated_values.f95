program tested
        implicit double precision(a-h,o-x,z),integer(I-N)
        implicit CHARACTER*5(y)
	  integer nbzone
      double precision, allocatable :: ak(:),l(:),r(:)
      double precision, allocatable :: kmin(:),kmax(:),kpas(:)
      double precision, allocatable :: lmin(:),lmax(:),lpas(:)
      double precision, allocatable :: rmin(:),rmax(:),rpas(:)
      double precision, allocatable :: z(:),x(:),am(:),bm(:)
      double precision, allocatable :: topo(:),bot(:)
	  integer, allocatable :: inbkpas(:),inbompas(:),inblpas(:),izone(:)
      integer, allocatable :: inbrpas(:),inum(:),inum2(:)
      Real, allocatable ::  ommin(:),om(:),ommax(:),ompas(:)
      open(unit=4,file='E_nb_zone.dat',action='read')
	  READ(4,*)nbzone,azbound
      allocate(ak(nbzone))
      allocate(l(nbzone))
      allocate(r(nbzone))
      allocate(kmin(nbzone))
      allocate(kmax(nbzone))
      allocate(kpas(nbzone))
      allocate(lmin(nbzone))
      allocate(lmax(nbzone))
      allocate(lpas(nbzone))
      allocate(rmin(nbzone))
      allocate(rmax(nbzone))
      allocate(rpas(nbzone))
      allocate(inbkpas(nbzone))
      allocate(inbompas(nbzone))
      allocate(inblpas(nbzone))
      allocate(inbrpas(nbzone))
      allocate(ommin(nbzone))
      allocate(om(nbzone))
      allocate(ommax(nbzone))
      allocate(ompas(nbzone))


      open(unit=1,file='E_tested_ranges.dat',action='read')
      open(unit=2,file='tested_values',action='write')
	  open(unit=22,file='E_zone.dat',action='write')


      call lecture_parametre(ilog,ipermh,ipermv,iec,irp,ith,dt,nitt,unitsim,al,&
     &az,ixy,imaille,itopo,reptop,repbot,icolonne,nmi,nci,nri,dx,dz,nclog,rho1,irho,g,amu,&
     &akrx,akx,akrz,akz,iom,omp,iss,sss,ia2,yunconfined,ivg,ans,asp,swres,itr,allg,alt,iriv,iqriv,qre,hbot,xberg,&
     &rug,pent,qriva,hriv,aklit,aklitv,tempriv,elit,akdrain,edrain,crconvp,crconvc,&
     &iteration,itsortie,unitsortie,icalvit,mcol,nmaille1,nmaille2,nmaille3,nmaille4,nmaille5,&
     &nmaille6,nmaille7,nmaille8,nmaille9,nmaille10,iparo,ibuilt)

      allocate(z(nmi))
      allocate(x(nmi))
      allocate(am(nmi))
      allocate(bm(nmi))
      allocate(inum(nmi))
      allocate(inum2(nmi))
      allocate(topo(nmi))
      allocate(bot(nmi))
      allocate(izone(nmi))

	  do i=1,nbzone
      READ(1,*)kmin(i),kmax(i),inbkpas(i)
	  enddo
	  do i=1,nbzone
      READ(1,*)ommin(i),ommax(i),inbompas(i)
	  enddo
	  do i=1,nbzone
      READ(1,*)lmin(i),lmax(i),inblpas(i)
	  enddo
	  do i=1,nbzone
      READ(1,*)rmin(i),rmax(i),inbrpas(i)
	  enddo

	  if (nbzone.eq.1) then
	  do i=1,nbzone
	  if(inbompas(i).gt.1) 	ompas(i)=dble((ommax(i)-ommin(i))/(inbompas(i)-1))
	  if(inbkpas(i).gt.1) 	kpas(i)=dble((kmax(i)-kmin(i))/(inbkpas(i)-1))
	  if(inblpas(i).gt.1) 	lpas(i)=dble((lmax(i)-lmin(i))/(inblpas(i)-1))
	  if(inbrpas(i).gt.1) 	rpas(i)=dble((rmax(i)-rmin(i))/(inbrpas(i)-1))
	  if(inbompas(i).le.1) ompas(i)=0
	  if(inbkpas(i).le.1) kpas(i)=0
	  if(inblpas(i).le.1) lpas(i)=0
	  if(inbrpas(i).le.1) rpas(i)=0
	  enddo

      do iin=0,inbompas(1)-1
		om(1)=ommin(1)+iin*ompas(1)
      	do ik=0,inbkpas(1)-1
	  	ak(1)=kmin(1)+dble(ik)*kpas(1)
      		do ij=0,inblpas(1)-1
	  		l(1)=lmin(1)+dble(ij)*lpas(1)
      			do ir=0,inbrpas(1)-1
	 			r(1)=rmin(1)+dble(ir)*rpas(1)
      write(2,20)ak(1),om(1),l(1),r(1)

      enddo
      	enddo
      		enddo
      			enddo
20    format(d9.2,f6.3,2d10.2)
	  else
	  do i=1,nbzone
	  if(inbompas(i).gt.1) 	ompas(i)=dble((ommax(i)-ommin(i))/(inbompas(i)-1))
	  if(inbkpas(i).gt.1) 	kpas(i)=dble((kmax(i)-kmin(i))/(inbkpas(i)-1))
	  if(inblpas(i).gt.1) 	lpas(i)=dble((lmax(i)-lmin(i))/(inblpas(i)-1))
	  if(inbrpas(i).gt.1) 	rpas(i)=dble((rmax(i)-rmin(i))/(inbrpas(i)-1))
	  if(inbompas(i).le.1) ompas(i)=0
	  if(inbkpas(i).le.1) kpas(i)=0
	  if(inblpas(i).le.1) lpas(i)=0
	  if(inbrpas(i).le.1) rpas(i)=0
	  enddo

      do iin=0,inbompas(1)-1
		om(1)=ommin(1)+iin*ompas(1)
      	do ik=0,inbkpas(1)-1
	  	ak(1)=kmin(1)+dble(ik)*kpas(1)
      		do ij=0,inblpas(1)-1
	  		l(1)=lmin(1)+dble(ij)*lpas(1)
      			do ir=0,inbrpas(1)-1
	 			r(1)=rmin(1)+dble(ir)*rpas(1)
      				do iin2=0,inbompas(2)-1
					om(2)=ommin(2)+iin2*ompas(2)
      					do ik2=0,inbkpas(2)-1
	  					ak(2)=kmin(2)+dble(ik2)*kpas(2)
      						do ij2=0,inblpas(2)-1
	  						l(2)=lmin(2)+dble(ij2)*lpas(2)
      							do ir2=0,inbrpas(2)-1
	  							r(2)=rmin(2)+dble(ir2)*rpas(2)
      write(2,30)ak(1),om(1),l(1),r(1),ak(2),om(2),l(2),r(2)

      enddo
      	enddo
      		enddo
      			enddo
      				enddo
      					enddo
      						enddo
      							enddo

30    format(d9.2,f6.3,2d10.2,d9.2,f6.3,2d10.2)
      endif


!!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C                                             C
!C                       MAILLAGE              C
!C                                             C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCC....CALCUL LIGNES COLONNES
        nc=nint(abs(al)/dx)
        nr=nint(abs(az)/dz)

!!!....INTERFROST TEST TH2
        if(ytest.eq."TH2") then
          dx=dble(0.1D+00/12D+00)
          dz=dble(0.1D+00/12D+00)
        endif
!!!....MAILLAGE MANUEL
        if (imaille.eq.1) then
!!!C....DEFINITION DU MAILLAGE BRUT
!!!cc....Premiere possibilite maillage a pas d espace dx variable
!!!cc....Progression logarithmique --> plus de precision au contact
!!!cc....reseau de surface
!!!cc....Progression log du pas d espace en X ILOG+1
        if(ilog.eq.1) then
        sum=0
        nc=nclog
        do ii=1,nc
        at=ii
        sum=sum+log(2.*at)
        enddo
        pas=abs(al)/sum
        do i=1,nc
        at=nc-i+1
        am(i)=dble(pas*log(2.*at))
        bm(i)=dble(dz)
        enddo
        else
        do i=1,nc*nr
        am(i)=dble(dx)
        bm(i)=dble(dz)
        enddo
        endif
        k=0
        do i=1,nr
        do j=1,nc
        k=k+1
        am(k)=dble(am(j))
        inum(k)=0
        inum2(k)=0

!!!cc....definition des coordonnees
!!!cc....modif  23/03/2011
!!!cc....log en z
        if(ilog.eq.1) then
        if(j.eq.1) then
        x(k)=(am(k)/2.D00)
        else
        x(k)=(x(k-1)+(am(k-1)+am(k))/2.D00)
        endif
        else
        x(k)=(dx/2.D00+(j-1)*dx)
        endif
        if(i.eq.1) then
        z(k)=(az-bm(k)/2.D00)
!!!cc....Modif 19-09-2014- point repere haut de colonne pour modele 1D
        if (itopo.eq.0.and.imaille.eq.1) then
        if (abs(reptop-repbot)+1D-08.gt.az.and.abs(reptop-repbot)-1D-08.lt.az) then
         z(k)=(az-bm(k)/2.D00+reptop-az)
        else
        print*, 'probleme dans votre constuction de modele'
        print*,'reptop-repbot=',reptop-repbot+1D-08,'et az=',az
        stop
        endif

        endif
        else
        z(k)=(z(k-nc)-(bm(k-nc)+bm(k))/2.D00)
        endif

        enddo
        enddo
!!!....RENUMEROTATION en FONCTION DE TOPO ET BOTTOM
        if (itopo.eq.1) then
        open(unit=10,file='E_geom.dat',form='formatted',status='old')
        do i=1,nc
        read(10,*) topo(i),bot(i)
        enddo
        else
        do i=1,nc
        topo(i)=(reptop)
        bot(i)=(repbot)
        enddo
        endif
        kr=0
        do i=1,nr
        do j=1,nc
        ii=(i-1)*nc+j
        if(z(ii).le.topo(j).and.z(ii).ge.bot(j)) then
        kr=kr+1
!!!....NOUVEAU NUMERO de MAILLE
        inum(kr)=ii
        inum2(ii)=kr
        endif
        enddo
        enddo
!!!....NM nombre de mailles reelles!!
        nm=kr
!!!....RECALCUL DE LA GEOMETRIE ET DU TABLEAU DE VOISINAGE
        do i=1,nm
        x(i)=x(inum(i))
        z(i)=z(inum(i))
        am(i)=am(inum(i))
        bm(i)=bm(inum(i))

!       write(11,*) i,x(i),z(i)

        enddo
        endif
		do i=1,nm
		if(nbzone.eq.1) then
		izone(i)=1
		else
		if (z(i).ge.-azbound) izone(i)=1
		if (z(i).lt.-azbound) izone(i)=2
		endif
		enddo
		do i=1,nm
		write(22,*)izone(i)
		nzone=max(nzone,izone(i))
		enddo

end

        subroutine lecture_parametre(ilog,ipermh,ipermv,iec,irp,ith,dt,nitt,unitsim,al,&
     &az,ixy,imaille,itopo,reptop,repbot,icolonne,nmi,nci,nri,dx,dz,nclog,rho1,irho,g,amu,&
     &akrx,akx,akrz,akz,iom,omp,iss,sss,ia2,yunconfined,ivg,ans,asp,swres,itr,allg,alt,iriv,iqriv,qre,hbot,xberg,&
     &rug,pent,qriva,hriv,aklit,aklitv,tempriv,elit,akdrain,edrain,crconvp,crconvc,&
     &iteration,itsortie,unitsortie,icalvit,mcol,nmaille1,nmaille2,nmaille3,nmaille4,nmaille5,&
     &nmaille6,nmaille7,nmaille8,nmaille9,nmaille10,iparo,ibuilt)
        implicit double precision(a-h,o-x,z),integer(I-N)
        implicit CHARACTER*5(y)
        open(3,file='E_parametre.dat',form='formatted',status='old')
        READ(3,'(4x,i1)')iec
        READ(3,'(3x,i1)')irp
        READ(3,'(4x,i1)')ith
        READ(3,'(3x,d9.0)')dt
        READ(3,'(5x,i10)')nitt
        READ(3,'(6x,d9.0)')unitsim
        READ(3,'(3x,d8.0)')al
        READ(3,'(3x,d9.0)')az
        READ(3,'(4x,i1)')ixy
        READ(3,'(8x,i1)')imaille
        READ(3,'(6x,i1)')itopo
        READ(3,'(7x,i1)')ibuilt
        READ(3,'(7x,d9.0)')reptop
        READ(3,'(7x,d9.0)')repbot
        READ(3,'(9x,i1)')icolonne
        READ(3,'(4x,i5)')nmi
        READ(3,'(4x,i5)')nci
        READ(3,'(4x,i5)')nri
        READ(3,'(3x,d8.0)')dx
        READ(3,'(3x,d8.0)')dz
        READ(3,'(6x,i5)')nclog
        READ(3,'(5x,i1)')ilog
        READ(3,'(5x,d8.0)')rho1
        READ(3,'(5x,i1)')irho
        READ(3,'(2x,d8.0)')g
        READ(3,'(4x,d8.0)')amu
        READ(3,'(7x,i1)')ipermh
        READ(3,'(7x,i1)')ipermv
        READ(3,'(5x,d8.0)')akrx
        READ(3,'(4x,d8.0)')akx
        READ(3,'(5x,d8.0)')akrz
        READ(3,'(4x,d8.0)')akz
        READ(3,'(4x,i1)')iom
        READ(3,'(4x,f5.3)')omp
        READ(3,'(4x,i1)')iss
        READ(3,'(3x,d8.0)')sss
        READ(3,'(4x,i1)')ia2
        READ(3,'(12x,A3)')yunconfined
        READ(3,'(4x,i1)')ivg
        READ(3,'(4x,d8.0)')ans
        READ(3,'(4x,d8.0)')asp
        READ(3,'(6x,f6.4)')swres
        READ(3,'(4x,i1)')itr
        READ(3,'(4x,d8.0)')allg
        READ(3,'(4x,d8.0)')alt
        READ(3,'(5x,i1)')iriv
        READ(3,'(6x,i1)')iqriv
        READ(3,'(4x,d8.0)')qre
        READ(3,'(5x,d8.0)')hbot
        READ(3,'(6x,d8.0)')xberg
        READ(3,'(4x,d8.0)')rug
        READ(3,'(5x,d8.0)')pent
        READ(3,'(6x,d8.0)')qriva
        READ(3,'(5x,d8.0)')hriv
        READ(3,'(6x,d8.0)')aklit
        READ(3,'(7x,d8.0)')aklitv
        READ(3,'(8x,d8.0)')tempriv
        READ(3,'(5x,d8.0)')elit
        READ(3,'(8x,d8.0)')akdrain
        READ(3,'(7x,d8.0)')edrain
        READ(3,'(8x,d8.0)')crconvp
        READ(3,'(8x,d8.0)')crconvc
        READ(3,'(10x,i8)')iteration
        READ(3,'(9x,i8)')itsortie
        READ(3,'(11x,d9.0)')unitsortie
        READ(3,'(8x,i1)')icalvit
        READ(3,'(5x,i5)')mcol
        READ(3,'(9x,i5)')nmaille1
        READ(3,'(9x,i5)')nmaille2
        READ(3,'(9x,i5)')nmaille3
        READ(3,'(9x,i5)')nmaille4
        READ(3,'(9x,i5)')nmaille5
        READ(3,'(9x,i5)')nmaille6
        READ(3,'(9x,i5)')nmaille7
        READ(3,'(9x,i5)')nmaille8
        READ(3,'(9x,i5)')nmaille9
        READ(3,'(10x,i5)')nmaille10
        READ(3,'(6x,i1)')iparo
        close(3)
        return
        end


