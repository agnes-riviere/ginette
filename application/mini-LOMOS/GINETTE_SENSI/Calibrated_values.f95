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
      open(unit=1,file='E_tested_ranges.dat',action='read')
      open(unit=2,file='tested_values',action='write')
	  open(unit=22,file='E_zone.dat',action='write')





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


		call pas_real(nbzone,ompas,ommax,ommin,inbompas)
		call pas_dble(nbzone,lpas,lmax,lmin,inblpas)
		call pas_dble(nbzone,rpas,rmax,rmin,inbrpas)
		call pas_log(nbzone,kpas,kmax,kmin,inbkpas)

	  if (nbzone.eq.1) then
      do iin=0,inbompas(1)-1
		om(1)=ommin(1)+iin*ompas(1)
      	do ik=0,inbkpas(1)-1
	  ak(1)=dble(10**(log10(kmax(1))-ik*kpas(1)))	  	
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

      do iin=0,inbompas(1)-1
		om(1)=ommin(1)+iin*ompas(1)
      	do ik=0,inbkpas(1)-1
	  	ak(1)=dble(10**(log10(kmax(1))-ik*kpas(1)))	
      		do ij=0,inblpas(1)-1
	  		l(1)=lmin(1)+dble(ij)*lpas(1)
      			do ir=0,inbrpas(1)-1
	 			r(1)=rmin(1)+dble(ir)*rpas(1)
      				do iin2=0,inbompas(2)-1
					om(2)=ommin(2)+iin2*ompas(2)
      					do ik2=0,inbkpas(2)-1
	  					ak(2)=dble(10**(log10(kmax(2))-ik2*kpas(2)))	
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



end

subroutine pas_real(nbzone,ompas,ommax,ommin,inbompas)
	  integer nbzone,i
      integer inbompas(nbzone)
      real ommin(nbzone),ommax(nbzone),ompas(nbzone)

	  do i=1,nbzone
	  if(inbompas(i).gt.1) 	then 
		ompas(i)=dble((ommax(i)-ommin(i))/(inbompas(i)-1))
		 else
	  ompas(i)=0
	  endif
	  enddo
return
end 


subroutine pas_dble(nbzone,ompas,ommax,ommin,inbompas)
	  integer nbzone,i
      integer inbompas(nbzone)
      double precision ommin(nbzone),ommax(nbzone),ompas(nbzone)
	  do i=1,nbzone
	  if(inbompas(i).gt.1) 	then 
		ompas(i)=dble((ommax(i)-ommin(i))/(inbompas(i)-1))
		 else
	  ompas(i)=0
	  endif
	  enddo
return
end 


subroutine pas_log(nbzone,kpas,kmax,kmin,inbkpas)
	  integer nbzone,i
      integer inbkpas(nbzone)
      double precision kmin(nbzone),kmax(nbzone),kpas(nbzone)
	  do i=1,nbzone
	  if(inbkpas(i).gt.1) 	then 
		kpas(i)=(-log10(kmin(i))+log10(kmax(i)))/(inbkpas(i)-1)
		 else
	  kpas(i)=0
	  endif
	  enddo
return
end 
