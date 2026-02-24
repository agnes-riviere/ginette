       program tested

      implicit integer (i)
      double precision k(8),l(8),r(8)
      double precision kmin(8),kmax(8),kpas(8)
      double precision lmin(8),lmax(8),lpas(8)
      double precision rmin(8),rmax(8),rpas(8)
	  integer inbkpas(8),inbnpas(8),inblpas(8)
      integer inbrpas(8)
      integer inbswrespas(8),inbaspas(8),inbanspas(8)
      Real  swresmin(8),swres(8),swresmax(8),swrespas(8)
      Real  nmin(8),n(8),nmax(8),npas(8)
      Real  ansmin(8),ans(8),ansmax(8),anspas(8)
      Real  asmin(8),as(8),asmax(8),aspas(8)
      open(unit=1,file='E_tested_ranges.dat',action='read')
      open(unit=2,file='tested_values',action='write')
	  do i=1,8
      READ(1,*)kmin(i),kmax(i),inbkpas(i)

	  enddo
	  do i=1,8
      READ(1,*)nmin(i),nmax(i),inbnpas(i)
	  enddo

		icompteur=0

	  do i=1,8
	  if(inbnpas(i).gt.1) 	npas(i)=dble((nmax(i)-nmin(i))/(inbnpas(i)-1))
	  if(inbkpas(i).gt.1) 	kpas(i)=(-log10(kmin(i))+log10(kmax(i)))/(inbkpas(i)-1)
	  if(inbnpas(i).le.1) npas(i)=0
	  if(inbkpas(i).le.1) kpas(i)=0
	  enddo
      do ik=0,inbkpas(1)-1

	  k(1)=dble(10**(log10(kmax(1))-ik))
     	do iin=0,inbnpas(1)-1
		n(1)=nmin(1)+iin*npas(1)
		do iin2=0,inbnpas(2)-1
		n(2)=nmin(2)+iin2*npas(2)
      do ik2=0,inbkpas(2)-1
	  k(2)=dble(10**(log10(kmax(2))-ik2))

      do iin3=0,inbnpas(3)-1
		n(3)=nmin(3)+iin3*npas(3)
      do ik3=0,inbkpas(3)-1
	  k(3)=dble(10**(log10(kmax(3))-ik3))

      do iin4=0,inbnpas(4)-1
		n(4)=nmin(4)+iin4*npas(4)
      do ik4=0,inbkpas(4)-1
	  k(4)=dble(10**(log10(kmax(4))-ik4))

      do iin5=0,inbnpas(5)-1
		n(5)=nmin(5)+iin5*npas(5)
      do ik5=0,inbkpas(5)-1
	  k(5)=dble(10**(log10(kmax(5))-ik5))



      do iin6=0,inbnpas(6)-1
		n(6)=nmin(6)+iin6*npas(6)
      do ik6=0,inbkpas(6)-1
	  k(6)=dble(10**(log10(kmax(6))-ik6))

      do iin7=0,inbnpas(7)-1
		n(7)=nmin(7)+iin7*npas(7)
      do ik7=0,inbkpas(7)-1
	  k(7)=dble(10**(log10(kmax(7))-ik7))

      do iin8=0,inbnpas(8)-1
		n(8)=nmin(8)+iin8*npas(8)
      do ik8=0,inbkpas(8)-1
	  k(8)=dble(10**(log10(kmax(8))-ik8))



		icompteur=icompteur+1
      write(2,20)k(1),n(1),&
     k(2),n(2),&
     k(3),n(3),&
     k(4),n(4),&
     k(5),n(5),&
	 k(6),n(6),k(7),n(7),&
	 k(8),n(8)
      enddo
      enddo

		enddo
      enddo

      enddo
      enddo

      enddo
      enddo

      enddo
      enddo

      enddo
      enddo

      enddo
      enddo

      enddo
      enddo

		print*,icompteur,"nombre de couples"
20    format(35d9.2)

       end
