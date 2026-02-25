        implicit double precision(a-h,o-x,z),integer(I-N)
character(100) fichier

!if (iargc() .ne. 2) then
!   write (6,*) ' usage: cree_fic fichier fichier_fic'
!   stop
!endif

!call getarg (1,fichier)
open (unit=1,file="Ch_R.txt")
open (unit=2,file="Ch_R.fic")
!call getarg (2,fichier)
open (unit=3,file="Ch_RD.txt")
open (unit=4,file="Ch_RD.fic")
open (unit=16,file="Ch_RG.fic")
!write(*,*)'fichier txt?'
!read(*,*)fichier
!open(1,file=fichier)
!write(*,*)'fichier .fic?'
!read(*,*)fichier
!open(2,file=fichier)

nsor = 0
read (1,*) date0,val
write (2,'(a,g12.4)') 'date',0
write (2,*) val
do
   read (1,*,end=1000) date,val
   write (2,'(a,g12.4)') 'date',anint(date)
   write (2,*) val
   nsor = nsor + 1
enddo
1000 continue
write (6,*) nsor,' sorties'

nsor=0
read (3,*) date0,val1
write (4,'(a,g12.4)') 'date',0
write (4,*) val1
write (16,'(a,g12.4)') 'date',0
write (16,*) val1
do
   read (3,*,end=1001) date,val1
   write (4,'(a,g12.4)') 'date',anint(date)
   write (4,*) val1
   write (16,'(a,g12.4)') 'date',anint(date)
   write (16,*) val1
   nsor = nsor + 1
enddo
1001 continue
write (6,*) nsor,' sorties'

end
