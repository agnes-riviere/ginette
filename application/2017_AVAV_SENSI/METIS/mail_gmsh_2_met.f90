implicit none
character(200) ligne
character(5) geom
integer(4) nmu,nbloc,ndim,i,typ,ns(8),ntags,tags(10),num,nbloc0,ndimrel,n_rattache,ncor,&
           i1,i2,i3,j,k,kprog,kper,nbloc_new,nmu_new
integer(4) typ_gmsh2met(6,3) ! la dimension 1 ne sert a rien
integer(4) elem_type_deb(2),elem_type_fin(2),n1,n2,n3,n4
data elem_type_deb /2,4/   ! numero du premier type d'element reconnu, a 2 et 3 dim
data elem_type_fin /3,6/   ! numero du dernier type
real(4) coord(3),pv
logical b_sup

real(4),allocatable :: x(:),y(:),z(:)
integer(4),allocatable :: type(:),nsom(:,:),region(:)
logical,allocatable :: b_rattache(:)

data typ_gmsh2met/6*0,3,2,1,0,0,0,0,0,0,7,5,6/
!
! tableau contenant le nombre de sommets de chaque type d'element
integer nstel
dimension nstel(20)
data nstel/4,3,2,2,8,6,4,4,3,2,    10*0/

if (iargc() .ne. 4) then
   write (6,*) ' Usage: mail_gmsh_2_met maillage_gmsh maillage_metis fichier_regions 2D/3D/2D-3D'
   stop
endif

call getarg (1,ligne)
open (unit=1,file=ligne,status='old')


call getarg (2,ligne)
open (unit=2,file=ligne)

call getarg (3,ligne)
open (unit=3,file=ligne)

call getarg (4,geom)
!
! nb de coordonnees par noeud
if (geom(1:2) .eq. '2D') then
   ndim = 2
   ndimrel = 2
   if (geom(3:5) .eq. '-3D') ndimrel = 3    ! nombre effectif de coordonnees par noeud
else
   ndim = 3
   ndimrel = 3
endif

! recherche du nombre de noeuds

do
   read (1,'(a)',end=1000) ligne
   if (ligne(1:6) .eq. '$Nodes') exit
enddo

read (1,*) nmu
write (6,*) nmu,' noeuds'

! recherche du nombre d'elements

do
   read (1,'(a)',end=1000) ligne
   if (ligne(1:9) .eq. '$Elements') exit
enddo

read (1,*) nbloc0
nbloc = 0
do i=1,nbloc0
   read (1,*) num,typ
   if (typ .ge. elem_type_deb(ndim-1) .and. typ .le. elem_type_fin(ndim-1)) nbloc = nbloc + 1
enddo
write (6,*) nbloc,' elements'

!
! on reprend

rewind 1

allocate (x(nmu))
allocate (y(nmu))
if (ndim .eq. 3) allocate (z(nmu))

do
   read (1,'(a)',end=1000) ligne
   if (ligne(1:6) .eq. '$Nodes') exit
enddo

read (1,*) nmu


!
! lecture coord

do i=1,nmu
   read (1,*) num,coord(1:ndimrel)
   x(i) = coord(1)
   if (ndim .eq. 2) then
      if (ndimrel .eq. 2) then
        y(i) = coord(2)
      else
        y(i) = coord(3) ! Pour EPTB, la coordonnee y ne sert a rien (option 2D-3D)
      endif
   else
      y(i) = coord(2)
      z(i) = coord(3)
   endif
enddo

! lecture des elements

do
   read (1,'(a)',end=1000) ligne
   if (ligne(1:9) .eq. '$Elements') exit
enddo

read (1,*) nbloc0

allocate (type(nbloc0))
if (ndim .eq. 2) then
   allocate (nsom(3,nbloc0))
else
    allocate (nsom(4,nbloc0))
endif

allocate (region(nbloc0))

!
! lecture

nbloc = 0
do i=1,nbloc0
   read (1,'(a)') ligne
   read (unit=ligne,fmt=*) num,typ
    !read (1,*) num,typ,ntags,tags(1:ntags),ns(1:nstel(typ_gmsh2met(typ,ndim)))
   if (typ .ge. elem_type_deb(ndim-1) .and. typ .le. elem_type_fin(ndim-1)) then
      read (unit=ligne,fmt=*) num,typ,ntags,tags(1:ntags),ns(1:nstel(typ_gmsh2met(typ,ndim)))
      nbloc = nbloc+1
      type(nbloc) = typ_gmsh2met(typ,ndim)
      region(nbloc) = tags(ntags)
      nsom (1:nstel(typ_gmsh2met(typ,ndim)),nbloc) = ns(1:nstel(typ_gmsh2met(typ,ndim)))
   endif
enddo

!
! verification du maillage
!
! 1) Recherche de noeuds isoles
!

   allocate (b_rattache(nmu))
   b_rattache(:) = .FALSE.
   do i=1,nbloc
      do j=1,nstel(type(i))
         b_rattache (nsom(j,i)) = .TRUE.
      enddo
   enddo

   n_rattache = 0
   do i=1,nmu
      if (b_rattache(i)) then
         n_rattache = n_rattache + 1
      else
         write (6,*) 'noeud',i,' non rattache'
      endif
   enddo
   write (6,*) nmu-n_rattache,' noeuds non rattaches'

!
! suppression des elements relies a des noeuds non rattaches

  if (nmu .ne. n_rattache) then
    write (6,*) ' Suppression des elements non rattaches ...'
    nbloc_new = nbloc
    kprog = 0  ! pour afficher la progression
    kper = (nmu-n_rattache) / 10
    if (kper .eq. 0) kper = 1
    do i=nbloc,1,-1
       b_sup = .FALSE.
       somtri:do j=1,nstel(type(i))
          if (.not. b_rattache(nsom(j,i))) then
              b_sup = .TRUE.
              exit somtri
          endif
       enddo somtri
       if (b_sup) then
          kprog = kprog + 1
          if (mod(kprog,kper) .eq. 0) write (6,*) kprog
          do j=i,nbloc_new-1
              type(j) = type(j+1)
              nsom(1:nstel(type(j)),j)=nsom(1:nstel(type(j+1)),j+1)
              region(j) = region(j+1)
          enddo
          nbloc_new = nbloc_new-1
       endif
    enddo

    write (6,*) ' ...termine',kprog,' faces supprimees, reste',nbloc_new,' elements'
    nbloc = nbloc_new

!
! elimination des noeuds non rattaches
    nmu_new = nmu
    kprog = 0  ! pour afficher la progression
    kper = (nmu-n_rattache) / 10
    if (kper .eq. 0) kper = 1
    write (6,*) ' Suppression des noeuds non rattaches ...'
    do i=nmu,1,-1
       if (.not. b_rattache(i)) then
          kprog = kprog + 1
          if (mod(kprog,kper) .eq. 0) write (6,*) kprog,' / ',nmu-n_rattache
          !write (6,*) 'Suppression du noeud',i
          !
          ! decalage des coordonnees
          do j=i,nmu_new-1
             x(j) = x(j+1)
             y(j) = y(j+1)
             if (ndim .eq. 3) z(j) = z(j+1)
          enddo
          !
          ! decalage des numeros dans nsom
          do j=1,nbloc
              do k=1,nstel(type(j))
                  if (nsom(k,j) .gt. i) nsom(k,j) = nsom(k,j)-1
              enddo
          enddo
          nmu_new = nmu_new - 1
      endif
    enddo
    write (6,*) '... termine,',nmu_new,' noeuds'
    nmu = nmu_new
  endif
!
! 2) orientation des elements

if (ndim .eq. 2) then
   ncor = 0
   do i=1,nbloc
      i1=nsom(1,i)
      i2=nsom(2,i)
      i3=nsom(3,i)
      pv = (x(i2)-x(i1))*(y(i3)-y(i1))-(y(i2)-y(i1))*(x(i3)-x(i1))
      if (pv .lt. 0) then
         ncor = ncor + 1
         nsom(2,i) = i3
         nsom(3,i) = i2
      endif
   enddo
   write (6,*) ncor,' elements reorientes'

else

  do i=1,nbloc
     if (type(i) .ne. 7) then
        write (6,*) ' Element',i,' non tetraedrique'
        stop
     endif
     !
     ! verification
     n1 = nsom(1,i)
     n2 = nsom(2,i)
     n3 = nsom(3,i)
     n4 = nsom(4,i)
     pv = ((y(n2)-y(n1))*(z(n3)-z(n1))-(z(n2)-z(n1))*(y(n3)-y(n1)))*(x(n4)-x(n1)) &
        - ((x(n2)-x(n1))*(z(n3)-z(n1))-(z(n2)-z(n1))*(x(n3)-x(n1)))*(y(n4)-y(n1)) &
        + ((x(n2)-x(n1))*(y(n3)-y(n1))-(y(n2)-y(n1))*(x(n3)-x(n1)))*(z(n4)-z(n1))
     !write (6,*) i,pv/6.
     if (pv .le. 0) then
        write (6,*) ' Element',i,' volume =',pv/6.
        nsom(1,i) = n1
        nsom(2,i) = n3
        nsom(3,i) = n2
        nsom(4,i) = n4
     endif
  enddo

endif
!
! ecriture des fichiers

write (2,*) nmu,nbloc
do i=1,nmu
    if (ndim .eq. 2) then
        write (2,*) x(i),y(i)
    else
        write (2,*) x(i),y(i),z(i)
    endif
enddo
do i=1,nbloc
   write (2,'(i2,6i8)')type(i),nsom(1:nstel(type(i)),i)
enddo

do i=1,nbloc
   write (3,'(i8)')region(i)
enddo

stop ' Termine'


1000 continue
     write (6,*) ' Erreur: fin de fichier'
     stop

end

