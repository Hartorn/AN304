  program extract_connec
	implicit none
	integer	:: node, element, dim, unit, info
	integer	:: i,j
	integer, dimension(:,:), allocatable	:: connectivity
	integer, dimension(:), allocatable	:: d 
	real*8, dimension(:,:), allocatable	:: xyz
	character*32	:: nom

	unit = 10
	open(unit, file='grid.mesh',status='old')
	read(unit,*) nom
	read(unit,*) info
	write(*,*) nom,info
	read(unit,*)
	read(unit,*) nom
	read(unit,*) dim
	write(*,*) nom,dim
	read(unit,*)
	read(unit,*) nom
	read(unit,*) node
	write(*,*) nom,node
	allocate(xyz(dim,node),d(node))
	do i=1,node
	   read(unit,*) xyz(1,i),xyz(2,i),d(i)
	enddo
	read(unit,*)
	read(unit,*) nom
	read(unit,*) element
	write(*,*) nom,element
	if( nom .ne. 'Triangles' ) then
		write(*,*) 'Error not Triangles'
		stop
	endif
	allocate(connectivity(4,element))
	do i =1, element
		read(unit,*) (connectivity(j,i),j=1,4)
	enddo
	!Pas besoin de plus pour partitionner et visualiser la parttition
	close(unit)
	unit=11
	open(unit,file='dualformetis.dat',status='unknown')
	write(unit,*) element, 1 !1:triangle
	do i = 1, element
		write(unit,*) (connectivity(j,i),j=1,3)
	enddo
	close(unit)


  deallocate(xyz,d,connectivity)
  stop
  end program extract_connec
