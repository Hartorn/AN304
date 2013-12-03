   program postprocess
	implicit none
	integer	:: i, n1, n2, n3, d
	integer	:: node,nelem
	integer :: meshunit,partitionunit,tecunit,verif
	integer, dimension(:,:), allocatable	:: connectivity
	integer, dimension(:), allocatable	:: domain 
	real*8, dimension(:,:), allocatable	:: xyz

	meshunit=10
	partitionunit=11
	tecunit=12

	open(meshunit,file='mesh',status='old')
	!open(partitionunit,file='dual.dpart',status='old') !for metis
	open(partitionunit,file='dual.map',status='old') !for scotch 
	open(tecunit,file='decomp.plt',status='unknown')
	
	!Read number of nodes and number of elements
	read(meshunit,*) node,nelem
	!Allocate the memory
	allocate(xyz(2,node));allocate(connectivity(4,nelem)) !one more for the connectivity to store the domain from the partition
	allocate(domain(node))

	do i=1,node
	   read(meshunit,*) xyz(1,i), xyz(2,i)
	enddo
	read(partitionunit,*) verif !for scotch only
	do i=1,nelem
	   read(meshunit,*) connectivity(1:3,i) !read the connectivity table
	   !for metis
	   !read(partitionunit,*) connectivity(4,i) !read the partition number of the element
	   !for scotch
	   read(partitionunit,*) verif,connectivity(4,i) !read the partition number of the element
	enddo

	!Fast way to visualize the partition, not good for Boundaries between 2 partitions
	do i = 1, nelem
	   n1 = connectivity(1,i)
	   n2 = connectivity(2,i)
	   n3 = connectivity(3,i)
	   d  = connectivity(4,i)
	   domain(n1) = d !associate a domain to the node, for boundaries node the las one win
	   domain(n2) = d
	   domain(n3) = d 
	enddo

	!store the solution in a tecplot format
	

	write(tecunit,*) 'TITLE      =  "2D FE Unstructured grid data" '
	write(tecunit,*) 'VARIABLES = "X", "Y", "Dom"'
	write(tecunit,*) 'ZONE T="P1", F=FEPOINT, N=',node,' E=',nelem,' ET=TRIANGLE'
   
   
	!write the nodes coordinates
	do i= 1,node
	   write(tecunit,'( 2(ES16.8),3x,I6 )') xyz(1,i),xyz(2,i),domain(i)
	enddo
   
	do i=1, nelem
	   write(tecunit,*) connectivity(1:3,i) 
	enddo
	close(partitionunit)
	close(meshunit)
	close(tecunit)

	deallocate(xyz,connectivity,domain)
	stop
   end program postprocess
