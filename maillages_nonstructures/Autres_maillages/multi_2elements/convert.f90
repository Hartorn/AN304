 program conv
	implicit none
	
	integer	:: i,j,inode,iface,nn,nf,n
	integer,parameter	:: nnode = 10024, nelem =19231, nface = 1638, neb = 819
	real*8	:: xy(3,nnode)
	integer	:: tri(3,nelem)
	integer :: face(2,nface)
	integer :: seg(2,neb)
	integer :: idx,type,ne,ne_b,dummy, bound
	integer :: tecfile, metisfile_dual

	open(11,file='nodes',status='old')
	read(11,*) type,nn,nf
	do i=1, nnode
	   read(11,*) inode
	   read(11,*) xy(1,i),xy(2,i)
	   xy(3,i) = 0
	enddo
	do i= 1,nface
	   read(11,*) iface,face(1,iface),face(2,iface)
	   xy(3,iface)=face(2,iface)
	enddo
	close(11)


	open(unit=12,file='grid',status='old')
	read(12,*) ne,ne_b
	do i=1, ne
	   read(12,*) idx,type
	   read(12,*) tri(1,idx),tri(2,idx),tri(3,idx)
	   read(12,*) dummy,dummy,dummy
	enddo
        do i=1, neb
	   read(12,*) idx,type,bound
	   read(12,*) seg(1,i), seg(2,i)
	   read(12,*) dummy,dummy
	enddo
	close(12)	

   ! write mesh for postprocessing
   open(13,file='mesh',status='unknown')
	write(13,*) nnode,nelem
	do i =1,nnode
           write(13,'( 2(ES16.8) )') xy(1,i),xy(2,i)
	enddo
	do i = 1, nelem
           write(13,*) tri(1,i), tri(2,i), tri(3,i)
	enddo
   close(13)
   ! write the mesh in format tecplot: mesh.plt to see it with visit

   tecfile = 10
   open(tecfile, file = 'file.plt', status='unknown', form = 'formatted')
   write(tecfile,*) 'TITLE      =  "2D FE Unstructured grid data" '
   write(tecfile,*) 'VARIABLES = "X", "Y"'
   write(tecfile,*) 'ZONE T="P1", F=FEPOINT, N=',nnode,' E=',nelem,' ET=TRIANGLE'


   !write the nodes coordinates
   do n= 1,nnode
      write(tecfile,'( 2(ES16.8) )') xy(1,n),xy(2,n)
   enddo

   do n=1, nelem
      write(tecfile,*) tri(1,n), tri(2,n), tri(3,n)
   enddo

   close(tecfile)
   metisfile_dual=11
   open(metisfile_dual, file = 'dualformetis.dat', status='unknown', form = 'formatted')
   write(metisfile_dual,*) nelem, 1 ! number of elements,type of elements:  1:triangle,2:tetrahedra,3:hexahedra,4:quadrilateral
   do i=1,nelem
      write(metisfile_dual,*) tri(1:3,i)
   enddo
   close(metisfile_dual)

   stop
 end program conv
