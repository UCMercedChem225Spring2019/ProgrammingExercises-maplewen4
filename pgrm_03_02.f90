      program pgrm_03_02
!
!     This program symmetrizes a square matrix and reports the results as a
!     linearized array (vector).
!
!     At execution time, the program expects 2 command line arguments: (1) nDim,
!     the leading dimension of the square matrix; and (2) an input file
!     containing the matrix elements in full storage form.
!
!
      implicit none
      integer,parameter::unitIn=10
      integer::i,j,iError,nDim,lenSym
      real,dimension(:),allocatable::symMatrix
      real,dimension(:,:),allocatable::sqMatrix
      character(len=256)::cmdlineArg
!
!
!     Begin by reading the leading dimension of the matrix and the input file
!     name from the command line. Then, open the file and read the input matrix,
!     sqMatrix
!
      call Get_Command_Argument(1,cmdlineArg)
      read(cmdlineArg,'(I)') nDim
      lenSym = (nDim*(nDim+1))/2
      allocate(sqMatrix(nDim,nDim),symMatrix(lenSym))
!
      call Get_Command_Argument(2,cmdlineArg)
      open(Unit=unitIn,File=TRIM(cmdlineArg),Status='OLD',IOStat=iError)
      if(iError.ne.0) then
        write(*,*)' Error opening input file.'
        STOP
      endIf
      do i = 1,nDim
        do j = 1,nDim
          read(unitIn,*) sqMatrix(i,j)
        endDo
      endDo
      close(Unit=unitIn)
!
!     Print the input matrix.
!
      write(*,*)' Input Matrix:'
      call Print_Matrix_Full_Real(sqMatrix,nDim,nDim)
      call Sq2SymMatrix(nDim,sqMatrix,symMatrix)
      write(*,*)
      write(*,*)' Output Matrix:'
      do i = 1,lenSym
        write(*,*) symMatrix(i)
      endDo
!
      end program pgrm_03_02
