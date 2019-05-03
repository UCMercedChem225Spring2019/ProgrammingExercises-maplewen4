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


      Subroutine Sq2SymMatrix(N,AMatIn,ArrayOut)
!
      Implicit None
      Integer,Intent(In)::N
      Real,Dimension((N*(N+1))/2),Intent(Out)::ArrayOut
      Real,Dimension(N,N),Intent(In)::AMatIn
!
      Integer::i,j,k
!
!
!
      k = 1
        Do i = 1, N
          Do j = 1, i
            ArrayOut(k) = AmatIn(i,j)
            k = k + 1
          endDo
        endDo
!
!
!
      Return
      End Subroutine Sq2SymMatrix

      Subroutine Print_Matrix_Full_Real(AMat,M,N)
!
!     This subroutine prints a real matrix that is fully dimension -
!     i.e.,
!     not stored in packed form. AMat is the matrix, which is
!     dimensioned
!     (M,N).
!
!     The output of this routine is sent to unit number 6 (set by the
!     local
!     parameter integer IOut).
!
!
!     Variable Declarations
!
      implicit none
      integer,intent(in)::M,N
      real,dimension(M,N),intent(in)::AMat
!
!     Local variables
      integer,parameter::IOut=6,NColumns=5
      integer::i,j,IFirst,ILast
!
 1000 Format(1x,A)
 2000 Format(5x,5(7x,I7))
 2010 Format(1x,I7,5F14.6)
!
      Do IFirst = 1,N,NColumns
        ILast = Min(IFirst+NColumns-1,N)
        write(IOut,2000) (i,i=IFirst,ILast)
        Do i = 1,M
          write(IOut,2010) i,(AMat(i,j),j=IFirst,ILast)
        endDo
      endDo
!
      Return
      End Subroutine Print_Matrix_Full_Real





