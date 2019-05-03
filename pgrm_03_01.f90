      program pgrm_03_01
!
!     This program computes the inverse square-root of a matrix.
!
!     At execution time, the program expects 2 command line arguments: (1) nDim;
!     and (2) the matrix being raised to the (1/2) power.
!
!
      implicit none
      integer,parameter::unitIn=10
      integer::i,iError,nDim,lenSym
      real,dimension(:),allocatable::inputSymMatrix
      real,dimension(:,:),allocatable::inputSqMatrix,invSqrtInputMatrix
      character(len=256)::cmdlineArg
!
!
!     Begin by reading the leading dimension of the matrix and the input file
!     name from the command line. Then, open the file and read the input matrix,
!     inputSymMatrix
!
      call Get_Command_Argument(1,cmdlineArg)
      read(cmdlineArg,'(I)') nDim
      lenSym = (nDim*(nDim+1))/2
      allocate(inputSymMatrix(lenSym),inputSqMatrix(nDim,nDim),  &
        invSqrtInputMatrix(nDim,nDim))
      call Get_Command_Argument(2,cmdlineArg)
      open(Unit=unitIn,File=TRIM(cmdlineArg),Status='OLD',IOStat=iError)
      if(iError.ne.0) then
        write(*,*)' Error opening input file.'
        STOP
      endIf
      do i = 1,lenSym
        read(unitIn,*) inputSymMatrix(i)
      endDo
      close(Unit=unitIn)
!
!     Form the square-root of inputSymMatrix. The result is loaded into square
!     (full storage) matrix invSqrtInputMatrix.
!
      write(*,*)' The matrix loaded (column) upper-triangle packed:'
      call SymmetricPacked2Matrix_UpperPacked(nDim,inputSymMatrix,  &
        inputSqMatrix)
      write(*,*)' Input Matrix:'
      call Print_Matrix_Full_Real(inputSqMatrix,nDim,nDim)
      call InvSQRT_SymMatrix(nDim,inputSymMatrix,invSqrtInputMatrix)
      write(*,*)' Inverse SQRT Matrix:'
      call Print_Matrix_Full_Real(invSqrtInputMatrix,nDim,nDim)
      write(*,*)' Matrix product that should be the identity.'
      call Print_Matrix_Full_Real(MatMul(MatMul(invSqrtInputMatrix,  &
        invSqrtInputMatrix),inputSqMatrix),nDim,nDim)
!
      end program pgrm_03_01
!
!
!
      Subroutine InvSQRT_SymMatrix(N,ArrayIn,AMatOut)
!
!
      Implicit None
      Integer,Intent(In)::N
      Integer::IError1,j,k
      Real,Dimension((N*(N+1))/2),Intent(In)::ArrayIn
      Real,Dimension(N,N)::EVecs,Temp1_Matrix,Temp2_Matrix,Temp3_Matrix
      Real,Dimension(N)::EVals
      Real,Dimension(3*N)::Temp_Vector
      Real,Dimension(N,N),Intent(Out)::AMatOut
!
!
      call SSPEV('V','U',N,ArrayIn,Evals,Evecs,N, Temp_Vector, IError1)
       If(IError1.ne.0) then
        Write(*,*)' Failure in DSPEV.'
        STOP
      endIf
!
!      Do j=1, N
!        EVals(j) = 1/SQRT(Evals(j))
!      endDo
!
!
      Do j=1, N
        Do k=1, N
          If(j.EQ.k) then
             Temp1_Matrix(j,k)=1/SQRT(EVals(j))
          ELSE If(j.NE.k) then
             Temp1_Matrix(j,k)=0.0
          endIf
        endDo
      endDo
!
      Do j=1, N
        Do k=1, N
          Temp2_Matrix(j,k) = EVecs(k,j)
        endDo
      endDo
!
      AMatout=(MatMul(MatMul(EVecs, Temp1_matrix),Temp2_Matrix))
!
      return
!
      End Subroutine InvSQRT_SymMatrix
!
!
      Subroutine SymmetricPacked2Matrix_UpperPacked(N,ArrayIn,AMatOut)
!
!     This subroutine accepts an array, ArrayIn, that is (N*(N+1))/2
!     long.
!     It then converts that form to the N-by-N matrix AMatOut taking
!     ArrayIn to be in upper-packed storage form. Note: The storage mode
!     also assumes the upper-packed storage is packed by columns.
!
      Implicit None
      Integer,Intent(In)::N
      Real,Dimension((N*(N+1))/2),Intent(In)::ArrayIn
      Real,Dimension(N,N),Intent(Out)::AMatOut
!
      Integer::i,j,k
!
!     Loop through the elements of AMatOut and fill them appropriately
!     from
!     Array_Input.
!
!
! *************************************************************************
! WRITE CODE HERE TO READ THE ARRAY ELEMENTS FROM THE INPUT FILE.
! *************************************************************************
!
        k = 1
        Do i = 1, N
          Do j = 1, i
            AmatOut(i,j) = ArrayIn(k)
            AmatOut(j,i) = ArrayIn(k)
            k = k + 1
          endDo
        endDo
!
!
!
      Return
      End Subroutine SymmetricPacked2Matrix_UpperPacked
!
!
!
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


