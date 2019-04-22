      program pgrm_03_03
!
!     This program computes components of the Hartree-Fock energy and the number
!     of electrons using contraction of various matrices loaded from
!     user-provided files (that presumably come from Gaussian calculations).
!
!     At run time, the program expects 5 command line arguments:
!       (1) the number of electrons;
!       (2) the number of atomic orbital basis functions;
!       (3) an input file containing core-Hamiltonian matrix elements (in
!           symmetric upper/column storage form);
!       (4) an input file containing the overlap matrix elements (in
!           symmetric upper/column storage form); and
!       (5) an input file containing Fock matrix elements (in symmetric
!           upper/column storage form).
!
!     At run time, the program outputs 6 items:
!       (1) the MO energies;
!       (2) the MO coefficients;
!       (3) the density matrix;
!       (4) the one-electron contribution to the total electronic energy;
!       (5) the two-electron contribution to the total electronic energy; and
!       (6) the tr(PS), which should equal the number of electrons.
!
!
!
      implicit none
      integer,parameter::unitIn=10
      integer::i,j,k,iError,nElectrons,nOcc,nBasis,lenSym
      real::oneElectronEnergy,twoElectronEnergy,tracePS
      real,dimension(:),allocatable::symFock,symCoreHamiltonian,symOverlap,  &
        moEnergies,tempSymMatrix
      real,dimension(:,:),allocatable::sqFock,fockTilde,sqCoreHamiltonian,  &
        InvSqrtOverlap,moCoefficients,densityMatrix,tempSqMatrix,  &
        SSPEV_Scratch,sqOverlap
      character(len=256)::cmdlineArg
!
!
!     Begin by reading the number of basis functions, allocating array memory,
!     and loading the symmetric Fock and overlap matrices from input files
!     provided on the command line.
!
      call Get_Command_Argument(1,cmdlineArg)
      read(cmdlineArg,'(I)') nElectrons
      nOcc = nElectrons/2
!
      call Get_Command_Argument(2,cmdlineArg)
      read(cmdlineArg,'(I)') nBasis
      lenSym = (nBasis*(nBasis+1))/2
      allocate(symFock(lenSym),symCoreHamiltonian(lenSym),  &
        symOverlap(lenSym),tempSymMatrix(lenSym))
      allocate(moEnergies(nBasis))
      allocate(sqFock(nBasis,nBasis),fockTilde(nBasis,nBasis),  &
        sqCoreHamiltonian(nBasis,nBasis),invSqrtOverlap(nBasis,nBasis),  &
        moCoefficients(nBasis,nBasis),densityMatrix(nBasis,nBasis),  &
        tempSqMatrix(nBasis,nBasis),SSPEV_Scratch(nBasis,3),sqOverlap(nBasis,nBasis))
!
      call Get_Command_Argument(3,cmdlineArg)
      open(Unit=unitIn,File=TRIM(cmdlineArg),Status='OLD',IOStat=iError)
      if(iError.ne.0) then
        write(*,*)' Error opening input file.'
        STOP
      endIf
      do i = 1,lenSym
        read(unitIn,*) symCoreHamiltonian(i)
      endDo
      close(Unit=unitIn)
!
      call Get_Command_Argument(4,cmdlineArg)
      open(Unit=unitIn,File=TRIM(cmdlineArg),Status='OLD',IOStat=iError)
      if(iError.ne.0) then
        write(*,*)' Error opening input file.'
        STOP
      endIf
      do i = 1,lenSym
        read(unitIn,*) symOverlap(i)
      endDo
      close(Unit=unitIn)
!
      call Get_Command_Argument(5,cmdlineArg)
      open(Unit=unitIn,File=TRIM(cmdlineArg),Status='OLD',IOStat=iError)
      if(iError.ne.0) then
        write(*,*)' Error opening input file.'
        STOP
      endIf
      do i = 1,lenSym
        read(unitIn,*) symFock(i)
      endDo
      close(Unit=unitIn)
!


!
!     Form the square-root of the overlap matrix.
!
      tempSymMatrix=symOverlap
      call InvSQRT_SymMatrix(nBasis,symOverlap,invSqrtOverlap)
      symOverlap=tempSymMatrix
!
!     Form fTilde and solve for the MO energies and coefficients.
!
      call SymmetricPacked2Matrix_UpperPacked(nBasis,symFock,sqFock)
      tempSqMatrix = MatMul(invSqrtOverlap,sqFock)
      fockTilde = MatMul(tempSqMatrix,invSqrtOverlap)
      call Sq2SymMatrix(nBasis,fockTilde,tempSymMatrix)
      call SSPEV('V','U',nBasis,tempSymMatrix,moEnergies,  &
        moCoefficients,nBasis,SSPEV_Scratch,iError)
      If(iError.ne.0) then
        write(*,*)' Failure in SSPEV.'
        STOP
      endIf
      write(*,*)' MO Energies:'
      call Print_Matrix_Full_Real(Reshape(moEnergies,[nBasis,1]),nBasis,1)
      tempSqMatrix = moCoefficients
      moCoefficients = MatMul(invSqrtOverlap,tempSqMatrix)
      write(*,*)' MO Coefficients:'
      call Print_Matrix_Full_Real(moCoefficients,nBasis,nBasis)
!
!
      Do i = 1, nBasis
         Do j = 1, nBasis
           densityMatrix(i,j) = 0
           Do k = 1, nOcc
             densityMatrix(i,j) =densityMatrix(i,j)+2*moCoefficients(i,k)*moCoefficients(j,k)
           endDo
         endDo
      endDo
      write(*,*)' density matrix:'
      call Print_Matrix_Full_Real(densityMatrix,nBasis,nBasis)
!
!
      call SymmetricPacked2Matrix_UpperPacked(nBasis,symCoreHamiltonian,sqCoreHamiltonian)
      tempSqMatrix = MatMul(densityMatrix,sqCoreHamiltonian)
      oneElectronEnergy = 0
      Do i = 1 , nBasis
         oneElectronEnergy = oneElectronEnergy + tempSqMatrix(i,i)
      endDo
      write(*,*)' One-electron energy contribution =  ', oneElectronEnergy
!      
      tempSqMatrix = MatMul(densityMatrix,(sqFock-sqCoreHamiltonian))
      oneElectronEnergy = 0
      Do i = 1 , nBasis
         twoElectronEnergy = twoElectronEnergy + 0.5*tempSqMatrix(i,i)
      endDo
      write(*,*)' Two-electron energy contribution =  ', twoElectronEnergy
!
      call SymmetricPacked2Matrix_UpperPacked(nBasis,symOverlap,sqOverlap)
      tempSqMatrix = MatMul(densityMatrix,sqOverlap)
      tracePS = 0
      Do i = 1 , nBasis
         tracePS = tracePS + tempSqMatrix(i,i)
      endDo
      write(*,*)' tr(PS) = <PS> =  ', tracePS
      end program pgrm_03_03
!
! 
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
      Subroutine Print_Matrix_Full_Real(AMat,M,N)
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
!
!
!



