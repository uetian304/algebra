module algebra
 implicit none
 contains
! =======================================
! FUNCTIONS:
!     det3()
!     dot_vv()
!     doubledot_22()
!     doubledot_42()
!     doubledot_44()
!     dyad_mm()
!     dyad_vv()
!     kdelta()
!     mat_subt()
!     norm()
!     sp_24()
!     sp_42()
!     sp_mv()
!     sp_vm()
!     vect_add()
!     vect_subt()
! SUBROUTINES:
!     FindInv()
!     idmat()
!     init_sqmat()
!     init_tensor4()
!     init_vect()
!     inverse3()
!     matrixderivative()
! =======================================

! Kronecker Delta  
 function kdelta(i,j)
 implicit none
 integer, intent(in)::i,j
 integer ::kdelta
 if (i==j) then
  kdelta=1
 else
  kdelta=0
 end if
 end function
 
 ! Derivative of matrix 'a' with respect to matrix 'b'
 ! Output in matrix 'c'
 !
 subroutine matrixderivative(a,b,c)
 implicit none
 integer ::i,j,k,l,m,n
 real*8, dimension(3,3),intent(in)::a,b
 real*8, dimension(3,3,3,3), intent(out)::c
 integer:: sa1,sa2,sb1,sb2
 sa1 = size(a,1) 
 sa2 = size(a,2)
 sb1 = size(b,1)
 sb2 = size(b,2)
 if (sa1/=3 .or. sa2 /= 3 .or. sb1 /= 3 .or. sb2 /= 3) then
        print *, "Subroutine matrixderivatives: Orders of the input matrices do not match."
        c=0.00D0
        return
 else
       do l=1,3
        do k=1,3
         do j=1,3
          do i=1,3
           c(i,j,k,l)=kdelta(j,k)*kdelta(i,l)
          end do
         end do
        end do
       end do
 end if
 end subroutine
  
 !Subroutine to find the inverse of a square matrix
 !Author : Louisda16th a.k.a Ashwith J. Rego
 !Reference : Algorithm has been well explained in:
 !http://math.uww.edu/~mcfarlat/inverse.htm
 !http://www.tutor.ms.unimelb.edu.au/matrix/matrix_inverse.html
  SUBROUTINE FindInv(matrix, inverse, n, errorflag)
  IMPLICIT NONE
  !Declarations
  INTEGER, INTENT(IN) :: n
  INTEGER, INTENT(OUT) :: errorflag  !Status. -1= error, 0= normal
  REAL*8, INTENT(IN), DIMENSION(n,n) :: matrix  !Input matrix
  REAL*8, INTENT(OUT), DIMENSION(n,n) :: inverse !Inverted matrix

  LOGICAL :: FLAG = .TRUE.
  INTEGER :: i, j, k 
  REAL*8 :: m
  REAL*8, DIMENSION(n,2*n) :: augmatrix !augmented matrix

  !Augment input matrix with an identity matrix
  DO i = 1, n
  DO j = 1, 2*n
  IF (j <= n ) THEN
  augmatrix(i,j) = matrix(i,j)
  ELSE IF ((i+n) == j) THEN
  augmatrix(i,j) = 1.D0
  Else
  augmatrix(i,j) = 0.D0
  ENDIF
  END DO
  END DO

  !Reduce augmented matrix to upper traingular form
  DO k =1, n-1
  IF (abs(augmatrix(k,k)) .LE. 1.0D-15) THEN
  FLAG = .FALSE.
  DO i = k+1, n
    IF (abs(augmatrix(i,k)) .GE. 1.0D-15) THEN
    DO j = 1,2*n
    augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)
    END DO
    FLAG = .TRUE.
    EXIT
    ENDIF

    IF (FLAG .EQV. .FALSE.) THEN
  PRINT*, "Matrix is non - invertible"
  inverse = 0.D0
  errorflag = -1
  return
    ENDIF
  END DO
  ENDIF
  DO j = k+1, n
  m = augmatrix(j,k)/augmatrix(k,k)
  DO i = k, 2*n
  augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)
  END DO
  END DO
  END DO

  !Test for invertibility
  DO i = 1, n
  IF (abs(augmatrix(i,i)) .LE. 1.0D-15) THEN
  PRINT*, "Matrix is non - invertible"
  inverse = 0.D0
  errorflag = -1
  return
  ENDIF
  END DO

  !Make diagonal elements as 1
  DO i = 1 , n
  m = augmatrix(i,i)
  DO j = i,(2*n)
   augmatrix(i,j) = (augmatrix(i,j)/m)
  END DO
  END DO

  !Reduced right side half of augmented matrix to identity matrix
  DO k = n-1, 1, -1
  DO i =1, k
  m = augmatrix(i,k+1)
  DO j = k, (2*n)
  augmatrix(i,j) = augmatrix(i,j)-augmatrix(k+1,j)*m
  END DO
  END DO
  END DO

  !store answer
  DO i =1, n
  DO j = 1, n
  inverse(i,j) = augmatrix(i,j+n)
  END DO
  END DO
  errorflag = 0
 END SUBROUTINE FINDinv
 !================================== 
 
 ! Scalar Porduct of Matrix and Vector
 function sp_mv(M, V)
 real*8, dimension(:, :), intent(in):: M
 real*8, dimension(:), intent(in):: V
 real*8, dimension(size(M,1)):: sp_mv
 integer :: i, j
 if (size(M,2) /= size(V,1)) then
      print *,"subroutine sp_mv: arrays lengths mismatch"
      sp_mv=0.00D0
 else
       do i=1,size(M,1)
         sp_mv(i)=0.0000D0
             do j=1,size(V,1)
             sp_mv(i)=sp_mv(i)+M(i,j)*V(j)
             end do
       end do
 end if
 end function sp_mv
 !=================================

! Scalar Porduct of Vector and Matrix
 function sp_vm(V, M)
 real*8, dimension(:), intent(in):: V
 real*8, dimension(:, :), intent(in):: M
 real*8, dimension(size(M,2)):: sp_vm
 integer :: i, j
 if(size(v,1).ne. size(m,1)) then
      print *,"subroutine sp_vm: arrays lengths mismatch"
      sp_vm=0.00D0
 else
       do i=1,size(M,2)
        sp_vm(i)=0.000000D0
        do j=1,size(V,1)
         sp_vm(i) = sp_vm(i) +V(j)*M(j,i)
        end do
       end do
  end if
 end function sp_vm 
 !=================================
 
 ! Dot Product of two Vectors (similar to dot_product)
 function dot_vv(A, B)
 real*8, dimension(:), intent(in) :: A, B
 real*8 :: dot_vv
 integer :: i
 dot_vv=0.00D0
 do i=1,size(A)
  dot_vv = dot_vv + A(i)*B(i)
 end do
 end function dot_vv
 !==================================

 ! Dyadic product of two Vectors
 function dyad_vv(v1,v2)
 real*8, intent(in), dimension(:)::v1,v2
 real*8, dimension(size(v1,1),size(v2,1))::dyad_vv
 integer i,j
 do i=1,size(v1,1)
       do j=1,size(v2,1)
            dyad_vv(i,j)=v1(i)*v2(j)
       end do 
 end do
 end function
 !================================== 
 
 ! Dyadic product of two Matrices
 function dyad_mm(a,b)
 implicit none
 real*8, dimension(3,3), intent(in)::a,b
 real*8, dimension(3,3,3,3)::dyad_mm
 integer :: i,j,k,l
 dyad_mm=0.00D0
 do i=1,3
  do j=1,3
   do k=1,3
    do l=1,3
     dyad_mm(l,k,j,i)=a(l,k)*b(j,i)
    end do
   end do
  end do
 end do
 
 end function
 !================================================
 
 ! 
 Subroutine inverse3(matrix, inverse)
 implicit none
 real*8,intent(in),dimension(3,3) :: matrix 
 real*8,intent(out),dimension(3,3):: inverse
 real*8 :: d,d1,d2,d3
 if((size(matrix,1) /= 3) .or. (size(matrix,2) /= 3)) then
  print*, "Subroutine inverse3 error: input matrix must be a 3x3 matrix"
  stop
 else
 d1=matrix(1,1)*(matrix(2,2)*matrix(3,3)-matrix(2,3)*matrix(3,2))
 d2=matrix(1,2)*(matrix(2,1)*matrix(3,3)-matrix(2,3)*matrix(3,1))
 d3=matrix(1,3)*(matrix(2,1)*matrix(3,2)-matrix(2,2)*matrix(3,1))
 d=d1-d2+d3
 inverse(1,1)=(matrix(2,2)*matrix(3,3)-matrix(2,3)*matrix(3,2))/d 
 inverse(2,1)=(matrix(2,1)*matrix(3,3)-matrix(2,3)*matrix(3,1))/d 
 inverse(3,1)=(matrix(2,1)*matrix(3,2)-matrix(2,2)*matrix(3,1))/d 
 inverse(1,2)=(matrix(1,2)*matrix(3,3)-matrix(1,3)*matrix(3,2))/d 
 inverse(2,2)=(matrix(1,1)*matrix(3,3)-matrix(1,3)*matrix(3,1))/d 
 inverse(3,2)=(matrix(1,1)*matrix(3,2)-matrix(1,2)*matrix(3,1))/d 
 inverse(1,3)=(matrix(1,2)*matrix(2,3)-matrix(1,3)*matrix(2,2))/d 
 inverse(2,3)=(matrix(1,1)*matrix(2,3)-matrix(1,3)*matrix(2,1))/d 
 inverse(3,3)=(matrix(1,1)*matrix(2,2)-matrix(1,2)*matrix(2,1))/d 
 endif
 end subroutine

 !================================================
 ! Function det3 calculates the determinant of a 3x3 matrix
 real*8 function det3(matrix)
 implicit none
 real*8,intent(in),dimension(3,3) :: matrix 
 real*8 :: d1,d2,d3
 
 if((size(matrix,1) .ne. 3) .and. (size(matrix,2) .ne. 3)) then
  print*, "error: Det3 - input matrix must be a 3x3 matrix"
  stop
 else 
  d1=matrix(1,1)*(matrix(2,2)*matrix(3,3)-matrix(2,3)*matrix(3,2))
  d2=matrix(1,2)*(matrix(2,1)*matrix(3,3)-matrix(2,3)*matrix(3,1))
  d3=matrix(1,3)*(matrix(2,1)*matrix(3,2)-matrix(2,2)*matrix(3,1))
  Det3=d1-d2+d3
 endif
 end function

 !================================================
 
 ! identity matrix of order 'order'
 subroutine idmat(eye,order)
 implicit none
 integer,intent(in):: order
 real*8,intent(out),dimension(order,order):: eye
 integer :: i,j
 if (order .gt. 1) then
      do j=1,order
            do i=1,order
            eye(i,j)=0.0000D0
            end do
            eye(j,j)=1.0000D0
      end do
 else
      stop "Subroutine idmat error: order must be greater than 1"
 end if
 end subroutine

!================================================
!Subroutine to set all values of 'matrix' of order nxn, equal to zero.
subroutine init_sqmat(matrix,n)
implicit none
integer, intent(in)::n
real*8, intent(out), dimension(n,n) :: matrix
integer i,j
      do i=1,n
      do j=1,n
      matrix(i,j)=0.00d0
      end do 
      end do 
return
end subroutine init_sqmat

!Subroutine to set all values of 'vector' of length n equal to zero.
subroutine init_vect(vector,n)
implicit none
integer, intent(in)::n
real*8,intent(out),dimension(n)::vector
integer i
	do i=1,n
	vector(i)=0.00d0
	end do 
return
end subroutine

!Function to calculate the double contraction of two tensors of rank 2
function doubledot_22(a,b)
implicit none
real*8, dimension(:,:), intent(in)::a,b
real*8 :: doubledot_22
integer :: i,j
doubledot_22 = 0.00d0
do i=1,3
 do j=1,3
  doubledot_22 = doubledot_22 + a(j,i)*b(j,i)
 end do
end do
end function
 
!Function to calculaate double dot product of two tensors, first of rank 4, other of rank 2.
function doubledot_42(a,b)
implicit none
real*8, dimension(3,3,3,3), intent(in)::a
real*8, dimension(3,3), intent(in)::b
real*8, dimension(3,3)::doubledot_42,x
integer :: i,j,k,l
doubledot_42=0.00d0
do j=1,3
 do i=1,3
  x=a(i,j,:,:)
  doubledot_42(i,j) = doubledot_42(i,j) + doubledot_22(x,b)
 end do
end do
end function

!Function to calculate double dot product of two tensors each of rank 4.
function doubledot_44(a,b)
implicit none
real*8, dimension(3,3,3,3), intent(in)::a,b
real*8, dimension(3,3,3,3):: doubledot_44
integer :: i,j,k,l
real*8 :: x
doubledot_44=0.00d0
do i=1,3
 do j=1,3
  do k=1,3
   do l=1,3
    x = doubledot_22(a(l,k,:,:),b(:,:,j,i))
    doubledot_44(l,k,j,i)=doubledot_44(l,k,j,i) + x
   end do
  end do
 end do
end do
end function

!Function to calculate dot product of two tensors, first of rank 2, other of rank 4.
function sp_24(a,b)
real*8, dimension(3,3,3,3), intent(in)::b
real*8, dimension(3,3), intent(in)::a
real*8, dimension(3,3,3,3)::sp_24
integer::k,l
sp_24=0.00d0
do k=1,3
 do l=1,3
  sp_24(:,:,l,k) = sp_24(:,:,l,k) + matmul(a,b(:,:,l,k))
 end do
end do
end function

!Function to calculate dot product of two tensors, first of rank 4, other of rank 2. 
function sp_42(a,b)
real*8, dimension(3,3,3,3), intent(in)::a
real*8, dimension(3,3), intent(in)::b
real*8, dimension(3,3,3,3)::sp_42
integer::k,l
sp_42=0.00d0
do k=1,3
 do l=1,3
  sp_42(l,k,:,:) = sp_42(l,k,:,:) + matmul(a(l,k,:,:),b)
 end do
end do
end function
 
!Subroutine to set all elements of a 4th rank tensor, equal to zero.
subroutine init_tensor4(tensor)
implicit none
real*8, dimension(3,3,3,3),intent(inout)::tensor
tensor=0.00d0
end subroutine
 
! Function returns a matrix resulting from subtraction of mat2 from mat1. 
function mat_subt(mat1,mat2)
implicit none
real*8, dimension(:,:), intent(in)::mat1,mat2
real*8, dimension(size(mat1,1),size(mat1,2)):: mat_subt
if (size(mat1,1)==size(mat2,1) .and. size(mat1,2)==size(mat2,2)) then 
   mat_subt = mat1-mat2
else
   print *,"function call: mat_subt:: matrices are not equal"
   mat_subt=0.00d0
end if
end function
 
!Function returns a vector after vector "b" is subtracted from vector "a" vectors.
function vect_subt(a,b)
implicit none
real*8, dimension(:), intent(in)::a,b
real*8, dimension(min(size(a,1),size(b,1)))::vect_subt
integer :: i
do i=1,min(size(a,1),size(b,1))
  vect_subt(i)=a(i)-b(i)
end do
end function 
 
!Function adds two vectors (of any length).
function vect_add(a,b)
implicit none
real*8, dimension(:), intent(in)::a,b
real*8, dimension(min(size(a,1),size(b,1)))::vect_add
integer :: i
do i=1,min(size(a,1),size(b,1))
 vect_add(i)=a(i)+b(i)
end do
end function

!Function returns the 2nd norm of a vector.
function norm(a)
implicit none 
real*8, dimension(:), intent(in)::a 
real*8:: norm,sn
integer::i
sn=0.00d0 
do i=1,size(a,1)
 sn=sn+a(i)**2
end do 
norm = sqrt(sn)
end function

end module
