!-------------------------------------------------------------------------!
! This code generates the values of absolute square of spherical harmonics!
! Y_lm and the square of real part of Y_lm.                               !
!-------------------------------------------------------------------------!

program SphericalHarmonics
   implicit none 
   
   real(8), parameter :: pi = 4.d0 * datan(1.0d0)
   real(8) :: c1, c2, c3, partC, d1, d2, d3
   integer(8) :: i, j, k, l, m, indicator, nr
   real(8), allocatable :: partE(:), partD(:), Y_ml(:), Y_ml2(:)
   real(8), allocatable :: Y_mlabs2(:), theta(:), phii(:)
   real(8), external :: factorial, dfactorial
   real(8) :: x, Pll, Pml, Pll1, f1, f2, Pllk 
   ! Pml is associated Legendre polynomial (ALP) corresponding to quantum 
   ! numbers m and l.
   ! Pll1 is ALP corresponding to m = l -1 and l = l
   !---------------------------------------------------------------------!
   10 print*, "Input the values of m and l"
   read(*,*) m, l

   if ( abs(m) .gt. l ) then
      print*, "*****************************************************"
      print*, "satisfy the relation: |m| < l "
      goto 10
   else 
      m = m
   end if

   if ( m .lt. 0 ) then
      m = m * (-1)
      indicator = 1           
      !----------------------------------------------------------------!
      ! It indicates that input m is negative and it has been converted!
      ! to positive for easier calculation.                            !
      !----------------------------------------------------------------!
   else
      m = m
      indicator = 0
   end if
   !----------------------------------------------------------------------!
   !            Re [ Y_lm(theta, phii) ] = partC * partD * partE          !   
   !    where                                                             !
   !             partD = p_lm(cos(theta))                                 !
   !             partE = cos(m*phii)                                      !
   !                                   _______________                    !
   !                                  /                                   !
   !                                 /  (2l+1)*(l-m)!                     !
   !                                /  _______________                    !
   !             partC     =       /                                      !
   !                           \  /     4*pi*(l+m)!                       !
   !                            \/                                        !
   !----------------------------------------------------------------------!

   !----------------------------------------------------------------------!
   !                  Calculation of partC                                !
   !----------------------------------------------------------------------!
     
   C1 = ( 2.0d0 * l + 1.0d0 ) / ( 4.d0 * pi ) 
   C2 = factorial(l-m)
   C3 = factorial(l+m)

   partC = sqrt(C1*C2/C3)
   
   nr = 16000      ! Number of random points to be generated for plotting

   allocate(partE(1:nr))
   allocate(partD(1:nr))
   allocate(Y_ml(1:nr))
   allocate(Y_ml2(1:nr))
   allocate(Y_mlabs2(1:nr))
   allocate(theta(1:nr))
   allocate(phii(1:nr))

   call srand(1)
   do i = 1, nr
      theta(i) = pi * dble(rand())
      phii(i)  = 2.0d0 * pi * dble(rand())
   
      !-------------------------------------------------------------------!
      !                      Calculation of part E                        !
      !-------------------------------------------------------------------!

      partE(i) = dcos( m * phii(i) )
      
      !-------------------------------------------------------------------!
      !                      Calculation of part D                        !
      !-------------------------------------------------------------------!
      
      if (dcos(theta(i)) .eq. 1.0d0) then
         if (m.ne.0) then
            Pml = 0.d0
         else
            Pml = 1.0d0
         end if
      else if (dcos(theta(i)).eq. -1.0d0) then
         if (m.ne.0) then
            Pml = 0.d0
         else
            Pml = (-1.0d0) ** l
         end if
      else
         Pll = (-1.0d0)**l * dfactorial(2*l-1) * (dsin(theta(i)))**l
   
         if (m.eq.l) then
            Pml = Pll 
         else if (m.eq.l-1) then
            Pll1 = (-1.0d0) * (dcos(theta(i)) / dsin(theta(i))) * Pll 
            pml = pll1
         else 
            Pll1 = (-1.0d0) * (dcos(theta(i)) / dsin(theta(i))) * Pll 
            do k = 2, l
               f1 = (-1.0d0) / (k*(2.0d0*l-k+1.0d0))
               f2 = 2.0d0 * (l-k+1.0d0) * (dcos(theta(i)) / dsin(theta(i)))
               Pllk = f1 * (Pll + f2*Pll1)
               Pll  = Pll1
               Pll1 = Pllk
               if (l-k.eq.m) then
                  exit
               end if
            end do ! k loop
            Pml = Pllk
         end if
      end if

      partD(i) = Pml
      Y_ml(i) = partC * partD(i) * partE(i)
      Y_ml2(i) = Y_ml(i) * Y_ml(i)
      Y_mlabs2(i) = Y_ml2(i) / partE(i) ** 2
   end do
   !----------------------------------------------------------------------!
   ! Now we print out the results to files

   open(20, file = "spherical-harmonics.dat")
   write(20,'(f0.10,2x,f0.10,2x,f0.10,2x,f0.10)') (phii(i),theta(i),&
                                             Y_ml2(i),Y_mlabs2(i),i=0,nr)
   !----------------------------------------------------------------------!
   ! For checking
   !----------------!

   print*, "For checking "
   print*, "--------------"
   print*, "theta =", theta(10), ", phi =", phii(10)
   print*, "Pml =", partD(10)
   print*, "cos(m*phi) =", partE(10)
   print*, "Part C =", partC 
   print*, "|Yml|**2 =", Y_mlabs2(10)
   print*, "Re[Yml] =", Y_ml(10)
   !----------------------------------------------------------------------!
   deallocate(partE)
   deallocate(partD)
   deallocate(Y_ml)
   deallocate(Y_ml2)
   deallocate(Y_mlabs2)
   deallocate(theta)
   deallocate(phii)
   !----------------------------------------------------------------------!
   print*, "The data have been written to spherical-harmonics.dat"
end program SphericalHarmonics    
!=========================================================================!
real(8) function factorial(m)
   implicit none
   integer(8) :: m, i

   if (m==0) then
      factorial = 1
   else if (m<0) then
      m = -1 * m
      factorial = 1
      do i = 1, m
         factorial = factorial * i
      end do
      factorial = -1 * factorial
   else
      factorial = 1
      do i = 1, m
         factorial = factorial * i
      end do
   end if
   factorial = dble(factorial)
end function factorial
!=========================================================================!
real(8) function dfactorial(x)
   implicit none
   integer(8) :: x, i
   ! we will only need double factorial of odd numbers and no negative 
   ! integers
   dfactorial = 1
   do i = 1, x, 2
      dfactorial = dfactorial * dble(i)
   end do
end function dfactorial
