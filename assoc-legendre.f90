!------------------------------------------------------------------------!
! This code calculates associated legendre polynomials(ALP).             !
! Pll represents ALP for m=l and l=l.                                    !
! Plli represents ALP for m=l-i and l=l                                  !
!------------------------------------------------------------------------!

program assoc_legendre
   implicit none
   real(8), parameter :: pi = 4.0d0 * datan(1.0d0)
   real(8) :: x, theta, Pll, Pml, Pll1, f1, f2, Pllk
   integer :: m, l, k
   real(8), external :: dfactorial
   !---------------------------------------------------------------------!
   10 print*, "Input the values of m, l."
   read(*,*) m, l
   
   if (m<0) then
      m = (-1) * m 
   else
      m = m
   end if

   if (m>l) then 
      print*, "Invalid values of m and l."
      goto 10
   else
      m = m
   end if

   15 print*, "Enter the value of x."
   read(*,*) x

   if (abs(x)>1.0d0) then
      print*, "x lies outside [-1,1]."
      goto 15
   else 
     m = m 
   end if
   !---------------------------------------------------------------------!
   theta = dacos(x)

   if (dcos(theta) .eq. 1.0d0) then
      if (m.ne.0) then
         Pml = 0.d0
      else
         Pml = 1.0d0
      end if
   else if (dcos(theta).eq. -1.0d0) then
      if (m.ne.0) then
         Pml = 0.d0
      else
         Pml = (-1.0d0) ** l
      end if
   else
      Pll = (-1.0d0)**l * dfactorial(2*l-1) * (dsin(theta))**l
      print*, "P for l=", l, "m=", l, "is", Pll
   
      if (m.eq.l) then
         Pml = Pll
      else if (m.eq.l-1) then
         Pll1 = (-1.0d0) * (dcos(theta) / dsin(theta)) * Pll
         pml = pll1
      else 
         Pll1 = (-1.0d0) * (dcos(theta) / dsin(theta)) * Pll
         print*, "P for l=", l, "m=", l-1, "is", Pll1
         do k = 2, l
            f1 = (-1.0d0) / (k*(2.0d0*l-k+1.0d0))
            f2 = 2.0d0 * (l-k+1.0d0) * (dcos(theta) / dsin(theta))
            Pllk = f1 * (Pll + f2*Pll1)
            print*, "P for l=", l, "m=", l-k, "is", Pllk
            Pll  = Pll1
            Pll1 = Pllk
            if (l-k.eq.m) then
               exit
            end if
         end do ! k loop
         Pml = Pllk
      end if
   end if
   print*, "Pml =", Pml
end program assoc_legendre 
!=========================================================================

real(8) function dfactorial(x)
   implicit none
   integer :: x, i

   dfactorial = 1 
   do i = 1, x, 2
      dfactorial = dfactorial * dble(i)
   end do
end function dfactorial
