program mycrm
   use input
   use colradfort
   use onion_module
   implicit none

   open (6, file='crm.out')

   call getinput

   call getadf04

   call alloc

   print *, mode

   select case (trim(mode))
      case ('astro')

         call colrad
      case ('levelscan')
         call levelscan
      case('tempdensscan')
         call tempDensScan
      case ('masscontour')
         call masscontour
      case ('lineplot')
         call lineplot
      case ('onion')
         call onion
      case default
         print *, ' Bad calculation mode requested. Check input. '
   end select

   close (6)
   call dealloc

end program