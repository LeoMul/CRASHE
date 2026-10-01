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
   if (mode .eq. 'astro') then

      shellVelocityOuterC = velocityExpansionC
      shellVelocityInnerC = 0.0_f64
      call getAtomicDensityLocal

      call colrad

   else if (mode .eq. 'levelscan') then
      call levelscan
   else if (mode .eq. 'masscontour') then
      call masscontour
   else if (mode .eq. 'lineplot') then
      call lineplot
   else if (mode .eq. 'onion') then
      call onion
   else
      print *, ' Bad calculation mode requested. Check input. '
   end if

   close (6)
   call dealloc

end program