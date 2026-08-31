c-------------------------------------------------------------------
      subroutine vf_export_all_walls(wall_quads_file_name)
c
c  export all walls, for view factor calculation.
c
      include 'SIZE'
      include 'TOTAL'
      include "view_factors/VIEW_FACTORS"

      integer e,eb,eg,i,j,k,ii,ix,iy,mid,kb
      integer i0,i1,j0,j1,k0,k1

      parameter (lblock=500)

      common /scrns_vf/ farea_a(nfe,lblock),
     & fnormal_a(3,nfe,lblock),
     & fvertice1_a(3,nfe,lblock),
     & fvertice2_a(3,nfe,lblock),
     & fvertice3_a(3,nfe,lblock), 
     & fvertice4_a(3,nfe,lblock),
     & wka(nfe*lblock),
     & wka2(3*nfe*lblock),
     & ibc(6,lblock),wk(6*lblock)

      character*1 s4(4)
      character*3 s3
      integer     i4
      equivalence(i4,s4)
      equivalence(s3,s4)

      character*80 wall_quads_file_name
      integer nwalls,fileid
      real sint2,sarea2,sn1(3)
	  
      integer nface,nlg,nemax
 
      nface = 2*ldim
      nlg = nelg(1)
 
      ! call blank(wall_quads_file_name,80)      

      if (nid.eq.0)  then 
      write(6,*) "VF: vf_export_all_walls "
      endif

      ! calculate total number of wall faces, and wall face areas, and wall face vertices
	  ! and wall face nornmal vector
	  ! 

      nwalls = 0

      do iel=1,nelt
      do ifc=1,2*ndim
		 
		 !if (cbc(ifc,iel,1).eq.'W  ') then
		 
		 ! export for all boundaries, not just wall
		 ! this is needed to enforce cloture, which is needed to normalized view factors
		 ! 
         if ((cbc(ifc,iel,1).ne.'E  ')
     & .and.(cbc(ifc,iel,1).ne.'   '))  then

           nwalls = nwalls + 1

		   call facind(i0,i1,j0,j1,k0,k1,nx1,ny1,nz1,ifc)

           do k=k0,k1
           do j=j0,j1
           do i=i0,i1
	
           call getSnormal(sn1,i,j,k,ifc,iel)	
	
           fnrmx(i,j,k,iel) =  -sn1(1)
		   fnrmy(i,j,k,iel) =  -sn1(2)
		   fnrmz(i,j,k,iel) =  -sn1(3)

           enddo
           enddo
           enddo

           if (i0.eq.i1) then
		   
             fvertice1(1,ifc,iel) = xm1(i0,j0,k0,iel)
             fvertice1(2,ifc,iel) = ym1(i0,j0,k0,iel)
             fvertice1(3,ifc,iel) = zm1(i0,j0,k0,iel)
			 
             fvertice2(1,ifc,iel) = xm1(i0,j0,k1,iel)
             fvertice2(2,ifc,iel) = ym1(i0,j0,k1,iel)
             fvertice2(3,ifc,iel) = zm1(i0,j0,k1,iel)

             fvertice3(1,ifc,iel) = xm1(i0,j1,k1,iel)
             fvertice3(2,ifc,iel) = ym1(i0,j1,k1,iel)
             fvertice3(3,ifc,iel) = zm1(i0,j1,k1,iel)
				 
             fvertice4(1,ifc,iel) = xm1(i0,j1,k0,iel)
             fvertice4(2,ifc,iel) = ym1(i0,j1,k0,iel)
             fvertice4(3,ifc,iel) = zm1(i0,j1,k0,iel)
				 
		   endif

           if (j0.eq.j1) then
		   
             fvertice1(1,ifc,iel) = xm1(i0,j0,k0,iel)
             fvertice1(2,ifc,iel) = ym1(i0,j0,k0,iel)
             fvertice1(3,ifc,iel) = zm1(i0,j0,k0,iel)
			                            
             fvertice2(1,ifc,iel) = xm1(i0,j0,k1,iel)
             fvertice2(2,ifc,iel) = ym1(i0,j0,k1,iel)
             fvertice2(3,ifc,iel) = zm1(i0,j0,k1,iel)
                                        
             fvertice3(1,ifc,iel) = xm1(i1,j0,k1,iel)
             fvertice3(2,ifc,iel) = ym1(i1,j0,k1,iel)
             fvertice3(3,ifc,iel) = zm1(i1,j0,k1,iel)
				                        
             fvertice4(1,ifc,iel) = xm1(i1,j0,k0,iel)
             fvertice4(2,ifc,iel) = ym1(i1,j0,k0,iel)
             fvertice4(3,ifc,iel) = zm1(i1,j0,k0,iel)
				 
		   endif


           if (k0.eq.k1) then
		   
             fvertice1(1,ifc,iel) = xm1(i0,j0,k0,iel)
             fvertice1(2,ifc,iel) = ym1(i0,j0,k0,iel)
             fvertice1(3,ifc,iel) = zm1(i0,j0,k0,iel)
			                            
             fvertice2(1,ifc,iel) = xm1(i0,j1,k0,iel)
             fvertice2(2,ifc,iel) = ym1(i0,j1,k0,iel)
             fvertice2(3,ifc,iel) = zm1(i0,j1,k0,iel)
                                        
             fvertice3(1,ifc,iel) = xm1(i1,j1,k0,iel)
             fvertice3(2,ifc,iel) = ym1(i1,j1,k0,iel)
             fvertice3(3,ifc,iel) = zm1(i1,j1,k0,iel)
				                        
             fvertice4(1,ifc,iel) = xm1(i1,j0,k0,iel)
             fvertice4(2,ifc,iel) = ym1(i1,j0,k0,iel)
             fvertice4(3,ifc,iel) = zm1(i1,j0,k0,iel)
				 
		   endif

         call surface_int(sint2,sarea2, fnrmx,iel,ifc)
         fnormal(1,ifc,iel) = sint2/sarea2

         call surface_int(sint2,sarea2, fnrmy,iel,ifc)
         fnormal(2,ifc,iel) = sint2/sarea2
		 
         call surface_int(sint2,sarea2, fnrmz,iel,ifc)
         fnormal(3,ifc,iel) = sint2/sarea2		 

         farea(ifc,iel) = sarea2

         endif

      enddo
      enddo

      nwalls = iglsum(nwalls,1)
      if (nid.eq.0)  then 
      write(6,*) "VF: Total number of wall faces: ", nwalls
      endif

      ! export all information above to file 

      if (nid.eq.0) then  ! open file at 1st rank

       !wall_quads_file_name = 'quads_file'

       write(6,*) "dumping quad information to ", wall_quads_file_name

       fileid = 377 + nid

       open(unit=fileid,file=wall_quads_file_name,status='unknown')
	   write(fileid,*) nwalls   ! total number of quads

       nwalls = 0

      endif  ! only open file at 1st rank
	   

c gather information per block into nid=0 processor
c and only write at nid=0 processor

      do eb=1,nlg,lblock
         nemax = min(eb+lblock-1,nlg)

         call izero(ibc, 6*lblock)
         call rzero(farea_a,nfe*lblock)
         call rzero(fnormal_a,3*nfe*lblock)
         call rzero(fvertice1_a,3*nfe*lblock)
         call rzero(fvertice2_a,3*nfe*lblock)
         call rzero(fvertice3_a,3*nfe*lblock)
         call rzero(fvertice4_a,3*nfe*lblock)

         call izero(wk, 6*lblock)
         call rzero(wka,nfe*lblock)
         call rzero(wka2,3*nfe*lblock)

         kb = 0
         do eg=eb,nemax    ! loop global element index for this block
            mid = gllnid(eg)   ! eg -> processors id
            e   = gllel (eg)   !  eg -> e, logal element index
            kb  = kb+1         ! local element account in this eb 
            if (mid.eq.nid) then ! if eg belong to this nid, then store info here.
               do i=1,nface
                  i4 = 0
                  call chcopy(s4,cbc(i,e,1),3)
                  ibc(i,kb) = i4

         if ((s3.ne.'E  ').and.(s3.ne.'   '))  then
                  !if (s3.eq.'W  ') then
                  ! store c_Cr_alloy to cra 

                  farea_a(i,kb) = farea(i,e)
				  
                  fnormal_a(1,i,kb) = fnormal(1,i,e)
                  fnormal_a(2,i,kb) = fnormal(2,i,e)
                  fnormal_a(3,i,kb) = fnormal(3,i,e)

                  fvertice1_a(1,i,kb) = fvertice1(1,i,e)
                  fvertice1_a(2,i,kb) = fvertice1(2,i,e)
                  fvertice1_a(3,i,kb) = fvertice1(3,i,e)

                  fvertice2_a(1,i,kb) = fvertice2(1,i,e)
                  fvertice2_a(2,i,kb) = fvertice2(2,i,e)
                  fvertice2_a(3,i,kb) = fvertice2(3,i,e)

                  fvertice3_a(1,i,kb) = fvertice3(1,i,e)
                  fvertice3_a(2,i,kb) = fvertice3(2,i,e)
                  fvertice3_a(3,i,kb) = fvertice3(3,i,e)

                  fvertice4_a(1,i,kb) = fvertice4(1,i,e)
                  fvertice4_a(2,i,kb) = fvertice4(2,i,e)
                  fvertice4_a(3,i,kb) = fvertice4(3,i,e)


                  endif
				  
               enddo
            endif
         enddo

         call igop(ibc,wk,'+  ', 6*lblock)                ! Sum across all processors
         call gop(farea_a,wka,'+  ',nfe*lblock)  ! Sum across all processors
         call gop(fnormal_a,wka,'+  ',3*nfe*lblock)  ! Sum across all processors
         call gop(fvertice1_a,wka,'+  ',3*nfe*lblock)  ! Sum across all processors
         call gop(fvertice2_a,wka,'+  ',3*nfe*lblock)  ! Sum across all processors
         call gop(fvertice3_a,wka,'+  ',3*nfe*lblock)  ! Sum across all processors
         call gop(fvertice4_a,wka,'+  ',3*nfe*lblock)  ! Sum across all processors


         ! write at nid=0
         if (nid.eq.0) then
            kb = 0
            do eg=eb,nemax
               kb  = kb+1

               do i=1,nface
                  i4 = ibc(i,kb)   ! equivalenced to s4 and s3
                  !write(6,*) s3
	         if ((s3.ne.'E  ').and.(s3.ne.'   '))  then
                  !if (s3.eq.'W  ') then  ! only dump for wall
                    !write(6,*) 'writing face data'
                    !write(510,*) 'writing face data',lx1,ly1

                nwalls = nwalls + 1

                write(fileid,*) nwalls, eg,i
                    
             write(fileid,*) fvertice1_a(1,i,kb),fvertice1_a(2,i,kb),
     & fvertice1_a(3,i,kb),fvertice2_a(1,i,kb),
     & fvertice2_a(2,i,kb),fvertice2_a(3,i,kb)

                write(fileid,*) fvertice3_a(1,i,kb),fvertice3_a(2,i,kb),
     & fvertice3_a(3,i,kb),fvertice4_a(1,i,kb),
     & fvertice4_a(2,i,kb),fvertice4_a(3,i,kb)

                write(fileid,*) fnormal_a(1,i,kb),fnormal_a(2,i,kb),
     & fnormal_a(3,i,kb),farea_a(i,kb)

                  endif    


               enddo
            enddo
         endif

      enddo
	  
      if (nid.eq.0) close(fileid)



      return
      end
!c------------------------------------------------------------------
!c-------------------------------------------------------------------
      subroutine vf_read_view_factors(view_factor_file)
c
c  read view factor 
c
      include 'SIZE'
      include 'TOTAL'
      include "view_factors/VIEW_FACTORS"
	  
      
	  
      character*80 view_factor_file
      character*80 dummyline

      integer nwalls,fileid,ivwalls,iw,ieg,ifc,nvwalls,iwall
      integer jw,jeg,jfc
      real   view_factor

      if (nid.eq.0) then 
       write(6,*) 'VF:Reading view factor file...'
      endif
	 
c reading is parallel
c so open file for all nid

      ierr  = 0

      !view_factor_file =  'view_factors'

      fileid = 511+nid

      if (nid.eq.0) then
      open(unit=fileid,file=view_factor_file,status='old',iostat=ierr)      
	  if (ierr.gt.0) call exitti('Cannot open view factor file !$',1) 
      endif
	
      if (nid.ne.0)open(unit=fileid,file=view_factor_file,status='old')      

      read(fileid,*) nwalls
      if (nid.eq.0) then
      write(6,*) 'VF:Reading view factor file: nwalls=',nwalls
      endif
	  
      if ((nwalls.gt.max_walls).and.(nid.eq.0)) then
      write(6,*) 'please increase max_walls to ', nwalls
      endif 

      call izero(egf_to_iwall,6*lelg)

      do iwall = 1,nwalls
	  
        read(fileid,*) iw,ieg,ifc,nvwalls

        egf_to_iwall(ifc,ieg) = iw
        iwall_to_eg(iw) = ieg
        iwall_to_f(iw) = ifc
		
        if ((iwall.ne.iw).and.(nid.eq.0)) then
        write(6,*) 'please check view factor file, something is wrong'
        write(6,*) 'iwall is ', iwall, ' but iw in vf file is ', iw
        write(6,*) 'ieg is ', ieg
        write(6,*) 'ifc is ', ifc
        write(6,*) 'nvwalls is ', nvwalls
        endif 

        if ((nvwalls.gt.max_visible_walls).and.(nid.eq.0)) then
        write(6,*) 'please increase max_visible_walls to ', nvwalls
        endif 

        if (gllnid(ieg).eq.nid) then        ! if this element belong to this nid 

          iel=gllel(ieg)
          n_visible_walls(ifc,iel) = nvwalls
          do ivwalls = 1,nvwalls
          read(fileid,*) jw,jeg,jfc,view_factor
           view_factor_value(ivwalls,ifc,iel) = view_factor
           jindex_to_jwall(ivwalls,ifc,iel)=jw
          enddo

        else

          do ivwalls = 1,nvwalls
          read(fileid,*) dummyline
          enddo

        endif

      enddo 

      close(fileid)

      if (nid.eq.0) then
      write(6,*) 'Done: reading view factor file'
      endif

      vf_crf_firstCalled =1

      return
      end
!c------------------------------------------------------------------
!c-------------------------------------------------------------------
      subroutine vf_export_partitioned_vf_files(vf_folder)
c
c  read view factor 
c
      include 'SIZE'
      include 'TOTAL'
      include "view_factors/VIEW_FACTORS"

      character*80 vf_folder
      character*256 cmdline
      character*80 vf_part_file_name
      integer fileid


      if (nid.eq.0) then 
       write(6,*) 'VF:exporting partitioned view factor file...'
      endif

      if (nid.eq.0) then
       ! delete folder
       cmdline = 'rm -rf ' // trim(vf_folder)
       call execute_command_line(trim(cmdline))   ! or system(trim(cmdline))

       ! recreate folder
       cmdline = 'mkdir ' // trim(vf_folder)
       call execute_command_line(trim(cmdline))   ! or system(trim(cmdline))
      endif

      call nekgsync()

      write(vf_part_file_name, '(A, A, I0)') 
     & trim(vf_folder), '/vfp', nid

      fileid = 367 + nid
      open(unit=fileid,file=vf_part_file_name,status='new')
      
	  nwall = 0
      do iel = 1,lelt
        eg = lglel(iel)
        do ifc= 1,6
            if ((cbc(ifc,iel,1).ne.'E  ')
     & .and.(cbc(ifc,iel,1).ne.'   '))  then
            nwall = nwall+1
            endif
        enddo
      enddo

      write(fileid,*)nwall

   71 format(i0,' ',i0,' ',i0,' ',i0)	
   72 format(i0,' ',i0,' ',i0,' ',f12.8)	

      do iel = 1,lelt
        eg = lglel(iel)
        do ifc= 1,6

            if ((cbc(ifc,iel,1).ne.'E  ')
     & .and.(cbc(ifc,iel,1).ne.'   '))  then
        
            iwall = egf_to_iwall(ifc,eg)
	        nvw = n_visible_walls(ifc,iel)
            ieg = int(eg)
            write(fileid,71) iwall,ieg,ifc,nvw
   
            do ivw = 1,nvw

            jwall = jindex_to_jwall(ivw,ifc,iel)
            jg = iwall_to_eg(jwall)
            jf = iwall_to_f(jwall)
            fij = view_factor_value(ivw,ifc,iel)
            write(fileid,72) jwall,jg,jf,fij

            enddo

            endif
        enddo
      enddo

      close(fileid)

      return
      end
!c------------------------------------------------------------------
!c-------------------------------------------------------------------
      subroutine vf_read_partitioned_vf_files(vf_folder)
c
c  read view factor 
c
      include 'SIZE'
      include 'TOTAL'
      include "view_factors/VIEW_FACTORS"

      common /scrns_vf4/ 
     & wka5(6*lelg),wka6(max_walls)


      character*80 vf_folder
      character*80 vf_part_file_name
      integer fileid


      if (nid.eq.0) then 
       write(6,*) 'VF:reading partitioned view factor file...'
      endif

      write(vf_part_file_name, '(A, A, I0)') 
     & trim(vf_folder), '/vfp', nid

      fileid = 775 + nid
      open(unit=fileid,file=vf_part_file_name,status='old')

      call izero(egf_to_iwall,6*lelg)
      call izero(iwall_to_eg,max_walls)
      call izero(iwall_to_f,max_walls)

      read(fileid,*) nwalls

      do iwall = 1,nwalls
	  
        read(fileid,*) iw,ieg,ifc,nvwalls

        egf_to_iwall(ifc,ieg) = iw
        iwall_to_eg(iw) = ieg
        iwall_to_f(iw) = ifc

        iel=gllel(ieg)
        n_visible_walls(ifc,iel) = nvwalls
         
		do ivwalls = 1,nvwalls
          read(fileid,*) jw,jeg,jfc,view_factor
           view_factor_value(ivwalls,ifc,iel) = view_factor
           jindex_to_jwall(ivwalls,ifc,iel)=jw
        enddo

      enddo 

      close(fileid)

      call igop(egf_to_iwall,wka5,'+  ',6*lelg)
      call igop(iwall_to_eg,wka6,'+  ',max_walls)
      call igop(iwall_to_f,wka6,'+  ',max_walls)

      vf_crf_firstCalled =1
	  
      return
      end
!c-----------------------------------------------------------------
!c-------------------------------------------------------------------
      subroutine vf_calculate_radiation_heat_flux(fac0,
     & open_face_option)
c
c calculate radiation heat flux using view factor method
c
      include 'SIZE'
      include 'TOTAL'
      include "view_factors/VIEW_FACTORS"

      common /scrns_vf2/ 
     & wka3(max_walls),wka4(max_walls)

      real eps,sigma,Cp,rho,U0,T0,f0,fac0,fac1

      real frad_min,frad_max,Jr_min,Jr_max,Gr_min,Gr_max

      integer i0,i1,j0,j1,k0,k1,iel,eg,ifc,i,j,k
      real frad1,sint2,sarea2
      real view_factor,Fsum
      integer iwall,jwall,jindex,open_face_option

       ntot = lx1*ly1*lz1*lelt
  
       if (nid.eq.0) then
       write(6,*) 'VF: starting vf_calculate_radiation_heat_flux'
       endif
  
c 1s step.
c mapping grid point temperature to  element faces
c ues element face averaged temperature

        call rzero(temp_ef,max_walls)
        call rzero(indicator_real_wall_ef,max_walls)

            do iel = 1,lelt
              eg = lglel(iel)
          
            do ifc= 1,6
          
               ! export is done for for all boundaries, not just wall
               ! this is needed to enforce cloture, which is needed to normalized view factors
            if ((cbc(ifc,iel,1).ne.'E  ')
     & .and.(cbc(ifc,iel,1).ne.'   '))  then
                     
            iwall = egf_to_iwall(ifc,eg)
                     
            call surface_int(sint2,sarea2,t(1,1,1,1,1),iel,ifc)
            temp_ef(iwall)= sint2/sarea2  
            indicator_real_wall_ef(iwall) = 0.0    ! 0.0 for other boundaries, like inlet and outlet
                     
            if(cbc(ifc,iel,1).eq.'W  ') then ! if this boundary is a real wall
            indicator_real_wall_ef(iwall) = 1.0
            endif

            endif
                
             enddo
             enddo

c 2nd step,  broad cast to all mpi ranks
        call gop(temp_ef,wka3,'+  ',max_walls) 
        call gop(indicator_real_wall_ef,wka3,'+  ',max_walls) 

       temp_ef_min = glmin(temp_ef,max_walls)
       temp_ef_max = glmax(temp_ef,max_walls)
       if (nid.eq.0) then
       write(6,*) 'VF: temp_ef_min:',temp_ef_min
       write(6,*) 'VF: temp_ef_max:',temp_ef_max
       endif


c check view factor sum
c maybe normalize here ? 
c
c for on boundary surfaces, its view factor sum should be 1
c for printed out information Fsum should always be 1.
c
          do iel = 1,lelt
              eg = lglel(iel)
              do ifc= 1,6

            if ((cbc(ifc,iel,1).ne.'E  ')
     & .and.(cbc(ifc,iel,1).ne.'   '))  then

             !if (cbc(ifc,iel,1).eq.'W  ') then
                     
              iwall = egf_to_iwall(ifc,eg)
              Fsum = 0.0
	
               do jindex = 1,n_visible_walls(ifc,iel)

                jwall = jindex_to_jwall(jindex,ifc,iel)               
			   view_factor = view_factor_value(jindex,ifc,iel)

               Fsum = Fsum+ view_factor

               enddo                     
           
          !if (nid.eq.0) write(6,*) 'VF: Fsum:',Fsum,cbc(ifc,iel,1)
					 
               endif

              enddo
           enddo



c 3rd step, calculate radiosity and irradiation, in order to clculate radiation heat flux
c
       ! T0 = 1000.0
       ! U0 = 9.21E-02
       ! rho = 0.353
       ! Cp = 1150
       ! f0 = T0*U0*rho*Cp
       ! eps = 0.85
       ! sigma = 5.67e-8 ! Stefan-Boltzmann constant
       ! fac0 = eps*sigma*T0**4.0/f0
       eps =vf_eps

       if (nid.eq.0) then
       write(6,*) 'VF: solving Jr and Gr iteratively.'
       write(6,*) 'VF: fac0:',fac0
       endif

        if (vf_crf_firstCalled.eq.1) then
        call rzero(Jr,max_walls)
        call rzero(Gr,max_walls)
		vf_crf_firstCalled = 0
        endif

        !
        !open_face_option = 1 ! ideal diffuse reflection for open face 
        !open_face_option = 2 ! radiation leaves open face but never come back, to match DOM data.


        niter = 5
	
        if (open_face_option.eq.1) then

        !niter = 10
        do iter = 1,niter
              
            call rzero(Jr,max_walls)

            do iel = 1,lelt
            eg = lglel(iel)
            do ifc= 1,6
          
            if ((cbc(ifc,iel,1).ne.'E  ')
     & .and.(cbc(ifc,iel,1).ne.'   '))  then
            !if (cbc(ifc,iel,1).eq.'W  ') then

         iwall = egf_to_iwall(ifc,eg)
         
		 ! using ideal diffuse reflection  for open surfaces (like inlet and outlet)
		 fac1 =  indicator_real_wall_ef(iwall) ! this is 0 for non-wall boundaries
         Jr(iwall) = fac1*fac0*temp_ef(iwall)**4.0
     & +(1.0-eps*fac1)*Gr(iwall)
                     
             endif
            enddo
            enddo

           call gop(Jr,wka3,'+  ',max_walls) 

           call rzero(Gr,max_walls)
                  
           do iel = 1,lelt
              eg = lglel(iel)
              do ifc= 1,6

            if ((cbc(ifc,iel,1).ne.'E  ')
     & .and.(cbc(ifc,iel,1).ne.'   '))  then

             !if (cbc(ifc,iel,1).eq.'W  ') then
                     
              iwall = egf_to_iwall(ifc,eg)
                     
               do jindex = 1,n_visible_walls(ifc,iel)
                jwall = jindex_to_jwall(jindex,ifc,iel)
                view_factor = view_factor_value(jindex,ifc,iel)

               Gr(iwall)  = Gr(iwall)+ view_factor*Jr(jwall)

               enddo                     
                     
               endif
              enddo
           enddo

       call gop(Gr,wka3,'+  ',max_walls) 


       if (nid.eq.0) write(6,*) 'VF: iter:',iter


       Jr_min = glmin(Jr,max_walls)
       Jr_max = glmax(Jr,max_walls)
       if (nid.eq.0) then
       write(6,*) 'VF: Jr_min:',Jr_min
       write(6,*) 'VF: Jr_max:',Jr_max
       endif


       Gr_min = glmin(Gr,max_walls)
       Gr_max = glmax(Gr,max_walls)
       if (nid.eq.0) then
       write(6,*) 'VF: Gr_min:',Gr_min
       write(6,*) 'VF: Gr_max:',Gr_max    
       endif
       enddo


       elseif (open_face_option.eq.2) then

        Jr_sky=0.0  ! for black hole
        !Jr_sky = fac0 !  assuming open faces are at T0

       ! get view factor for sky for all wall
	   
         do iel = 1,lelt
            eg = lglel(iel)
            do ifc= 1,6
     
            if (cbc(ifc,iel,1).eq.'W  ') then

                iwall = egf_to_iwall(ifc,eg)
                F_sky(iwall) = 1.0

               do jindex = 1,n_visible_walls(ifc,iel)
                jwall = jindex_to_jwall(jindex,ifc,iel)
                view_factor = view_factor_value(jindex,ifc,iel)
                fac1 =  indicator_real_wall_ef(jwall)                   ! this is 0 for non-wall boundaries
                F_sky(iwall) =  F_sky(iwall) -  view_factor*fac1 !  F_sky = 1- sum*(F_wall)
               enddo
			   
             endif
            enddo
            enddo
	   

        !niter = 20
        do iter = 1,niter
              
            call rzero(Jr,max_walls)

            do iel = 1,lelt
            eg = lglel(iel)
            do ifc= 1,6
          
            if (cbc(ifc,iel,1).eq.'W  ') then

         iwall = egf_to_iwall(ifc,eg)
         
		 ! using ideal diffuse reflection  for open surfaces (like inlet and outlet)
		 fac1 =  indicator_real_wall_ef(iwall) ! this is 0 for non-wall boundaries
         Jr(iwall) = fac0*temp_ef(iwall)**4.0
     & +(1.0-eps)*Gr(iwall)
                     
             endif
            enddo
            enddo

           call gop(Jr,wka3,'+  ',max_walls) 

           call rzero(Gr,max_walls)
                  
           do iel = 1,lelt
              eg = lglel(iel)
              do ifc= 1,6


             if (cbc(ifc,iel,1).eq.'W  ') then
                     
              iwall = egf_to_iwall(ifc,eg)
                     
               do jindex = 1,n_visible_walls(ifc,iel)
                jwall = jindex_to_jwall(jindex,ifc,iel)
                view_factor = view_factor_value(jindex,ifc,iel)
                fac1 =  indicator_real_wall_ef(jwall)
               Gr(iwall)  = Gr(iwall)+ view_factor*Jr(jwall)*fac1

               enddo                     
                     
               Gr(iwall)  = Gr(iwall) + F_sky(iwall)*Jr_sky
					 
               endif
              enddo
           enddo

       call gop(Gr,wka3,'+  ',max_walls) 


       if (nid.eq.0) write(6,*) 'VF: iter:',iter


       Jr_min = glmin(Jr,max_walls)
       Jr_max = glmax(Jr,max_walls)
       if (nid.eq.0) then
       write(6,*) 'VF: Jr_min:',Jr_min
       write(6,*) 'VF: Jr_max:',Jr_max
       endif


       Gr_min = glmin(Gr,max_walls)
       Gr_max = glmax(Gr,max_walls)
       if (nid.eq.0) then
       write(6,*) 'VF: Gr_min:',Gr_min
       write(6,*) 'VF: Gr_max:',Gr_max    
       endif
       enddo

       endif

       if (nid.eq.0) then
       write(6,*) 'VF: Done solving Jr and Gr.'
       endif

c 3rd step, calulate radiation heat flux for each element face
c with nondimensionalization
C  only apply radiation heat flux on wall.
c
        call rzero(frad_ef,max_walls)

        do iel = 1,lelt
           eg = lglel(iel)
           do ifc= 1,6
            if ((cbc(ifc,iel,1).ne.'E  ')
     & .and.(cbc(ifc,iel,1).ne.'   '))  then
            !if (cbc(ifc,iel,1).eq.'W  ') then
            iwall = egf_to_iwall(ifc,eg)
            frad_ef(iwall)  =  Jr(iwall) - Gr(iwall)   ! + for radiation heat flux leaving solid to fluid       
            endif

!            if ((cbc(ifc,iel,1).eq.'v  ')
!     & .or.(cbc(ifc,iel,1).eq.'O  '))  then
!            !if (cbc(ifc,iel,1).eq.'W  ') then
!            iwall = egf_to_iwall(ifc,eg)
!               
!       write(6,*) 'VF:frad_ef nonwall:',frad_ef(iwall),cbc(ifc,iel,1)
!
!            endif
			
           enddo
        enddo

       call gop(frad_ef,wka3,'+  ',max_walls) 

       frad_min = glmin(frad_ef,max_walls)
       frad_max = glmax(frad_ef,max_walls)
       if (nid.eq.0) then
       write(6,*) 'VF: frad_ef_min:',frad_min
       write(6,*) 'VF: frad_ef_max:',frad_max
       endif

c
c 4th step, map radiation heat flux from element face level to grid point level
c only apply on wall
c
        call rzero(frad(1,1,1,1),ntot)

        do iel = 1,lelt
          eg = lglel(iel)
          do ifc= 1,6
          
          if (cbc(ifc,iel,1).eq.'W  ') then
                     
             iwall = egf_to_iwall(ifc,eg)
             frad1 = frad_ef(iwall)
                      
           call facind(i0,i1,j0,j1,k0,k1,nx1,ny1,nz1,ifc)

           do k=k0,k1
           do j=j0,j1
           do i=i0,i1
             frad(i,j,k,iel) =  frad1
           enddo  
           enddo 
           enddo 

           endif
       
        enddo
       enddo

      call col2(frad(1,1,1,1),bm1,ntot)
      call dssum(frad(1,1,1,1),nx1,ny1,nz1)
      call col2(frad(1,1,1,1),binvm1,ntot)

      frad_min = glmin(frad,ntot)
      frad_max = glmax(frad,ntot)
      if (nid.eq.0) then
      write(6,*) 'VF: frad_min:',frad_min
      write(6,*) 'VF: frad_max:',frad_max    
      endif

      return
      end
!c------------------------------------------------------------------
!c-----------------------------------------------------------------------
      subroutine vf_print_radiation_heat_flux(ss)
c      implicit none
      include 'SIZE'
      include 'TOTAL'

      include "view_factors/VIEW_FACTORS"

      integer ss,e,iside
      real sint2,sarea2,totalAA,qrwt,totT

      totalAA = 0.0
      qrwt = 0.0
	  
      totT = 0.0
      do e=1,nelv
      do iside=1,2*ndim 
         if (boundaryID(iside,e).eq.ss) then
           sint2 = 0.0
           sarea2 = 0.0
           call surface_int(sint2,sarea2,frad,e,iside)
           totalAA = totalAA  + sarea2
           qrwt = qrwt + sint2
           call surface_int(sint2,sarea2,t(1,1,1,1,1),e,iside)
           totT = totT  + sint2
         endif
      enddo
      enddo

      !if (nid.eq.0) write(6,*) 'flag4-4'
      totT = glsum(totT,1)
      totalAA= glsum(totalAA,1)
      qrwt= glsum(qrwt,1)

      if (nid.eq.0) write(6,*) ' qrwt:', qrwt
      if (nid.eq.0) write(6,*) 'totalAA:',totalAA

      qrwt = qrwt /totalAA
      if (nid.eq.0) then
      write(6,*) 'sideset ',ss,' average radiation heat flux:', qrwt
      endif
      totT = totT /totalAA
      if (nid.eq.0) then
      write(6,*) 'sideset ',ss,' average temperature:', totT
      endif


      return
      end
c-----------------------------------------------------------------------