REAL FUNCTION posteriorChiSqTest(err,forest,rep,decidflag,RealityErr,numConstraints,constraints,countD)
! This is the betaV4 posterior chi-sq test
! changes: 
! - predefined cut-off array
! - only tests for constraints listed in constraints file

implicit none
  	integer, intent(in) :: numConstraints,decidflag,rep

        real, intent(in) :: RealityErr
  	real, intent(in) :: err(numConstraints,2)
  	real, intent(in) :: constraints(numConstraints+1)
	real :: cutoff(numConstraints),countD(numConstraints)
        character(len=*), intent(in) :: forest
 	integer :: i
        integer :: flag2
	
        cutoff=(/1.00546346458884,1.00692690089577,1.38317684535951,1.59871791721053,&
        		&1.47136430769351,1.44385683792429,1.44385683792429,100.0,100.0,100.0,&
        		&1.94486008493371,1.04071385522242,1.04071385522242,1.00704707310333,&
        		&1.44385683792429,1.45700207905303,1.44385683792429,100.0,1.00427857706379,&
        		&1.44385683792429,1.44385683792429,1.44385683792429,1.44385683792429,&
        		&1.44385683792429,1.44385683792429,1.44385683792429,1.44385683792429,&
        		&2.70554345409542,100.0/)
	
        flag2=1
!	if (decidflag.eq.1)then
		do i=1,numConstraints
			if (constraints(i).eq.1)then
				if (err(i,2).ge.cutoff(i))then
					flag2=0	! is the error is outside acceptable range, flag to reject
				endif
			endif
		enddo
!	endif
	if (RealityErr*constraints(30).gt.0.1)  flag2=0	! realityErr is constraint 29

	

posteriorChiSqTest=flag2
return
END


	
	
