! The FoBAAR model
! tjk, 10 Apr 2013
! MC polished 2015

module obsdrivers

! for use in printing execution time....
real, dimension(2) :: tarray
real :: result
 
! set the home directory 
character(len=1), parameter:: home = '.'	

! set as a run identifier for running different configurations
integer,parameter :: rep =1              

! the following parameters control the optimization
integer, parameter :: r = 1		                  ! number of wandering iterations (from 1 to 10000. not absolutely necessary)
integer, parameter :: q = 2                       ! number of optimization iterations (from 2 to 100000. usually about 10000 is sufficient)       
integer, parameter :: initial =1                  ! starting value 0: new optimization 1: start from previous best
integer, parameter :: explore =1 	          ! maximum number of parameter sets in explore stage (program exits if reached)

! ---------- ---------- ---------- ---------- ---------- ----------

! Here we set the length of the simulation
! Part of the time series is used for optimization
! while the other part can be used only for testing the model
character(len=*), parameter:: forest = 'HA'
integer, parameter :: NOAA=1
integer, parameter :: stday =1+365*0
integer, parameter :: stdaycal =stday+(365*14)+3
integer, parameter :: ndaycal= stday+(365*16)+3
integer, parameter :: nday=ndaycal
integer, parameter :: nyears = nday/365
integer,parameter :: subDaily = 24. 


! ---------- ---------- ---------- ---------- ---------- ----------
! ---------- 	Set the parameters for the optimization (do not change)

real, parameter :: acceptanceRate=0.15  
integer, parameter :: exploreLength =(explore/acceptanceRate)/1.3 	! maximum number of iterations in search stage 
integer, parameter :: numOuterRuns=1							! both obsolete
integer, parameter :: numInnerRuns=1
integer, parameter :: nparm = 50+10	! Modified for coupling DAMM (add 10 parameters)							! parameters dynamic based on the number of soil layers
 
! ---------- ---------- ---------- ---------- ---------- ----------
! ---------- 	Set the values of files to be read in and location of data within files. 
! ----------  Each number identifies the column location of the data in the data file
! ----------  Do not change unless the input files are changed

integer, parameter :: numConstraints = 28 + 1 ! min added for swc

integer,parameter :: numColumnsMet=11
integer,parameter :: numColumnsFlux=19
integer,parameter :: numColumnsBio=16

integer,parameter :: indexMETProjectYear=1
integer,parameter :: indexMETProjectDay=2
integer,parameter :: indexMETDayHour=3
integer,parameter :: indexMETtemp=4    ! C
integer,parameter :: indexMETpar=5     ! ppfd 
integer,parameter :: indexMETvpd=6     ! kpa
integer,parameter :: indexMETswc10=7     ! fraction
integer,parameter :: indexMETswc50=8     ! fraction
integer,parameter :: indexMETsoilT=9   ! C
integer,parameter :: indexMETco2=10    ! ppm
integer,parameter :: indexMETprecip=11 ! mm

integer,parameter :: indexFLUXnee=4    ! umol/m2/s
integer,parameter :: indexFLUXneeE=5   ! umol/m2/s
integer,parameter :: indexFLUXneegf=6  ! umol/m2/s
integer,parameter :: indexFLUXle=7     ! W/m2
integer,parameter :: indexFLUXleE=8    ! W/m2
integer,parameter :: indexFLUXlegf=9   ! W/m2
integer,parameter :: indexFLUXgpp=10   ! umol/m2/s
integer,parameter :: indexFLUXgppgf=11 ! umol/m2/s
integer,parameter :: indexFLUXre=12    ! umol/m2/s
integer,parameter :: indexFLUXregf=13  ! umol/m2/s
integer,parameter :: indexFLUXrs=14    ! umol/m2/s
integer,parameter :: indexFLUXrsE=15   ! umol/m2/s
integer,parameter :: indexFLUXAutoRs=16    ! umol/m2/s
integer,parameter :: indexFLUXAutoRsE=17   ! umol/m2/s
integer,parameter :: indexFLUXAutoTRs=18    ! umol/m2/s
integer,parameter :: indexFLUXAutoTRsE=19   ! umol/m2/s

integer,parameter :: indexBioLAI=4
integer,parameter :: indexBioLAIe=5
integer,parameter :: indexBioLitterF=6
integer,parameter :: indexBioLitterFe=7
integer,parameter :: indexBioPhenology=8
integer,parameter :: indexBioPhenologyE=9
integer,parameter :: indexBioCw=10
integer,parameter :: indexBioCwE=11
integer,parameter :: indexBioCwInc=12
integer,parameter :: indexBioCwIncE=13
integer,parameter :: indexBioSoilTotC=14
integer,parameter :: indexBioSoilTotCe=15

integer,parameter :: numSubDailyVariablesOut=35
integer,parameter :: numDailyVariablesOut=70

! ---------- ---------- ---------- ---------- ---------- ----------
! 	Extras
integer:: iter,lastrun,innerRun,outerRun  
real :: airT,ma_AirT,soilt,rad,rh,par,nit,day,ca,lat,yearday,projectday,projectyear
real :: neemeas(3),neemeasDailyDay(3),neemeasDailyNight(3),neemeasday(2),neemeasnight(2),neemeas1,neemeas2,err1,err2
real :: neemodDailyDay(1),neemodDailyNight(1)
real :: litterfallmeas(2),bbdate,laimeas(2)
real :: cwmeas,cwmeasE,cwmeasInc,cwmeasIncE,cwInc,cwPreviousYear,cwPrevious
real :: soilTotCmeas,soilTotCmeasE
real :: rsoilmeas(2), AutoRsoilmeas(2), rsoilmeasDaily(2), AutoRsoilmeasDaily(3)
real :: floatingIter,floatingAcc, floatingAns, floatingStepflag
real :: swc, swhc, swcdaily, wateredge
real :: swcmeas,swcmeasdaily(2), swcmeas10(1), swcmeas50(1), swcmeasdaily10, swcmeasdaily50

! set how many soil pools the model should used (depreciated? do not change)
integer,parameter :: numSoilPools = 3 	! 1,2 (fast, slow), 3 (fast, intermediate, slow)

! ---------- ---------- ---------- ---------- ---------- ----------
end module obsdrivers 
! ---------- ---------- ---------- ---------- ---------- ----------
! ---------- ---------- ---------- ---------- ---------- ----------


! ---------- ---------- ---------- ---------- ---------- ----------
! Begin the integrated FoBAAR MDF code
! ---------- ---------- ---------- ---------- ---------- ----------

Program FoBAAR
use obsdrivers
 
implicit none

! Initialize the many variables that the model will use.
! Any new variable introduced to the code must be initialized here.
real :: G,GDaily,GDaily2,GDailyPrevious,GsubDaily,PhotoSynth
real :: Gc,GcDaily,GcDAily2,GcSubDaily,ETsubDaily,ETDaily(2),ETmeasDaily(3),ETmeas(2)
real :: NeeDaily,NeeDailyMeas,posteriorChiSqTest
real :: PPFDsun,PPFDshd,radDaily,precipDaily
real :: VPD
real :: LAI,AssignLai,LAISun
real :: leafout,leafin,leafinmeas(2),leafoutmeas(2)
real :: Trate,TrateDaily,TrateS,TrateDailyS,TrateW,TrateF,TrateRoot,Tadj,Tadj2,gdd(2),rnorm

! autotrophic respirations
real :: Ra,Rroot,Dresp
real :: RaDaily,RrootDaily,DrespDaily
real :: iRa,iRroot,iDresp
! heterotrophic respirations
real :: Rhh,RhLit,Rh1,Rh2,Rh3
real :: RhDaily,RhLitDaily,Rh1Daily,Rh2Daily,Rh3Daily,RhhDaily,RhhDaily2
real :: iRh,iRhLit,iRh1,iRh2,iRh3
! aboveground and belowground respirations
real :: RBG
real :: RBGDaily, RBGDaily2
real :: iRBG 
! ecosystem respiration
real :: Re
real :: ReDaily
real :: iRe

real :: NEE,iGPP,iNEE(4),iNEEmeas(4)
real :: iRsoil(3),RsoilModDaily(3),iRsoilmeas
real :: TrenchedRsoilmeas(2), TrenchedRsoilmeasDaily(2), RhsoilDaily
real :: annualNEE(nyears),annualGPP(nyears),annualRa(nyears),annualRh(nyears)
real :: cwMeasFirstLast(2),cwModFirstLast(2),de,RealityErr,dayflag(nday,subDaily)
real :: P(nparm),incP(nparm),stepSize(nparm),absoluteBoundsP(nparm,2),boundsP(nparm,2)
real :: parameters(numInnerRuns,(q*(numOuterRuns)+(exploreLength)),nparm)
real :: error(numInnerRuns,(q*(numOuterRuns)+exploreLength),3)
real :: lastvalueCwMeas,lastvalueCwMod
real :: constraints(numConstraints+1)

real :: accparm(explore,numInnerRuns,nparm+1),oldP(nparm),bestP(nparm),deltap(nparm)
real :: LMA,max_fol,multtl,multtf,ralabfrom,ralabto, ralabtoDaily,ralabfromDaily
real :: Af,Aw,Ar,Lf(2),Lw,iLw,Lr,iLr,Atolab,Afromlab,npp,iAf,iAw,iAr
real :: Cf,Cw,Cr,Clab,Clit,Csom,CsomPools(numSoilPools)
real :: annealstart,annealtemp,cooling(1),besterr,allbest
real :: err(numConstraints,2),toterr,TotalErrorV4,TotalErrorV4bayes,countD(numConstraints),qerr(numConstraints)
real :: ran,expo,fAparDaily
real :: Xfang,Vcmax,EaVcmax,EdVcmax,EaJmax,EdJmax,SJmax,Rd,Rdt,gs,VQ10    	
real :: Dlit,D1,D2,iD1,iD2,iDlit	

real :: innerRunParameters(2,numInnerRuns,nparm)
real :: innerRunIncP(2,numInnerRuns,nparm)
real :: innerRunBoundsP(2,numInnerRuns,nparm,2)
real :: bestInnerRunTotError

real :: bSoilMeas(3)
real :: prange(nparm,2)
real :: innerRunOutSubDaily(numInnerRuns,nday,subDaily,numSubDailyVariablesOut)
real :: BestSubDaily(numInnerRuns,nday,subDaily,numSubDailyVariablesOut)
real :: innerRunOutDaily(numInnerRuns,nday,numDailyVariablesOut)
real :: pred(numInnerRuns,nday,numDailyVariablesOut)
real :: BestDaily(numInnerRuns,nday,numDailyVariablesOut)
real :: posteriorFluxComponents(explore,nyears*4)         ! the posterior flux components are NEE GPP Ra Rh
real :: posteriorErrorTerms(explore,numConstraints)         ! the posterior error terms
real :: cl(nday,numDailyVariablesOut,2)
real :: clLow(numInnerRuns,nday,numDailyVariablesOut)
real :: clHigh(numInnerRuns,nday,numDailyVariablesOut)
real :: Raratio_high(365),Raratio_low(365)

real :: cl_sub(nday,subDaily,numSubDailyVariablesOut,2)
real :: clLow_sub(numInnerRuns,nday,subDaily,numSubDailyVariablesOut)
real :: clHigh_sub(numInnerRuns,nday,subDaily,numSubDailyVariablesOut)


real :: metData(nday,subDaily,numColumnsMet)
real :: fluxData(nday,subDaily,numColumnsFlux)
real :: bioData(nday,numColumnsBio)
real :: soilTemp(nday,subDaily)
real :: maxt(nday),mint(nday)
real :: period, hour							! the number of time intervals during a day

! Modification for coupling DAMM
real :: alpha(8), Ea(4), KMcp, KMO2
real :: BD,PD,Dliq,Dgas,O2,theta,Enzpool,Enz
! Modification end

real :: sf

real :: drainage,runoff,soilWaterContent                                      ! soil water variables
real :: waterStressP

! temporary holding variables
real :: xx,tmp,tmp2,tmp3,tmp4					
real :: longTerm                        ! this will be set to 0 for normal run, 1 for longTerm run (>50 yrs)

real:: DailyAverageData
real:: DailyAverageModel
real:: CountN
real:: DailyAverageTrenchData
real:: DailyAverageTrenchModel



double PRECISION :: toterrAll

integer :: acc(numInnerRuns),ans,stepflag,laststep,acceptflag,jfunct,decidflag,DecidOrEvg
integer :: i,ii,j,jj,k,firstexplore,startIteration,endIteration
integer :: bestInnerRun
integer :: countParamSets,countTotalIterations,countExplore,flag
integer :: nparams
integer :: year

integer :: seed_size
integer,allocatable :: seed(:)

character(len=150):: filename,filename1
character(len=15):: filenum(5)
character(len=15):: repChar

call random_seed() 		 		! initialize with system generated seed
call random_seed(size=seed_size) 	! find out size of seed
allocate(seed(seed_size))
call random_seed(get=seed) 	 	! get system generated seed

! depending on the number of soil pools, set the number of parameters.
if(numSoilPools.eq.1)then
	nparams=40
else if (numSoilPools.eq.2) then
	nparams=42
else
	nparams=50+10 ! Modification for coupling DAMM (adding 10 parameters)
endif

! define the time-integrating period for rates, etc.
period=24./subDaily

! check whether this is a longTerm or standard run
longTerm = 0.0
if (nyears.gt.50)then
    longTerm=1.0
endif

! set whether the current site is deciduous or evergreen
DecidOrEvg = 0 ! evergreen
decidflag = DecidOrEvg
jfunct = 2

write(*,*) 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Read the driver data from file

call readObservedData(bioData,maxt,mint,&
     &fluxData,metData,Dayflag,soilTemp,longTerm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! print information about the current simulation to screen
write(repChar,*)rep
print *, "Running ...."
print *, 'Site Name:'
print *, forest
print *, 'Decid (1)/Evergreen (0)'
print *, decidflag
print *, 'Rep:'
print *, rep


! Load the data that tells the model which constraints to use in current simulation
open(unit=1,file=&
	&home//'/Constraints/'//&
	&'constraints_'//trim(adjustl(repChar))//'.csv',status='old')
	read(1,*)	! skip the header info
	read(1,*)(constraints(i),i=1,numConstraints+1)
close(1)

! Modified for coupling DAMM
Dliq=3.17
Dgas=1.67
BD=0.80
PD=2.52
! Modification end


countTotalIterations=1
OUTER: DO outerRun = initial,numOuterRuns

	! INPUT PARAMETERS
	! only on first pass...
	if(initial.eq.1)then
		if (decidflag.eq.1) then
			open(unit=33,file=home//'/Results/initialRun_D_'//trim(adjustl(repChar))//'.csv',status='old')	
		else
			open(unit=33,file=home//'/Results/initialRun_E_'//trim(adjustl(repChar))//'.csv',status='old')		
		endif
	
		read(33,*)(P(i),i=1,nparams)
		close(33)
	endif
 
	! initialize the parameter bounds
	if(outerRun.eq.initial)then		! read in only on first pass
	
		open(unit=27,file=home//'/initial_P_rangesV5rep'//trim(adjustl(repChar))//'.csv',status='old')
                do i=1,nparams
                        read(27,*) absoluteBoundsP(i,1),absoluteBoundsP(i,2),boundsP(i,1),boundsP(i,2)
                end do	
		read(27,*) lat
		lat=lat*3.14/180.0		!convert to radians
		close(27)
	endif

	! reset incP at 5% of bounds. This is the initial step size for the markov chain walk.
	incP(:)=0.075*(boundsP(:,2)-boundsP(:,1))

	if (outerRun.eq.0) then
		P(:)=(boundsP(:,1)+boundsP(:,2))/2.
		
		! 40% of initial parameter range. This value is arbitrary, but 40% gives fast convergence.
		incP(:)=0.4*(boundsP(:,2)-boundsP(:,1)) 
	endif

	prange(:,1)=P(:)
	prange(:,2)=P(:)

	firstexplore = 0
	besterr=9999999999.
	bestInnerRunTotError=9999999999.
	allbest=besterr
	bestInnerRun=0

	InnerRunIncP=0
	InnerRunBoundsP=0
	InnerRunParameters=0
	
	accparm=0
	acc=1
	oldp=P
	bestp=P
	err=0

	! define the initial anneal temperature
	annealstart=100
	annealtemp=annealstart

	! define the speed at which the temperature cools
	cooling(1)=.001**(1/real(q-r))	

	! initialize variables
	qerr=1    
	cl(:,:,1)=-10E6
	cl(:,:,2)=10E6
	clLow=-10E6
	clHigh=10E6

        cl_sub(:,:,:,1)=-10E6
        cl_sub(:,:,:,2)=10E6
        clLow_sub=-10E6
        clHigh_sub=10E6


	Raratio_high=-10e6
        Raratio_low=10e6

	countParamSets=1 
	stepflag = 0
	laststep=1
	acceptflag = 0
	startIteration=1
	stepSize=0.05
	tmp4=0


	Do innerRun=1,numInnerRuns
		countExplore=0
		countParamSets=1 
		
		if ((outerRun.gt.0)) then
			startIteration=r
			endIteration=q
		else
			endIteration=r
		endif
	
		! ---------------------------------------------------------
		! ---------------------------------------------------------
		! ---------- This is the start of the optimization loop
		! ---------- Everything within this loop is performed with a different set of model parameters
	
		mcloop: do iter=startIteration,endIteration+exploreLength
		
			! this block of code generates the set of parameters for the current iteration 
			if ((iter.gt.1).and.(iter.le.q)) then
				if (iter.gt.r) then	! cool anneal temperature
					annealtemp = annealtemp*cooling(1)
				endif
				
				do j=1,nparams	! move from current parameter set
					if (stepflag.eq.0) then
						deltap(j) = rnorm()*incp(j)
					endif
					
					if (stepflag.eq.1) then
						deltap(j) = deltap(j)
					endif
		
					p(j) = oldp(j)+deltap(j)
										
					! if the parameters fall outside the absolute limit, bounce them back
					if (p(j).lt.boundsp(j,1)) then 
						p(j) = boundsp(j,1)+(boundsp(j,1)-p(j))/2
					endif
					if (p(j).gt.boundsp(j,2)) then
						p(j) = boundsp(j,2)+(boundsp(j,2)-p(j))/2
					endif
!						
					if(boundsp(j,1).eq.boundsp(j,2)) then
					    p(j) = boundsp(j,1)
					endif
				end do	
			
				if ((iter.eq.r).or.(iter.eq.endIteration)) then 
					P = bestP
					stepflag = -1
				endif
			endif
			
			! Post-optimization, the parameters can explore beyond the prior parameter ranges
			! this block of code generates the set of parameters for the post-optimization exploration of the parameter space
			if(iter.gt.q) then		! in parameter space explore stage
				do j = 1,nparams	! move parameter
					deltap(j) = rnorm()*incp(j)
        			        p(j) = oldp(j)+deltap(j)
				
					if((iter-q).lt.5000) then
						tmp4 = 0.15	! initial fast expansion of range
					else
						tmp4 = 0.1
					endif
							
					! if the parameters fall more than 10% outside the range, bounce them back
					if (p(j).lt.(prange(j,1)-(tmp4*abs(prange(j,1))))) then 
						p(j) = prange(j,1)+(prange(j,1)-p(j))/2
					endif
					if (p(j).gt.prange(j,2)+(tmp4*abs(prange(j,2)))) then
						p(j) = prange(j,2)+(prange(j,2)-p(j))/2
					endif

					if(boundsp(j,1).eq.boundsp(j,2)) then
						p(j) = boundsp(j,1)
					endif

					! set absolute bounds
					if (p(j).lt.(absoluteBoundsP(j,1))) then
						p(j) = absoluteBoundsP(j,1)
					endif
					if (p(j).gt.(absoluteBoundsP(j,2))) then
						p(j) = absoluteBoundsP(j,2)
					endif		
				end do		
			endif

! Modified for coupling DAMM
alpha(1)=10**p(51)
alpha(2)=10**p(52)
alpha(3)=10**p(53)
alpha(4)=10**p(54)
Ea(1)=p(55)
!Ea(2)=p(56)
!Ea(3)=p(57)
!Ea(4)=p(58)
Enzpool=p(58)
KMcp=p(59)
KMO2=p(60)
! Modification end

			! the following code initializes variables to parameter values 
            ! and zeros other variables 
			if(decidflag.eq.1) then
			    Cf = 0
			else
			    Cf = p(12)
			endif

                        !write(*,*) Cf
		

			Cr = P(19)
			Cw = P(20)
			Clit = P(21)
			Clab = P(23)
			swc=p(16)                                ! initialize soil water content to be at holding capacity
            swhc=p(16)
           
			CsomPools(1) = P(22) 
			CsomPools(2) = P(42) 
			CsomPools(3) = p(43) 
			
 ! pseudo (assumed!) constants 
			Xfang=0.5	
			Vcmax= P(30)
			EaVcmax= 76459
			EdVcmax=220121	
			EaJmax= 44370
			EdJmax= 213908
			SJmax= 710 
			Rd= P(36)*0.01
			VQ10= P(37)  
				
			G = 0

			iGPP = 0
			iNEE = 0
			iNEEmeas=0
			iRa = 0
            iRroot = 0
  			iDresp=0
			iRh = 0
			iRhLit=0
			iRh1 = 0
			iRh2 = 0
			iRh3 = 0
            iRe = 0
			iD1=0
			iD2=0
			iDlit=0
			iRsoil = 0
			iRsoilmeas=0
			iLw=0
                        iLr=0
                        iAw=0
                        iAr=0
                        iAf=0

			atolab=0
			afromlab=0
			Lf=0
			ca = 358
			err = 0
			RealityErr=0
			toterr=0
			leafout=0
			leafin=0
			cwMeasFirstLast=-10E6
			cwModFirstLast=-10E6
			lastvalueCwMeas=0
			lastvalueCwMod=0
			cwPreviousYear=0
                        cwPrevious=0
                        cwInc=0
			year=0

                        GDailyPrevious=0
			
			countD=0	! occurance of observations set to zero for each loop

			! -----------------------------------------------------------------
			! -----------------------------------------------------------------
			! ------ We have now assigned the parameters and initialized variables
			! ------ So we run the model for the current parameter set				
			
			! loop through each day of the simulation
			timeloop: do i=stday,nday         
				
				! check the current day of year	
				yearday=metData(i,1,indexMETProjectDay)
				
				projectyear=metData(i,1,indexMETProjectYear)	
				projectday = metData(i,1,indexMETProjectDay)	
				
				! reset some variables at the start of each year
				if((yearday.ge.1.0).and.(yearday.lt.2)) then			
        		                year = year+1
					leafin=0
					leafout=0
					leafinmeas=-999
					leafoutmeas=-999
					ma_AirT=0
					cwPreviousYear=cw
				endif
                              
                                sf = P(38)
 			
				! calculate moving average airTemperature for leaf drop estimate
				! start from mid-point of each year (only used for leaf fall date estimates) 
				if (yearday.gt.150) then
					ma_AirT = (sum(mint(i-10:i)))*0.1
				endif

                ! grab the wood biomass measurements	
				if(bioData(i,indexBioCw).gt.0) then
					cwmeas = bioData(i,indexBioCw)
                    cwmeasE = bioData(i,indexBioCwE)
				else
					cwmeas = -999
				endif
	
               
                ! grab the wood increment measurements	
				if(bioData(i,indexBioCwInc).gt.0) then
					cwmeasInc=bioData(i,indexBioCwInc)
					cwmeasIncE=bioData(i,indexBioCwIncE)
				else
					cwmeasInc=-999
				endif
				
                ! grab the soil c measurements
				if(bioData(i,indexBioSoilTotC).gt.0) then
					soilTotCmeas =bioData(i,indexBioSoilTotC)
                    soilTotCmeasE =bioData(i,indexBioSoilTotCe)
				else
					soilTotCmeas = -999
				endif


				litterfallmeas(1) = bioData(i,indexBioLitterF)
				if(litterfallmeas(1).gt.0) then
					litterfallmeas(2) = bioData(i,indexBioLitterFe)
				else	
					litterfallmeas(2) = 0
				endif
			
				! grab bud-burst and leaf fall measurements
				if((bioData(i,indexBioPhenology).eq.1))then
					leafoutmeas(1) = yearday
					leafoutmeas(2) = bioData(i,indexBioPhenologyE)
				endif
				if((bioData(i,indexBioPhenology).eq.2))then
					leafinmeas(1) = yearday
					leafinmeas(2) = bioData(i,indexBioPhenologyE)
				endif
				

				! at new year, start incremental values
				if((i.eq.stday).or.((yearday.ge.1.0).and.(yearday.lt.2))) then
			        ! if at end of year, record annual values
			        if(year.gt.1)then
				        annualNEE(year-1)=iNEE(1) ! this is the annual NEE
				        annualGPP(year-1)=iGPP
				        annualRa(year-1)=iRa
				        annualRh(year-1)=iRh
			        endif
			        
					iGPP = 0
					iNEE = 0
					iNEEmeas=0
                                        iRa = 0
                                        iRroot = 0
                                        iDresp=0
                                        iRh = 0
                                        iRhLit=0
                                        iRh1 = 0
                                        iRh2 = 0
                                        iRh3 = 0
                                        iRBG = 0
                                        iRe = 0
					iD1=0
					iD2=0
					iDlit=0
					iRsoil = 0
					iRsoilmeas=0
					iLw=0
                                        iLr=0
                                        iAw=0
                                        iAf=0
                                        iAr=0
					Lf(2)=0
				endif
	 				
				! calculate the growing degree days (gdd)
				! max_fol, multtf, multtl values. These control phenology later in the code
				if (decidflag.eq.1) then
					if(yearday.le.p(25)) then
						gdd=0.
						max_fol=1.
					endif
				  
					!time switch defaults
					multtf=1.		! turnover of foliage on
					multtl=0.		! turnover of labile C off
					
				 	gdd(1)=gdd(1)+0.5*max(maxt(i)+mint(i),0.0)		!growing degree day heat sum from day p(25)
			
					if (gdd(1).ge.p(12)) then   !winter end (if we have got to the end of winter temperature)
						if (max_fol.eq.1) then	!spring
							if (leafout.eq.0) then
								leafout=yearday
							endif
							multtl=1.
							multtf=0.
						else			        !summer
							multtl=0.
							multtf=0.
						endif
					endif
			
					if (yearday.ge.200) then	
						max_fol=0.
						multtl=0.
					endif 

					if(((yearday.ge.200).and.(ma_AirT.lt.p(13))).or.&
					&((yearday.ge.200).and.(leafin.gt.0))) then		!drop leaves
						multtf=1.
						if(leafin.eq.0) then
							leafin=yearday+9
						endif
					endif
				endif
					
				if (decidflag.eq.0) then
					multtl = 1
					multtf = 1
					p(14) = 1
					if (yearday.le.p(25)) then
						gdd(1)=0.
					endif
					gdd(1)=gdd(1)+0.5*max(maxt(i)+mint(i),0.0)
				endif
				if(decidflag.eq.1)then	
					tadj = 1
					tadj2 =1
				else
						tadj =max(0.0,min(1.0,(p(33)/maxt(i))))
						tadj2 =min(1.0,max(0.0,maxt(i)/40))
				endif
				
				if ((yearday.gt.1).and.(yearday.lt.2)) then 
					gdd(2)=0
				endif
				
				gdd(2)=min(p(27),max(gdd(2)+(maxt(i)+mint(i))/2.,0.0))
				
				LMA = P(24)
				LAI=max(0.001,Cf/LMA)

				laimeas(1)=bioData(i,indexBioLAI)		
				laimeas(2)=bioData(i,indexBioLAIe)
				if((decidflag.eq.0).and.(laimeas(1).lt.0))then
					laimeas(1)=-999
					laimeas(2)=-999
				endif
				GDaily=0
				GDaily2=0
				GcDaily=0
				GcDaily2=0
				ETDaily=0
				ETmeasDaily=0
				DrespDaily=0
				RaDaily = 0
                RrootDaily = 0
                DrespDaily = 0
                RhDaily = 0
				RhLitDaily=0
				Rh1Daily= 0
				Rh2Daily=0
				Rh3Daily=0
				RBGDaily = 0
                RBGDaily2 = 0
				RhhDaily = 0
                RhhDaily2 = 0
                ReDaily = 0
				RsoilModDaily=0
				ralabtoDaily=0
				ralabfromDaily=0
				NeeDaily=0
				NeeDailyMeas=0
				radDaily=0
				fAparDaily=0
				precipDaily=0
				swcdaily = 0
                RhsoilDaily=0
				
				DailyAverageData=0
				DailyAverageModel=0
				CountN=0
				DailyAverageTrenchData=0
				DailyAverageTrenchModel=0
			
				
				
				Ca = metData(i,1,indexMETco2)  		! ppm
				if (Ca.lt.100) then
					Ca = 358	! may be no data in Ca driver data
				endif

				neemeasDailyNight=0
				neemeasDailyDay=0
				neemodDailyNight=0
				neemodDailyDay=0

				swcmeas = 0
				swcmeasdaily = 0
				swcmeas10 = 0
				swcmeasdaily10 = 0
				swcmeas50 = 0
				swcmeasdaily50 = 0
				
				rsoilmeasDaily = 0
				AutoRsoilmeasDaily = 0
                TrenchedRsoilmeasDaily = 0
				hour = 0
				
				hourloop: do j = 1,subDaily
					hour = hour + 1.0
					GsubDaily = 0
					Dresp = 0
					GcSubDaily = 0
					ETsubDaily = 0
					
					! assign met data for current hour (all MET data sould be filled apriori)
					rad = metData(i,j,indexMETpar)   		  ! this is PAR (g_filled) in e-6mol/m2/s
					airT = metData(i,j,indexMETtemp)    	  ! air temperature, degrees C
					VPD = metData(i,j,indexMETvpd) 
					soilT = metData(i,j,indexMETsoilT)
					
                                        ! Modification for reading soil water
                                        ! content
                                        swcmeas10(1)=metData(i,j,indexMETswc10)
                                        swcmeas50(1)=metData(i,j,indexMETswc50)
                                        swcmeasdaily10=swcmeasdaily10+swcmeas10(1)/period
                                        swcmeasdaily50=swcmeasdaily50+swcmeas50(1)/period
                                        ! Modification end
					
					! daily sums
					radDaily=radDaily+rad
                    precipDaily=precipDaily+metData(i,j,indexMetPrecip)
						  				
					! *************************************************************************************************
					! Carbon fluxes
					! *************************************************************************************************
					LAIsun=0
					PPFDsun=0
					PPFDshd=0
					
					if(swcmeas50(1).gt.p(16))then
						waterStressP=1
					else
						waterStressP=(swcmeas50(1)/p(16))**p(46) 
				        endif
					
					! this section of the code runs the photosynthetic subroutine
					! run this code only if it is currently daytime and leaves are on the trees
       				if(lai.gt.0.001) then
						if(rad.gt.0.5) then
							! assign sunlit and shaded fractions of total LAI
							LAIsun = AssignLAI(lai,xfang,yearday,lat,j,&
											   &subDaily,rad,PPFDsun,PPFDshd)
							fAparDaily=fAparDaily+((PPFDsun+PPFDshd))

							! sun leaf photosynthesis
							G=PhotoSynth(airT,PPFDsun,VPD,Ca,Nit,LAI,Vcmax*WaterStressP,&
										EaVcmax,EdVcmax,EaJmax,EdJmax,SJmax,&
										&Rd,Rdt,VQ10,p(31),p(32),p(34),Gc)*&
										&(0.0432396*period)*tadj*gdd(2)/p(27)
							G=max(G,0.0)
							Gc=max(Gc,0.0)
							Gc = Gc * 3600 * period
                            GsubDaily = (G*LAIsun)
							GcSubDaily = (Gc*LAIsun)
							
							!Dresp=Rdt*0.0432396*period
							Dresp = 0

							GDaily2=GDaily2+G
							GcDaily2=GcDaily2+Gc
								
							! shade leaf photosynthesis
							G=PhotoSynth(airT,PPFDshd,VPD,Ca,Nit,LAI,Vcmax*WaterStressP,&
										&EaVcmax,EdVcmax,EaJmax,EdJmax,SJmax,&
										&Rd,Rdt,VQ10,p(31),p(32),p(34),Gc)*&
										&(0.0432396*period)*tadj*gdd(2)/p(27)
							G=max(G,0.0)
                            Gc=max(Gc,0.0)
                            Gc = Gc * 3600 * period
							GsubDaily = GsubDaily+ (G*(LAI-LAIsun))
							GcSubDaily = GcSubDaily+ (Gc*(LAI-LAIsun))
							
							!Dresp=Dresp+Rdt*0.0432396*period								
							ETsubDaily=GcSubDaily
							
							GDaily = GDaily + GsubDaily
							DrespDaily = DrespDaily+Dresp
							GcDaily = GcDaily+GcSubDaily
							
							ETDaily(1) = ETDaily(1)+ETSubDaily
							ETmeasDaily(1) = ETmeasDaily(1)+&
										&fluxData(i,j,indexFLUXlegf)/2.5e6*period*3600
							if (fluxData(i,j,indexFLUXle) .gt. -999) then
								ETmeasDaily(2) = ETmeasDaily(2)+&
										&fluxData(i,j,indexFLUXle)/2.5e6*period*3600
                                ETmeasDaily(3) = ETmeasDaily(3)+&
										&fluxData(i,j,indexFLUXleE)/2.5e6*period*3600
								ETdaily(2)=ETdaily(2)+ETsubDaily
							endif
						endif
					endif
               
			
                                        ! Root respiration	
                                        TrateRoot=0.5*exp(p(29)*soilT)
                                        Rroot=(10**p(18))*Cr*TrateRoot*period
					
					! respiration when carbon is moved to or from the labile carbon pools
					ralabfrom = 0
					ralabto = 0
					ralabfromDaily = ralabfromDaily + ralabfrom
					ralabtoDaily = ralabtoDaily + ralabto
									
                                        ! calculate autotrophic respiration
                                        Ra = Dresp + ralabfrom + ralabto + Rroot + P(2)*(GsubDaily-Dresp)

                                      
                   !calculate heterotrophic respiration
                            !Trate=0.5*exp(p(10)*airT)
                            !RhLit = (10**p(8))*Clit*Trate*period
                                        
					        !TrateS=0.5*exp(p(28)*soilT)
					        !Rh1=(10**p(9))*CsomPools(1)*TrateS*tadj2*period	
				            !Rh2=(10**p(26))*CsomPools(2)*TrateS*tadj2*period       
					        !Rh3=(10**p(44))*CsomPools(3)*TrateS*period	
							
					! Modified for coupling DAMM
                    theta=swcmeas10(1)
                    O2 = Dgas * 0.209 * ((1-BD/PD-theta)**(4/3))
                    Enz=Enzpool*Dliq*(theta**3)

					RhLit = alpha(1)*exp(-Ea(1)/(8.3e-3*(soilT+273.15))) * &
					&Clit * (Enz/(Enz+KMcp)) * (O2/(KMO2+O2))*period
								
							
                    Rh1 = alpha(2)*exp(-Ea(1)/(8.3e-3*(soilT+273.15))) * &
                    &CsomPools(1) * (Enz/(Enz+KMcp)) * (O2/(KMO2+O2))*period
								
																
                    Rh2 = alpha(3)*exp(-Ea(1)/(8.3e-3*(soilT+273.15))) * &
                    &CsomPools(2) * (Enz/(Enz+KMcp)) * (O2/(KMO2+O2))*period
								
							
								
                    Rh3 = alpha(4)*exp(-Ea(1)/(8.3e-3*(soilT+273.15))) * &
                    &CsomPools(3) * (Enz/(Enz+KMcp)) * (O2/(KMO2+O2))*period
								
					! Modification end
                                       

                                        Rhh = RhLit + Rh1 + Rh2 + Rh3
                
                                        ! calculate ecosystem respiration
                                        Re = Ra + Rhh
					
					! add these hourly values to the daily totals
                    RrootDaily = RrootDaily + Rroot
					RaDaily = RaDaily+Ra
					RhLitDaily = RhLitDaily+RhLit
					Rh1Daily = Rh1Daily+Rh1
					Rh2Daily = Rh2Daily+Rh2
					Rh3Daily = Rh3Daily+Rh3
                    RBGDaily = RBGDaily + RBG
         			RhhDaily = RhhDaily + Rhh
					ReDaily = ReDaily + Re

					
					! calculate the net ecosystem exchange by summing respiration and assimilation
					nee = Ra+Rhh-GsubDaily
					neeDaily = neeDaily+nee
					
					! assign the observed flux data for the current hour to neemeas
					if ((fluxData(i,j,indexFLUXnee).gt.-999)) then	! check quality control flags	
						neemeas(1)=fluxData(i,j,indexFLUXnee)*0.0432396*period		! e-6mol/m2/s convert to gC m-2 period-1i
                                                neemeas(2)=fluxData(i,j,indexFLUXneeE)*0.0432396*period
					else
						neemeas(1)=-999
						neemeas(2)=1
					endif

                                        if (dayflag(i,j).eq.0) then
                                                neemeasnight(1)=neemeas(1)
                                                neemeasnight(2)=neemeas(2)/1.4
                                                neemeasday(1)=-999
                                                neemeasday(2)=1
                                        else
                                                neemeasday(1)=neemeas(1)
                                                neemeasday(2)=neemeas(2)/0.6
                                                neemeasnight(1)=-999
                                                neemeasnight(2)=1
                                        endif

                                        if (fluxData(i,j,indexFLUXle) .gt. -999) then
                                                ETmeas(1) = fluxData(i,j,indexFLUXle)/2.5e6*period*3600
                                                ETmeas(2) = fluxData(i,j,indexFLUXleE)/2.5e6*period*3600
                                        else
                                                ETmeas(1)=-999
                                                ETmeas(2)=1
                                        endif
					
					! sum up the observed hourly NEE to daily, monthly and annual
					if(fluxData(i,j,indexFLUXneegf).gt.-999)then
                                           neeDailyMeas=neeDailyMeas+(fluxData(i,j,indexFLUXneegf))*(0.0432396*period)
					else
                                            neeDailyMeas=0
                                        endif

					iNEE(1) = iNEE(1) + nee
					! note NEE is not fully gap filled for Howland, so here includes -9999s
					! iNEEmeas(1) therefore does not make sense, so depreciated for now 
					iNEEmeas(1)=0 !iNEEmeas(1)+(fluxData(i,j,indexFLUXneegf))*(0.0432396*period)
					
					if(fluxData(i,j,indexFLUXnee).gt.-999) then	! annual NEE constraint non gap filled 
						iNEE(2) = iNEE(2) + nee
						iNEEmeas(2)=iNEEmeas(2)+(fluxData(i,j,indexFLUXnee))*(0.0432396*period)
					endif
					
					if(longTerm.eq.0) then
						if((mod(yearday,30.).ne.0).and.(fluxData(i,j,indexFLUXnee).gt.-999))then		! 30 day Nee constraint
							iNee(3) = iNee(3)+nee
							iNEEmeas(3)=iNEEmeas(3)+(fluxData(i,j,indexFLUXnee))*(0.0432396*period)
						endif
		                        endif
										
					iRsoil(1) = iRsoil(1) + RBG

					
					! assign the observed autochamber flux data for the current hour to AutoRsoilmeas
					if ((fluxData(i,j,indexFLUXAutoRs).gt.-999)) then       ! check quality control flags
                    	AutoRsoilmeas(1)=fluxData(i,j,indexFLUXAutoRs)*0.0432396*period*sf ! e-6mol/m2/s convert to gC m-2 period-1
                        AutoRsoilmeas(2)=fluxData(i,j,indexFLUXAutoRsE)*0.0432396*period*sf
                    else
						AutoRsoilmeas(1)=-999
                        AutoRsoilmeas(2)=1
                    endif

                    ! if there is an autochamber measurement, record RBG and the meas
					if ((AutoRsoilmeas(1).gt.-999) .and. (TrenchedRsoilmeas(1).gt.-999)) then
    					RBGDaily2=RBGDaily2+RBG
    					RhhDaily2=RhhDaily2+Rhh
    					CountN=CountN+1
					end if


                                        ! assign the observed Trenched
                                        ! autochamber flux
                                        ! data for the current hour to
                                        ! TrenchedRsoilmeas
                    if((fluxData(i,j,indexFLUXAutoTRs).gt.-999).and.(fluxData(i,j,indexFLUXAutoRs).gt.-999)) then       ! check quality control flags
                    	TrenchedRsoilmeas(1)=fluxData(i,j,indexFLUXAutoTRs)*0.0432396*period*sf!  e-6mol/m2/s convert to gC m-2 period-1
                        TrenchedRsoilmeas(2)=fluxData(i,j,indexFLUXAutoTRsE)*0.0432396*period*sf!
                    else
                    	TrenchedRsoilmeas(1)=-999
                        TrenchedRsoilmeas(2)=1
                    endif
					
					neemeas1=-999
                                        neemeas2=-999
                                        err1=-999
                                        err2=-999
                    if (neemeasday(1).gt.-999) then
                    	if ((((neemeasday(1)-nee)/(neemeasday(2)))**2) .lt. 9999) then
                        	neemeas1=neemeasday(1)
                            err1=(((neemeasday(1)-nee)/(neemeasday(2)))**2);
                            err(1,1) = err(1,1)+(((neemeasday(1)-nee)/(neemeasday(2)))**2)
                            countD(1)=countD(1)+1
                        endif
                    endif
                    
                    if ((neemeasnight(1).gt.-0.1) .and. (neemeasnight(1).lt.0.4)) then
                                                if ((((neemeasnight(1)-nee)/(neemeasnight(2)))**2) .lt. 9999) then
                                                        neemeas2=neemeasnight(1)
                                                        err2=(((neemeasnight(1)-nee)/(neemeasnight(2)))**2)
                                                        err(2,1) = err(2,1)+(((neemeasnight(1)-nee)/(neemeasnight(2)))**2)
                                                        countD(2)=countD(2)+1
                                                endif
                                        endif

                                        if (ETmeas(1).gt.-999) then
                                                err(19,1) = err(19,1)+(((ETmeas(1)-ETsubDaily)/(ETmeas(2)))**2)
                                                countD(19)=countD(19)+1
                                        endif
                                        
                                        ! trenched rsoil flux errors
                                        if (TrenchedRsoilmeas(1).gt.0) then
                                                err(13,1) = err(13,1)+((TrenchedRsoilmeas(1)-Rhh)/(TrenchedRsoilmeas(2)))**2
                                                countD(13)=countD(13)+1
                                        endif

                            
                                        ! Autochamber rsoil flux errors
                                        if (AutoRsoilmeas(1).gt.0) then
                                                err(14,1) = err(14,1)+((AutoRsoilmeas(1)-RBG)/(AutoRsoilmeas(2)))**2
                                                countD(14)=countD(14)+1
                                        endif
					
					! anything with an 'i' before it is an annual cumulative variable 
					iGPP = iGPP + GsubDaily
					iRa = iRa + Ra
                    iRroot = iRroot + Rroot
					iDresp=iDresp+Dresp
					iRh = iRh + Rhh
					iRh1=iRh1+Rh1
					iRh2=iRh2+Rh2
					iRh3=iRh3+Rh3
					iRhLit=iRhLit+RhLit
					iRe = iRe + Re
					

					if (iter.eq.endIteration) then  ! save the hourly data for output after inner loop
						innerRunOutSubDaily(innerRun,i,j,:) = &
							&(/metData(i,1,1),yearday,hour,dayflag(i,j),GsubDaily,neemeas(1),& !6
							&AutoRsoilmeas(1),TrenchedRsoilmeas(1),nee,Re,Ra,RBG,Rhh,& !13
							&neemeas1,neemeas2,neemeasday(2),neemeasnight(2),neemeas(2),& !18
							&litterfallmeas(1),Lf(2),cwmeas,cw,AutoRsoilmeas(2),TrenchedRsoilmeas(2),leafinmeas(1),leafin,& !26							
							&laimeas(1),lai,ETmeas(1),ETsubDaily,Rroot,RhLit,Rh1,Rh2,Rh3/) !35 
					endif

				end do hourloop

				GDailyPrevious=GDaily-DrespDaily
				
				! calculate the temperature rates to be used for decomposition and respiration
				TrateDaily=0.5*exp(p(11)*(0.5*(maxt(i)+mint(i))))
				TrateDailyS=0.5*exp(p(41)*(0.5*(maxval(soilTemp(i,1:24))+minval(soilTemp(i,1:24)))))
						
				! Calculate the transfer of Litter from the litter pool to the fast soil pool	
				Dlit = (10**p(1))*Clit*TrateDaily
				

				! Set transfer rate between soil pools
             			if(numSoilPools.ge.2)then	
					D1 = (p(35))*Rh1Daily !TrateDailyS	
					if(numSoilPools.eq.3)then
					    D2 = (p(45))*Rh2Daily !TrateDailyS	
					endif
				endif	
				iD1=iD1+D1
				iD2=iD2+D2
				iDlit=iDlit+Dlit
						
				! Daily allocation
						npp = GDaily-DrespDaily
                                                if(GDaily.gt.DRespDaily)then
                                                        npp=npp-P(2)*(GDaily-DRespDaily)	
						endif
						
						if ((multtf.gt.0).and.(decidflag.gt.0))then
							Atolab = (1.-p(14))*(10**p(5))*Cf	
						else
							Atolab =0
						endif

						if(multtl.gt.0)then	

						Afromlab=(10**p(15))*Clab*decidflag ! simplifying the leaf allocation routine here.
							if(npp.gt.0)then	
								Af= (npp*p(3)*multtl)+Afromlab	
								npp = npp-(npp*p(3)*multtl)		! allocate to foliage
							else
								Af=Afromlab
							endif
						else
							Af=0
							Afromlab=0
						endif
						
                                                if (npp.gt.0)then					
							Ar= npp*p(4)								! allocate to roots
							Aw=(1-p(4))*npp							! allocate to wood
						else
							Ar=0
							Aw=0
						endif
						
				! litterfall...leaf, wood, roots
				if(multtf.gt.0)then
					!if(lai.le.0.5)then	! leaves just drop if there are few left
							if(Cf .le. 20) then
						Lf(1) =Cf
					else
						Lf(1) = (10**p(5))*Cf*p(14)*multtf
					endif
				endif
						
				Lw = (10**p(6))*Cw
				Lr = (10**p(7))*Cr


                                iLw=iLw+Lw
                                iLr=iLr+Lr
                                iAr=iAr+Ar
                                iAw=iAw+Aw
                                iAf=iAf+Af
				
				!Daily Pools:
				Cf = Cf + Af - Lf(1) - Atolab - ralabtoDaily
				if(Cf.lt.0)then
					Cf=0
				endif	

                                Cw = Cw+Aw-Lw
				if(Cw.lt.0)then
					Cw=0
				endif						
				Cr =Cr+ Ar -Lr-RrootDaily	
				if(Cr.lt.0)then
					Cr=0
				endif						
													
				Clit =Clit + Lf(1) + Lr - RhLitDaily - Dlit
				if(Clit.lt.0)then
					Clit=0
				endif						
								
				CsomPools(1) = CsomPools(1)+ Dlit -D1- Rh1Daily +Lw	
				if(CsomPools(1).lt.0)then
					CsomPools(1)=0
				endif	
					
				CsomPools(2) = CsomPools(2)+D1-D2 - Rh2Daily 	
				if(CsomPools(2).lt.0)then
					CsomPools(2)=0
				endif	
					
				CsomPools(3) = CsomPools(3)+D2 - Rh3Daily 	
				if(CsomPools(3).lt.0)then
					CsomPools(3)=0
				endif
	
				Clab=Clab+Atolab-Afromlab-ralabfromDaily	
				if(Clab.lt.0)then
					Clab=0
				endif						

				Lf(2)=Lf(2)+Lf(1)
				
					
				!Evaluate Daily:
				Csom=sum(CsomPools)
					        ! daily soil water content
					        swc=soilWaterContent(swc,precipDaily,ETdaily(1),swhc,p(17),drainage,runoff) 
					
                                    
					if ((i.ge.stdaycal).and.(i.le.ndaycal)) then				        
				  
					if (laimeas(1).ge.0) then
						err(3,1) = err(3,1)+((laimeas(1)-lai)/laimeas(2))**2
						countD(3)=countD(3)+1
					endif
					
					if (cwmeas.gt.0)  then
						! record first and last to catch trend
						if (cwMeasFirstLast(1).lt.0)then
							cwMeasFirstLast(1)=cwmeas
							cwModFirstLast(1)=cw
						endif
						cwMeasFirstLast(2)=cwmeas
						cwModFirstLast(2)=cw

						! get error for total amount
						err(10,1) = err(10,1)+((cwmeas-cw)/(cwmeasE))**2
						countD(10)=countD(10)+1
							
   						lastvalueCwMeas=cwmeas
						lastvalueCwMod=cw
					endif
						
					! Error for wood carbon increment
                                        if (cwmeasInc.gt.0)  then
					
						if (projectyear .ge. 2008) then
							cwInc = cw - cwPrevious
						else
							cwInc=cw-cwPreviousYear				
						endif
						! get error for total amount
						err(4,1) = err(4,1)+((cwmeasInc -cwInc)/(cwmeasIncE))**2
						countD(4)=countD(4)+1
					endif
					
				
                   ! check tehe Csom sums match total Csom (prevents consistent biases)
                        if (soilTotCmeas.gt.0)then
                        soilTotCmeas=bioData(i,indexBioSoilTotC)
                        soilTotCmeasE=bioData(i,indexBioSoilTotCe)
	                    err(28,1)=err(28,1)+((soilTotCmeas-Csom)/(soilTotCmeasE))**2
						CountD(28)=CountD(28)+1
					    endif
						
						
					! calculate the daily swc errors
                                        !if (swcmeasdaily(1).gt.0) then
                                        !        err(29,1) = err(29,1) +((swcmeasdaily(1) - swcdaily)/swcmeasdaily(2)) ** 2
                                        !        countD(29) = countD(29) + 1
                                        !endif						
						
			
					if (litterfallmeas(1).gt.0) then
							err(5,1) = err(5,1)+((litterfallmeas(1)-Lf(2))/(1.5*litterfallmeas(2)))**2
							countD(5)=countD(5)+1
						endif
						
					!if(decidflag.eq.1)then
						if (((yearday.eq.365.).and.(leafoutmeas(1).gt.0)).or.((i.eq.ndaycal).and.(leafoutmeas(1).gt.0))) then
							err(11,1) = err(11,1)+((leafoutmeas(1)-leafout)/(leafoutmeas(2)))**2
							countD(11)=countD(11)+1
						endif
	
					!	if (((yearday.eq.365.).and.(leafinmeas(1).gt.0)).or.((i.eq.ndaycal).and.((leafinmeas(1).gt.0)))) then
							err(12,1) = err(12,1)+((leafinmeas(1)-leafin)/(leafinmeas(2)))**2
							countD(12)=countD(12)+1
						endif
					!endif
					
					! annual non gap filled  soil Respiration measurements
					if (((yearday.eq.365.).and.(iRsoilmeas.gt.10)).or.((i.eq.ndaycal).and.(iRsoilmeas.gt.10))) then	
						err(8,1) = err(8,1)+((iRsoilmeas-iRsoil(2))/(0.25*(iRsoilmeas)))**2
						CountD(8)=CountD(8)+1
					endif
					
						if (((mod(yearday,30.).eq.0)).and.(iNeemeas(3).ne.0.)) then	! monthly non.gap.filled measurements
							! some sites have very small NEE
							err(17,1) = err(17,1)+((iNeemeas(3)-iNee(3))/(10+0.15*abs(iNeemeas(3))))**2
							CountD(17)=CountD(17)+1
						endif
						
						
					if ((yearday.eq.365.).or.(i.eq.ndaycal)) then
						! carbon in roots relative to estimated initial value							
                                                err(7,1) =err(7,1)+((Cr-P(19))/(P(19)*0.2))**2
						countD(7)=countD(7)+1
							
						! carbon in litter relative to estimated initial value
						err(6,1) =err(6,1)+((Clit-P(21))/(P(21)*0.2))**2
						countD(6)=countD(6)+1
							
							! carbon in litter turnover time
							err(23,1) = err(23,1)+(((Clit/(iRhLit))-4)/3)**2
							countD(23)=countD(23)+1
	
						err(15,1)=err(15,1)+(max(0.5-iRBG/iRe,0.0)+max(iRBG/iRe-0.8,0.0))*100
                        countD(15)=countD(15)+1

		
							! annual non gap filled  (estimated as ~0.25*iNeemeas(2) by Barr et al. for NOAA synth)
							if(NOAA.ne.1)then
								! some sites (PFa) have very low NEE, which can generate large cost when scaling by error	
								err(16,1) = err(16,1)+((iNeemeas(1)-iNee(1))/(40+0.2*(iNeemeas(1))))**2
							else
								if(iNeemeas(2).ne.0)then	! some sites with very small NEE
										! PFa monthly is so small can't set large intercept of 40
										err(16,1) = err(16,1)+((iNeemeas(2)-iNee(2))/(40+0.2*abs(iNeemeas(2))))**2
								else
									CountD(16)=CountD(16)-1
								endif
							endif
							CountD(16)=CountD(16)+1
	
						! Cf reality error check (should be zero at end of year)
						if(decidflag.eq.1)then
							err(18,1) = err(18,1)+Cf
							CountD(18)=CountD(18)+1
						endif
							
							! SOM turnover time errors (these turnover times are mass weighted to get the average for each pool)
							! most values taken from Julia Gaudinski et al. 2001 (?)
								if(numSoilPools.eq.3)then
									! carbon in SOM1 turnover time (this is the Microbial pool)
									err(24,1) = err(24,1)+(((CsomPools(1)/(iRh1+iD1))-1.5)/1.2)**2 !hard coded 1.5 measurement of speed (years) / error of measurement (years)
									countD(24)=countD(24)+1
		
									! carbon in SOM2 turnover time (This is the fast pool)
									err(25,1) = err(25,1)+(((CsomPools(2)/(iRh2+iD2))-58)/30)**2
									countD(25)=countD(25)+1
									
									! carbon in SOM3 turnover time (This is the slow pool)
									err(26,1) = err(26,1)+(((CsomPools(3)/iRh3)-1354)/450)**2
									countD(26)=countD(26)+1
									
								else if(numSoilPools.eq.2)then
									! carbon in SOM1 turnover time (this is the O horizon)
									err(24,1) = err(24,1)+(((CsomPools(1)/iRh1+iD1)-58)/30)**2
									countD(24)=countD(24)+1
		
									! carbon in SOM2 turnover time (this is the A and B horizon)
									err(25,1) = err(25,1)+(((CsomPools(2)/iRh2+iD2)-545)/250)**2
									countD(25)=countD(25)+1
									
								else        ! carbon in SOM1 turnover time
									err(24,1) = err(24,1)+(((CsomPools(1)/iRh1)-473)/235)**2
									countD(24)=countD(24)+1
								endif
							 
							 
						! Lwood error 
						err(20,1)=err(20,1)+((40-iLw)/(0.2*40))**2
						CountD(20)=CountD(20)+1
					endif
						
						
				!endif 
				
				! test reality over whole period (this changed from V4, where reality was tested at end ndaycal
				if(i.eq.nday)then	
					! Carbon Pool Reality check - no pools should be emptying
					if(NOAA.eq.0)then		
						RealityErr =RealityErr+ (&
								&max((P(19))-Cr,0.0)+&		! pools should not end below their starting value
								&max((P(20))-Cw,0.0)+&
								&max((P(21))-Clit,0.0)+&
								&max((P(22))-(CsomPools(1)),0.0)+&
								&max((P(42))-(CsomPools(2)),0.0)+&
								&max((P(43))-(CsomPools(3)),0.0)+&
								&max((P(23))-Clab,0.0)) 
											
						RealityErr = RealityErr - &				! or grow by more than 50% over simulation period
								&min(0.,(2.5*P(19))-Cr)-&
								&min(0.,(1.5*P(20))-Cw)-&
								&min(0.,(2*P(21))-Clit)-&
								&min(0.,(1.5*P(22))-CsomPools(1))-&
								&min(0.,(1.5*P(42))-CsomPools(2))-&
								&min(0.,(1.5*P(43))-CsomPools(3))-&
								&min(0.,(1.5*P(23))-Clab)
					endif
					
					if(NOAA.eq.1)then
						if(decidflag.eq.1)then
							RealityErr =RealityErr+ (&
									&max((0.9*P(19))-Cr,0.0)+&		! pools should not end below their starting value
									&max((P(20))-Cw,0.0)+&
									&max((P(21))-Clit,0.0)+&
									&max((0.95*P(22))-(CsomPools(1)),0.0)+&
									&max((0.95*P(42))-(CsomPools(2)),0.0)+&
									&max((P(43))-(CsomPools(3)),0.0)+&
									&max((P(23))-Clab,0.0)) 
											
							RealityErr = RealityErr - &				! or grow by more than 50% over simulation period
									&min(0.,(2.5*P(19))-Cr)-&
									&min(0.,(2.5*P(20))-Cw)-&
									&min(0.,(2*P(21))-Clit)-&
									&min(0.,(1.5*P(22))-CsomPools(1))-&
									&min(0.,(1.5*P(42))-CsomPools(2))-&
									&min(0.,(1.5*P(43))-CsomPools(3))-&
									&min(0.,(1.1*P(23))-Clab)
					    else
							RealityErr =RealityErr+ (&
									&max((P(19))-Cr,0.0)+&		! pools should not end below their starting value
									&max((P(20))-Cw,0.0)+&
									&max((P(21))-Clit,0.0)+&
									&max((0.95*P(22))-(CsomPools(1)),0.0)+&
									&max((0.95*P(42))-(CsomPools(2)),0.0)+&
									&max((P(43))-(CsomPools(3)),0.0)+&
									&max((0.9*P(12))-Cf,0.0)) 	! 0.7 to allow for some interannual variability
									
							RealityErr = RealityErr - &				! or grow by more than 50% over simulation period
									&min(0.,(1.25*P(19))-Cr)-&
									&min(0.,(1.40*P(20))-Cw)-&
									&min(0.,(1.25*P(21))-Clit)-&
									&min(0.,(1.5*P(22))-CsomPools(1))-&
									&min(0.,(1.25*P(42))-CsomPools(2))-&
									&min(0.,(1.25*P(43))-CsomPools(3))-&
									&min(0.,(1.4*P(12))-Cf)
							        
						endif
					endif
				endif
					
					
					! prevent SOM pools having transient dynamics by testing each year.
					if ((yearday.eq.365.).or.(i.eq.ndaycal)) then
					
								if(numSoilPools.eq.3)then ! 3 pool model
									! 1. The microbial soil carbon pool as 1-3% of the total pool
									err(21,1)=err(21,1)+((0.02*Csom-CsomPools(1))/(0.01*Csom))**2
									CountD(21)=CountD(21)+1
									
									! 2. The slow SOC pool (equivalent to the fast pool when using two pool model)
									err(22,1)=err(22,1)+((p(42)+(20*(countD(22)-8))-CsomPools(2))/(p(42)*0.2))**2 !allowing 20g maximum change of carbon in this pool
									CountD(22)=CountD(22)+1
									! 3. The passive SOC pool should not change greatly over time
									err(27,1)=err(27,1)+((p(43)-CsomPools(3))/(p(43)*0.2))**2
									CountD(27)=CountD(27)+1
									
								else if(numSoilPools.eq.2)then ! 2 pool model
									err(21,1)=err(21,1)+((p(22)-CsomPools(1))/(p(22)*0.2))**2
									CountD(21)=CountD(21)+1
									err(22,1)=err(22,1)+((p(42)+(50*(countD(22)-8))-CsomPools(2))/(p(42)*0.2))**2
									CountD(22)=CountD(22)+1
								else ! 1 pool model
									err(21,1)=err(21,1)+((p(43)-CsomPools(1))/(p(43)*0.2))**2
									CountD(21)=CountD(21)+1
								endif
					endif 
		
				if (cwmeas.gt.-999)  then
					lastvalueCwMeas=cwmeas
					lastvalueCwMod=cw
				endif
						
						DailyAverageData=24*(AutoRsoilmeasDaily(1)/CountN)
						DailyAverageModel=24*(RBGDaily2/CountN)
						DailyAverageTrenchData=24*(TrenchedRsoilmeasDaily(1)/CountN)
						DailyAverageTrenchModel=24*(RhhDaily2/CountN)
						
						
						
				if (iter.ge.endIteration) then  ! save for output after inner loop
					fAparDaily=fAparDaily/radDaily
					pred(InnerRun,i,:) = (/metData(i,1,1),& !1
								yearday,&					!2
								GDaily,&					!3
								DrespDaily,&				!4
								RaDaily,&					!5
								iDlit,&					    !6
								iD1,&					    !7
								iD2,&					    !8
								RrootDaily,&				!9
								RhLitDaily,&				!10
								Rh1Daily,&					!11
								Rh2Daily,&					!12
								Rh3Daily,&					!13
								RhDaily,&                   !14
								lai,&                       !15
								DailyAverageModel,&         !16
								NeeDaily,&					!17
								RhsoilDaily,&				!18
								RsoilModDaily(2),&			!19
								Af,&              	        !20
								AutoRsoilmeasDaily(1),&		!21
								DailyAverageTrenchData,&	!22
								DailyAverageData,&		    !23
								TrenchedRsoilmeasDaily(1),&	!24
								sf,&						!25
								iRBG/iRe,&                  !26
								swcmeasdaily10,&			!27
								swcmeasdaily50,&			!28
								D1,&						!29
								D2,&						!30
								iGPP,&						!31
								iDresp,&					!32
								RhhDaily2,&					!33
								Cf,&						!34
								ETmeas(1),&					!35
								iRa+iRh,&					!36
								iRroot,&					!37
								iRhLit,&					!38
								iRh1,&						!39
								iRh2,&						!40
								iRh3,&						!41
								iRh,&						!42
								npp,&            			!43
								iRBG,&			            !44
								iNee(1),&					!45
								iNee(2),&					!46
								laimeas(1),&				!47
								Cw,&						!48
								Cr,&						!49
								Clab,&						!50
								Clit,&						!51
								Csom,&						!52
								CsomPools(1),&				!53
								CsomPools(2),&				!54
								CsomPools(3),&				!55
								lai,&						!56
								neeDailyMeas,&				!57
								WaterStressP,&				!58
								RBGDaily2,&					!59
								DailyAverageTrenchModel,&	!60
								RhhDaily,&				    !61
								RBGDaily,&				    !62
								iNEEmeas(2),&				!63
								iLw,&				        !64
								lr,&					    !65
								Rhh,&					    !66
								RBG,&				        !67
								AutoRsoilmeas(1),&			!68
								TrenchedRsoilmeas(1),&		!69
								Lf(1)/)						!70
												        	
												        	
													
								

					if(iter.eq.endIteration) then
						innerRunOutDaily=pred
					endif
							
				endif
			
			
				if (litterfallmeas(1).gt.-999) then
					Lf(2)=0
				endif
					
				if((mod(yearday,30.).eq.0)) then		! 30 day Nee constraint reset
					iNee(3) = 0
					iNEEmeas(3)=0
				endif

			 
			END DO timeloop
		
			annualNEE(year)=iNEE(1) ! Save from final year
			annualGPP(year)=iGPP
			annualRa(year)=iRa
			annualRh(year)=iRh
			
			! error 9 is long term increment, no average error, just calculate here.
			err(9,1) =(((cwMeasFirstLast(2)-cwMeasFirstLast(1))&	! Cw  increment
						&-(cwModFirstLast(2)-&
						&cwModFirstLast(1)))/&
						&(0.1*(cwMeasFirstLast(2)-cwMeasFirstLast(1))))**2
			! calculate the total error for the current parameter set
                        ! and data constraints applied

                        err(1,1) = err(1,1)*5
                        err(2,1) = err(2,1)*5

                        toterr=TotalErrorV4bayes(err,RealityErr,rep,numConstraints,constraints,&
                                        &numSoilPools,toterrAll,countD,nparams,boundsP,P)

			if (iter.le.q) then
				de = toterr-besterr
				if (de<0) then 
				      ans=1
				else
				      ans=0
				end if
									
				!Apply the metropolis function....
				!call metrop(de,annealtemp,ans)
				call random_number(ran)
				expo=-de/annealtemp
				if (expo.gt.88) expo=88
				if (ran<exp(expo)) then	
				      ans = 2
				end if
				! end metropolis
								
				if ((toterr.le.allbest).and.(toterr.ge.0)) then
					if ((acceptflag.eq.0))then
						allbest=toterr
						bestp=p
						laststep=iter
					endif
				endif
				
				if (ans.eq.0) then
					stepflag = 0
					if ((iter.gt.(r))) then
						incP=incP*(0.999)
					endif
				endif						
				if ((ans.ge.1).and.(toterr.ge.0)) then
					besterr=toterr
					acc(innerRun)=acc(innerRun)+1
					if ((iter.gt.(r))) then
						oldp=p
						incP=incP*(1.002) !	fror 23% acceptance rate
					endif
					
					if ((ans.eq.1).and.(stepflag.ne.-1)) then
						stepflag = 1 
					else 
						stepflag = 0
					endif
				endif
			endif 

			parameters(1,countTotalIterations,:) = P(:)
			
			error(1,countTotalIterations,1) = toterr
			error(1,countTotalIterations,2) = besterr
			error(1,countTotalIterations,3) = allbest	
		
			countTotalIterations=countTotalIterations+1

			! this prints habitualy to screen / can't do in parallel
                        if (((outerRun.ge.0).or.(iter.eq.q)).and.&
                                &((iter.le.50).or.(mod(int(iter),1000).eq.0).or.&
			        &(iter.eq.r).or.(iter.eq.q-1)).and.(iter.le.q)) then

                        write(*,"(1i5.0,22f10.3,2f10.2)") iter,&
                        		&err(1,1),err(2,1),err(3,1),err(4,1),err(5,1),&
					&err(6,1),err(7,1),err(13,1),err(14,1),&
					&err(15,1),err(16,1),err(19,1),err(23,1),& 
					&err(24,1),err(25,1),err(26,1),err(28,1),&
					&REalityErr,toterr,allbest
					endif
                     		
			if ((iter.eq.q).or.(iter.eq.r))then
				write(*,*) "iter.  neeDay  neeNight   Lai  CwIncr  Lf  &
						&LitC  RootC   TrenchR AutoR  iRBG/iRe iNEE &
						&ET litC_TO   SOM1C_TO  SOM2C_TO  SOM3C_TO &  
						&TOT_SOC RyltyE  totE   allBest"
						
			endif
	
			floatingIter = iter+0.01
			floatingAcc = acc(innerRun)+0.01
			floatingAns = ans+0.01
			floatingStepflag = Stepflag+0.01

			! here we save the parameters of the first and second exploration to file ....
			! and then print all innerRun parameters after the inner loop has finished
						
			if (iter.eq.r) then
				innerRunParameters(1,innerRun,:) = bestP(:)
				innerRunIncP(1,innerRun,:) = incP(:)
				innerRunBoundsP(1,innerRun,:,1) =boundsP(:,1)
				innerRunBoundsP(1,innerRun,:,2) = boundsP(:,2)
			end if
			if (iter.eq.endIteration) then
				innerRunParameters(2,innerRun,:) = bestP(:)
				innerRunIncP(2,innerRun,:) = incP(:)
				innerRunBoundsP(2,innerRun,:,1) =boundsP(:,1)
				innerRunBoundsP(2,innerRun,:,2) = boundsP(:,2)
				acc(innerRun) = 0
			        stepflag = 0
				
                                prange(:,1)=P(:)
				prange(:,2)=P(:)
				
								
								if(outerRun.eq.numOuterRuns)then
									qerr(:)=max(0.0001,err(:,1))
									! for data streams with just one meas, normalization does not apply
									if (qerr(6).lt.1)qerr(6)=1
									if (qerr(7).lt.1)qerr(7)=1
                                    if (qerr(8).lt.1)qerr(8)=1		! added for cumRsoil constraint to keep consistent with data
									if (qerr(9).lt.1)qerr(9)=1
									if (qerr(10).lt.1)qerr(10)=1
                                    if (qerr(15).lt.1)qerr(15)=1
                                    if (qerr(16).lt.1)qerr(16)=1	! added for annual NEE constraint (few data points, can get if constrained to iNEE alone)
                                    if (qerr(20).lt.1)qerr(20)=1
									if (qerr(21).lt.1)qerr(21)=1
									if (qerr(22).lt.1)qerr(22)=1
									if (qerr(23).lt.1)qerr(23)=1
									if (qerr(24).lt.1)qerr(24)=1
									if (qerr(25).lt.1)qerr(25)=1
									if (qerr(26).lt.1)qerr(26)=1
									if (qerr(27).lt.1)qerr(27)=1
									if (qerr(28).lt.1)qerr(28)=1
									
									err(:,2)=err(:,1)/qerr(:)
								endif
							endif
	
        		if (iter.ge.endIteration) then
				! reset for Chi/sqr testing....
				err(:,2)=err(:,1)/qerr(:)
			endif

                        ! entering posterior parameters
			if (iter.ge.endIteration) then
			        if ((acc(innerRun).ge.explore).or.(outerRun.lt.numOuterRuns).or.(iter.ge.(endIteration+exploreLength))) then 
                                        ! if we have gotten to the end, and found the parameters, and then we exit...
					exit							
                                        ! if we don't accept explore parameters within 10000 runs then we leave anyway. 
				endif
				
				flag=0
				countExplore=countExplore+1
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! Now in SEARCH STAGE
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				
                                ! check if the current parameter set is acceptable based on chi-squared test
				flag=posteriorChiSqTest(err,forest,rep,decidflag,RealityErr,numConstraints,constraints,countD)
                                if (flag.eq.1)then
					acc(innerRun) = acc(innerRun) + 1
					accparm(acc(innerRun),innerRun,1:nparams)=p(1:nparams)
                                        accparm(acc(innerRun),innerRun,nparams+1)=toterrAll

					! adaptive step size algorithm for search stage
					! my own invention - so step wisely
					! will adjust step sizes dynamically to get a 25% (or whatever specified) acceptance rate
					! this lets the parameter range expand beyond the priors in the posterior exploration 
				
					tmp=acc(innerRun)
					tmp2=countExplore
					tmp=tmp/tmp2
					
					if (tmp.gt.acceptanceRate)then	! if more than x% of parameters are getting accepted, increase step size
						stepSize=stepSize*1.01
					else	! otherwise, decrease step size to increase acceptance rate
						stepSize=stepSize*0.99
					endif

					!auto scaling of prange
					prange(:,1)=min(prange(:,1),p(:))
					prange(:,2)=max(prange(:,2),p(:))

					! auto scaling in IncP
					incP=max((prange(:,2)-prange(:,1))*stepSize(:),abs(bestP(:))*0.00001)
					
					if ((acc(innerRun).le.100.0).and.(firstexplore.lt.1)) then
						acc(innerRun) = 1
						firstexplore=firstexplore+1
					endif
						
					if(mod(acc(innerRun),100).eq.0)then	
						write(*,*) "*****************Increasing Range************************"
						write(*,*)acc(innerRun), iter, (prange(:,2)-prange(:,1)),toterrAll
			
					endif
					
					cl(:,:,1)=max(cl(:,:,1),Pred(InnerRun,:,:))
					cl(:,:,2)=min(cl(:,:,2),Pred(InnerRun,:,:))

					cl_sub(:,:,:,1)=max(cl_sub(:,:,:,1),innerRunOutSubDaily(InnerRun,:,:,:))
					cl_sub(:,:,:,2)=min(cl_sub(:,:,:,2),innerRunOutSubDaily(InnerRun,:,:,:))

					
					posteriorFluxComponents(acc(innerRun),1:nyears)=annualNEE
					posteriorFluxComponents(acc(innerRun),(nyears+1):nyears*2)=annualGPP
					posteriorFluxComponents(acc(innerRun),(nyears*2+1):nyears*3)=annualRa
					posteriorFluxComponents(acc(innerRun),(nyears*3+1):nyears*4)=annualRh
					
						do i =1,numConstraints
							! get average error for each data stream
				       		!if(err(i,1).gt.0)err(i,1)=err(i,1)/countD(i) 
				       		posteriorErrorTerms(acc(innerRun),i)=err(i,1) 
							
						enddo
					
					clHigh(innerRun,:,:)=cl(:,:,1)
					clLow(innerRun,:,:) = cl(:,:,2)

					clHigh_sub(innerRun,:,:,:)=cl_sub(:,:,:,1)
					clLow_sub(innerRun,:,:,:) = cl_sub(:,:,:,2)					
					
					oldp = p
				
				else
				! parameter set not accepted - decrease step size
					incP=incP*0.99
						
					tmp=acc(innerRun)
					tmp2=countExplore
					tmp=tmp/tmp2
					
					if (tmp.lt.0.225)then
						stepSize=stepSize*0.99
					endif
				endif
				
			end if	
		END DO mcloop
	End DO

	
	!***********************************************************************
	!! here need to choose best parameter set and use for next outer run
	!***********************************************************************
	bestInnerRun=1

	close(30)				
	bestinnerRun=1
	bestp=innerRunParameters(2,bestInnerRun,:)
	P = bestP
	incP = innerRunIncP(2,bestInnerRun,:)
	boundsP(:,1)=innerRunBoundsP(2,bestInnerRun,:,1) 
	boundsP(:,2)=innerRunBoundsP(2,bestInnerRun,:,2) 

	print *, "Finished ...."
	print *, 'Site Name: ', forest
	print *, 'Decid (1)/Evergreen (0): ', decidflag
	print *, 'Rep: ', rep

	! print final output corresponding to best parameter set
	if(outerRun.eq.numOuterRuns)then

		write (filenum (1), *) decidflag
		write (filenum (2), *) rep
		write (filenum (3), *) jfunct
		write (filenum (4), *) numSoilPools
	
	
		filename =home//"/Results/"//"_outputDaily_"//&
			&trim(adjustl(filenum(1)))//"_"//&
			&trim(adjustl(filenum(2)))//"_"//trim(adjustl(filenum(3)))//"_"//trim(adjustl(filenum(4)))//".csv"
		do i=1,nday
			if (i.eq.1) then !write headers on first day step of last run
				open(unit=26,file=filename,status='unknown')	! outputs
			endif
			write(26,"(70f10.3)")(innerRunOutDaily(bestInnerRun,i,j), j=1,numDailyVariablesOut)
			
			if (i.eq.nday) then
				close(26)
			endif
		End Do			
	
		! print final output corresponding to best parameter set
		filename = home//"/Results/"//"_outputSubDaily_"//&
			&trim(adjustl(filenum(1)))//"_"//trim(adjustl(filenum(2)))//"_"//&
			&trim(adjustl(filenum(3)))//"_"//trim(adjustl(filenum(4)))//".csv"
		do i=1,nday
			if (i.eq.1) then !write headers on first day step of last run
				open(unit=26,file=filename,status='unknown')	! outputs
			endif
			do j =1,subDaily		
				write(26,"(35f10.3)")(innerRunOutSubDaily(bestInnerRun,i,j,k), k=1,numSubDailyVariablesOut)
			end do	
			if (i.eq.nday) then
				close(26)
			endif
		End Do			
	
		! print parameters used in each iteration of the MC3 simulation	
		filename1 = home//"/Results/"//"parameterEvolution_"//&
			&trim(adjustl(filenum(1)))//"_"//trim(adjustl(filenum(2)))//"_"&
			&//trim(adjustl(filenum(3)))//"_"//trim(adjustl(filenum(4)))//".csv"
		do i=1,(q*(numOuterRuns+1)+exploreLength)
			if (i.eq.1) then !write headers on first day step of last run
				open(unit=26,file=filename1,status='unknown')	! parameters
			endif
		End Do		
	
		close(26)
		close(27)
		close(28)
		close(29)
		close(30)
		close(31)

		filename = home//"/Results/"//"AccParams_"//trim(adjustl(filenum(1)))//"_"//&
			&trim(adjustl(filenum(2)))//"_"//trim(adjustl(filenum(3)))//"_"//trim(adjustl(filenum(4)))//".txt"
		open(unit=30,file=filename)
		!write(30,*) "P1 P2 P3 P4 P5 P6 P7 P8 P9 P10 P11 P12 P13 P14 P15 P16 P17 P18 P19 P20 P21 P22 &
		!	&P23 P24 P25 P26 P27 P28 P29 P30 P31 P32 P33 P34 P35 P36 P37 P38 P39 P40 P41 P42"
		do i=1,acc(bestInnerRun)
			write(30,"(60f15.3)") accparm(i,bestInnerRun,1:nparams)
		end do
		close(30)

        ! print accepted toterr on it's own
		filename = home//"/Results/"//"TotalError_"//trim(adjustl(filenum(1)))//"_"//&
			&trim(adjustl(filenum(2)))//"_"//trim(adjustl(filenum(3)))//"_"//trim(adjustl(filenum(4)))//".txt"
		open(unit=30,file=filename)
		write(30,*) "TotErrAll"
		do i=1,acc(bestInnerRun)
			write(30,*) accparm(i,bestInnerRun,nparams+1)
		end do
		close(30)

		filename = home//"/Results/"//"FoBAAR_UCL_"//trim(adjustl(filenum(1)))//"_"//&
			&trim(adjustl(filenum(2)))//"_"//trim(adjustl(filenum(3)))//"_"//trim(adjustl(filenum(4)))//".txt"
		open(unit=32,file=filename)	
		write (32,*) "Random header"
		do j=1,nday
			write (32,"(70f15.3)") clHigh(bestInnerRun,j,:)
		enddo
		close(32)

		
		filename = home//"/Results/"//"FoBAAR_UCL_SubDaily"//trim(adjustl(filenum(1)))//"_"//&
			&trim(adjustl(filenum(2)))//"_"//trim(adjustl(filenum(3)))//"_"//trim(adjustl(filenum(4)))//".txt"
		open(unit=32,file=filename)	
		write (32,*) "Random header"
		do j=1,nday
		        do jj=1,subDaily
			        write (32,"(35f15.3)") clHigh_sub(bestInnerRun,j,jj,:)
		        enddo
		enddo
		close(32)		
		
		filename = home//"/Results/"//"FoBAAR_LCL_"//trim(adjustl(filenum(1)))//"_"//&
			&trim(adjustl(filenum(2)))//"_"//trim(adjustl(filenum(3)))//"_"//trim(adjustl(filenum(4)))//".txt"
		open(unit=32,file=filename)
		write (32,*) "Random header"
		do j=1,nday
			write (32,"(70f15.3)") clLow(bestInnerRun,j,:)
		enddo
		close(32)

		
		filename = home//"/Results/"//"FoBAAR_LCL_SubDaily"//trim(adjustl(filenum(1)))//"_"//&
			&trim(adjustl(filenum(2)))//"_"//trim(adjustl(filenum(3)))//"_"//trim(adjustl(filenum(4)))//".txt"
		open(unit=32,file=filename)	
		write (32,*) "Random header"
		do j=1,nday
		        do jj=1,subDaily
			        write (32,"(35f15.3)") clLow_sub(bestInnerRun,j,jj,:)
		        enddo
		enddo
		close(32)		
		
		filename = home//"/Results/"//"BestParams_"//trim(adjustl(filenum(1)))//"_"//&
	        &trim(adjustl(filenum(2)))//"_"//trim(adjustl(filenum(3)))//"_"//trim(adjustl(filenum(4)))//".txt"
		open(unit=30,file=filename)
		write (30,*) (bestP(i),i=1,nparams)
		close(30)
	
		filename = home//"/Results/"//"PostFluxComponents_"//trim(adjustl(filenum(1)))//"_"//&
	        &trim(adjustl(filenum(2)))//"_"//trim(adjustl(filenum(3)))//"_"//trim(adjustl(filenum(4)))//".txt"
		open(unit=30,file=filename)
		do j=1,explore
			write (30,"(5f15.3)") (posteriorFluxComponents(j,i),i=1,nyears*4)
		enddo
		close(30)	
		
		filename = home//"/Results/"//"PostErrorTerms_"//trim(adjustl(filenum(1)))//"_"//&
	        &trim(adjustl(filenum(2)))//"_"//trim(adjustl(filenum(3)))//"_"//trim(adjustl(filenum(4)))//".txt"
		open(unit=30,file=filename)
		do j=1,explore
			write (30,"(29f15.3)") (posteriorErrorTerms(j,i),i=1,numConstraints)
		enddo
		close(30)
		
	
	endif

	close(33)

	if (decidflag.eq.1) then
		open(unit=33,file=home//'/Results/'//'initialRun_D_'//trim(adjustl(repChar))//'.csv',&
			&status='unknown')		
	else
		open(unit=33,file=home//'/Results/'//'initialRun_E_'//trim(adjustl(repChar))//'.csv',&
			&status='unknown')		
	endif

	write (33,*) (bestP(i),i=1,nparams)
	close(33)

ENDDO outer

close(26)

!print execution time.
call ETIME(tarray, result) 
 print *, "Execution time (seconds)...."
              print *, tarray(1)	!user time in seconds
              print *, tarray(2) !system time in seconds     


END


!*************************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*************************************************************************


FUNCTION rnorm() RESULT( fn_val )

!   Generate a random normal deviate using the polar method.
!   Reference: Marsaglia,G. & Bray,T.A. 'A convenient method for generating
!              normal variables', Siam Rev., vol.6, 260-264, 1964.

IMPLICIT NONE
REAL  :: fn_val

! Local variables

REAL            :: u, sumx
REAL, SAVE      :: v, sln
LOGICAL, SAVE   :: second = .FALSE.
REAL, PARAMETER :: one = 1.0, vsmall = TINY( one )

IF (second) THEN
! If second, use the second random number generated on last call

  second = .false.
  fn_val = v*sln

ELSE
! First call; generate a pair of random normals

  second = .true.
  DO
    CALL RANDOM_NUMBER( u )
    CALL RANDOM_NUMBER( v )
    u = SCALE( u, 1 ) - one
    v = SCALE( v, 1 ) - one
    sumx = u*u + v*v + vsmall         ! vsmall added to prevent LOG(zero) / zero
    IF(sumx < one) EXIT
  END DO
  sln = SQRT(- SCALE( LOG(sumx), 1 ) / sumx)
  fn_val = u*sln
END IF
RETURN
END FUNCTION rnorm


! --------------------------------------------------------------------------------------------------


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!! Read the driver data from file
subroutine readObservedData(bioData,maxt,mint,fluxData,metData,Dayflag,soilTemp,longTerm)
! this routine will read in the driver data and calculate the soil temperature

use obsdrivers

implicit none

real :: fluxData(nday,subDaily,numColumnsFlux)
real :: metData(nday,subDaily,numColumnsMet)
real :: bioData(nday,numColumnsBio)
real :: maxt(nday),mint(nday)
real, intent(in) :: longTerm
real :: dayflag(nday,subDaily)
real :: soilTemp(nday,subDaily)

integer :: i,j,k

! DRIVERS

! Read in the met data
open(unit=26,file=home//'/Data/'//'/HOmetdata.csv',status='old')
if (longTerm.eq.0)then
	! Read in the flux data
	open(unit=25,file=home//'/Data/'//'/HOfluxdata.csv',status='old')	
	! Read in the biometric data
	open(unit=29,file=home//'/Data/'//'/HObiodata.csv',status='old')
endif

bioData = -999
maxt = -999
mint = 999
  
! READ IN OBSERVATIONS AND DRIVERS

DO i=1,nday
	DO j = 1,subDaily
		Read(26,*)(metData(i,j,k),k=1,numColumnsMet) 
		if(metData(i,j,indexMETpar).gt.0)Dayflag(i,j)=1		! if filled PAR > 0, set dayflag to day 

                if (longTerm.eq.0.0)then
			Read(25,*)(fluxData(i,j,k),k=1,numColumnsFlux)
		endif		
		! calculate daily max and min temperatures
		if ((metData(i,j,indexMETtemp).lt.mint(i)).and.(metData(i,j,indexMETtemp).gt.-100)) then
			mint(i) = metData(i,j,indexMETtemp)
		endif
		
		if (metData(i,j,indexMETtemp).gt.maxt(i).and.(metData(i,j,indexMETtemp).gt.-100)) then
			maxt(i) = metData(i,j,indexMETtemp)
		endif
				
	END DO
                if (longTerm.eq.0)then
	
			Read(29,*)(bioData(i,k),k=1,numColumnsBio)
		endif		
	
END DO

close(26)

	if (longTerm.eq.0)then
		close(25)
		close(29)
	endif

end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




