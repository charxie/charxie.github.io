   	program ET
   	
c********************************************************************
c * Predictor-corrector method                                      *
c * Floquet space method					    *
c********************************************************************   	
   	
	implicit real*8 (a-h,o-z)
	parameter (n=3,nflo=n*81,ncut=100)
	parameter (neigd=2000)
	dimension perio(n),omega(n),ampli(n),delta(n)
	dimension perhop(n,n),omehop(n,n),amphop(n,n),delhop(n,n)
	dimension hamf(n,n,0:ncut)
	logical restart,static,floquet,dynamics,model
	logical scanfreq,scanener
        character yesorno*3
	common /eigv/ omcm(neigd),zr(neigd,neigd),zi(neigd,neigd),
     :                ne,nblw

	data VD1/1.0d0/
	data VBB/1.0d0/
	data VNA/1.0d0/
	data ED /0.0d0/
	data EA /0.0d0/
	data EB /5.0d0/
	data tau/0.001d0/
	data nstep/100000/
	data iseed/1111/
	data dperiod/0.1/
	data denergy/0.2/

c	call rnum(iseed)
	
	twopi=8.d0*datan(1.d0)
	model=.true.
	scanener=.false.
	scanfreq=.not.scanener
	restart=.false.
	static=.false.
	floquet=.false.
	dynamics=.not.floquet
	ncase=1
 
885	write(*,*) ' restart ? (y/n) '
	read(*,'(a)') yesorno
	if(yesorno.eq.'y  ') then
	   restart=.true.
	elseif(yesorno.eq.'n  ') then
	   restart=.false.
	else
	   write(*,*) ' answer yes or no '
	   goto 885
	endif
	
	open(14,file='omax')
	open(16,file='quae')
	
	if(model) call sta(tau,n,ED,EA,EB,VD1,VBB,VNA)	

881	write(*,*) ' stop ?  (y/n)'
	read(*,'(a)') yesorno
	if(yesorno.eq.'y  ') then
	   stop
	elseif(yesorno.eq.'n  ') then
	   goto 880
	else
	   write(*,*) ' answer yes or no '
	   goto 881
	endif

880     continue

        incromega=1

930	if(model) call defdynpar(tau,n,perio,ampli,delta,omega,
     :	               perhop,amphop,delhop,omehop,    	
     :                 ncase,incromega,dperiod)

	nstep1=2.0*incromega*dperiod/tau
	tau1=tau

940	if(model.and.floquet) call htime(tau1,n,nstep1,
     :                         ED,EA,EB,VD1,VBB,VNA,
     :                         ampli,omega,delta,
     :                         amphop,omehop,delhop)
     
        nstep2=200
     
        if(.not.model.and.floquet) then
            call hmod(n,nstep2)
            nstep1=nstep2
	endif

	if(floquet) call fourier(n,tau1,nstep1,ncut,hamf)
	
	if(floquet) call flo(n,nflo,ncut,hamf,nstep1,tau1,omax)

	if(model.and.dynamics) call dyn(tau,n,nstep,restart,static,
     :                        ED,EA,EB,VD1,VBB,VNA,
     :	                      ampli,omega,delta,
     :                        amphop,omehop,delhop,
     :                        omax)

        if(scanfreq)write(* ,'(2f10.5)') incromega*dperiod,omax
        if(scanfreq)write(14,'(2f10.5)') incromega*dperiod,omax
        if(scanener)write(* ,'(2f10.5)') EB,omax
        if(scanener)write(14,'(2f10.5)') EB,omax
        if(floquet.and.scanfreq)  
     :     write(16,'(100e10.3)') incromega*dperiod,
     :                            (omcm(i),i=nflo/2-10,nflo/2+10)


	if(.not.model) stop

 	incromega=incromega+1
 	
	if(scanener.and.incromega.le.100) then
	   EB=EB+denergy
	   goto 940
	endif

 	if(scanfreq.and.incromega.le.100) goto 930

 	stop
	end