 	  program article_fkopp 
	  implicit none 

      double precision alfas,W,rp2,m,S2,yi,yf
c==========PT====================================
c for the LHC , Psi

      alfas=0.25d0
      W=14000d0

c     particle properties     
      m=3.096d0
      rp2=0.81d0
     
      S2=0.029d0
c     rapitity limits of integration
      yi=-6.5d0
      yf=6.5d0
     
c    subroutine for the calculation of dsigma/dpt
        call dsigmadpttt(alfas,W,m,rp2,S2)

c    subroutine for the calculation of integral in pt     
       call integraldadmul(m,alfas,rp2,W,s2,yi,yf)

c============= Y ===================================
c    subroutine for the calculation of dsigma/dy     

c      call dsigdy_xcxb()
      
c    subroutine for the integration in y
c      call integration_xcxb()
      end program article_fkopp 

  	 include 'quadpack.f'

	
         subroutine integration_xcxb()
         integer par, det
         double precision mass1,rp21,alpha_s1,ecm,
     &   a1,b1,s2med1 
         print*, "Particles: xc0 (1) or xb0(2) "
         read(*,*)  par 
         Print*, "Detectors: Tevatron(1) or LHC(2)"
         read(*,*) det 
    
        print*, "----Calling integration routine----"
        call integration(par, det)
        print*, "----integration routine ended----------"
        print*, "                          "
        print*, "                          "
        print*, "                          "
	     
	     RETURN
	     end
	     
	     subroutine dsigdy_xcxb()
         integer par, det
         print*, "Particles: xc0 (1) or xb0(2) "
         read(*,*)  par 
         Print*, "Detectors: Tevatron(1) or LHC(2)"
         read(*,*) det 
    
        print*, "----Calling rapidity routine-----"
        call rapitity(par, det)
        print*, "----rapidity routine ended----"
        print*, "                          "
        print*, "                          "

	     
	     RETURN
	     end
	     
	     
	     
	     subroutine dsigmadpttt(alphas,ecm,m,rp2,s2)
	     implicit none 
         double precision y,yf,alphas,rp2,ecm,pt
         double precision dsigmadpt,res,m,s2
         integer p
         y=0.0d0
         yf=y
 	 print*, "==================================="
	 print*, "     Pt[GeV]   " ,  "  NB/GeV"
         do p=6,100,1
         pt=0.1*p
         res=dsigmadpt(y,yf,m,alphas,rp2,ecm,pt)*s2
	   
         print*, pt , res 
         end do 
         print*, "==============================="

         return 
         end 
            
	     
	     
	    subroutine integraldadmul(m1,alphas1,rp21,ecm1,s2,a1,b1)
	    implicit none
	    EXTERNAL WF,F
        double precision WF,F,EPS,RESULT
	    integer N
	    double precision A(2),B(2),a1,b1
	    double precision m1,alphas1,rp21,ecm1
	    double precision m,alphas,rp2,ecm,s2

	    common/var/m,alphas,rp2,ecm
	    
c	    N-> numero de dimensioes
c	    A(1)= limite inferior de integracao - variavel x
c	    B(1)= limite superior de integracao - variavel x
c           EPS => precisao

        m=m1 
        alphas=alphas1
        rp2=rp21
        ecm=ecm1


	    EPS=1.0D-6
	    
	    N=2
	    A(1)=a1
	    A(2)=0.0d0
	    
	    B(1)=b1
        B(2)=100.0d0
        CALL DADMUL(WF,N,A,B,EPS,RESULT)
	    print*, "==================================="
	    WRITE(*,*) "Total sigma" , RESULT*s2  , "NB"
	    print*, "==================================="
	    return
	    end 

	    
	    
      


        subroutine integration(particle, detec)
        external SGS1,F1,aux_int
        double precision dsigmady,y,alfas_s,alfas_ll,mass,b1,alpha_s1
        double precision rp2,s2med,SGS1,F1,aux_int,resint,w,rg,alpha_s
        double precision mass1,rp21,w1,alfas_s1,s2med1,rg1,a1,ecm
        integer i,particle, detec
        common/varint/mass,rp2,w,alfas_s,s2med,rg
        common/var_alfas/alpha_s
           
        
        

        if (particle .eq. 1) then
c	xc0
        mass=3.414d0
        rp2=0.075d0/1.45d0
        else 
c	xb0
        mass=9.859d0
        rp2=1.42d0
        endif


        if(detec .eq. 1)then
c	Tevatron
        w=1960.0d0
        a1=4.0d0
        
        s2med=0.058d0 
        rg=1.0d0
        
        else 
c	LHC
        w=14000.0d0
        s2med=0.028d0 
	    rg=1.0d0
	    a1=8.0d0
        end if


          if(mass .eq. 3.414d0 .or. mass.eq. 3.096d0)then
	     alpha_s1=0.335d0
	     else
	     alpha_s1=0.25d0
	     end if
        
        
        alpha_s=alpha_s1
        alfas_s=alpha_s
        
         
        
        
        
 
        resint=SGS1(-a1,a1,1d-3,aux_int)
        print*, resint , "nb"
        return 
        end 






        double precision function aux_int(y)
        implicit none 
        double precision mass,rp2,y,s2med,rg
        double precision w,alfas_s,dsigmady
        common/varint/mass,rp2,w,alfas_s,s2med,rg

        aux_int=dsigmady(mass,rp2,y,w,alfas_s,
     &rg)*s2med
        return 
        end     

        subroutine rapitity(particle, detec)
        implicit none 
        double precision dsigmady,y,alfas_s,alpha_s1,mass,ecm
        double precision rp2,s2med,a,alpha_s,ans,rg,w
        integer particle,detec, i
        common/var_alfas/alpha_s1

        open (unit=10,file="gbw_xc0_y_lhc.dat",
     &status="unknown")

      


        if (particle .eq. 1) then
c	xc0
        mass=3.414d0
        rp2=0.075d0/1.45d0
        else 
c	xb0
        mass=9.859d0
        rp2=1.42d0
        endif


        if(detec .eq. 1)then
c	Tevatron
        w=1960.0d0

        
        s2med=0.058d0 
        rg=1.0d0
        
        else 
c	LHC
        w=14000.0d0
        s2med=0.028d0 
	    rg=1.0d0
        end if


          if(mass .eq. 3.414d0 .or. mass.eq. 3.096d0)then
	     alpha_s=0.335d0
	     else
	     alpha_s=0.25d0
	     end if
        
        
        alpha_s1=alpha_s

        do i=1,900,10
        y=0.01*i
        ans=dsigmady(mass,rp2,y,w,alpha_s1,
     &rg)*s2med
        
        if(ans .gt. 1e-10)then
c        if(y .lt. 0.1d0)then
  
        print*,y, ans  
        write(10,*) y , ans
        else
        goto 50
        end if
        end do
   50 continue 
        close(10)
        return 
        end

        double precision function dsigmady(M1,R1,y1,s12,alfas_s1,RG1)
        implicit none 
        external SGS1,fauxdpe
        double precision M1,R1,y1,s12
        double precision b,pi,alfas,SGS1,fauxdpe
        double precision const,res,RG1,alfas_s1
        double precision RG,alfas_s
        double precision m,y,s
       
        integer limit
        parameter ( limit = 1000 )
        integer lenw
        parameter ( lenw = limit * 4 )
        double precision abserr
        double precision epsabs
        double precision epsrel
        double precision testf
        integer ier
        integer iwork(limit)
        integer key
        integer last
        integer neval
        double precision result
        double precision true,alpha_s
        double precision work(lenw)
	   common/var1/m,y
	   common/var2/s,alfas
       common/var3/RG
       common/var_alfas/alpha_s

       alfas=alfas_s1
       alpha_s=alfas
	   pi=3.1415926535d0
       RG=RG1
	   m=M1
	   y=y1
	   s=s12**2.0d0	
	   b=4.0d0


	   const=(1.0*pi**4.0*M1*R1*0.389d0*alfas**2)/ 
     &(b**2.0)



        call dqag ( fauxdpe,0.0d0,1.0d0,1d-3,1d-5,5, result, abserr, 
     &  neval, ier, limit, lenw, last, iwork, work )
	
	     res=const*(result**2.0)
	     dsigmady=res*1d6*RG**4.0
	     return 
	     end

	




	     double precision function fauxdpe(qta)
	     implicit none 
	     double precision qt2,x1,x2,y,m,s,qt,qta
	     double precision fg_gbw,cut,alfas,alpha_s
	     double precision a1,b1,c1,jac,t,fgtmd
	     common/var1/m,y
	     common/var2/s,alfas
         common/var_alfas/alpha_s

         alfas=alpha_s
	     qt=(1.0-qta)/qta
         jac=(1.0/qta**2.0)
           
         qt2=qt*qt
	     x1=(m/sqrt(s))*exp(y)
	     x2=(m/sqrt(s))*exp(-y)
	     a1=(1.0/qt2)**2.0

         b1=fg_gbw(x1,qt2,m)*fg_gbw(x2,qt2,m)

	     c1=(2.0*qt)/((m**2.0 + 1.0*qt2)**2.0)
	     fauxdpe=a1*b1*c1*jac
          
	     return 
	     end 


	   

	     double precision function alfas_ll(q2)
	     implicit none 
	     double precision q2,pi,nf,b0,b1,t
	     double precision lambdaqcd2,c1,c2,ms2
	     pi=3.1415926535d0
	     lambdaqcd2=0.25d0**2.0d0
c	     lambdaqcd2=0.35d0**2.0d0
	     nf=4.0d0
	     t=log(q2/lambdaqcd2)
	     b0=(33.0d0 - 2*nf)/(12.0*pi)
	     c1=1.0d0/(b0*t)
	     alfas_ll=c1
	     return 
	     end 


          double precision function fg_gbw(x,kt2,m)
          implicit none 
	     double precision x,kt2,sigma0,x0,lamb,pi,ro2,m
	     double precision alpha_s
	     common/var_alfas/alpha_s
	     pi=3.141592d0
c	     sigma0=29.12/0.389d0
c	     x0=0.41d-4
c	     lamb=0.277d0
	     
	     
         sigma0=27.32/0.389d0
	     x0=4.2d-5
	     lamb=0.248d0
	     ro2=(x/x0)**(lamb)
	     
	     if(m .eq. 3.414d0 .or. m.eq. 3.096d0)then
	     alpha_s=0.335d0
	     else
	     alpha_s=0.25d0
	     end if

c         alpha_s ja incluso em termos da escala de massa 
c         multiplicamos a updf por kt2 para torna-la adimensional
c         definicao-> fg_gbw= (3.0*sigma0*ro2*kt2**2.0*exp(-ro2*kt2))/(4.0*pi**2.0*alfas)

	     fg_gbw=((3.0*sigma0*ro2*kt2**2.0*exp(-ro2*kt2))/
     &(4.0*pi**2.0*alpha_s))*(1.0-x)**5.0
	     return 
	     end 
	     
	
      
      
      
       DOUBLE PRECISION FUNCTION WF(N,XX)
       IMPLICIT DOUBLE PRECISION (A-H,J-Z)    
	   INTEGER N
	   DOUBLE PRECISION XX,yf,pt,rp2,ecm,alphas,m
	   DIMENSION XX(2)
	   common/var/m,alphas,rp2,ecm
       N=2	
	   yf=XX(1)
	   pt=XX(2)

	   WF=dsigmadpt(0.0d0,yf,m,alphas,rp2,ecm,pt)
	   RETURN 
	   END
	     
	     
	     
	     
	    double precision function dsigmadpt(y,yf,m1,alphas,rp2,ecm,pt)
	    implicit none 
	    external aux_intdp
	    double precision x1,x2,mt,pt,m,ec2,ae,a3,ay,b,pi,n0     
	    double precision alphas ,aux_intdp,m1,y,yf,rp2,ecm
        integer limit
        parameter ( limit = 1000 )
        integer lenw
        parameter ( lenw = limit * 4 )
        double precision abserr
        double precision epsabs
        double precision epsrel
        double precision testf
        integer ier
        integer iwork(limit)
        integer key
        integer last
        integer neval
        double precision result
        double precision true,alfas1
        double precision work(lenw)
        common/var4/x1,x2,m
        
        m=m1
        mt=sqrt(m**2.0 + pt**2)
        x1=(mt/ecm)*exp(y) + (pt/ecm)*exp(yf)
        x2=(mt/ecm)*exp(-y) + (pt/ecm)*exp(-yf)
        
        a3=(pt/ecm)*exp(yf)
        ay=(mt/ecm)*exp(y)
        
        b=4.0d0
        ec2=4.0d0/9.0d0
        ae=1.0d0/137.0d0
        pi=3.141592d0
         
        call dqag ( aux_intdp,0.0d0,1.0d0,1d-3,1d-5,5, result, abserr, 
     &  neval, ier, limit, lenw, last, iwork, work )
	
        n0=(2.0*pi**2.0*alphas**2.0*ae*ec2*a3**2.0*ay**2.0*m1*rp2)/
     &(3.0*b**2.0*x1**2*(a3*mt**2.0 + ay*pt**2.0)**2.0*mt**4.0)   
	
	    dsigmadpt=2.0*pt*pi*n0*result**2.0*0.389d6
	    return 
	    end
	     
	     
	     
	     
	     
	     double precision function aux_intdp(qta)
	     implicit none 
	     double precision qta,qt,pi,a1,b1,jac,fg_gbw
	     double precision x1,x2,m,qt2
	     common/var4/x1,x2,m
          
         pi=3.141592d0
         qt=(1.0-qta)/qta
         jac=(1.0/qta**2.0)
         qt2=qt*qt
	     a1=(1.0/qt2)**2.0
         b1=fg_gbw(x1,qt2,m)*fg_gbw(x2,qt2,m)
         
	     aux_intdp=2.0*qt*jac*a1*b1

         return 
         end


         
         
         
         
CCCCCCCCCCCCCCC Rotina de integracao DADMUL (CERN MathLib) CCCCCCCCCCCCC


      SUBROUTINE DADMUL(F,N,A,B,EPS,RESULT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      double precision RELERR
      INTEGER N,NFNEVL,MINPTS/1.D3/,MAXPTS/1.D6/,IWK/11.D5/
      integer IFAIL
      DIMENSION WK(1100000)
      LOGICAL LDV
      DIMENSION A(*),B(*)
      DIMENSION CTR(15),WTH(15),WTHL(15),Z(15)
      DIMENSION W(2:15,5),WP(2:15,3)

      PARAMETER (R1 = 1, HF = R1/2)

      PARAMETER (XL2 =  0.35856 85828 00318 073D0)
      PARAMETER (XL4 =  0.94868 32980 50513 796D0)
      PARAMETER (XL5 =  0.68824 72016 11685 289D0)

      PARAMETER (W2 =  980*R1/6561, W4 = 200*R1/19683)
      PARAMETER (WP2 =  245*R1/486, WP4 = 25*R1/729)

      DATA (W(N,1),W(N,3),N=2,15)
     1     /-0.193872885230909911D+00,  0.518213686937966768D-01,
     2     -0.555606360818980835D+00,  0.314992633236803330D-01,
     3     -0.876695625666819078D+00,  0.111771579535639891D-01,
     4     -0.115714067977442459D+01, -0.914494741655235473D-02,
     5     -0.139694152314179743D+01, -0.294670527866686986D-01,
     6     -0.159609815576893754D+01, -0.497891581567850424D-01,
     7     -0.175461057765584494D+01, -0.701112635269013768D-01,
     8     -0.187247878880251983D+01, -0.904333688970177241D-01,
     9     -0.194970278920896201D+01, -0.110755474267134071D+00,
     A     -0.198628257887517146D+01, -0.131077579637250419D+00,
     B     -0.198221815780114818D+01, -0.151399685007366752D+00,
     C     -0.193750952598689219D+01, -0.171721790377483099D+00,
     D     -0.185215668343240347D+01, -0.192043895747599447D+00,
     E     -0.172615963013768225D+01, -0.212366001117715794D+00/

      DATA (W(N,5),W(N+1,5),N=2,14,2)
     1     / 0.871183254585174982D-01,  0.435591627292587508D-01,
     2     0.217795813646293754D-01,  0.108897906823146873D-01,
     3     0.544489534115734364D-02,  0.272244767057867193D-02,
     4     0.136122383528933596D-02,  0.680611917644667955D-03,
     5     0.340305958822333977D-03,  0.170152979411166995D-03,
     6     0.850764897055834977D-04,  0.425382448527917472D-04,
     7     0.212691224263958736D-04,  0.106345612131979372D-04/

      DATA (WP(N,1),WP(N,3),N=2,15)
     1     /-0.133196159122085045D+01,  0.445816186556927292D-01,
     2     -0.229218106995884763D+01, -0.240054869684499309D-01,
     3     -0.311522633744855959D+01, -0.925925925925925875D-01,
     4     -0.380109739368998611D+01, -0.161179698216735251D+00,
     5     -0.434979423868312742D+01, -0.229766803840877915D+00,
     6     -0.476131687242798352D+01, -0.298353909465020564D+00,
     7     -0.503566529492455417D+01, -0.366941015089163228D+00,
     8     -0.517283950617283939D+01, -0.435528120713305891D+00,
     9     -0.517283950617283939D+01, -0.504115226337448555D+00,
     A     -0.503566529492455417D+01, -0.572702331961591218D+00,
     B     -0.476131687242798352D+01, -0.641289437585733882D+00,
     C     -0.434979423868312742D+01, -0.709876543209876532D+00,
     D     -0.380109739368998611D+01, -0.778463648834019195D+00,
     E     -0.311522633744855959D+01, -0.847050754458161859D+00/

      RESULT=0
      ABSERR=0
      IFAIL=3
      IF(N .LT. 2 .OR. N .GT. 15) RETURN
      IF(MINPTS .GT. MAXPTS) RETURN

      IFNCLS=0
      LDV=.FALSE.
      TWONDM=2**N
      IRGNST=2*N+3
      IRLCLS=2**N+2*N*(N+1)+1
      ISBRGN=IRGNST
      ISBRGS=IRGNST
      IF(MAXPTS .LT. IRLCLS) RETURN
      DO 10 J = 1,N
         CTR(J)=(B(J)+A(J))*HF
 10   WTH(J)=(B(J)-A(J))*HF

 20   RGNVOL=TWONDM
      DO 30 J = 1,N
         RGNVOL=RGNVOL*WTH(J)
 30   Z(J)=CTR(J)
      SUM1=F(N,Z)

      DIFMAX=0
      SUM2=0
      SUM3=0
      DO 40 J = 1,N
         Z(J)=CTR(J)-XL2*WTH(J)
      F2=F(N,Z)
      Z(J)=CTR(J)+XL2*WTH(J)
      F2=F2+F(N,Z)
      WTHL(J)=XL4*WTH(J)
      Z(J)=CTR(J)-WTHL(J)
      F3=F(N,Z)
      Z(J)=CTR(J)+WTHL(J)
      F3=F3+F(N,Z)
      SUM2=SUM2+F2
      SUM3=SUM3+F3
      DIF=ABS(7*F2-F3-12*SUM1)
      DIFMAX=MAX(DIF,DIFMAX)
      IF(DIFMAX .EQ. DIF) IDVAXN=J
 40   Z(J)=CTR(J)

      SUM4=0
      DO 70 J = 2,N
         J1=J-1
      DO 60 K = J,N
         DO 50 L = 1,2
         WTHL(J1)=-WTHL(J1)
      Z(J1)=CTR(J1)+WTHL(J1)
      DO 50 M = 1,2
         WTHL(K)=-WTHL(K)
      Z(K)=CTR(K)+WTHL(K)
 50   SUM4=SUM4+F(N,Z)
 60   Z(K)=CTR(K)
 70   Z(J1)=CTR(J1)

      SUM5=0
      DO 80 J = 1,N
         WTHL(J)=-XL5*WTH(J)
 80   Z(J)=CTR(J)+WTHL(J)
 90   SUM5=SUM5+F(N,Z)
      DO 100 J = 1,N
         WTHL(J)=-WTHL(J)
      Z(J)=CTR(J)+WTHL(J)
      IF(WTHL(J) .GT. 0) GO TO 90
 100  CONTINUE

      RGNCMP=RGNVOL*(WP(N,1)*SUM1+WP2*SUM2+WP(N,3)*SUM3+WP4*SUM4)
      RGNVAL=W(N,1)*SUM1+W2*SUM2+W(N,3)*SUM3+W4*SUM4+W(N,5)*SUM5
      RGNVAL=RGNVOL*RGNVAL
      RGNERR=ABS(RGNVAL-RGNCMP)
      RESULT=RESULT+RGNVAL
      ABSERR=ABSERR+RGNERR
      IFNCLS=IFNCLS+IRLCLS

      IF(LDV) THEN
 110     ISBTMP=2*ISBRGN
       IF(ISBTMP .GT. ISBRGS) GO TO 160
       IF(ISBTMP .LT. ISBRGS) THEN
          ISBTPP=ISBTMP+IRGNST
        IF(WK(ISBTMP) .LT. WK(ISBTPP)) ISBTMP=ISBTPP
      ENDIF
       IF(RGNERR .GE. WK(ISBTMP)) GO TO 160
       DO 130 K = 0,IRGNST-1
 130      WK(ISBRGN-K)=WK(ISBTMP-K)
       ISBRGN=ISBTMP
       GO TO 110
      ENDIF
 140  ISBTMP=(ISBRGN/(2*IRGNST))*IRGNST
      IF(ISBTMP .GE. IRGNST .AND. RGNERR .GT. WK(ISBTMP)) THEN
         DO 150 K = 0,IRGNST-1
 150      WK(ISBRGN-K)=WK(ISBTMP-K)
       ISBRGN=ISBTMP
       GO TO 140
      ENDIF

 160  WK(ISBRGN)=RGNERR
      WK(ISBRGN-1)=RGNVAL
      WK(ISBRGN-2)=IDVAXN
      DO 170 J = 1,N
         ISBTMP=ISBRGN-2*J-2
      WK(ISBTMP+1)=CTR(J)
 170  WK(ISBTMP)=WTH(J)
      IF(LDV) THEN
         LDV=.FALSE.
       CTR(IDVAX0)=CTR(IDVAX0)+2*WTH(IDVAX0)
       ISBRGS=ISBRGS+IRGNST
       ISBRGN=ISBRGS
       GO TO 20
      ENDIF
      RELERR=ABSERR/ABS(RESULT)
      IF(ISBRGS+IRGNST .GT. IWK) IFAIL=2
      IF(IFNCLS+2*IRLCLS .GT. MAXPTS) IFAIL=1
      IF(RELERR .LT. EPS .AND. IFNCLS .GE. MINPTS) IFAIL=0
      IF(IFAIL .EQ. 3) THEN
         LDV=.TRUE.
       ISBRGN=IRGNST
       ABSERR=ABSERR-WK(ISBRGN)
       RESULT=RESULT-WK(ISBRGN-1)
       IDVAX0=WK(ISBRGN-2)
       DO 190 J = 1,N
          ISBTMP=ISBRGN-2*J-2
       CTR(J)=WK(ISBTMP+1)
 190   WTH(J)=WK(ISBTMP)
       WTH(IDVAX0)=HF*WTH(IDVAX0)
       CTR(IDVAX0)=CTR(IDVAX0)-WTH(IDVAX0)
       GO TO 20
      ENDIF
      NFNEVL=IFNCLS
      RETURN
      END


!       ***********************************************************
!              		     Integração SGS
!       ***********************************************************

        DOUBLE PRECISION FUNCTION SGS1(A,B,EPS,F1 )
        DOUBLE PRECISION A,B,S,U,V,SF,C,SL,SG
        DOUBLE PRECISION SP,SA,F1,SGS9,EPS,abb
        EXTERNAL F1
        S = 0.d0
        U = A
    1   V = B
        IF ( U .LT. B ) THEN
        SF = SGS9 ( F1,U,V )
    2   C  = (U+V)/2
        SL = SGS9 ( F1,U,C )
        SG = SGS9 ( F1,C,V )
        SP = SL+SG
        abb=abs(sf)
        if(abb.ne.0.0) go to 5
        abb=1.
    5   SA = ABS(SP-SF)/(abb*EPS)
        IF (  SA.GE.1.0 ) THEN
        V  = C
        SF = SL
        GOTO 2
        END IF
        U = V
        S = S+SP
        GOTO 1
        END IF
        SGS1=S
        RETURN
        END

        DOUBLE PRECISION FUNCTION SGS9 ( F1,A,B )
        DOUBLE PRECISION A,B,H,S,C,X,F1
        EXTERNAL F1
        C = (A+B)/2
        H = (B-A)/2
        X = .96028985E0*H
        S = .10122853E0*(F1(C+X)+F1(C-X))
        X = .79666647E0*H
        S = S + .22238103E0*(F1(C+X)+F1(C-X))
        X = .52553240E0*H
        S = S + .31370664E0*(F1(C+X)+F1(C-X))
        X = .18343464E0*H
        S = S + .36268378E0*(F1(C+X)+F1(C-X))
        SGS9 = S * H
        RETURN
        END


