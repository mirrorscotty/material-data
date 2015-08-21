C       Achanta's effective diffusvity program
        REAL K,N,P
        TC = 71.0
        TK = 273.15 + TC
         P = 0.77
        D0 = 149.7e-12
        DV = 1.0e-10*((TK/328.15)**1.85)*(2.0/P)
        A1 = 7.967
        B1 = 1668.0
        C1 = 228.0
        D1 = A1 - B1/(C1 + TC)
        PV0= ALOG(D1)/760.0
        EPS= 0.27
        RHOS = 1.21
        RHORAT = 1.21
        P1 = 24.0e3
        P2 = 25.0
        K  = 0.176-1.748e-3*TC
        N  = 0.182+6.946e-3*TC
        A  = PV0/(4.56e-3*TK)
        X  = 0.04
 40     X   = X + 0.01
        IF (X.GT.0.30)GOTO 400
        EB = P1*EXP(-P2*X)
	G  = 1.987
	DL = D0*EXP(-EB/(10*G*TK))
        PAR1= (K/X)**(1.0/N)+1.0
        PAR2= A*(PAR1-1.0)/((PAR1**2)*N*X)
        RHO1= A/PAR1
        PAR3= (353.4/TK) - 0.611*RHO1
        PAR4= (353.4/TK)*PAR2/(PAR3**2)
        H   = PAR3*DV/(1.0D3*(1-EPS)*(RHO1/PAR3))
        PAR5= H*PAR4/RHORAT
        DEFF= DL + PAR5
        WRITE(*,*)DEFF,X
        GOTO 40
 400	CONTINUE
	STOP
	END

