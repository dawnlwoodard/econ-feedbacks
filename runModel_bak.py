import modelFuncsTest2 as mf
import numpy as np
import sys
import matplotlib.pylab as py

reload(mf)
#**************************** Parse arguments for run ********************#
if (len(sys.argv) == 3):
        scenario = sys.argv[1]
        nyears = int(sys.argv[2])
else:
        scenario = input("Enter the run: ")
        nyears = input("Enter the number of years: ")



#**************************** Set parameter values  ********************#
initialCappm = 280.0               # initial atmospheric mass in ppm (1ppm = 2.1 Pg C)
initialCaPgC = initialCappm*2.1  # initial atmospheric mass in Pg C
oceanarea    = 3.61e14           # area of the ocean surface
kgas         = 0.062             # gas exchange rate mol per m^2 per ppm per year
dic          = 0.00191           # dissolved inorganic carbon (DIC, mol/L)
taumtod      = 9.0               # exchange between mixed layer and deep ocean
depthm       = 75.0              # depth of mixed layer
deptho       = 3800.0            # depth of full ocean
Revellef     = 11.1              # Revelle factor 
deltat       = 1.0               # time step of Euler's integration
mwc          = 12.01             # molecular weight of carbon
taub         = 20.0              # turnover time of the terrestrial biosphere
kb           = 1.0/taub          # decomposition rate constant for terrestrial biosphere
NPPinitial   = 50.0              # Pg C/yr Net primary production
alpha        = 0.0062            # K/ppm is the transient climate sensitivity
deltat       = 1.0               # time step of Euler's integration
beta         = 0.6              # sensitivity parameter for NPP to CO2
Q10          = 2.0               # sensitivity of respiration to temperature
deltaF       = 5.35*np.log(1140/280)

cat = 8.0
tautempfac = 2.5

aland = -0.0022648
bland = 0.04661
cland = 0.789

#K0 = -0.4194842
K0 = -0.4194842
K1 = 0.9330434
K2 = -0.1868783
#K3 = 0.02036102
#K3 = 0.0187
K3 = 0.01516102


rcpMidyearConc = np.genfromtxt("RCP85_MIDYR_CONC.DAT",skip_header=35)
rcpCO2 = rcpMidyearConc[3:,3]
rcpCO2 = rcpCO2[36:336]

def divZero(a,b):
	if (len(a)!=len(b)):
		return
	else:
		out = np.zeros(len(a))
		for i in range(len(a)):
			if b[i]!=0:
				out[i] = float(a[i])/b[i]
			else:
				out[i] = 0
		return out
#K0 = -0.4194842
#K1 = 0.9330434
#K2 = -0.1868783

#x = [ -8.45182171e-04,-2.36273779e-02,1.89598788e-02,1.17399917e+00,6.66046984e-01,4.89989197e+01,2.00558292e+01,1.39181246e+01,8.53741862e+01,3.64501171e+00,1.08250727e+01,1.49786807e-02,2.50000000e+02]

x = np.genfromtxt('optimizedParams.txt',delimiter=',')
print x
aland,bland,K0,K1,K2,K3,Q10,beta,NPPinitial,taub,Revellef,depthm,tautempfac,cat,kgas,initialCappm=x

#K3 = K3+0.004
#aland,K3,Q10,beta,Revellef,tautempfac=x

#params = {"aland":aland,"K3":K3,"Q10":Q10,"beta":beta,"NPPinitial":NPPinitial,"taub":taub,"Revellef":Revellef,"depthm":depthm,"tautempfac":tautempfac,"cat":cat,"kgas":kgas,"initialCappm":initialCappm]


if (scenario=='debug'):
        dCo,dCl,dCa,dT,tair,gdp,old_gdp,KayaFF,Ff,fonet,fbnet,Ca,date,pop,pop_adj,KayaE,KayaCE,ee,ED,EE = mf.runModel('debug','85','NFC','base',False,nyears,aland, bland, K0,K1,K2,K3,Q10,beta,NPPinitial,taub,Revellef,depthm,tautempfac,cat,kgas,initialCappm)
	print ' \n***************** Running NFC ********************'	
	print 'RCP8.5 (1850-2100): dT,dTm,dCl,dCo,dCa',tair[299]-tair[61],np.mean(tair[281:299])-np.mean(tair[186:205]),np.sum(fbnet[50:]),np.sum(fonet[50:]),Ca[-1]/2.12
        print 'Historical (1850-2005): dT,dCl,dCo,dCa',np.mean(tair[186:205])-np.mean(tair[51:79]),np.sum(fbnet[50:205]),np.sum(fonet[50:205]),Ca[205]/2.12
	print '***************** End NFC ********************\n '	

        dCo,dCl,dCa,dT,tair,gdp,old_gdp,KayaFF,Ff,fonet,fbnet,Ca,date,pop,pop_adj,KayaE,KayaCE,ee,ED,EE = mf.runModel('debug','85','NFC','min',False,nyears,aland, bland, K0,K1,K2,K3,Q10,beta,NPPinitial,taub,Revellef,depthm,tautempfac,cat,kgas,initialCappm)
	print ' \n***************** Running NFC ********************'	
	print 'RCP8.5 (1850-2100): dT,dTm,dCl,dCo,dCa',tair[299]-tair[61],np.mean(tair[281:299])-np.mean(tair[186:205]),np.sum(fbnet[50:]),np.sum(fonet[50:]),Ca[-1]/2.12
        print 'Historical (1850-2005): dT,dCl,dCo,dCa',np.mean(tair[186:205])-np.mean(tair[51:79]),np.sum(fbnet[50:205]),np.sum(fonet[50:205]),Ca[205]/2.12
	print '***************** End NFC ********************\n '	
        
	dCo,dCl,dCa,dT,tair,gdp,old_gdp,KayaFF,Ff,fonet,fbnet,Ca,date,pop,pop_adj,KayaE,KayaCE,ee,ED,EE = mf.runModel('debug','85','NFC','max',False,nyears,aland, bland, K0,K1,K2,K3,Q10,beta,NPPinitial,taub,Revellef,depthm,tautempfac,cat,kgas,initialCappm)
	print ' \n***************** Running NFC ********************'	
	print 'RCP8.5 (1850-2100): dT,dTm,dCl,dCo,dCa',tair[299]-tair[61],np.mean(tair[281:299])-np.mean(tair[186:205]),np.sum(fbnet[50:]),np.sum(fonet[50:]),Ca[-1]/2.12
        print 'Historical (1850-2005): dT,dCl,dCo,dCa',np.mean(tair[186:205])-np.mean(tair[51:79]),np.sum(fbnet[50:205]),np.sum(fonet[50:205]),Ca[205]/2.12
	print '***************** End NFC ********************\n '	

	
	Cadiff = Ca[1:]-Ca[:-1]
        rcpdiff = rcpCO2[1:]-rcpCO2[:-1]
        py.close(19)
        py.figure(19,figsize=(8,6),dpi=150)
        #py.plot(date[50:211], -fonet[50:211], color="dodgerblue", linestyle='-',linewidth=3.0,label='Modelled')
        #py.plot(date[50:211], oceanHoff, color = "dodgerblue",linestyle='--',label='Data')
        py.plot(date, rcpCO2*2.12, color = "gray",linestyle='--',label='Data')
        py.plot(date, Ca*2.12, color="gray", linestyle='-',linewidth=3.0,label='Modelled')
        #py.plot(co2GCP2Year, co2GCP2Ff,color = "gray",linestyle='--')
        #py.plot(co2GCP2Year, co2GCP2Luc,color = "goldenrod",linestyle='--')
        #py.plot(co2GCP2Year, co2GCP2Atm,color = "mediumslateblue",linestyle='--')
        py.legend(loc='upper left')
        py.ylabel('Atmospheric Carbon (PgC)')
        py.xlabel('Time (years)')
        py.savefig('atmosphereConcTuning.png')

        py.close(19)
        py.figure(19,figsize=(8,6),dpi=150)
        py.plot(date[1:], Cadiff*2.12, color="gray", linestyle='-',linewidth=3.0,label='Modelled')
        py.plot(date[1:], rcpdiff*2.12, color = "gray",linestyle='--',label='Data')
        py.plot(date,Ff,color='k',label='Fossil Fuel')
        #py.plot(co2GCP2Year, co2GCP2Ff,color = "gray",linestyle='--')
        #py.plot(co2GCP2Year, co2GCP2Luc,color = "goldenrod",linestyle='--')
        #py.plot(co2GCP2Year, co2GCP2Atm,color = "mediumslateblue",linestyle='--')
        py.legend(loc='upper left')
        py.ylabel('Atmospheric Carbon Flux (PgC/yr)')
        py.xlabel('Time (years)')
        py.savefig('atmosphereFluxTuning.png')


	#dCo,dCl,dCa,dT,tair,gdp,old_gdp,ff,Ff,fonet,fbnet,ca,date = mainFun('85','DFC','base',False)
        #print 'dCo,dCl,dT,dCa,Ff,KayaFF =',dCo,dCl,dT,dCa,np.sum(Ff),np.sum(KayaFF)


elif (scenario=='1pct'):
	nyears = 141
	for j in ['BGC','NFC']:
			dCo,dCl,dCa,dT,tair,Ff,fonet,fbnet,date = mf.runModel(scenario,'none',j,'base',True,nyears,aland, bland, K0,K1,K2,K3,Q10,beta,NPPinitial,taub,Revellef,depthm,tautempfac,cat,kgas,initialCappm)
			dCO2m = dCa/2.12
			if j == 'BGC':
				C_buc = dCa
				Ff_bgc = np.sum(Ff)
				betaO=dCo/dCO2m
				betaL=dCl/dCO2m
				betaA = -(betaO+betaL)
				print 'betaA,betaL,betaO',betaA,betaL,betaO
			elif j == 'NFC':
				print 'dCo,dCl,dT,dCa =',dCo,dCl,dT,dCa
				gamO = (dCo-betaO*dCa)/dT
				gamL = (dCl-betaL*dCa)/dT
				gamA = -(gamO+gamL)
				C_uc = dCa
				Ff_nfc = np.sum(Ff)
				gain = (Ff_bgc - Ff_nfc)/Ff_bgc
				gain2 = (gamA*dT)/(dCa*(2.12+betaL + betaO))
				gain1 = 1 - C_buc/C_uc
				print 'gain1,gain2,gamA,gamL,gamO = ',gain,gain2,gamA,gamL,gamO
                        	print 'dT,dCl,dCo,dCa',dT,dCl,dCo,dCa



elif (scenario=='all'):
        output = open('outputTable2.dat','w')
        outputBGC = open('BGCoutputTable2.dat','w')
        outcsv = open('outputTable2.txt','w')
        outcsvBGC = open('BGCoutputTable2.txt','w')
        outcsvBGC2 = open('BGCoutputTable.txt','w')
        outBGCParams = open('BGCparams.csv','w')
        outParams = open('outParams.csv','w')
        outcsvECON = open('econOutputTable.csv','w')
	otherEconOutput = open('econAllVals.csv','w')
        output.write('{:<3} {} {:<5} {:<4} {:<6} {:<6} {:<5} {:<5} {:<7} {:<9} {:<6} {:<6} {:<7}\n'.format('RCP','Run','Level','dT','dCa','dFF','dCo','dCl','Gain','GamA','GamO','GamL','GamFF'))
        tempAll=[]
	gainCalc = {}
        for k in ['85']:
                #for j in ['BGC','NFC']:
                        #for l in ['min','base','max']:
                for j in ['BGC','LFC','OFC','NFC','PFC','GFC','EFC','DFC','AFC','FC']:
			gainCalc[j] = {}
			for l in ['max','base','min']:
				gainCalc[j][l] = {}
                                #for m in [False,True]:
                                for m in [True,False]:
					gainCalc[j][l][m] = {}
                                        #if (j!='FC' and j!='AFC' and j!='GFC' and j!='EFC' and j!='DFC' and l!='base'):
                                        if (l!='base' and (j in ['LFC','OFC']) ):
                                                continue
                                        print k,j,l,m
					dCo,dCl,dCa,dT,Tair,gdp,old_gdp,KayaFF,Ff,fonet,fbnet,Ca,date,pop,pop_adj,KayaE,KayaCE,ee,ED,EE = mf.runModel(scenario,k,j,l,m,nyears,aland, bland, K0,K1,K2,K3,Q10,beta,NPPinitial,taub,Revellef,depthm,tautempfac,cat,kgas,initialCappm)
					oldKayaE = ee/(old_gdp*pop)
					oldKayaCE = Ff/ee
                                        CO2m = dCa/2.12
					gainCalc[j][l][m] = Ca-Ca[0]
					print 'dCo,dCl,dT,dCa,Ff,KayaFF =',dCo,dCl,dT,dCa,np.sum(Ff),np.sum(KayaFF)
                                        if j == 'BGC':
                                                C_buc = dCa
                                                betaO=dCo/CO2m
                                                betaL=dCl/CO2m
                                                #betaFF = np.sum(Ff)/CO2m
                                                #print betaFF
                                                betaA2 = dCa/CO2m
                                                betaA = -(betaO+betaL)
                                                print 'betaA,betaA2',betaA,betaA2
						gainTS = (Ca-Ca[0])*0
                                                if not m:
							np.savetxt("ModelData/gainCumulative"+k+j+l+".csv",gainTS,delimiter=",")
                                        #               output.write('{:>3} {:>3} {:>5} {:0<4.3} {:0<6.5} {:0<6.5} {:0<5.4} {:0<5.4} {:-<7} {:-<9} {:-<6} {:->6} {:->7}\n'.format(k,j,l,dT,dCa,np.sum(Ff),dCo,dCl,'-','-','-','-','-'))
                                        #               outcsv.write('{:>3},{:>3},{:>5},{:0<4.3},{:0<6.5},{:0<6.5},{:0<5.4},{:0<5.4},{},{},{},{},{}\n'.format(k,j,l,dT,dCa,np.sum(Ff),dCo,dCl,0,0,0,0,0))
                                                if m:
							np.savetxt("ModelData/gainCumulativeBGC"+k+j+l+".csv",gainTS,delimiter=",")
							outParams.write('{:>3},{:>3},{:>5},{:0<4.5},{:0<6.5},{:0<6.5},{:0<5.4},{:0<5.4},{},{},{},{},{}\n'.format(k,j,l,dT,dCa,np.sum(KayaFF),dCo,dCl,0,0,0,0,0))
                                                        outputBGC.write('{:>3} {:>3} {:>5} {:0<4.3} {:0<6.5} {:0<6.5} {:0<5.4} {:0<5.4} {:-<7} {:-<9} {:-<6} {:->6} {:->7}\n'.format(k,j,l,dT,dCa,np.sum(Ff),dCo,dCl,'-','-','-','-','-'))
                                                        #outcsvBGC.write('{:>3};{:>3};{:>5};{:.2f};{:.0f};{:.0f};{:.0f};{:.0f};{};{};{};{};{}\n'.format(k,j,l,dT,dCa,np.sum(Ff),dCo,dCl,'---','---','---','---','---'))
                                                        outcsv.write('{:>3};{:>3};{:>5};{:.2f};{:.0f};{:.0f};{:.0f};{:.0f};{};{};{};{};{}\n'.format(k,j,l,dT,dCa,np.sum(Ff),dCo,dCl,'---','---','---','---','---'))
                                                        outcsvBGC2.write('{:>3};{:>3};{:>5};{:0<4.3};{:0<6.5};{:0<6.5};{:0<5.4};{:0<5.4};{};{};{};{};{}\n'.format(k,j,l,dT,dCa,np.sum(Ff),dCo,dCl,'---','---','---','---','---'))
                                                        outBGCParams.write('{:>3},{:>3},{:>5},{:0<4.5},{:0<6.5},{:0<6.5},{:0<5.4},{:0<5.4},{},{},{},{},{}\n'.format(k,j,l,dT,dCa,np.sum(Ff),dCo,dCl,0,0,0,0,0))
							if (l=='base'):
								#outcsv.write('{:>3},;{:>3},;{:>5},;{:0<4.3},;{:0<6.5},;{:0<6.5},;{:0<5.4},;{:0<5.4},;{: .4f},;{: 9.4f},;{: 6.2f},;{: 6.2f},;{: 7.2f},\n'.format(k,j,l,dT,dCa,np.sum(KayaFF),dCo,dCl,gain,gamA,gamO,gamL,gamFF))
								outcsvBGC.write('{:>3},;{:>3},;{:>5},;{:.2f},;{:.0f},;{:.0f},;{:.0f},;{:.0f},;{},;{},;{},;{},;{},\n'.format(k,j,l,dT,dCa,np.sum(KayaFF),dCo,dCl,'---','---','---','---','---'))
							elif (l=='min'):
								outcsvBGC.write('{:>3});{:>3});{:>5});{:.2f});{:.0f});{:.0f});{:.0f});{:.0f});{});{});{});{});{})\n'.format(k,j,l,dT,dCa,np.sum(KayaFF),dCo,dCl,'---','---','---','---','---'))
							elif (l=='max'):
								outcsvBGC.write('({:>3},;({:>3},;({:>5},;({:.2f},;({:.0f},;({:.0f},;({:.0f},;({:.0f},;({},;({},;({},;({},;({},\n'.format(k,j,l,dT,dCa,np.sum(KayaFF),dCo,dCl,'---','---','---','---','---'))
                                                        #outBGCParams.write('{:>3},{:>3},{:>5},{:0<4.5},{:0<6.5},{:0<6.5},{:0<5.4},{:0<5.4},{},{},{},{},{}\n'.format(k,j,l,dT,dCa,np.sum(Ff),dCo,dCl,0,0,0,0,0))
                                        elif j in ['NFC','LFC','OFC']:
						#oldKayaE = KayaE
						#oldKayaCE = KayaCE
                                                gamO = (dCo-betaO*CO2m)/dT
                                                gamL = (dCl-betaL*CO2m)/dT
                                                #gamFF = (np.sum(KayaFF)-betaFF*CO2m)/dT
                                                #gamFF = (np.sum(KayaFF))/dT
                                                gamFF = (np.sum(KayaFF)-np.sum(Ff))/dT
                                                gamA = gamFF-(gamO+gamL)
                                                gamA2 = (dCa-betaA*CO2m)/dT
                                                C_uc = dCa
						C_test = Ca[-1]-Ca[0]
						Call_uc = Ca-Ca[0]
						print 'Carbon BGC, NFC',C_uc,C_buc
                                                gain = 1 - C_buc/C_uc
						gainTest = C_buc/C_test
						print 'Gain Test = ',gainTest
						#gainTS = 1 - divZero(np.cumsum(Call_buc),np.cumsum(Call_uc))
						#gainTS = 1 - divZero(Call_buc,Call_uc)
						gainTS = 1 - divZero(gainCalc['BGC'][l][m],gainCalc[j][l][m])
						gainTS[0] = 0
						gain = gainTS[-1]
						#gainAll = 1 - divZero(Call_buc,Call_uc)
						#gainAll[0] = 0
						#gainTS = np.cumsum(gainAll) 
						#gainTS[0] = 0
                                                print 'gain1,gamA,gamL,gamO,gamFF = ',gain,gamA,gamL,gamO,gamFF
						print 'RCP8.5 (1850-2100): dT,dTm,dCl,dCo,dCa',Tair[299]-Tair[61],np.mean(Tair[281:299])-np.mean(Tair[186:205]),np.sum(fbnet[50:]),np.sum(fonet[50:]),Ca[-1]/2.12

        					print 'Historical (1850-2005): dT,dCl,dCo,dCa',np.mean(Tair[186:205])-np.mean(Tair[51:79]),np.sum(fbnet[50:205]),np.sum(fonet[50:205]),Ca[205]/2.12
                                                if not m:
							np.savetxt("ModelData/gainCumulative"+k+j+l+".csv",gainTS,delimiter=",")
							output.write('{:>3} {:>3} {:>5} {:0<4.3} {:0<6.5} {:0<6.5} {:0<5.4} {:0<5.4} {: .4f} {: 9.4f} {: 6.2f} {: 6.2f} {: 7.2f}\n'.format(k,j,l,dT,dCa,np.sum(KayaFF),dCo,dCl,gain,gamA,gamO,gamL,gamFF))
							otherEconOutput.write('{:>3},{:>3},{:>5},{: .4f},{: .4f},{: .4f},{: .4f},{: .4f},{: .4f}\n'.format(k,j,l,dT,pop_adj[-1],gdp[-1],ED[-1],EE[-1],np.sum(KayaFF)))
							#otherEconOutput.write('{:>3} {:>3} {:>5} {:0<4.3} {:0<6.5} {:0<6.5} {:0<5.4} {:0<5.4} {: 6.2f}\n'.format(k,j,l,dT,pop_adj[-1],gdp[-1],ED[-1],EE[-1],np.sum(KayaFF)))
	
							outParams.write('{:>3},{:>3},{:>5},{:0<4.5},{:0<6.5},{:0<6.5},{:0<5.4},{:0<5.4},{: .2f},{: 9.2f},{: 6.2f},{: 6.2f},{: 7.2f}\n'.format(k,j,l,dT,dCa,np.sum(KayaFF),dCo,dCl,gain,gamA,gamO,gamL,gamFF))
							#outcsv2.write('{:>3},{:>3},{:>5},{:0<4.3},{:0<6.5},{:0<6.5},{:0<5.4},{:0<5.4},{: .4f},{: 9.4f},{: 6.2f},{: 6.2f},{: 7.2f},\n'.format(k,j,l,dT,dCa,np.sum(KayaFF),dCo,dCl,gain,gamA,gamO,gamL,gamFF))
							#outcsv.write('{:>3},;{:>3},;{:>5},;{:0<4.3},;{:0<6.5},;{:0<6.5},;{:0<5.4},;{:0<5.4},;{: .4f},;{: 9.4f},;{: 6.2f},;{: 6.2f},;{: 7.2f},\n'.format(k,j,l,dT,dCa,np.sum(KayaFF),dCo,dCl,gain,gamA,gamO,gamL,gamFF))
							#outcsv.write('{:>3},;{:>3},;{:>5},;{:.2f},;{:.0f},;{:.0f},;{:.0f},;{:.0f},;{:.3f},;{:.2f},;{:.2f},;{:.2f},;{:.2f},\n'.format(k,j,l,dT,dCa,np.sum(KayaFF),dCo,dCl,gain,gamA,gamO,gamL,gamFF))
							outcsv.write('{:>3};{:>3};{:>5};{:.2f};{:.0f};{:.0f};{:.0f};{:.0f};{:.3f};{:.2f};{:.2f};{:.2f};{:.2f}\n'.format(k,j,l,dT,dCa,np.sum(KayaFF),dCo,dCl,gain,gamA,gamO,gamL,gamFF))
							#elif (l=='min'):
								#outcsv.write('{:>3});{:>3});{:>5});{:.2f});{:.0f});{:.0f});{:.0f});{:.0f});{:.3f});{:.2f});{:.2f});{:.2f});{:.2f})\n'.format(k,j,l,dT,dCa,np.sum(KayaFF),dCo,dCl,gain,gamA,gamO,gamL,gamFF))
							#elif (l=='max'):
								#outcsv.write('({:>3},;({:>3},;({:>5},;({:.2f},;({:.0f},;({:.0f},;({:.0f},;({:.0f},;({:.3f},;({:.2f},;({:.2f},;({:.2f},;({:.2f},\n'.format(k,j,l,dT,dCa,np.sum(KayaFF),dCo,dCl,gain,gamA,gamO,gamL,gamFF))
                                                else:
							np.savetxt("ModelData/gainCumulativeBGC"+k+j+l+".csv",gainTS,delimiter=",")
                                                        outputBGC.write('{:>3} {:>3} {:>5} {:0<4.3} {:0<6.5} {:0<6.5} {:0<5.4} {:0<5.4} {: .4f} {: 9.4f} {: 6.2f} {: 6.2f} {: 7.2f}\n'.format(k,j,l,dT,dCa,np.sum(KayaFF),dCo,dCl,gain,gamA,gamO,gamL,gamFF))
							#outcsv2.write('{:>3},{:>3},{:>5},{:0<4.3},{:0<6.5},{:0<6.5},{:0<5.4},{:0<5.4},{: .4f},{: 9.4f},{: 6.2f},{: 6.2f},{: 7.2f},\n'.format(k,j,l,dT,dCa,np.sum(KayaFF),dCo,dCl,gain,gamA,gamO,gamL,gamFF))
							outBGCParams.write('{:>3},{:>3},{:>5},{:0<4.5},{:0<6.5},{:0<6.5},{:0<5.4},{:0<5.4},{: .4f},{: 9.4f},{: 6.2f},{: 6.2f},{: 7.2f}\n'.format(k,j,l,dT,dCa,np.sum(KayaFF),dCo,dCl,gain,gamA,gamO,gamL,gamFF))
                                                        if (l=='base'):
								outcsvBGC.write('{:>3},;{:>3},;{:>5},;{:0<4.3},;{:0<6.5},;{:0<6.5},;{:0<5.4},;{:0<5.4},;{: .4f},;{: 9.4f},;{: 6.2f},;{: 6.2f},;{: 7.2f},\n'.format(k,j,l,dT,dCa,np.sum(KayaFF),dCo,dCl,gain,gamA,gamO,gamL,gamFF))
							elif (l=='min'):
								outcsvBGC.write('{:>3});{:>3});{:>5});{:0<4.3});{:0<6.5});{:0<6.5});{:0<5.4});{:0<5.4});{: .4f});{: 9.4f});{: 6.2f});{: 6.2f});{: 7.2f})\n'.format(k,j,l,dT,dCa,np.sum(KayaFF),dCo,dCl,gain,gamA,gamO,gamL,gamFF))
							elif (l=='max'):
								outcsvBGC.write('({:>3},;({:>3},;({:>5},;({:0<4.3},;({:0<6.5},;({:0<6.5},;({:0<5.4},;({:0<5.4},;({: .4f},;({: 9.4f},;({: 6.2f},;({: 6.2f},;({: 7.2f},\n'.format(k,j,l,dT,dCa,np.sum(KayaFF),dCo,dCl,gain,gamA,gamO,gamL,gamFF))
					
					else:
                                                C_fc = dCa
						Call_fc = Ca-Ca[0]
                                                gamO = (dCo-betaO*CO2m)/dT
                                                gamL = (dCl-betaL*CO2m)/dT
                                                #gamFF = -(np.sum(KayaFF)-betaFF*CO2m)/dT
                                                #gamFF = (np.sum(KayaFF))/dT
                                                gamFF = (np.sum(KayaFF)-np.sum(Ff))/dT
                                                gamA = gamFF-(gamO+gamL)
                                                gamA2 = (dCa-betaA*CO2m)/dT
                                                #print 'alph numbers',Tair[-1],Tair[0],Ca[-1],Ca[0]
                                                alph = (Tair[-1]-Tair[0])/(Ca[-1]*(1/2.13)-Ca[0]*(1/2.13)) # don't forget to convert Ca to ppm!
                                                #gain = -1.9*(gamL+gamO)/(1+betaL+betaO)
                                                gain = 1 - C_buc/C_fc
						#gainTS = 1 - divZero(np.cumsum(Call_buc),np.cumsum(Call_fc))
						#gainTS = 1 - divZero(Call_buc,Call_fc)
						gainTS = 1 - divZero(gainCalc['BGC'][l][m],gainCalc[j][l][m])
						gainTS[0] = 0
						gain = gainTS[-1]
						#gainAll = 1 - divZero(Call_buc,Call_uc)
						#gainAll[0] = 0
						#gainTS = np.cumsum(gainAll) 
                                                print 'gain1,gamA,gamL,gamO,gamFF = ',gain,gamA,gamL,gamO,gamFF
                                                if m:
							np.savetxt("ModelData/gainCumulativeBGC"+k+j+l+".csv",gainTS,delimiter=",")
							#print gai
                                                        outputBGC.write('{:>3} {:>3} {:>5} {:0<4.3} {:0<6.5} {:0<6.5} {:0<5.4} {:0<5.4} {: .4f} {: 9.4f} {: 6.2f} {: 6.2f} {: 7.2f}\n'.format(k,j,l,dT,dCa,np.sum(KayaFF),dCo,dCl,gain,gamA,gamO,gamL,gamFF))
                                                        outcsvBGC2.write('{:>3},{:>3},{:>5},{:0<4.3},{:0<6.5},{:0<6.5},{:0<5.4},{:0<5.4},{: .4f},{: 9.4f},{: 6.2f},{: 6.2f},{: 7.2f}\n'.format(k,j,l,dT,dCa,np.sum(KayaFF),dCo,dCl,gain,gamA,gamO,gamL,gamFF))
                                                        outBGCParams.write('{:>3},{:>3},{:>5},{:0<4.5},{:0<6.5},{:0<6.5},{:0<5.4},{:0<5.4},{: .4f},{: 9.4f},{: 6.2f},{: 6.2f},{: 7.2f}\n'.format(k,j,l,dT,dCa,np.sum(KayaFF),dCo,dCl,gain,gamA,gamO,gamL,gamFF))
							if (l=='base'):
								#outcsv.write('{:>3},;{:>3},;{:>5},;{:0<4.3},;{:0<6.5},;{:0<6.5},;{:0<5.4},;{:0<5.4},;{: .4f},;{: 9.4f},;{: 6.2f},;{: 6.2f},;{: 7.2f},\n'.format(k,j,l,dT,dCa,np.sum(KayaFF),dCo,dCl,gain,gamA,gamO,gamL,gamFF))
								outcsvBGC.write('{:>3},;{:>3},;{:>5},;{:.2f},;{:.0f},;{:.0f},;{:.0f},;{:.0f},;{:.3f},;{:.2f},;{:.2f},;{:.2f},;{:.2f},\n'.format(k,j,l,dT,dCa,np.sum(KayaFF),dCo,dCl,gain,gamA,gamO,gamL,gamFF))
							elif (l=='min'):
								outcsvBGC.write('{:>3});{:>3});{:>5});{:.2f});{:.0f});{:.0f});{:.0f});{:.0f});{:.3f});{:.2f});{:.2f});{:.2f});{:.2f})\n'.format(k,j,l,dT,dCa,np.sum(KayaFF),dCo,dCl,gain,gamA,gamO,gamL,gamFF))
							elif (l=='max'):
								outcsvBGC.write('({:>3},;({:>3},;({:>5},;({:.2f},;({:.0f},;({:.0f},;({:.0f},;({:.0f},;({:.3f},;({:.2f},;({:.2f},;({:.2f},;({:.2f},\n'.format(k,j,l,dT,dCa,np.sum(KayaFF),dCo,dCl,gain,gamA,gamO,gamL,gamFF))
							if (j=='FC'):
								outcsvECON.write('{},{},{},{},{},{}\n'.format(j,l,m,'','',''))
								outcsvECON.write('{:0<.2},{:0<.2},{:0<.2},{:0<.2},{:0<.2},{:0<.2}\n'.format(pop[230],pop_adj[230],100*-1*(pop[230]-pop_adj[230])/pop[230],pop[-1],pop_adj[-1],100*-1*(pop[-1]-pop_adj[-1])/pop[-1]))
								outcsvECON.write('{:0<.2},{:0<.2},{:0<.2},{:0<.2},{:0<.2},{:0<.2}\n'.format(old_gdp[230],gdp[230],100*-1*(old_gdp[230]-gdp[230])/old_gdp[230],old_gdp[-1],gdp[-1],100*-1*(old_gdp[-1]-gdp[-1])/old_gdp[-1]))
								outcsvECON.write('{:0<.2},{:0<.2},{:0<.2},{:0<.2},{:0<.2},{:0<.2}\n'.format(oldKayaE[230],KayaE[230],100*-1*(oldKayaE[230]-KayaE[230])/oldKayaE[230],oldKayaE[-1],KayaE[-1],100*-1*(oldKayaE[-1]-KayaE[-1])/oldKayaE[-1]))
								outcsvECON.write('{:0<.2},{:0<.2},{:0<.2},{:0<.2},{:0<.2},{:0<.2}\n'.format(oldKayaCE[230],KayaCE[230],100*-1*(oldKayaCE[230]-KayaCE[230])/oldKayaCE[230],oldKayaCE[-1],KayaCE[-1],100*-1*(oldKayaCE[-1]-KayaCE[-1])/oldKayaCE[-1]))
								outcsvECON.write('{:0<.2},{:0<.2},{:0<.2},{:0<.2},{:0<.2},{:0<.2}\n'.format(Ff[230],KayaFF[230],100*-1*(Ff[230]-KayaFF[230])/Ff[230],Ff[-1],KayaFF[-1],100*-1*(Ff[-1]-KayaFF[-1])/Ff[-1]))
								outcsvECON.write('{},{},{},{},{},{}\n'.format('','','','','',''))
							"""
							if (l=='base'):
								outcsvBGC.write('{:>3},;{:>3},;{:>5},;{:0<4.3},;{:0<6.5},;{:0<6.5},;{:0<5.4},;{:0<5.4},;{: .4f},;{: 9.4f},;{: 6.2f},;{: 6.2f},;{: 7.2f},\n'.format(k,j,l,dT,dCa,np.sum(KayaFF),dCo,dCl,gain,gamA,gamO,gamL,gamFF))
							elif (l=='min'):
								outcsvBGC.write('{:>3});{:>3});{:>5});{:0<4.3});{:0<6.5});{:0<6.5});{:0<5.4});{:0<5.4});{: .4f});{: 9.4f});{: 6.2f});{: 6.2f});{: 7.2f})\n'.format(k,j,l,dT,dCa,np.sum(KayaFF),dCo,dCl,gain,gamA,gamO,gamL,gamFF))
							elif (l=='max'):
								outcsvBGC.write('({:>3},;({:>3},;({:>5},;({:0<4.3},;({:0<6.5},;({:0<6.5},;({:0<5.4},;({:0<5.4},;({: .4f},;({: 9.4f},;({: 6.2f},;({: 6.2f},;({: 7.2f},\n'.format(k,j,l,dT,dCa,np.sum(KayaFF),dCo,dCl,gain,gamA,gamO,gamL,gamFF))
							"""
                                                elif not m:
							np.savetxt("ModelData/gainCumulative"+k+j+l+".csv",gainTS,delimiter=",")
                                                        outputBGC.write('{:>3} {:>3} {:>5} {:0<4.3} {:0<6.5} {:0<6.5} {:0<5.4} {:0<5.4} {: .4f} {: 9.4f} {: 6.2f} {: 6.2f} {: 7.2f}\n'.format(k,j,l,dT,dCa,np.sum(KayaFF),dCo,dCl,gain,gamA,gamO,gamL,gamFF))
                                                        outParams.write('{:>3},{:>3},{:>5},{:0<4.5},{:0<6.5},{:0<6.5},{:0<5.4},{:0<5.4},{: .4f},{: 9.4f},{: 6.2f},{: 6.2f},{: 7.2f}\n'.format(k,j,l,dT,dCa,np.sum(KayaFF),dCo,dCl,gain,gamA,gamO,gamL,gamFF))
							otherEconOutput.write('{:>3},{:>3},{:>5},{: .4f},{: .4f},{: .4f},{: .4f},{: .4f},{: .4f}\n'.format(k,j,l,dT,pop_adj[-1],gdp[-1],ED[-1],EE[-1],np.sum(KayaFF)))
							"""
							if (l=='base' and j=='FC'):
								outcsvECON.write('{:0<.2},{:0<.2},{:0<.2},{:0<.2},{:0<.2},{:0<.2}\n'.format(pop[230],pop_adj[230],100*-1*(pop[230]-pop_adj[230])/pop[230],pop[-1],pop_adj[-1],100*-1*(pop[-1]-pop_adj[-1])/pop[-1]))
								outcsvECON.write('{:0<.2},{:0<.2},{:0<.2},{:0<.2},{:0<.2},{:0<.2}\n'.format(old_gdp[230],gdp[230],100*-1*(old_gdp[230]-gdp[230])/old_gdp[230],old_gdp[-1],gdp[-1],100*-1*(old_gdp[-1]-gdp[-1])/old_gdp[-1]))
								outcsvECON.write('{:0<.2},{:0<.2},{:0<.2},{:0<.2},{:0<.2},{:0<.2}\n'.format(oldKayaE[230],KayaE[230],100*-1*(oldKayaE[230]-KayaE[230])/oldKayaE[230],oldKayaE[-1],KayaE[-1],100*-1*(oldKayaE[-1]-KayaE[-1])/oldKayaE[-1]))
								outcsvECON.write('{:0<.2},{:0<.2},{:0<.2},{:0<.2},{:0<.2},{:0<.2}\n'.format(oldKayaCE[230],KayaCE[230],100*-1*(oldKayaCE[230]-KayaCE[230])/oldKayaCE[230],oldKayaCE[-1],KayaCE[-1],100*-1*(oldKayaCE[-1]-KayaCE[-1])/oldKayaCE[-1]))
								outcsvECON.write('{:0<.2},{:0<.2},{:0<.2},{:0<.2},{:0<.2},{:0<.2}\n'.format(Ff[230],KayaFF[230],100*-1*(Ff[230]-KayaFF[230])/Ff[230],Ff[-1],KayaFF[-1],100*-1*(Ff[-1]-KayaFF[-1])/Ff[-1]))
							"""
                                                        if (l=='base'):
								#outcsv.write('{:>3},;{:>3},;{:>5},;{:0<4.3},;{:0<6.5},;{:0<6.5},;{:0<5.4},;{:0<5.4},;{: .4f},;{: 9.4f},;{: 6.2f},;{: 6.2f},;{: 7.2f},\n'.format(k,j,l,dT,dCa,np.sum(KayaFF),dCo,dCl,gain,gamA,gamO,gamL,gamFF))
								outcsv.write('{:>3},;{:>3},;{:>5},;{:.2f},;{:.0f},;{:.0f},;{:.0f},;{:.0f},;{:.3f},;{:.2f},;{:.2f},;{:.2f},;{:.2f},\n'.format(k,j,l,dT,dCa,np.sum(KayaFF),dCo,dCl,gain,gamA,gamO,gamL,gamFF))
							elif (l=='min'):
								#outcsv.write('{:>3});{:>3});{:>5});{:0<4.3});{:0<6.5});{:0<6.5});{:0<5.4});{:0<5.4});{: .4f});{: 9.4f});{: 6.2f});{: 6.2f});{: 7.2f})\n'.format(k,j,l,dT,dCa,np.sum(KayaFF),dCo,dCl,gain,gamA,gamO,gamL,gamFF))
								outcsv.write('{:>3});{:>3});{:>5});{:.2f});{:.0f});{:.0f});{:.0f});{:.0f});{:.3f});{:.2f});{:.2f});{:.2f});{:.2f})\n'.format(k,j,l,dT,dCa,np.sum(KayaFF),dCo,dCl,gain,gamA,gamO,gamL,gamFF))
							elif (l=='max'):
								#outcsv.write('({:>3},;({:>3},;({:>5},;({:0<4.3},;({:0<6.5},;({:0<6.5},;({:0<5.4},;({:0<5.4},;({: .4f},;({: 9.4f},;({: 6.2f},;({: 6.2f},;({: 7.2f},\n'.format(k,j,l,dT,dCa,np.sum(KayaFF),dCo,dCl,gain,gamA,gamO,gamL,gamFF))
								outcsv.write('({:>3},;({:>3},;({:>5},;({:.2f},;({:.0f},;({:.0f},;({:.0f},;({:.0f},;({:.3f},;({:.2f},;({:.2f},;({:.2f},;({:.2f},\n'.format(k,j,l,dT,dCa,np.sum(KayaFF),dCo,dCl,gain,gamA,gamO,gamL,gamFF))
                                                        #output.write('{:>3} {:>3} {:>5} {:0<4.3} {:0<6.5} {:0<6.5} {:0<5.4} {:0<5.4} {: .4f} {: 9.4f} {: 6.2f} {: 6.2f} {: 7.2f}\n'.format(k,j,l,dT,dCa,np.sum(KayaFF),dCo,dCl,gain,gamA,gamO,gamL,gamFF))
                                                        #outcsv.write('{:>3},{:>3},{:>5},{:0<4.3},{:0<6.5},{:0<6.5},{:0<5.4},{:0<5.4},{: .4f},{: 9.4f},{: 6.2f},{: 6.2f},{: 7.2f}\n'.format(k,j,l,dT,dCa,np.sum(KayaFF),dCo,dCl,gain,gamA,gamO,gamL,gamFF))

        output.close()
        outcsv.close()
        outputBGC.close()
        outcsv.close()
	outcsvECON.close()
	otherEconOutput.close()
	outBGCParams.close()
