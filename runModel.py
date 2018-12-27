import modelFunctions as mf
import numpy as np
import sys
import matplotlib.pylab as py
import matplotlib.pyplot as plt
import inspect

reload(mf)
#**************************** Parse arguments for run ********************#
scenario = 'all' #run through all model scenarios. 
nyears = 300 #run scenario from 1800-2100


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

#Parameters from the parameterization of NPP versus temperature
aland = -0.0022648
bland = 0.04661
cland = 0.789

#Parameters for the impulse function
K0 = -0.4194842
K1 = 0.9330434
K2 = -0.1868783
K3 = 0.01516102


rcpMidyearConc = np.genfromtxt("InputData/RCP85_MIDYR_CONC.DAT",skip_header=35)
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

#read in initial parameters from values output from fitted model
x = np.genfromtxt('InputData/optimizedParams.txt',delimiter=',')
aland,bland,K0,K1,K2,K3,Q10,beta,NPPinitial,taub,Revellef,depthm,tautempfac,cat,kgas,initialCappm=x



if (scenario=='1pct'):
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
				print inspect.getframeinfo(inspect.currentframe()).lineno, 'betaA,betaL,betaO',betaA,betaL,betaO
			elif j == 'NFC':
				print inspect.getframeinfo(inspect.currentframe()).lineno, 'dCo,dCl,dT,dCa =',dCo,dCl,dT,dCa
				gamO = (dCo-betaO*dCa)/dT
				gamL = (dCl-betaL*dCa)/dT
				gamA = -(gamO+gamL)
				C_uc = dCa
				Ff_nfc = np.sum(Ff)
				gain = (Ff_bgc - Ff_nfc)/Ff_bgc
				gain2 = (gamA*dT)/(dCa*(2.12+betaL + betaO))
				gain1 = 1 - C_buc/C_uc
				print inspect.getframeinfo(inspect.currentframe()).lineno, 'gain1,gain2,gamA,gamL,gamO = ',gain,gain2,gamA,gamL,gamO
                        	print inspect.getframeinfo(inspect.currentframe()).lineno, 'dT,dCl,dCo,dCa',dT,dCl,dCo,dCa



elif (scenario=='all'):
        output = open('outputTable2.dat','w')
        outputBGC = open('BGCoutputTable2.dat','w')
        outcsv = open('outputTable2.txt','w')
        outcsvBGC = open('BGCoutputTable2.txt','w')
        outcsvBGC2 = open('BGCoutputTable.txt','w')
        outBGCParams = open('BGCparams.csv','w')
        outParams = open('outParams.csv','w')
        outcsvECON = open('econOutputTable.csv','w')
        econ_21stCentury = open('econ_21stCentury.csv','w')
	otherEconOutput = open('econAllVals.csv','w')
        output.write('{:<3} {} {:<5} {:<4} {:<6} {:<6} {:<5} {:<5} {:<7} {:<9} {:<6} {:<6} {:<7}\n'.format('RCP','Run','Level','dT','dCa','dFF','dCo','dCl','Gain','GamA','GamO','GamL','GamFF'))
        tempAll=[]
	gainCalc = {}
        for k in ['85']:
                for j in ['BGC','LFC','OFC','NFC','PFC','GFC','EFC','DFC','AFC','FC']:
			gainCalc[j] = {}
			for l in ['max','base','min']:
				gainCalc[j][l] = {}
                                for m in [True,False]:
					gainCalc[j][l][m] = {}
                                        if (l!='base' and (j in ['LFC','OFC']) ):
                                                continue
                                        print inspect.getframeinfo(inspect.currentframe()).lineno, k,j,l,m
					dCo,dCl,dCa,dT,Tair,gdp,old_gdp,KayaFF,Ff,fonet,fbnet,Ca,date,pop,pop_adj,KayaE,KayaCE,ee,ED,EE = mf.runModel(scenario,k,j,l,m,nyears,aland, bland, K0,K1,K2,K3,Q10,beta,NPPinitial,taub,Revellef,depthm,tautempfac,cat,kgas,initialCappm)
					KayaE = ED/(old_gdp*pop)
					oldKayaE = ee/(old_gdp*pop)
					oldKayaCE = Ff/ee
                                        CO2m = dCa/2.12
					gainCalc[j][l][m] = Ca-Ca[0]
					print inspect.getframeinfo(inspect.currentframe()).lineno, 'dCo,dCl,dT,dCa,Ff,KayaFF =',dCo,dCl,dT,dCa,np.sum(Ff),np.sum(KayaFF)
                                        if j == 'BGC':
                                                C_buc = dCa
                                                betaO=dCo/CO2m
                                                betaL=dCl/CO2m
                                                betaA2 = dCa/CO2m
                                                betaA = -(betaO+betaL)
                                                print inspect.getframeinfo(inspect.currentframe()).lineno, 'betaA,betaA2',betaA,betaA2
						gainTS = (Ca-Ca[0])*0
                                                if not m:
							np.savetxt("OutputData/gainCumulative"+k+j+l+".csv",gainTS,delimiter=",")
                                                if m:
							np.savetxt("OutputData/gainCumulativeBGC"+k+j+l+".csv",gainTS,delimiter=",")
							outParams.write('{:>3},{:>3},{:>5},{:0<4.5},{:0<6.5},{:0<6.5},{:0<5.4},{:0<5.4},{},{},{},{},{}\n'.format(k,j,l,dT,dCa,np.sum(KayaFF),dCo,dCl,0,0,0,0,0))
                                                        outputBGC.write('{:>3} {:>3} {:>5} {:0<4.3} {:0<6.5} {:0<6.5} {:0<5.4} {:0<5.4} {:-<7} {:-<9} {:-<6} {:->6} {:->7}\n'.format(k,j,l,dT,dCa,np.sum(Ff),dCo,dCl,'-','-','-','-','-'))
                                                        outcsv.write('{:>3};{:>3};{:>5};{:.2f};{:.0f};{:.0f};{:.0f};{:.0f};{};{};{};{};{}\n'.format(k,j,l,dT,dCa,np.sum(Ff),dCo,dCl,'---','---','---','---','---'))
                                                        outcsvBGC2.write('{:>3};{:>3};{:>5};{:0<4.3};{:0<6.5};{:0<6.5};{:0<5.4};{:0<5.4};{};{};{};{};{}\n'.format(k,j,l,dT,dCa,np.sum(Ff),dCo,dCl,'---','---','---','---','---'))
                                                        outBGCParams.write('{:>3},{:>3},{:>5},{:0<4.5},{:0<6.5},{:0<6.5},{:0<5.4},{:0<5.4},{},{},{},{},{}\n'.format(k,j,l,dT,dCa,np.sum(Ff),dCo,dCl,0,0,0,0,0))
							if (l=='base'):
								outcsvBGC.write('{:>3},;{:>3},;{:>5},;{:.2f},;{:.0f},;{:.0f},;{:.0f},;{:.0f},;{},;{},;{},;{},;{},\n'.format(k,j,l,dT,dCa,np.sum(KayaFF[199:299]),dCo,dCl,'---','---','---','---','---'))
							elif (l=='min'):
								outcsvBGC.write('{:>3});{:>3});{:>5});{:.2f});{:.0f});{:.0f});{:.0f});{:.0f});{});{});{});{});{})\n'.format(k,j,l,dT,dCa,np.sum(KayaFF[199:299]),dCo,dCl,'---','---','---','---','---'))
							elif (l=='max'):
								outcsvBGC.write('({:>3},;({:>3},;({:>5},;({:.2f},;({:.0f},;({:.0f},;({:.0f},;({:.0f},;({},;({},;({},;({},;({},\n'.format(k,j,l,dT,dCa,np.sum(KayaFF[199:299]),dCo,dCl,'---','---','---','---','---'))
                                        elif j in ['NFC','LFC','OFC']:
                                                gamO = (dCo-betaO*CO2m)/dT
                                                gamL = (dCl-betaL*CO2m)/dT
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
						gainTS = 1 - divZero(gainCalc['BGC'][l][m],gainCalc[j][l][m])
						gainTS[0] = 0
						gain = gainTS[-1]
                                                print 'gain1,gamA,gamL,gamO,gamFF = ',gain,gamA,gamL,gamO,gamFF
						print 'RCP8.5 (1850-2100): dT,dTm,dCl,dCo,dCa',Tair[299]-Tair[61],np.mean(Tair[281:299])-np.mean(Tair[186:205]),np.sum(fbnet[50:]),np.sum(fonet[50:]),Ca[-1]/2.12

        					print 'Historical (1850-2005): dT,dCl,dCo,dCa',np.mean(Tair[186:205])-np.mean(Tair[51:79]),np.sum(fbnet[50:205]),np.sum(fonet[50:205]),Ca[205]/2.12
                                                if not m:
							np.savetxt("OutputData/gainCumulative"+k+j+l+".csv",gainTS,delimiter=",")
							output.write('{:>3} {:>3} {:>5} {:0<4.3} {:0<6.5} {:0<6.5} {:0<5.4} {:0<5.4} {: .4f} {: 9.4f} {: 6.2f} {: 6.2f} {: 7.2f}\n'.format(k,j,l,dT,dCa,np.sum(KayaFF),dCo,dCl,gain,gamA,gamO,gamL,gamFF))
							otherEconOutput.write('{:>3},{:>3},{:>5},{: .4f},{: .4f},{: .4f},{: .4f},{: .4f},{: .4f}\n'.format(k,j,l,dT,pop_adj[-1],gdp[-1],ED[-1],EE[-1],np.sum(KayaFF)))
	
							outParams.write('{:>3},{:>3},{:>5},{:0<4.5},{:0<6.5},{:0<6.5},{:0<5.4},{:0<5.4},{: .2f},{: 9.2f},{: 6.2f},{: 6.2f},{: 7.2f}\n'.format(k,j,l,dT,dCa,np.sum(KayaFF),dCo,dCl,gain,gamA,gamO,gamL,gamFF))
							outcsv.write('{:>3};{:>3};{:>5};{:.2f};{:.0f};{:.0f};{:.0f};{:.0f};{:.3f};{:.2f};{:.2f};{:.2f};{:.2f}\n'.format(k,j,l,dT,dCa,np.sum(KayaFF),dCo,dCl,gain,gamA,gamO,gamL,gamFF))
                                                else:
							np.savetxt("OutputData/gainCumulativeBGC"+k+j+l+".csv",gainTS,delimiter=",")
                                                        outputBGC.write('{:>3} {:>3} {:>5} {:0<4.3} {:0<6.5} {:0<6.5} {:0<5.4} {:0<5.4} {: .4f} {: 9.4f} {: 6.2f} {: 6.2f} {: 7.2f}\n'.format(k,j,l,dT,dCa,np.sum(KayaFF),dCo,dCl,gain,gamA,gamO,gamL,gamFF))
							outBGCParams.write('{:>3},{:>3},{:>5},{:0<4.5},{:0<6.5},{:0<6.5},{:0<5.4},{:0<5.4},{: .4f},{: 9.4f},{: 6.2f},{: 6.2f},{: 7.2f}\n'.format(k,j,l,dT,dCa,np.sum(KayaFF),dCo,dCl,gain,gamA,gamO,gamL,gamFF))
                                                        if (l=='base'):
								outcsvBGC.write('{:>3},;{:>3},;{:>5},;{:0<4.3},;{:0<6.5},;{:0<6.5},;{:0<5.4},;{:0<5.4},;{: .4f},;{: 9.4f},;{: 6.2f},;{: 6.2f},;{: 7.2f},\n'.format(k,j,l,dT,dCa,np.sum(KayaFF[199:299]),dCo,dCl,gain,gamA,gamO,gamL,gamFF))
							elif (l=='min'):
								outcsvBGC.write('{:>3});{:>3});{:>5});{:0<4.3});{:0<6.5});{:0<6.5});{:0<5.4});{:0<5.4});{: .4f});{: 9.4f});{: 6.2f});{: 6.2f});{: 7.2f})\n'.format(k,j,l,dT,dCa,np.sum(KayaFF[199:299]),dCo,dCl,gain,gamA,gamO,gamL,gamFF))
							elif (l=='max'):
								outcsvBGC.write('({:>3},;({:>3},;({:>5},;({:0<4.3},;({:0<6.5},;({:0<6.5},;({:0<5.4},;({:0<5.4},;({: .4f},;({: 9.4f},;({: 6.2f},;({: 6.2f},;({: 7.2f},\n'.format(k,j,l,dT,dCa,np.sum(KayaFF[199:299]),dCo,dCl,gain,gamA,gamO,gamL,gamFF))
					
					else:
                                                C_fc = dCa
						Call_fc = Ca-Ca[0]
                                                gamO = (dCo-betaO*CO2m)/dT
                                                gamL = (dCl-betaL*CO2m)/dT
                                                gamFF = (np.sum(KayaFF)-np.sum(Ff))/dT
                                                gamA = gamFF-(gamO+gamL)
                                                gamA2 = (dCa-betaA*CO2m)/dT
                                                alph = (Tair[-1]-Tair[0])/(Ca[-1]*(1/2.13)-Ca[0]*(1/2.13)) # don't forget to convert Ca to ppm!
                                                gain = 1 - C_buc/C_fc
						gainTS = 1 - divZero(gainCalc['BGC'][l][m],gainCalc[j][l][m])
						gainTS[0] = 0
						gain = gainTS[-1]
                                                print inspect.getframeinfo(inspect.currentframe()).lineno, 'gain1,gamA,gamL,gamO,gamFF = ',gain,gamA,gamL,gamO,gamFF
                                                if m:
							np.savetxt("OutputData/gainCumulativeBGC"+k+j+l+".csv",gainTS,delimiter=",")
							if (l=='base'):
								outcsvBGC.write('{:>3},;{:>3},;{:>5},;{:.2f},;{:.0f},;{:.0f},;{:.0f},;{:.0f},;{:.3f},;{:.2f},;{:.2f},;{:.2f},;{:.2f},\n'.format(k,j,l,dT,dCa,np.sum(KayaFF[199:299]),dCo,dCl,gain,gamA,gamO,gamL,gamFF))
							elif (l=='min'):
								outcsvBGC.write('{:>3});{:>3});{:>5});{:.2f});{:.0f});{:.0f});{:.0f});{:.0f});{:.3f});{:.2f});{:.2f});{:.2f});{:.2f})\n'.format(k,j,l,dT,dCa,np.sum(KayaFF[199:299]),dCo,dCl,gain,gamA,gamO,gamL,gamFF))
							elif (l=='max'):
								outcsvBGC.write('({:>3},;({:>3},;({:>5},;({:.2f},;({:.0f},;({:.0f},;({:.0f},;({:.0f},;({:.3f},;({:.2f},;({:.2f},;({:.2f},;({:.2f},\n'.format(k,j,l,dT,dCa,np.sum(KayaFF[199:299]),dCo,dCl,gain,gamA,gamO,gamL,gamFF))
							if (j=='FC'):
								outcsvECON.write('{},{},{},{},{},{}\n'.format(j,l,m,'','',''))
								outcsvECON.write('{:0<.2},{:0<.2},{:0<.2},{:0<.2},{:0<.2},{:0<.2}\n'.format(pop[230],pop_adj[230],100*-1*(pop[230]-pop_adj[230])/pop[230],pop[-1],pop_adj[-1],100*-1*((pop[-1]-pop[-100])-(pop_adj[-1]-pop_adj[-100]))/(pop[-1]-pop[-100])))
								outcsvECON.write('{:0<.2},{:0<.2},{:0<.2},{:0<.2},{:0<.2},{:0<.2}\n'.format(old_gdp[230],gdp[230],100*-1*(old_gdp[230]-gdp[230])/old_gdp[230],old_gdp[-1]-old_gdp[-100],gdp[-1]-gdp[-100],100*-1*((old_gdp[-1]-old_gdp[-100])-(gdp[-1]-gdp[-100]))/(old_gdp[-1]-old_gdp[-100])))
								outcsvECON.write('{:0<.2},{:0<.2},{:0<.2},{:0<.2},{:0<.2},{:0<.2}\n'.format(oldKayaE[230],KayaE[230],100*-1*(oldKayaE[230]-KayaE[230])/oldKayaE[230],oldKayaE[-1],KayaE[-1],100*-1*((oldKayaE[-1]-oldKayaE[-100])-(KayaE[-1]-KayaE[-100]))/(oldKayaE[-1]-oldKayaE[-100])))
								outcsvECON.write('{:0<.2},{:0<.2},{:0<.2},{:0<.2},{:0<.2},{:0<.2}\n'.format(oldKayaCE[230],KayaCE[230],100*-1*(oldKayaCE[230]-KayaCE[230])/oldKayaCE[230],oldKayaCE[-1]-oldKayaCE[-100],KayaCE[-1]-oldKayaCE[-100],100*-1*((oldKayaCE[-1]-oldKayaCE[-100])-(KayaCE[-1]-KayaCE[-100]))/(oldKayaCE[-1]-oldKayaCE[-100])))
								outcsvECON.write('{:0<.2},{:0<.2},{:0<.2},{:0<.2},{:0<.2},{:0<.2}\n'.format(Ff[230],KayaFF[230],100*-1*(Ff[230]-KayaFF[230])/Ff[230],Ff[-1],KayaFF[-1],100*-1*(Ff[-1]-KayaFF[-1])/Ff[-1]))
								econ_21stCentury.write('{},{},{},{},{},{}\n'.format('','','','','',''))
								
								econ_21stCentury.write('{},{},{},{},{},{}\n'.format(j,l,m,'','',''))
								econ_21stCentury.write('{},{},{},{},{},{}\n'.format('','Base (0)','Base (2000)','FC (2000)','Base (2100)','FC (2100)'))
								econ_21stCentury.write('{},{:0<.6},{:0<.6},{:0<.6},{:0<.6},{:0<.6}\n'.format('pop',pop[0],pop[200],pop_adj[200],pop[-1],pop_adj[-1]))
								econ_21stCentury.write('{},{:0<.6},{:0<.6},{:0<.6},{:0<.6},{:0<.6}\n'.format('GDP',old_gdp[0],old_gdp[200],gdp[200],old_gdp[-1],gdp[-1]))
								econ_21stCentury.write('{},{:0<.6},{:0<.6},{:0<.6},{:0<.6},{:0<.6}\n'.format('KayaE',oldKayaE[0],oldKayaE[200],KayaE[200],oldKayaE[-1],KayaE[-1]))
								econ_21stCentury.write('{},{:0<.6},{:0<.6},{:0<.6},{:0<.6},{:0<.6}\n'.format('KayaCE',oldKayaCE[0],oldKayaCE[200],KayaCE[200],oldKayaCE[-1],KayaCE[-1]))
								econ_21stCentury.write('{},{:0<.6},{:0<.6},{:0<.6},{:0<.6},{:0<.6}\n'.format('Ff',Ff[0],Ff[200],KayaFF[200],Ff[-1],KayaFF[-1]))
								econ_21stCentury.write('{},{},{},{},{},{}\n'.format('','','','','',''))
								
						
                                                elif not m:
							np.savetxt("OutputData/gainCumulative"+k+j+l+".csv",gainTS,delimiter=",")
                                                        outputBGC.write('{:>3} {:>3} {:>5} {:0<4.3} {:0<6.5} {:0<6.5} {:0<5.4} {:0<5.4} {: .4f} {: 9.4f} {: 6.2f} {: 6.2f} {: 7.2f}\n'.format(k,j,l,dT,dCa,np.sum(KayaFF),dCo,dCl,gain,gamA,gamO,gamL,gamFF))
                                                        outParams.write('{:>3},{:>3},{:>5},{:0<4.5},{:0<6.5},{:0<6.5},{:0<5.4},{:0<5.4},{: .4f},{: 9.4f},{: 6.2f},{: 6.2f},{: 7.2f}\n'.format(k,j,l,dT,dCa,np.sum(KayaFF),dCo,dCl,gain,gamA,gamO,gamL,gamFF))
							otherEconOutput.write('{:>3},{:>3},{:>5},{: .4f},{: .4f},{: .4f},{: .4f},{: .4f},{: .4f}\n'.format(k,j,l,dT,pop_adj[-1],gdp[-1],ED[-1],EE[-1],np.sum(KayaFF)))
                                                        if (l=='base'):
								outcsv.write('{:>3},;{:>3},;{:>5},;{:.2f},;{:.0f},;{:.0f},;{:.0f},;{:.0f},;{:.3f},;{:.2f},;{:.2f},;{:.2f},;{:.2f},\n'.format(k,j,l,dT,dCa,np.sum(KayaFF),dCo,dCl,gain,gamA,gamO,gamL,gamFF))
							elif (l=='min'):
								outcsv.write('{:>3});{:>3});{:>5});{:.2f});{:.0f});{:.0f});{:.0f});{:.0f});{:.3f});{:.2f});{:.2f});{:.2f});{:.2f})\n'.format(k,j,l,dT,dCa,np.sum(KayaFF),dCo,dCl,gain,gamA,gamO,gamL,gamFF))
							elif (l=='max'):
								outcsv.write('({:>3},;({:>3},;({:>5},;({:.2f},;({:.0f},;({:.0f},;({:.0f},;({:.0f},;({:.3f},;({:.2f},;({:.2f},;({:.2f},;({:.2f},\n'.format(k,j,l,dT,dCa,np.sum(KayaFF),dCo,dCl,gain,gamA,gamO,gamL,gamFF))

        output.close()
        outcsv.close()
        outputBGC.close()
        outcsv.close()
	outcsvECON.close()
	otherEconOutput.close()
	outBGCParams.close()
	econ_21stCentury.close()
