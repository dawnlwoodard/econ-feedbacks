import numpy as np
import time
import matplotlib.pylab as py
import matplotlib.pyplot as plt
import sys 
import csv

def runModel(scenario,rcp,coupling,level,bgcBase,nyears,aland, bland, K0,K1,K2,K3,Q10,beta,NPPinitial,taub,Revellef,depthm,tautempfac,cat,kgas,initialCappm,dic=0.00191,deptho=3800.0,deltat=1,alpha=0.0062,mwc = 12.01):
	if level=='min':
		tautempfac=tautempfac*(4.0/3.0)
		Q10 = 1.3
		beta = 0.6
		taub = 10
		depthm = 75
		Revellef = 11
	if level == 'max':
		tautempfac=tautempfac*(2.0/3.0)

		Q10 = 1.0
		beta = 0.7
		taub = 20
		Revellef = 12.0
		depthm = 125
			
	print '\n****************\n',coupling,level,'\n***************\n'			
			
	deltaF	     = 5.35*np.log(1140./280)
	oceanarea    = 3.61e14		 # area of the ocean surface
	kb = 1.0/taub
	header=1
	initialCaPgC = initialCappm*2.1  # initial atmospheric mass in Pg C
	tb = 1#.8
	taub_o = taub
	changeDepth = True
	diffusionJR = True
	kfrac = 1.0/20
	if diffusionJR:
		Kdeep_o = 2500
		nboxes = 25
		deptho = 3800.0
		hboxes = [(deptho-depthm)/float(nboxes)]*nboxes
        else:
                nboxes = 1
                hboxes=[(deptho-depthm)]
	Ca	     = np.zeros(nyears) # Carbon pool in atmosphere
	Fao	     = np.zeros(nyears) # flux from atmosphere to ocean
	Foa	     = np.zeros(nyears) # flux from ocean to atmosphere
	Fmd	     = np.zeros(nyears) # fluxe from mixed layer to deep ocean
	Fdm	     = np.zeros(nyears) # flux from deep ocean to mixed layer
	Com	     = np.zeros(nyears) # Carbon pool in mixed year
	Cod	     = np.zeros((nboxes,nyears)) # Carbon pool in deep ocean
	Tair	     = np.zeros(nyears) # Air temperature
	Fonet	     = np.zeros(nyears) # net ocean flux (positive is to atm; Pg C/yr)
	Goa	     = np.zeros(nyears) # gross ocean flux (positive is to atm Pg C/yr)
	Gao	     = np.zeros(nyears) # net ocean flux (positive is to atm Pg C/yr)
	Gmd	     = np.zeros(nyears) # net ocean flux (positive is to atm Pg C/yr)
	Gdm	     = np.zeros(nyears) # net ocean flux (positive is to atm Pg C/yr)


	Cb	 = np.zeros(nyears)	# Carbon pool in land biosphere
	NPP	 = np.zeros(nyears)	# Flux into land biosphere
	Rh	 = np.zeros(nyears)	# Flux out of the land biosphere
	Fb	 = np.zeros(nyears)	# Net land flux (Rh-NPP) Positive is to atm.
	RF	 = np.zeros(nyears)	# Radiative forcing
	dTair	 = np.zeros(nyears)	#
	taub     = np.zeros(nyears)     #

	GDP      = np.zeros(nyears)
	pop_adj  = np.zeros(nyears)
	KayaC    = np.zeros(nyears)
	KayaFF   = np.zeros(nyears)
	KayaE    = np.zeros(nyears)
	KayaCE   = np.zeros(nyears)
	KayaG    = np.zeros(nyears)
	EE       = np.zeros(nyears)
	ED       = np.zeros(nyears)
	taubLUCfrac2       = np.zeros(nyears)
	opec_adj         = np.zeros(nyears)
	pe_adj           = np.zeros(nyears)
	pe_econly                = np.zeros(nyears)
	pe_transonly             = np.zeros(nyears)
	cpec_adj         = np.zeros(nyears)
	ngpec_adj        = np.zeros(nyears)
	transpe_adj     = np.zeros(nyears)
	ed_tempEffects = {}
	ed_tempEffects['demand'] = np.zeros(nyears)
	ed_tempEffects['temp'] = np.zeros(nyears)


	t1 = range(1,nyears+1);
	t1 = np.array(t1)
	t = np.log(t1)
	y1 = K0 + K1*t + K2*t**2 + K3*t**3
	y2 = np.exp(y1)
	tempinc = np.zeros(nyears)
	for i in range(1,nyears):
		tempinc[i] = y2[i]-y2[i-1]
	tempinc[0] = y2[0]
	tempinc = tempinc/deltaF


	#*********************************** Define Energy Relationships **********************#

	pe = np.genfromtxt("InputData/energy85_all.csv",delimiter=",")
	luc = np.genfromtxt("InputData/co2Data.csv",delimiter=",")
	lucYr = luc[:,0]
	co2Old = luc[36:,2]
	luc = luc[:,1]
	luc = luc[35:]
	

	f_cpec = 0.0001314341*pe+0.01309835794
	f_ngpec = 1.641838e-10*pe**3-6.539296e-7*pe**2+0.0008115*pe-0.167908
	f_opec = -1.3539426e-5*pe+0.042859
	f_ngpec[f_ngpec<0] = 0
	f_opec[f_opec<0] = 0
	f_cpec[f_cpec<0] = 0
	f_tpe = 0.0002397144*pe+0.1673878785
	transpe = f_tpe*pe

	cpec = f_cpec*pe #coal primary energy for electricity

	ngpec = f_ngpec*pe

	opec = f_opec*pe
	ffec = opec+ngpec+cpec
	pe_nonff = pe - ffec - transpe
	pDeaths=0


	#READ IN DATA
	def readData(rcp):
		old_gdp = np.genfromtxt("InputData/gdp"+rcp+"_all.csv",delimiter=",")
		old_gdp = np.array(old_gdp).astype(float)
		old_gdp[0] = old_gdp[1]
		dold_gdp = (old_gdp[1:]-old_gdp[:-1])/old_gdp[:-1]
		dold_gdp = np.append(dold_gdp[0],dold_gdp)

		burke_temp = np.genfromtxt("InputData/burke_temp_base.csv",delimiter=",")
		temp_min = np.genfromtxt("InputData/burke_temp_max.csv",delimiter=",")
		temp_max = np.genfromtxt("InputData/burke_temp_min.csv",delimiter=",")


		eFracs = {}
		with open("InputData/eFracs"+rcp+".csv") as infile:
			reader = csv.reader(infile)
			eFracs = {row[0]:row[1:] for row in reader}
		for key in eFracs:
			eFracs[key] = [np.float(k) for k in eFracs[key]]

		eLosses = {}
		with open("InputData/eLosses.csv") as infile:
			reader = csv.reader(infile)
			eLosses = {row[0]:row[1] for row in reader}
		eKeys = eLosses.keys()
		for thing in eKeys:
			eLosses[thing] = float(eLosses[thing])
			eLosses[thing+"_adj"] = 0
		finE = np.genfromtxt("InputData/finalE85.csv",delimiter=",")
		ee = np.genfromtxt("InputData/energy"+rcp+"_all.csv",delimiter=",")
		ff = np.genfromtxt("InputData/Fossilhistandrcp"+rcp+".txt",skip_header=1)
		pop = np.genfromtxt("InputData/pop"+rcp+"_all.csv",delimiter=",")
		return ff[:,0],old_gdp[:-1],burke_temp,temp_min,temp_max,ff[:,1],pop,ee,eFracs,eLosses,finE

	def fitGDP(tcurve,xvals):
		#fit burke tempeature data with a*x^2 + b*x + c model for base, min, and max
		fit = np.polyfit(tcurve[:,0],tcurve[:,1],2)
		temp_fit = [fit[2] + i*fit[1] + fit[0]*i**2 for i in xvals]
		return fit,temp_fit

	def getGDP(fit,Tc,T0,old_gdp_i): #inputs are curve fit for temp-response, current temp, initial temp, old_gdp from GCAM at current timestep

		tdiff = Tc-T0
		dnew_gdp = fit[2]+tdiff*fit[1]+fit[0]*tdiff**2
		dnew_gdp = dnew_gdp*0.01
		new_gdp = old_gdp_i+dnew_gdp*old_gdp_i
		return new_gdp

	def getPop(fit,old_pop_i,old_gdp_i): #inputs are curve fit for temp-response, current temp, initial temp, old_gdp from GCAM at current timestep

		dnew_pop = fit[0]*old_gdp_i+fit[1]
		dnew_pop = dnew_pop*0.01
		new_pop = old_pop_i+dnew_pop*old_pop_i
		return new_pop

	def runOceanDiffu(Ca_o,Com_o,Cod_o,Tair_o,Tair0):
                Fod_n = np.zeros(nboxes)
                Cod_n = np.zeros(nboxes)
                Kdeep = Kdeep_o*(1.0-(Tair_o-Tair0)*tautempfac/100.0)
                minKd = Kdeep*kfrac
                if (changeDepth):
                        KdeepList = [Kdeep - i*(Kdeep-minKd)/float(nboxes) for i in range(nboxes)]
                else:
                        KdeepList = [Kdeep]*nboxes
                Fao_n   = kgas*(Ca_o/2.12)
                Foa_n   = kgas*(Ca[0]/2.12)*(1.0 + Revellef*(Com_o-Com[0])/Com[0])
                for i in range(nboxes):
                        if i==0:
                                k_mj = float(KdeepList[i])/(depthm*(hboxes[i]+depthm)/2.0)
                                k_jm = float(KdeepList[i])/(hboxes[i]*(hboxes[i]+depthm)/2.0)
                                k_jk = float(KdeepList[i])/hboxes[i]**2
                                k_kj = k_jk
                                Fod_n[i] = k_mj*Com_o+k_kj*Cod_o[i+1]-k_jk*Cod_o[i]-k_jm*Cod_o[i]
                        elif i==(nboxes-1):
                                k_ij = float(KdeepList[i])/(hboxes[i-1]*(hboxes[i]+hboxes[i-1])/2.0)
                                k_ji = float(KdeepList[i])/(hboxes[i]*(hboxes[i-1]+hboxes[i])/2.0)
                                Fod_n[i] = k_ij*Cod_o[i-1] - k_ji*Cod_o[i]
                        else:
                                k_ij = float(KdeepList[i])/(hboxes[i-1]*(hboxes[i]+hboxes[i-1])/2.0)
                                k_ji = float(KdeepList[i])/(hboxes[i]*(hboxes[i-1]+hboxes[i])/2.0)
                                k_kj = float(KdeepList[i])/(hboxes[i+1]*(hboxes[i]+hboxes[i+1])/2.0)
                                k_jk = float(KdeepList[i])/(hboxes[i]*(hboxes[i+1]+hboxes[i])/2.0)
                                Fod_n[i] = k_ij*Cod_o[i-1]+k_kj*Cod_o[i+1]-k_jk*Cod_o[i]-k_ji*Cod_o[i]
                        Cod_n[i] = Cod_o[i]+(Fod_n[i])*deltat

                Com_n   = Com_o + (Fao_n - Foa_n - k_mj*Com_o+k_jm*Cod_o[0])*deltat
                
		Gao_n   = Fao_n*mwc*oceanarea/1.0e15
                Goa_n   = Foa_n*mwc*oceanarea/1.0e15
                Fonet_n = Goa_n-Gao_n
                Gmd_n   = (k_mj*Com_o)*mwc*oceanarea/1.0e15
                Gdm_n   = (k_jm*Cod_o[0])*mwc*oceanarea/1.0e15
                return Fonet_n,Com_n,Cod_n



	#Ocean Subroutine
	def runOcean(Ca_o,Com_o,Cod_o,Tair_o,Tair0):
		taumtod = cat*(1.0+(Tair_o-Tair0)*tautempfac/100.0)
		Fao_n   = kgas*(Ca_o/2.1)
		Foa_n   = kgas*(Ca[0]/2.1)*(1.0 + Revellef*(Com_o-Com[0])/Com[0])
		Fmd_n   = Com_o*(1.0/taumtod)
		Fdm_n   = Cod_o*(1.0/taumtod)* (depthm/(deptho-depthm))
		Com_n   = Com_o + (Fao_n - Foa_n + Fdm_n - Fmd_n)*deltat
		Cod_n   = Cod_o + (Fmd_n - Fdm_n)*deltat
		Gao_n   = Fao_n*mwc*oceanarea/1.0e15
		Goa_n   = Foa_n*mwc*oceanarea/1.0e15
		Fonet_n = Goa_n-Gao_n
		return Fonet_n,Com_n,Cod_n

	#Land Subroutine
	def runLand(Ca_o,Tair_o,Cb_o,luc_o,taub_i):
		if (Ca_o > 0) and (Ca[0] > 0):
                        NPP_n  = NPP[0]*(1 + beta*np.log(Ca_o/Ca[0]))*(1+aland*(Tair_o-Tair[0])**2+bland*(Tair_o-Tair[0]))
                else:
                        NPP_n  = NPP[0]*(1 + beta*0)*(1+aland*(Tair_o-Tair[0])**2+bland*(Tair_o-Tair[0]))
		Rh_n	 = (1.0/taub_i)*Q10**((Tair_o -Tair[0])/10.0)*Cb_o
		Cb_n	 = Cb_o + (NPP_n - Rh_n - luc_o)*deltat
		Fb_n	 = Rh_n - NPP_n + luc_o 
		return Fb_n,NPP_n,Rh_n,Cb_n

		
	def mainFun(scenario,rcp,coupling,level,bgcBase):
		# atmosphere initial conditions
		Ca[0]	     = initialCaPgC		    # initial boundary condition
		Ca[0]	     = co2Old[0]*2.12
		# ocean initial conditions
		Fao[0]	     = kgas*(Ca[0]/2.1)		    # initial boundary condition
		Foa[0]	     = Fao[0]			    # initial boundary condition
		Com[0]	     = depthm*dic*1000.0	    # initial boundary condition
		for i in range(nboxes):
			Cod[i,0] = (hboxes[i])*dic*1000.0
		Fmd[0]	     = Com[0]*(1.0/cat)	    # initial boundary condition
		Fdm[0]	     = Fmd[0]			    # initial boundary condition
		taub[0] = taub_o
		Tair[0]      = 15.0			    # initial boundary condition
		Gao[0]	     = Fao[0]*mwc*oceanarea/1.0e15
		Goa[0]	     = Fao[0]*mwc*oceanarea/1.0e15  
		Fonet[0]     = Goa[0]-Gao[0]
		# land initial conditions
		NPP[0]	    = NPPinitial		    # initial boundary condition
		Rh[0]	    = NPPinitial		    # initial boundary condition
		Fb[0]	    = 0.0			    # initial boundary condition
		Cb[0]	    = NPPinitial*taub[0]		    # initial boundary condition
		# air temperature initial conditions
		Tair[0]     = 15	    # initial boundary condition
		dTair	 = np.zeros(nyears)	#
		pDeaths=0
		
		#************************** Run Model ***********************#
		if (scenario=='1pct'):
			Ff = np.zeros(nyears)
			Ca[0] = 285*2.12
			date = range(0,nyears)
			date = [i+(2100-nyears) for i in date]
			tic = time.clock()
			for i in range(1,nyears):
				if (diffusionJR):
					Fonet[i],Com[i],Cod[:,i] = runOceanDiffu(Ca[i-1],Com[i-1],Cod[:,i-1],Tair[i-1],Tair[0])
				else:
					Fonet[i],Com[i],Cod[:,i] = runOcean(Ca[i-1],Com[i-1],Cod[:,i-1],Tair[i-1],Tair[0])
					
				Fb[i],NPP[i],Rh[i],Cb[i] = runLand(Ca[i-1],Tair[i-1],Cb[i-1],0,taub[0])

				Ca[i] = Ca[i-1]*1.01
				Ff[i] = (Ca[i]-Ca[i-1])/deltat - Fonet[i] - Fb[i]
				if (coupling=='BGC'):
					Tair[i] = Tair[0]
				else:
					RF[i] = 5.35*np.log(Ca[i-1]/Ca[0])
					for j in range(1,i+1):
						dTair[i] = dTair[i]+RF[i-j]*tempinc[j-1]
					Tair[i] = Tair[0] + dTair[i]
				

			#************************* Display Results ********************#
			CO2m = Ca/2.1
			return -1*np.sum(Fonet),-1*np.sum(Fb),Ca[-1]-Ca[0],Tair[-1]-Tair[0],Tair,Ff,Fonet,Fb,date
		else:
			print '\n\n****************\n',coupling,level,tautempfac,Q10,'\n***************\n\n'			
			date, old_gdp, burke_temp, temp_min,temp_max,Ff, pop,ee,eFracs,eLosses,finE = readData(rcp)

			old_gdp = old_gdp[36:]
			pop = pop[36:]
			pop = pop*1e6
			pGrowth = pop[1:] - pop[:-1]
			pop_adj = np.copy(pop)
			Ff = Ff[36:]
			date = range(0,nyears)
			date = [i+(2100-nyears) for i in date]
			GDP[0] = old_gdp[0]
			EE[0] = ee[0]
			ED[0] = ee[0]

			ec_frac = 0.1333353474*np.log(pe) - 0.4276365585
			pec = ec_frac*pe
			nonec = pe-pec



			KayaC = Ff/old_gdp
			KayaFF[0] = Ff[0]
			#fit data from Burke et al to use for temperature response of GDP
			if (level=='base'):
				fit, temp_fit = fitGDP(burke_temp,np.arange(1,5,0.1))
			elif (level=='max'):
				fit, temp_fit = fitGDP(temp_max,np.arange(1,5,0.1))
			else: #level == 'min'
				fit, temp_fit = fitGDP(temp_min,np.arange(1,5,0.1))

			#energy demand relationship
			ed_base = 0.04
			ed_high = 0.08
			ed_low = 0.005


			#************************** Run Model ***********************#
			f = open('taubLUCData'+rcp+coupling+level+'.csv','w')
			f.write(rcp+coupling+level+'\n')
	 		for i in range(1,nyears):
			  taubLUCfrac = (Cb[0]-tb*np.sum(luc[:i]))/Cb[0]
			  taubLUCfrac2[i] = (Cb[0]-tb*np.sum(luc[:i]))/Cb[0]
			
			  f.write(str(i)+','+str(taubLUCfrac)+'\n')
			  taub[i] = taub[0]*(Cb[0]-tb*np.sum(luc[:i]))/Cb[0]
			  
			  if (coupling in ["FC","NFC"]):
				  Fonet[i],Com[i],Cod[:,i] = runOceanDiffu(Ca[i-1],Com[i-1],Cod[:,i-1],Tair[i-1],Tair[0])
				  Fb[i],NPP[i],Rh[i],Cb[i] = runLand(Ca[i-1],Tair[i-1],Cb[i-1],luc[i],taub[i])
			  elif (coupling in ["LFC"]):
				  Fonet[i],Com[i],Cod[:,i] = runOceanDiffu(Ca[i-1],Com[i-1],Cod[:,i-1],Tair[0],Tair[0])
				  Fb[i],NPP[i],Rh[i],Cb[i] = runLand(Ca[i-1],Tair[i-1],Cb[i-1],luc[i],taub[i])

			  elif (coupling in ["OFC"]):
				  Fonet[i],Com[i],Cod[:,i] = runOceanDiffu(Ca[i-1],Com[i-1],Cod[:,i-1],Tair[i-1],Tair[0])
				  Fb[i],NPP[i],Rh[i],Cb[i] = runLand(Ca[i-1],Tair[0],Cb[i-1],luc[i],taub[i])

			  elif (coupling in ["BGC"]):
				  Fonet[i],Com[i],Cod[:,i] = runOceanDiffu(Ca[i-1],Com[i-1],Cod[:,i-1],Tair[0],Tair[0])
				  Fb[i],NPP[i],Rh[i],Cb[i] = runLand(Ca[i-1],Tair[0],Cb[i-1],luc[i],taub[i])
			  
			  else: #coupling is "FC","BGC","NFC"
				  if bgcBase:
					Fb[i],NPP[i],Rh[i],Cb[i] = runLand(Ca[i-1],Tair[0],Cb[i-1],luc[i],taub[i])
				  	Fonet[i],Com[i],Cod[:,i] = runOceanDiffu(Ca[i-1],Com[i-1],Cod[:,i-1],Tair[0],Tair[0])
				  else:
					Fb[i],NPP[i],Rh[i],Cb[i] = runLand(Ca[i-1],Tair[i-1],Cb[i-1],luc[i],taub[i])
				  	Fonet[i],Com[i],Cod[:,i] = runOceanDiffu(Ca[i-1],Com[i-1],Cod[:,i-1],Tair[i-1],Tair[0])
			  RF[i] = 5.35*np.log(Ca[i-1]/Ca[0])
			  for j in range(1,i+1):
				dTair[i] = dTair[i]+RF[i-j]*tempinc[j-1]
			  Tair[i] = Tair[0] + dTair[i]
			  
			  if (i<195):
				newDeaths = 0
			  else:	
				newDeaths = 199283*(Tair[i]-Tair[195])
			  if (level=='max'):
			  	if ((Tair[i]-Tair[0]) > 0.8):
					GDP[i] = getGDP(fit,Tair[i],Tair[0],old_gdp[i])
				else:
					GDP[i] = old_gdp[i]
				pDeaths = pDeaths + 1.5*newDeaths
				pop_adj[i] = pop[i]-pDeaths
			  elif (level=='base'):
			  	if ((Tair[i]-Tair[0]) > 0.8):
					GDP[i] = getGDP(fit,Tair[i],Tair[0],old_gdp[i])
				else:
					GDP[i] = old_gdp[i]
				pDeaths = pDeaths + 1*newDeaths
				pop_adj[i] = pop[i]-pDeaths
			  elif (level=='min'):
				GDP[i] = old_gdp[i]-old_gdp[i]*0.00236*(Tair[i]-Tair[0])**2
				pDeaths = pDeaths + 0.5*newDeaths
				pop_adj[i] = pop[i]-pDeaths
			  tdiff = Tair[i] - Tair[0]
			  tdiff_mdrn = Tair[i] - Tair[200]
			  eExtra = 3.36*1e3*(-0.49*(tdiff-0.613)**3+3.57*tdiff**2-4.344*tdiff+0.014)*1
			  ed_tempEffects['demand'][i] = (-0.472*(tdiff_mdrn)**3+4.34*tdiff_mdrn**2-4.65*tdiff_mdrn)*1.56
			  ed_tempEffects['temp'][i] = tdiff_mdrn 
			  if (level == 'base'):
				if i>=200:
					ED[i] = ee[i]+(-0.472*(tdiff_mdrn)**3+4.34*tdiff_mdrn**2-4.65*tdiff_mdrn)*1.56
				else:
					ED[i] = ee[i]
				transpe_adj[i] = transpe[i]*.60*(1+0.003*(Tair[i]-Tair[0]))+transpe[i]*0.22*(1+0.01*(Tair[i]-Tair[0]))+transpe[i]*0.08+transpe[i]*0.04*(1+0.00373*(Tair[i]-Tair[0]))+transpe[i]*0.02*(1-0.0043*(Tair[i]-Tair[0]))+transpe[i]*0.04
				opec_adj[i] = opec[i]*(1-0.0075*(Tair[i]-Tair[0]))*(1-0.0087*(Tair[i]-Tair[0]))
				cpec_adj[i] = cpec[i]*(1-0.0105*(Tair[i]-Tair[0]))*(1-0.0087*(Tair[i]-Tair[0]))
				ngpec_adj[i] = ngpec[i]*(1-0.0075*(Tair[i]-Tair[0]))*(1-0.0087*(Tair[i]-Tair[0]))

			  elif (level == 'max'):
				if i>=200:
					ED[i] = ee[i]+(tdiff_mdrn*3.55+(-0.472*(tdiff_mdrn)**3+4.34*tdiff_mdrn**2-4.65*tdiff_mdrn)*1.56)
				else:
					ED[i] = ee[i]
				transpe_adj[i] = transpe[i]*(1+0.008*(Tair[i]-Tair[0]))
				opec_adj[i] = opec[i]*(1-0.003*(Tair[i]-Tair[0]))*(1-0.0033*(Tair[i]-Tair[0]))
				cpec_adj[i] = cpec[i]*(1-0.006*(Tair[i]-Tair[0]))*(1-0.0033*(Tair[i]-Tair[0]))
				ngpec_adj[i] = ngpec[i]*(1-0.003*(Tair[i]-Tair[0]))*(1-0.0033*(Tair[i]-Tair[0]))
			  else:
				if i>=200:
					ED[i] = ee[i]+(-tdiff_mdrn*3.55+(-0.472*(tdiff_mdrn)**3+4.34*tdiff_mdrn**2-4.65*tdiff_mdrn)*1.56)
				else:
					ED[i] = ee[i]
				transpe_adj[i] = transpe[i]
				opec_adj[i] = opec[i]*(1-0.01*(Tair[i]-Tair[0]))*(1-0.014*(Tair[i]-Tair[0]))
				cpec_adj[i] = cpec[i]*(1-0.015*(Tair[i]-Tair[0]))*(1-0.014*(Tair[i]-Tair[0]))
				ngpec_adj[i] = ngpec[i]*(1-0.01*(Tair[i]-Tair[0]))*(1-0.014*(Tair[i]-Tair[0]))
			  pe_adj[i]=opec_adj[i] + cpec_adj[i] + ngpec_adj[i] + transpe_adj[i] + pe_nonff[i]
			  pe_econly[i]=opec_adj[i] + cpec_adj[i] + ngpec_adj[i] + transpe[i] + pe_nonff[i]
			  pe_transonly[i]=opec[i] + cpec[i] + ngpec[i] + transpe_adj[i] + pe_nonff[i]

			  EE[i] = pe_adj[i]
			  if (coupling=="GFC"):
				KayaE[i] = (ee[i]/(old_gdp[i]*pop[i]))
				KayaCE[i] = Ff[i]/ee[i]
				KayaFF[i] = pop[i]*GDP[i]*KayaE[i]*KayaCE[i] #GDP on, E off
			  elif (coupling=="EFC"):
				KayaE[i] = (ee[i]/(old_gdp[i]*pop[i]))
				KayaCE[i] = Ff[i]/EE[i]
				KayaFF[i] = pop[i]*old_gdp[i]*KayaE[i]*KayaCE[i] #GDP off, E on
			  elif (coupling=='DFC'):
				KayaE[i] = (ED[i]/(old_gdp[i]*pop[i]))
				KayaCE[i] = Ff[i]/ee[i]
				KayaFF[i] = pop[i]*old_gdp[i]*KayaE[i]*KayaCE[i] #GDP off, E off, ED on

			  elif (coupling=='PFC'):
				KayaE[i] = (ee[i]/(old_gdp[i]*pop[i]))
				KayaCE[i] = Ff[i]/ee[i]
				KayaG[i] = old_gdp[i]
				KayaFF[i] = pop_adj[i]*KayaG[i]*KayaE[i]*KayaCE[i]
			  elif (coupling in ['AFC','FC']): 
				KayaE[i] = (ee[i]/(old_gdp[i]*pop[i]))
				KayaCE[i] = Ff[i]/EE[i]
				KayaG[i] = GDP[i]
				KayaFF[i] = pop[i]*KayaG[i]*KayaE[i]*KayaCE[i] 
			  elif (coupling in ['NFC','BGC','LFC','OFC']):
				KayaE[i] = ee[i]/(old_gdp[i]*pop[i])
				KayaCE[i] = Ff[i]/ee[i]
				KayaFF[i] = Ff[i]
			  Ca[i] = Ca[i-1] + (KayaFF[i] + Fonet[i] + Fb[i])*deltat

			#************************* Save Results ********************#
			f.close()
			xe = range(2000,2100)

			CO2m = Ca/2.1
			if bgcBase:
				np.savetxt('OutputData/edTempEffects_Temp_'+coupling+level+'.csv',ed_tempEffects['temp'],delimiter=',')
				np.savetxt('OutputData/edTempEffects_energy_'+coupling+level+'.csv',ed_tempEffects['demand'],delimiter=',')
				np.savetxt('OutputData/bgcBaseTair'+rcp+coupling+level+'.csv',Tair,delimiter=',')
				np.savetxt('OutputData/bgcBaseCO2m'+rcp+coupling+level+'.csv',CO2m,delimiter=',')
				np.savetxt('OutputData/bgcBaseKayaFF'+rcp+coupling+level+'.csv',KayaFF,delimiter=',')
				np.savetxt('OutputData/bgcBaseKayaCE'+rcp+coupling+level+'.csv',KayaCE,delimiter=',')
				np.savetxt('OutputData/bgcBaseKayaE'+rcp+coupling+level+'.csv',KayaE,delimiter=',')
				np.savetxt('OutputData/bgcBaseEE'+rcp+coupling+level+'.csv',EE,delimiter=',')
				np.savetxt('OutputData/bgcBaseED'+rcp+coupling+level+'.csv',ED,delimiter=',')
				np.savetxt('OutputData/bgcBasepe_transonly'+rcp+coupling+level+'.csv',pe_transonly,delimiter=',')
				np.savetxt('OutputData/bgcBasepe_econly'+rcp+coupling+level+'.csv',pe_econly,delimiter=',')
				np.savetxt('OutputData/bgcBaseCa'+rcp+coupling+level+'.csv',Ca,delimiter=',')
				np.savetxt('OutputData/bgcBaseFonet'+rcp+coupling+level+'.csv',Fonet,delimiter=',')
				np.savetxt('OutputData/bgcBaseFb'+rcp+coupling+level+'.csv',Fb,delimiter=',')
				np.savetxt('OutputData/bgcBaseGDP'+rcp+coupling+level+'.csv',GDP,delimiter=',')
				np.savetxt('OutputData/bgcBasepop'+rcp+'.csv',pop,delimiter=',')
				np.savetxt('OutputData/bgcBasePop'+rcp+coupling+level+'.csv',pop_adj,delimiter=',')
				np.savetxt('OutputData/bgcBaseFf'+rcp+'.csv',Ff,delimiter=',')
				np.savetxt('OutputData/bgcBaseold_gdp'+rcp+'.csv',old_gdp,delimiter=',')
				np.savetxt('OutputData/bgcBaseee'+rcp+'.csv',ee,delimiter=',')
				np.savetxt('OutputData/bgcBasedate.csv',date,delimiter=',')
				np.savetxt('OutputData/bgcBaseLUCFrac'+rcp+coupling+level+'.csv',taubLUCfrac2,delimiter=',')
			if not bgcBase:
				np.savetxt('OutputData/Tair'+rcp+coupling+level+'.csv',Tair,delimiter=',')
				np.savetxt('OutputData/CO2m'+rcp+coupling+level+'.csv',CO2m,delimiter=',')
				np.savetxt('OutputData/KayaFF'+rcp+coupling+level+'.csv',KayaFF,delimiter=',')
				np.savetxt('OutputData/KayaCE'+rcp+coupling+level+'.csv',KayaCE,delimiter=',')
				np.savetxt('OutputData/KayaE'+rcp+coupling+level+'.csv',KayaE,delimiter=',')
				np.savetxt('OutputData/EE'+rcp+coupling+level+'.csv',EE,delimiter=',')
				np.savetxt('OutputData/ED'+rcp+coupling+level+'.csv',ED,delimiter=',')
				np.savetxt('OutputData/pe_transonly'+rcp+coupling+level+'.csv',pe_transonly,delimiter=',')
				np.savetxt('OutputData/pe_econly'+rcp+coupling+level+'.csv',pe_econly,delimiter=',')
				np.savetxt('OutputData/Ca'+rcp+coupling+level+'.csv',Ca,delimiter=',')
				np.savetxt('OutputData/Fonet'+rcp+coupling+level+'.csv',Fonet,delimiter=',')
				np.savetxt('OutputData/Fb'+rcp+coupling+level+'.csv',Fb,delimiter=',')
				np.savetxt('OutputData/GDP'+rcp+coupling+level+'.csv',GDP,delimiter=',')
				np.savetxt('OutputData/pop'+rcp+'.csv',pop,delimiter=',')
				np.savetxt('OutputData/Pop'+rcp+coupling+level+'.csv',pop_adj,delimiter=',')
				np.savetxt('OutputData/Ff'+rcp+'.csv',Ff,delimiter=',')
				np.savetxt('OutputData/old_gdp'+rcp+'.csv',old_gdp,delimiter=',')
				np.savetxt('OutputData/ee'+rcp+'.csv',ee,delimiter=',')
				np.savetxt('OutputData/date.csv',date,delimiter=',')
				np.savetxt('OutputData/LUCFrac'+rcp+coupling+level+'.csv',taubLUCfrac2,delimiter=',')
			return -1*np.sum(Fonet),-1*np.sum(Fb),Ca[-1]-Ca[0],Tair[-1]-Tair[0],Tair,GDP,old_gdp,KayaFF,Ff,Fonet,Fb,Ca,date,pop,pop_adj,KayaE,KayaCE,ee,ED,EE

	return mainFun(scenario,rcp,coupling,level,bgcBase)
	
