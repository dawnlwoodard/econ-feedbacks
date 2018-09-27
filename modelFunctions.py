import numpy as np
import time
import matplotlib.pylab as py
import matplotlib.pyplot as plt
import sys 
import csv

def runModel(scenario,rcp,coupling,level,bgcBase,nyears,aland, bland, K0,K1,K2,K3,Q10,beta,NPPinitial,taub,Revellef,depthm,tautempfac,cat,kgas,initialCappm,dic=0.00191,deptho=3800.0,deltat=1,alpha=0.0062,mwc = 12.01):
	#Q10=0
	#if scenario in ['NFC','BGC','FC']:
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
#	else:
#		tautempfac=tautempfac*(3.0/3.0)
#		Q10 = 1.1
#		beta = 0.65
			
	print '\n\n****************\n',coupling,level,tautempfac,Q10,'\n***************\n\n'			
	#print level,tautempfac,Q10			
			
	deltaF	     = 5.35*np.log(1140./280)
	#nyears = 141
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
#	taumtod = 9.0
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
	#key energy suppply relationships

	pe = np.genfromtxt("/Volumes/GoogleDrive/My Drive/Research/BoxModel/model2/energy85_all.csv",delimiter=",")
	luc = np.genfromtxt("co2Data.csv",delimiter=",")
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
		old_gdp = np.genfromtxt("/Volumes/GoogleDrive/My Drive/Research/BoxModel/model2/gdp"+rcp+"_all.csv",delimiter=",")
		old_gdp = np.array(old_gdp).astype(float)
		old_gdp[0] = old_gdp[1]
		dold_gdp = (old_gdp[1:]-old_gdp[:-1])/old_gdp[:-1]
		dold_gdp = np.append(dold_gdp[0],dold_gdp)

		#temp_data = np.genfromtxt("../model2/rcp85tas.csv",delimiter=",")
		#data from KNMI (Netherlands meteorological institute) climate explorer
		#temp_data = np.genfromtxt("tas_multimodelmean.dat")
		#year = temp_data[:,0]
		#temp_data = temp_data[:,1]-273.15 #convert to degrees celsius

		#tdiff = [x-temp_data[0] for x in temp_data]

		#tempHist = np.mean(temp_data[120:150])

		#mean_dgdp= np.mean(dold_gdp[120:150]) 
		burke_temp = np.genfromtxt("InputData/burke_temp_base.csv",delimiter=",")
		temp_min = np.genfromtxt("InputData/burke_temp_max.csv",delimiter=",")
		temp_max = np.genfromtxt("InputData/burke_temp_min.csv",delimiter=",")

		dgdp_data = np.genfromtxt("InputData/burke_base.csv",delimiter=",")
		dgdp_data = dgdp_data[1:,:]
		dgdp_high = np.genfromtxt("InputData/burke_max2.csv",delimiter=",")
		dgdp_low = np.genfromtxt("InputData/burke_min2.csv",delimiter=",")

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
		ee = np.genfromtxt("/Volumes/GoogleDrive/My Drive/Research/BoxModel/model2/energy"+rcp+"_all.csv",delimiter=",")
		ff = np.genfromtxt("/Volumes/GoogleDrive/My Drive/Research/BoxModel/model2/Fossilhistandrcp"+rcp+".txt",skip_header=1)
		pop = np.genfromtxt("/Volumes/GoogleDrive/My Drive/Research/BoxModel/model2/pop"+rcp+"_all.csv",delimiter=",")
		return ff[:,0],old_gdp[:-1],burke_temp,temp_min,temp_max,ff[:,1],pop,ee,eFracs,eLosses,finE

	def fitGDP(tcurve,xvals):
		#fit burke tempeature data with a*x^2 + b*x + c model for base, min, and max
		fit = np.polyfit(tcurve[:,0],tcurve[:,1],2)
		#x = np.arange(0,5,0.1)
		temp_fit = [fit[2] + i*fit[1] + fit[0]*i**2 for i in xvals]
		return fit,temp_fit

	def getGDP(fit,Tc,T0,old_gdp_i): #inputs are curve fit for temp-response, current temp, initial temp, old_gdp from GCAM at current timestep

		tdiff = Tc-T0
		dnew_gdp = fit[2]+tdiff*fit[1]+fit[0]*tdiff**2
		dnew_gdp = dnew_gdp*0.01
		new_gdp = old_gdp_i+dnew_gdp*old_gdp_i
		return new_gdp

	def getPop(fit,old_pop_i,old_gdp_i): #inputs are curve fit for temp-response, current temp, initial temp, old_gdp from GCAM at current timestep

		#dnew_pop = fit[0]*np.exp(old_gdp_i*fit[1])+fit[2]
		dnew_pop = fit[0]*old_gdp_i+fit[1]
		dnew_pop = dnew_pop*0.01
		new_pop = old_pop_i+dnew_pop*old_pop_i
		return new_pop

	def runOceanDiffu(Ca_o,Com_o,Cod_o,Tair_o,Tair0):
                #print '***in ocean module, nboxes=',nboxes
                Fod_n = np.zeros(nboxes)
                Cod_n = np.zeros(nboxes)
                Kdeep = Kdeep_o*(1.0-(Tair_o-Tair0)*tautempfac/100.0)
                minKd = Kdeep*kfrac
                if (changeDepth):
                        KdeepList = [Kdeep - i*(Kdeep-minKd)/float(nboxes) for i in range(nboxes)]
                else:
                        KdeepList = [Kdeep]*nboxes
                #Kdeep = Kdeep_o
                Fao_n   = kgas*(Ca_o/2.12)
                Foa_n   = kgas*(Ca[0]/2.12)*(1.0 + Revellef*(Com_o-Com[0])/Com[0])
        #       print 'kgas,Ca_o,Ca[0],Fao,Foa,Com_o,Com[0]',kgas,Ca_o,Ca[0],Fao_n,Foa_n,Com_o,Com[0]
                #Fmd_n = Kdeep*Com_o/((depthm+hboxes)/2.0)
                #Fdm_n = Kdeep*Cod_o/((depthm+hboxes)/2.0)
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
                        #Cod_n[i] = (Cod_o[i] + deltat*0.5*(Fod_n[i] 
                #Com_n   = Com_o + (Fao_n*(1e15)/(mwc*oceanarea) - Foa_n*(1e15)/(mwc*oceanarea) - k_mj*Com_o+k_jm*Cod_o[0])*deltat

                Com_n   = Com_o + (Fao_n - Foa_n - k_mj*Com_o+k_jm*Cod_o[0])*deltat
                #Com_n = (Com_o + 
                # the three steps below convert from mol/m2/yr to Pg C/yr for the ocean fluxes
                Gao_n   = Fao_n*mwc*oceanarea/1.0e15
                Goa_n   = Foa_n*mwc*oceanarea/1.0e15
                Fonet_n = Goa_n-Gao_n
                Gmd_n   = (k_mj*Com_o)*mwc*oceanarea/1.0e15
                Gdm_n   = (k_jm*Cod_o[0])*mwc*oceanarea/1.0e15
                #Fonet_n = Foa_n-Fao_n
                #print 'ocean-atm fluxes',Fao_n,Foa_n,Fmd_n,Fdm_n,Com_o,Com[0],Revellef*(Com_o-Com[0])/Com[0]
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
		# the three steps below convert from mol/m2/yr to Pg C/yr for the ocean fluxes
		Gao_n   = Fao_n*mwc*oceanarea/1.0e15
		Goa_n   = Foa_n*mwc*oceanarea/1.0e15
		Fonet_n = Goa_n-Gao_n
		return Fonet_n,Com_n,Cod_n

	#Land Subroutine
	def runLand(Ca_o,Tair_o,Cb_o,luc_o,taub_i):
		#print '**********************************************************'
                #print 'i,Ca_o,Tair_o,Cb_o,luc_o,taub_i',i,Ca_o,Tair_o,Cb_o,luc_o,taub_i
                #print '**********************************************************'
		if (Ca_o > 0) and (Ca[0] > 0):
                        NPP_n  = NPP[0]*(1 + beta*np.log(Ca_o/Ca[0]))*(1+aland*(Tair_o-Tair[0])**2+bland*(Tair_o-Tair[0]))
                else:
                        NPP_n  = NPP[0]*(1 + beta*0)*(1+aland*(Tair_o-Tair[0])**2+bland*(Tair_o-Tair[0]))
		#NPP_n  = NPP[0]*(1 + beta*np.log(Ca_o/Ca[0]))*(1+aland*(Tair_o-Tair[0])**2+bland*(Tair_o-Tair[0]))
		#NPP_n  = NPP[0]*(1 + beta*np.log(Ca_o/Ca[0]))*(1+mland*(Tair_o-Tair[0]))
		Rh_n	 = (1.0/taub_i)*Q10**((Tair_o -Tair[0])/10.0)*Cb_o
		Cb_n	 = Cb_o + (NPP_n - Rh_n - luc_o)*deltat
		Fb_n	 = Rh_n - NPP_n + luc_o 
		return Fb_n,NPP_n,Rh_n,Cb_n

		
	def mainFun(scenario,rcp,coupling,level,bgcBase):
		# atmosphere initial conditions
		Ca[0]	     = initialCaPgC		    # initial boundary condition
		Ca[0]	     = co2Old[0]*2.12
		#print '\n \n *********** Initial C (PgC)',Ca[0],co2Old[0],'****************\n\n'
		# ocean initial conditions
		Fao[0]	     = kgas*(Ca[0]/2.1)		    # initial boundary condition
		Foa[0]	     = Fao[0]			    # initial boundary condition
		Com[0]	     = depthm*dic*1000.0	    # initial boundary condition
		#Cod[0]	     = (deptho-depthm)*dic*1000.0   # initial boundary condition
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
			#Ca[0] = Ca[0]*1.01
			date = range(0,nyears)
			date = [i+(2100-nyears) for i in date]
			tic = time.clock()
			for i in range(1,nyears):
				#the steps below describe the ocean model flux and pool components
				#Fonet[i],Com[i],Cod[i] = runOcean(Ca[i-1],Com[i-1],Cod[:,i-1],Tair[i-1],Tair[0])
				if (diffusionJR):
					Fonet[i],Com[i],Cod[:,i] = runOceanDiffu(Ca[i-1],Com[i-1],Cod[:,i-1],Tair[i-1],Tair[0])
				else:
					Fonet[i],Com[i],Cod[:,i] = runOcean(Ca[i-1],Com[i-1],Cod[:,i-1],Tair[i-1],Tair[0])
					
				# the following lines update the fluxes and pools in the terrestrial biosphere
				Fb[i],NPP[i],Rh[i],Cb[i] = runLand(Ca[i-1],Tair[i-1],Cb[i-1],0,taub[0])
				# estimate the global mean surface air temperature 

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
			#print ("The total anthro. ocean carbon inventory is: %.1f " % (-1*np.sum(Fonet)))
			#print ("The total anthro. land carbon inventory is: %.1f " % (-1*np.sum(Fb)))
			#print ("The total atmos carbon inventory is: %.1f " % (CO2m[-1]))
			#print ("The total anthro. atmos carbon inventory is: %.1f " % (Ca[-1]-Ca[0]))
			#print ("The total emitted fossil fuel carbon is: %.1f " % np.sum(Ff))
			#print ("The temperature change from pre-industrial is: %.5f " % (Tair[-1]-Tair[0]))
			
			#np.savetxt('ModelData/Tair'+scenario+coupling+level+'.csv',Tair,delimiter=',')
			#np.savetxt('ModelData/CO2m'+scenario+coupling+level+'.csv',CO2m,delimiter=',')
			#np.savetxt('ModelData/Ff'+scenario+'.csv',Ff,delimiter=',')
			#np.savetxt('ModelData/Ca'+scenario+coupling+level+'.csv',Ca,delimiter=',')
			#np.savetxt('ModelData/Fonet'+scenario+coupling+level+'.csv',Fonet,delimiter=',')
			#np.savetxt('ModelData/Fb'+scenario+coupling+level+'.csv',Fb,delimiter=',')
			#np.savetxt('ModelData/date.csv',date,delimiter=',')
			
			#print ("The total emitted fossil fuel carbon is: %.1f " % np.sum(Ff))
			return -1*np.sum(Fonet),-1*np.sum(Fb),Ca[-1]-Ca[0],Tair[-1]-Tair[0],Tair,Ff,Fonet,Fb,date
		else:
#read in socioeconomic data
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
				#print '\n\n\n *********** GDP parameters ****************\n\n',fit,temp_fit,'*******************************\n\n\n'
			elif (level=='max'):
				fit, temp_fit = fitGDP(temp_max,np.arange(1,5,0.1))
			else: #level == 'min'
				fit, temp_fit = fitGDP(temp_min,np.arange(1,5,0.1))

			#energy demand relationship
			ed_base = 0.04
			ed_high = 0.08
			ed_low = 0.005

			#print '\n\n****************\n',level,tautempfac,Q10,'\n***************\n\n'			

			#************************** Run Model ***********************#
			f = open('taubLUCData'+rcp+coupling+level+'.csv','w')
			f.write(rcp+coupling+level+'\n')
	 		for i in range(1,nyears):
			  #print '*******************************\n',scenario,i
			  taubLUCfrac = (Cb[0]-tb*np.sum(luc[:i]))/Cb[0]
			  taubLUCfrac2[i] = (Cb[0]-tb*np.sum(luc[:i]))/Cb[0]
			
			  f.write(str(i)+','+str(taubLUCfrac)+'\n')
			  taub[i] = taub[0]*(Cb[0]-tb*np.sum(luc[:i]))/Cb[0]
			  
			  if (coupling in ["FC","NFC"]):
				  Fonet[i],Com[i],Cod[:,i] = runOceanDiffu(Ca[i-1],Com[i-1],Cod[:,i-1],Tair[i-1],Tair[0])
				  # the following lines update the fluxes and pools in the terrestrial biosphere
				  Fb[i],NPP[i],Rh[i],Cb[i] = runLand(Ca[i-1],Tair[i-1],Cb[i-1],luc[i],taub[i])
			  elif (coupling in ["LFC"]):
				  Fonet[i],Com[i],Cod[:,i] = runOceanDiffu(Ca[i-1],Com[i-1],Cod[:,i-1],Tair[0],Tair[0])
				  # the following lines update the fluxes and pools in the terrestrial biosphere
				  Fb[i],NPP[i],Rh[i],Cb[i] = runLand(Ca[i-1],Tair[i-1],Cb[i-1],luc[i],taub[i])

			  elif (coupling in ["OFC"]):
				  Fonet[i],Com[i],Cod[:,i] = runOceanDiffu(Ca[i-1],Com[i-1],Cod[:,i-1],Tair[i-1],Tair[0])
				  Fb[i],NPP[i],Rh[i],Cb[i] = runLand(Ca[i-1],Tair[0],Cb[i-1],luc[i],taub[i])

			  elif (coupling in ["BGC"]):
				  Fonet[i],Com[i],Cod[:,i] = runOceanDiffu(Ca[i-1],Com[i-1],Cod[:,i-1],Tair[0],Tair[0])
				  Fb[i],NPP[i],Rh[i],Cb[i] = runLand(Ca[i-1],Tair[0],Cb[i-1],luc[i],taub[i])
			  
			  else: #coupling is "FC","BGC","NFC"
				  #the steps below describe the ocean model flux and pool components
				  # the following lines update the fluxes and pools in the terrestrial biosphere
				  if bgcBase:
					Fb[i],NPP[i],Rh[i],Cb[i] = runLand(Ca[i-1],Tair[0],Cb[i-1],luc[i],taub[i])
				  	Fonet[i],Com[i],Cod[:,i] = runOceanDiffu(Ca[i-1],Com[i-1],Cod[:,i-1],Tair[0],Tair[0])
				  else:
					Fb[i],NPP[i],Rh[i],Cb[i] = runLand(Ca[i-1],Tair[i-1],Cb[i-1],luc[i],taub[i])
				  	Fonet[i],Com[i],Cod[:,i] = runOceanDiffu(Ca[i-1],Com[i-1],Cod[:,i-1],Tair[i-1],Tair[0])
			  #KayaC[i] = Ff[i]/old_gdp[i]
			  RF[i] = 5.35*np.log(Ca[i-1]/Ca[0])
			  for j in range(1,i+1):
				dTair[i] = dTair[i]+RF[i-j]*tempinc[j-1]
			  Tair[i] = Tair[0] + dTair[i]
			  
			  if (i<195):
				newDeaths = 0
			  else:	
				newDeaths = 199283*(Tair[i]-Tair[195])
				#newDeaths = 250000
			  if (level=='max'):
			  	if ((Tair[i]-Tair[0]) > 0.8):
					GDP[i] = getGDP(fit,Tair[i],Tair[0],old_gdp[i])
				#       if (coupling in ['AFC','FC']): #add in Energy-GDP component
				#               GDP[i] = GDP[i]+GDP[i]*0.6*((EE[i-1]-EE[i-2])/EE[i-1])
				else:
					GDP[i] = old_gdp[i]
				pDeaths = pDeaths + 1.5*newDeaths
				pop_adj[i] = pop[i]-pDeaths
			  	#pop_adj[i] = getPop((2.2686,-0.0001176,0.9593),pop[i],GDP[i])
			  	#pop_adj[i] = getPop((-5.4237e-5,2.019+1),pop[i],GDP[i])
			  elif (level=='base'):
			  	if ((Tair[i]-Tair[0]) > 0.8):
					GDP[i] = getGDP(fit,Tair[i],Tair[0],old_gdp[i])
				else:
					GDP[i] = old_gdp[i]
				#GDP[i] = old_gdp[i]-old_gdp[i]*0.00236*(Tair[i]-Tair[0])**2
				pDeaths = pDeaths + 1*newDeaths
				pop_adj[i] = pop[i]-pDeaths
			  	#pop_adj[i] = getPop((2.01274341,-3.3493644e-4,6.04879572e-1),pop[i],GDP[i])
			  	#pop_adj[i] = getPop((1.8559,-0.0003388,0.6027),pop[i],GDP[i]-GDP[0])
			  	#pop_adj[i] = getPop((-5.4237e-5,2.019),pop[i],GDP[i])
			  elif (level=='min'):
				GDP[i] = old_gdp[i]-old_gdp[i]*0.00236*(Tair[i]-Tair[0])**2
				#GDP[i] = old_gdp[i]-old_gdp[i]*0.00016*(Tair[i]-Tair[0])**2
				pDeaths = pDeaths + 0.5*newDeaths
				pop_adj[i] = pop[i]-pDeaths
			  	#pop_adj[i] = getPop((1.2367,-0.0006706,0.06776),pop[i],GDP[i]-GDP[0])
			  	#pop_adj[i] = getPop((-5.4237e-5,2.019-1),pop[i],GDP[i])
			  #pop_adj[i] = pop_adj[i-1]+pGrowth[i-1]-199283*(Tair[i]-Tair[195])
			  #if (coupling=='FC'):
			  #	print 'pop_adj,old pop,gdp',pop_adj[i],pop[i],GDP[i]
			  #pop_adj[i] = pop[i] - 256667*(Tair[i]-Tair[0])
			  tdiff = Tair[i] - Tair[0]
			  eExtra = 3.36*1e3*(-0.49*(tdiff-0.613)**3+3.57*tdiff**2-4.344*tdiff+0.014)*1
			  if (level == 'base'):
				#ED[i] = pec[i]+pec[i]*ed_base*(Tair[i]-Tair[0])+nonec[i]
				ED[i] = ee[i]+3.36*(-0.49*(tdiff-0.613)**3+3.57*tdiff**2-4.344*tdiff+0.014)*1
				#transpe_adj[i] = transpe[i]
				transpe_adj[i] = transpe[i]*.60*(1+0.003*(Tair[i]-Tair[0]))+transpe[i]*0.22*(1+0.01*(Tair[i]-Tair[0]))+transpe[i]*0.08+transpe[i]*0.04*(1+0.00373*(Tair[i]-Tair[0]))+transpe[i]*0.02*(1-0.0043*(Tair[i]-Tair[0]))+transpe[i]*0.04
			  	#pe_adj[i]=opec_adj[i] + cpec_adj[i] + ngpec_adj[i] + transpe_adj[i] + pe_nonff[i]
				opec_adj[i] = opec[i]*(1-0.0075*(Tair[i]-Tair[0]))*(1-0.0087*(Tair[i]-Tair[0]))
				cpec_adj[i] = cpec[i]*(1-0.0105*(Tair[i]-Tair[0]))*(1-0.0087*(Tair[i]-Tair[0]))
				ngpec_adj[i] = ngpec[i]*(1-0.0075*(Tair[i]-Tair[0]))*(1-0.0087*(Tair[i]-Tair[0]))

			  elif (level == 'max'):
				#ED[i] = pec[i]+pec[i]*ed_high*(Tair[i]-Tair[0])+nonec[i]
				ED[i] = ee[i]+3.36*(-0.49*(tdiff-0.613)**3+3.57*tdiff**2-4.344*tdiff+0.014)*.5
				#transpe_adj[i] = transpe[i]
				transpe_adj[i] = transpe[i]*(1+0.008*(Tair[i]-Tair[0]))
			  	#pe_adj[i]=opec_adj[i] + cpec_adj[i] + ngpec_adj[i] + transpe_adj[i] + pe_nonff[i]
				opec_adj[i] = opec[i]*(1-0.003*(Tair[i]-Tair[0]))*(1-0.0033*(Tair[i]-Tair[0]))
				cpec_adj[i] = cpec[i]*(1-0.006*(Tair[i]-Tair[0]))*(1-0.0033*(Tair[i]-Tair[0]))
				ngpec_adj[i] = ngpec[i]*(1-0.003*(Tair[i]-Tair[0]))*(1-0.0033*(Tair[i]-Tair[0]))
			  	#pe_adj[i]=(opec_adj[i] + cpec_adj[i] + ngpec_adj[i] + transpe_adj[i] + pe_nonff[i])*1.5
			  	#opec_adj[i] = opec[i]*(1-0.01*(Tair[i]-Tair[0]))*(1-0.014*(Tair[i]-Tair[0]))
			  	#cpec_adj[i] = cpec[i]*(1-0.015*(Tair[i]-Tair[0]))*(1-0.014*(Tair[i]-Tair[0]))
			  	#ngpec_adj[i] = ngpec[i]*(1-0.01*(Tair[i]-Tair[0]))*(1-0.014*(Tair[i]-Tair[0]))
			  else:
				ED[i] = ee[i]+3.36*(-0.49*(tdiff-0.613)**3+3.57*tdiff**2-4.344*tdiff+0.014)*1.5
				#ED[i] = pec[i]+pec[i]*ed_low*(Tair[i]-Tair[0])+nonec[i]
				transpe_adj[i] = transpe[i]
				opec_adj[i] = opec[i]*(1-0.01*(Tair[i]-Tair[0]))*(1-0.014*(Tair[i]-Tair[0]))
				cpec_adj[i] = cpec[i]*(1-0.015*(Tair[i]-Tair[0]))*(1-0.014*(Tair[i]-Tair[0]))
				ngpec_adj[i] = ngpec[i]*(1-0.01*(Tair[i]-Tair[0]))*(1-0.014*(Tair[i]-Tair[0]))
			  	#pe_adj[i]=(opec_adj[i] + cpec_adj[i] + ngpec_adj[i] + transpe_adj[i] + pe_nonff[i])*0.5
			  	#opec_adj[i] = opec[i]*(1-0.005*(Tair[i]-Tair[0]))*(1-0.0033*(Tair[i]-Tair[0]))
			  	#cpec_adj[i] = cpec[i]*(1-0.015*(Tair[i]-Tair[0]))*(1-0.0033*(Tair[i]-Tair[0]))
			  	#ngpec_adj[i] = ngpec[i]*(1-0.005*(Tair[i]-Tair[0]))*(1-0.0033*(Tair[i]-Tair[0]))
			  pe_adj[i]=opec_adj[i] + cpec_adj[i] + ngpec_adj[i] + transpe_adj[i] + pe_nonff[i]
			  pe_econly[i]=opec_adj[i] + cpec_adj[i] + ngpec_adj[i] + transpe[i] + pe_nonff[i]
			  pe_transonly[i]=opec[i] + cpec[i] + ngpec[i] + transpe_adj[i] + pe_nonff[i]

			  EE[i] = pe_adj[i]
			  #else:
			#       pe_econly[i] = ee[i]
			#       pe_transonly[i] = ee[i]
			#       EE[i] = ee[i]
			#       ED[i] = ee[i]
			  #tofix
			  if (coupling=="GFC"):
				KayaE[i] = (ee[i]/(old_gdp[i]*pop[i]))
				KayaCE[i] = Ff[i]/ee[i]
				KayaFF[i] = pop[i]*GDP[i]*KayaE[i]*KayaCE[i] #GDP on, E off
			  elif (coupling=="EFC"):
				KayaE[i] = (ee[i]/(old_gdp[i]*pop[i]))
				KayaCE[i] = Ff[i]/EE[i]
				KayaFF[i] = pop[i]*old_gdp[i]*KayaE[i]*KayaCE[i] #GDP off, E on
			  elif (coupling=='DFC'):
				KayaE[i] = (ED[i]/(old_gdp[i]*pop[i])
				KayaCE[i] = Ff[i]/ee[i]
				KayaFF[i] = pop[i]*old_gdp[i]*KayaE[i]*KayaCE[i] #GDP off, E off, ED on

			  elif (coupling=='PFC'):
				KayaE[i] = (ee[i]/(old_gdp[i]*pop[i]))
				KayaCE[i] = Ff[i]/ee[i]
				#KayaG[i] = GDP[i]/pop[i]
				KayaG[i] = old_gdp[i]/pop[i]
				KayaFF[i] = pop_adj[i]*KayaG[i]*KayaE[i]*KayaCE[i]
			  elif (coupling in ['AFC','FC']): #coupling == AFC
				KayaE[i] = (ee[i]/(old_gdp[i]*pop[i]))
				KayaCE[i] = Ff[i]/EE[i]
				KayaG[i] = GDP[i]/pop[i]
				KayaFF[i] = pop[i]*KayaG[i]*KayaE[i]*KayaCE[i] #GDP on, E on
			  elif (coupling in ['NFC','BGC','LFC','OFC']):
				KayaE[i] = ee[i]/(old_gdp[i]*pop[i])
				KayaCE[i] = Ff[i]/ee[i]
				KayaFF[i] = Ff[i]
			  #if (coupling == 'BGC'):
			  #	Tair[i] = Tair[0]
			  #print 'Ca[i-1],KayaFF[i],Ff[i],Fonet[i],Fb[i],deltat = ',Ca[i-1],KayaFF[i],Ff[i],Fonet[i],Fb[i],deltat
			  #print '**********************************************************'
                          #print 'run,i,Ca[i-1],Ff[i],Fonet[i],Fb[i]',scenario,i,Ca[i-1],KayaFF[i],Fonet[i],Fb[i]
                          #print '**********************************************************'
			  Ca[i] = Ca[i-1] + (KayaFF[i] + Fonet[i] + Fb[i])*deltat

			#************************* Display Results ********************#
			f.close()
			xe = range(2000,2100)
			"""
			if (coupling=='PFC'):
				print 'Population diff',pop[-10:-1]-pop_adj[-10:-1]
			print ("The total anthro. ocean carbon inventory is: %.1f " % (-1*np.sum(Fonet)))
			print ("The total anthro. land carbon inventory is: %.1f " % (-1*np.sum(Fb)))
			print ("The total anthro. atmos carbon inventory is: %.1f " % (Ca[-1]-Ca[0]))
			#print ("The total emitted fossil fuel carbon is: %.1f " % np.sum(Ff))
			if (coupling in ["AFC","FC","DFC","EFC","GFC","PFC"]):
				print ("The total emitted Kaya fossil fuel carbon is: %.1f " % np.sum(KayaFF))
				print "GDP ", old_gdp[-70],GDP[-70],old_gdp[-1],GDP[-1]
				print "Energy supply ",ee[-70],EE[-70],ee[-1],EE[-1]
				print "Demand", ED[-70],ED[-1]
				print "Population", pop[-70],pop_adj[-70],pop[-1],pop_adj[-1]
			print ("The temperature change from pre-industrial is: %.1f " % (Tair[-1]-Tair[0]))
			"""	
			#print '\n\n****************\n',level,tautempfac,Q10,'\n***************\n\n'			

			CO2m = Ca/2.1
			#print 'the losses in coal, ng, oil, and dist are: ',eLosses['coal_adj'],eLosses['ng_adj'],eLosses['oil_adj'],eLosses['dist_adj']
			if bgcBase:
				print 'saving files'
				np.savetxt('ModelData/bgcBaseTair'+rcp+coupling+level+'.csv',Tair,delimiter=',')
				np.savetxt('ModelData/bgcBaseCO2m'+rcp+coupling+level+'.csv',CO2m,delimiter=',')
				np.savetxt('ModelData/bgcBaseKayaFF'+rcp+coupling+level+'.csv',KayaFF,delimiter=',')
				np.savetxt('ModelData/bgcBaseKayaCE'+rcp+coupling+level+'.csv',KayaCE,delimiter=',')
				np.savetxt('ModelData/bgcBaseKayaE'+rcp+coupling+level+'.csv',KayaE,delimiter=',')
				np.savetxt('ModelData/bgcBaseEE'+rcp+coupling+level+'.csv',EE,delimiter=',')
				np.savetxt('ModelData/bgcBaseED'+rcp+coupling+level+'.csv',ED,delimiter=',')
				np.savetxt('ModelData/bgcBasepe_transonly'+rcp+coupling+level+'.csv',pe_transonly,delimiter=',')
				np.savetxt('ModelData/bgcBasepe_econly'+rcp+coupling+level+'.csv',pe_econly,delimiter=',')
				np.savetxt('ModelData/bgcBaseCa'+rcp+coupling+level+'.csv',Ca,delimiter=',')
				np.savetxt('ModelData/bgcBaseFonet'+rcp+coupling+level+'.csv',Fonet,delimiter=',')
				np.savetxt('ModelData/bgcBaseFb'+rcp+coupling+level+'.csv',Fb,delimiter=',')
				np.savetxt('ModelData/bgcBaseGDP'+rcp+coupling+level+'.csv',GDP,delimiter=',')
				np.savetxt('ModelData/bgcBasepop'+rcp+'.csv',pop,delimiter=',')
				np.savetxt('ModelData/bgcBasePop'+rcp+coupling+level+'.csv',pop_adj,delimiter=',')
				np.savetxt('ModelData/bgcBaseFf'+rcp+'.csv',Ff,delimiter=',')
				np.savetxt('ModelData/bgcBaseold_gdp'+rcp+'.csv',old_gdp,delimiter=',')
				np.savetxt('ModelData/bgcBaseee'+rcp+'.csv',ee,delimiter=',')
				np.savetxt('ModelData/bgcBasedate.csv',date,delimiter=',')
				np.savetxt('ModelData/bgcBaseLUCFrac'+rcp+coupling+level+'.csv',taubLUCfrac2,delimiter=',')
			if not bgcBase:
				print 'saving files'
				np.savetxt('ModelData/Tair'+rcp+coupling+level+'.csv',Tair,delimiter=',')
				np.savetxt('ModelData/CO2m'+rcp+coupling+level+'.csv',CO2m,delimiter=',')
				np.savetxt('ModelData/KayaFF'+rcp+coupling+level+'.csv',KayaFF,delimiter=',')
				np.savetxt('ModelData/KayaCE'+rcp+coupling+level+'.csv',KayaCE,delimiter=',')
				np.savetxt('ModelData/KayaE'+rcp+coupling+level+'.csv',KayaE,delimiter=',')
				np.savetxt('ModelData/EE'+rcp+coupling+level+'.csv',EE,delimiter=',')
				np.savetxt('ModelData/ED'+rcp+coupling+level+'.csv',ED,delimiter=',')
				np.savetxt('ModelData/pe_transonly'+rcp+coupling+level+'.csv',pe_transonly,delimiter=',')
				np.savetxt('ModelData/pe_econly'+rcp+coupling+level+'.csv',pe_econly,delimiter=',')
				np.savetxt('ModelData/Ca'+rcp+coupling+level+'.csv',Ca,delimiter=',')
				np.savetxt('ModelData/Fonet'+rcp+coupling+level+'.csv',Fonet,delimiter=',')
				np.savetxt('ModelData/Fb'+rcp+coupling+level+'.csv',Fb,delimiter=',')
				np.savetxt('ModelData/GDP'+rcp+coupling+level+'.csv',GDP,delimiter=',')
				np.savetxt('ModelData/pop'+rcp+'.csv',pop,delimiter=',')
				np.savetxt('ModelData/Pop'+rcp+coupling+level+'.csv',pop_adj,delimiter=',')
				np.savetxt('ModelData/Ff'+rcp+'.csv',Ff,delimiter=',')
				np.savetxt('ModelData/old_gdp'+rcp+'.csv',old_gdp,delimiter=',')
				np.savetxt('ModelData/ee'+rcp+'.csv',ee,delimiter=',')
				np.savetxt('ModelData/date.csv',date,delimiter=',')
				np.savetxt('ModelData/LUCFrac'+rcp+coupling+level+'.csv',taubLUCfrac2,delimiter=',')
			#print ("The total emitted fossil fuel carbon is: %.1f " % np.sum(Ff))
			return -1*np.sum(Fonet),-1*np.sum(Fb),Ca[-1]-Ca[0],Tair[-1]-Tair[0],Tair,GDP,old_gdp,KayaFF,Ff,Fonet,Fb,Ca,date,pop,pop_adj,KayaE,KayaCE,ee,ED,EE

	return mainFun(scenario,rcp,coupling,level,bgcBase)
	
