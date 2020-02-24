# Same emissions component as DICE2016, except a hard cap is enforced fo ruse in hotelling rule exploration
@defcomp capped_emissions begin
	GSIG	= Variable(index=[time])	#Change in sigma (cumulative improvement of energy efficiency)
	SIGMA	= Variable(index=[time])	#CO2-equivalent-emissions output ratio
    EIND    = Variable(index=[time])    #Industrial emissions (GtCO2 per year)
	ETREE	= Variable(index=[time])	#Emissions from deforestation
    E       = Variable(index=[time])    #Total CO2 emissions (GtCO2 per year)
	CUMETREE= Variable(index=[time])	#Cumulative from land
    CCA     = Variable(index=[time])    #Cumulative industrial emissions
	CCATOT	= Variable(index=[time])	#Cumulative total carbon emissions

    sig0	= Parameter()				#Carbon intensity 2010-2015 (kgCO2 per output 2010 USD)
	gsigma1 = Parameter()				#Initial growth of sigma (per year)
	dsig	= Parameter()				#Decline rate of decarbonization (per period)
	eland0	= Parameter()				#Carbon emissions from land 2015 (GtCO2 per year)
	deland	= Parameter()				#Decline rate of land emissions (per period)
	e0		= Parameter()				#Industrial emissions 2015 (GtCO2 per year)
    # MIU     = Parameter(index=[time])   #Emission control rate GHGs
    YGROSS  = Parameter(index=[time])   #Gross world product GROSS of abatement and damages (trillions 2010 USD per year)
    cca0    = Parameter()               #Initial cumulative industrial emissions
    cumetree0=Parameter()				#Initial emissions from deforestation (see GAMS code)
    
    # Added for this capped version:
    fosslim = Parameter()   # Fossil fuel limit on cumulative industrial emissions (GtC)
    initial_MIU = Parameter(index=[time])   # Emission control rate GHGs
    final_MIU   = Variable(index=[time])
    cap_hit::Bool = Variable()

    function run_timestep(p, v, d, t)

        if is_first(t)
            v.cap_hit = false
        end 

		#Define function for GSIG
		if is_first(t)
			v.GSIG[t] = p.gsigma1
		else
			v.GSIG[t] = v.GSIG[t-1] * ((1 + p.dsig)^5)
		end
		
		#Define function for SIGMA
		if is_first(t)
			v.SIGMA[t] = p.sig0
		else
			v.SIGMA[t] = v.SIGMA[t-1] * exp(v.GSIG[t-1] * 5)
		end
		
        #Define function for EIND
        # v.EIND[t] = v.SIGMA[t] * p.YGROSS[t] * (1 - p.MIU[t])
        # added for this capped version:
        if v.cap_hit
            v.final_MIU[t] = 1
        else
            v.final_MIU[t] = p.initial_MIU[t]
        end
        v.EIND[t] = v.SIGMA[t] * p.YGROSS[t] * (1 - v.final_MIU[t]) # change for this capped version        
            
		#Define function for ETREE
		if is_first(t)
			v.ETREE[t] = p.eland0
		else
			v.ETREE[t] = v.ETREE[t - 1] * (1 - p.deland)
		end
		
		#Define function for CUMETREE
		if is_first(t)
			v.CUMETREE[t] = p.cumetree0
		else
			v.CUMETREE[t] = v.CUMETREE[t-1] + v.ETREE[t-1]
		end
		
        #Define function for CCA
        if is_first(t)
            v.CCA[t] = p.cca0
        else
            v.CCA[t] = v.CCA[t-1] + v.EIND[t-1] * 5/3.666
        end

        # added for this capped version:
        next_CCA = v.CCA[t] + (v.EIND[t] * 5/3.66)
        if (!v.cap_hit) && (next_CCA > p.fosslim)

            v.EIND[t] = v.EIND[t] - (next_CCA - p.fosslim) * 3.66/5
            v.final_MIU[t] =  1 - v.EIND[t] / (v.SIGMA[t] * p.YGROSS[t])

            v.cap_hit = true
        end
            
        #Define function for E
        v.E[t] = v.EIND[t] + v.ETREE[t]

		#Define function for CCATOT
		v.CCATOT[t] = v.CCA[t] + v.CUMETREE[t]			
    end
end