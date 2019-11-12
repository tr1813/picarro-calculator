import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from math import exp,log,sqrt


def readRaw(filename):
	
	df = pd.read_csv(filename)

	renamed = dict([(i,i.strip()) for i in df]) # strip white space from column names
	df.rename(columns=renamed, inplace=True)	# and rename all headers of the dataFrame with new 'clean names'
	df["Identifier 1"] = [i.strip() for i in df["Identifier 1"].values]
	df["Identifier 2"] = [i.strip() for i in df["Identifier 2"].values]
	df["Line"] = [int(i) for i in df["Line"]]
	df.set_index(["Identifier 1","Identifier 2","Inj Nr"], inplace=True)
	return df

def checkEmpty(df):

	if np.any(np.isnan(df["H2O_Mean"].values)) == 'True':  # checks a boolean array which evaluates to True if
															# if one cell is empty
		
		print("Warning: there is at least one empty line!")
	else: 
		print("No empty cells. Proceeding...")

def checkVolume(df):

	conds = [df.H2O_Mean.values < 17000 , df.H2O_Mean.values > 23000]	# create the conditions for ignoring a certain analysis line
	choices = [1, 2]													# what the error code value will become

	out = np.select(conds, choices, default=0)							# mapping the error code on the conditions

	df["Error Code"]=out												# updating the dataFrame

	above = df.where(df["Error Code"] == 2)
	below = df.where(df["Error Code"] == 1)
	warning = pd.concat([above,below]).dropna()
	lines =[i[2] for i in warning.index.values]
	for i in lines:
		print("Warning: H20 value outside bounds on line {}...".format(i))

def runSummary(df):

	run_summary = df[["Line",
				   "H2O_Mean",
				   "Good",
				   "Error Code",
				   "d(18_16)Mean",
				   "d(D_H)Mean"]]
	
	run_summary = run_summary.sort_index(axis=0)
	return run_summary

def plotSummary(df):

	fig,(ax1,ax2,ax3) = plt.subplots(3,1,figsize= (6,8),sharex= True)
	
	ax1.plot(df["Line"],df["H2O_Mean"],'o',color = 'black',markerfacecolor='w')
	ax1.legend(loc='best')
	ax1.set_ylim(15000,30000)
	ax1.set_ylabel("H2O in ppmV")
	ax1.axhline(20000,color = 'red')
	ax1.fill_between(np.arange(-2,133,1),17000,23000,color='red',alpha=0.5)
	ax1.set_xlim(0,132)

	ax2.plot(df["Line"],df["d(18_16)Mean"],'o',color = 'black',markerfacecolor='w')
	ax2.legend(loc='best')
	ax2.set_ylim(-50,10)
	ax2.set_ylabel("d18O in permil")

	ax3.plot(df["Line"],df["d(D_H)Mean"],'o',color = 'black',markerfacecolor='w')
	ax3.legend(loc='best')
	ax3.set_ylim(-300,50)
	ax3.set_ylabel("d2H in permil")
	ax3.set_xlabel("Analysis line")


	for ax in (ax1,ax2,ax3):
		ax.grid()
		ax.tick_params(direction='in',top=True,right=True)

	plt.tight_layout()
	plt.show()

def IsotopeSelect(df,isotope):
	if isotope == "O":
		dat =  df[["Line","d(18_16)Mean","Ignore"]]
	else:
		dat =  df[["Line","d(D_H)Mean","Ignore"]]

	print(type(dat))
	return dat

def get_initial_coeffs(i):
	coeffs=[]
	c=0
	k = get_k(i)
	while c <10:
			
		coeffs.append(round(i*k**c,7))			
		c+=1
	return coeffs
	
def get_k(i):
	k = exp(-log(i)/9)
	return k

def initMemCoeffs():

	injection = np.arange(0,10,1)
	coeffs =get_initial_coeffs(0.849029498880193) 
	mem_coeffs =dict([(i+1,j) for i,j in zip(injection,coeffs)])


	init_mem_corr={
		"combined_sd":"None",
		"date":'None',
		"mem_coeffs":mem_coeffs
	}
	return init_mem_corr

def ignoreFirstSeven(df):
	df = df.where(df["Line"]>6).dropna()
	return df

def getSampleNames(df):
	return df.index.levels[0].values  

def getMeanSDs(df,isotope = "O"):

	if isotope == "O":
		col1="d(18_16)Mean"
		col2="d(18_16)mem_corrected"
		col3="d(18_16)drift_corrected"
		col4="d(18_16)vsmow_corrected"
		general_label= "d18O"
	else:
		col1="d(D_H)Mean"
		col2="d(D_H)mem_corrected"
		col3="d(D_H)drift_corrected"
		col4="d(D_H)vsmow_corrected"

		general_label= "d2H"

	val=[]
	
	for i in getSampleNames(df):
		if (i == "HAUS1") or (i=="HAUS2"):
			dat = df.loc[[i,'Standard']]
			ISOmean = dat.iloc[0:10][col1]
			ISOmem_corr = dat.iloc[0:10][col2]
			ISOdrift_corr = dat.iloc[0:10][col3]
			ISOvsmow_corr = dat.iloc[0:10][col4]
			ISO,ISO_mem,ISO_drift,ISO_smow = ISOmean.mean(),ISOmem_corr.mean(),ISOdrift_corr.mean(),ISOvsmow_corr.mean()
			sd,sd_mem,sd_drift,sd_smow = ISOmean.std(),ISOmem_corr.std(),ISOdrift_corr.std(),ISOvsmow_corr.std()
			val.append((i,"Standard",ISO,sd,ISO_mem,sd_mem,ISO_drift,sd_drift,ISO_smow,sd_smow))
			
		else:
			if i == "TAP":
				dat = df.loc[i]
				ISO = dat.loc["Conditioning"][col1]
				ISOmean = dat.loc["Conditioning"][col1]
				ISOmem_corr = dat.loc["Conditioning"][col2]
				ISOdrift_corr = dat.loc["Conditioning"][col3]
				ISOvsmow_corr = dat.loc["Conditioning"][col4]
				ISO,ISO_mem,ISO_drift,ISO_smow = ISOmean.mean(),ISOmem_corr.mean(),ISOdrift_corr.mean(),ISOvsmow_corr.mean()
				sd,sd_mem,sd_drift,sd_smow = ISOmean.std(),ISOmem_corr.std(),ISOdrift_corr.std(),ISOvsmow_corr.std()
				val.append((i,"Conditioning",ISO,sd,ISO_mem,sd_mem,ISO_drift,sd_drift,ISO_smow,sd_smow))

				
				ISO = dat.loc["Standard"][col1]
				ISOmean = dat.loc["Standard"][col1]
				ISOmem_corr = dat.loc["Standard"][col2]
				ISOdrift_corr = dat.loc["Standard"][col3]
				ISOvsmow_corr = dat.loc["Standard"][col4]
				ISO,ISO_mem,ISO_drift,ISO_smow = ISOmean.mean(),ISOmem_corr.mean(),ISOdrift_corr.mean(),ISOvsmow_corr.mean()
				sd,sd_mem,sd_drift,sd_smow = ISOmean.std(),ISOmem_corr.std(),ISOdrift_corr.std(),ISOvsmow_corr.std()
				val.append((i,"Standard",ISO,sd,ISO_mem,sd_mem,ISO_drift,sd_drift,ISO_smow,sd_smow))
				
				
				js = np.arange(0,3,1)
				for j in js:
					ISO = dat.loc['Control'].iloc[j*4:(j+1)]
					ISOmean = dat.loc["Control"].iloc[j*4:(j+1)*4][col1]
					ISOmem_corr = dat.loc["Control"].iloc[j*4:(j+1)*4][col2]
					ISOdrift_corr = dat.loc["Control"].iloc[j*4:(j+1)*4][col3]
					ISOvsmow_corr = dat.loc["Control"].iloc[j*4:(j+1)*4][col4]
					ISO,ISO_mem,ISO_drift,ISO_smow = ISOmean.mean(),ISOmem_corr.mean(),ISOdrift_corr.mean(),ISOvsmow_corr.mean()
					sd,sd_mem,sd_drift,sd_smow = ISOmean.std(),ISOmem_corr.std(),ISOdrift_corr.std(),ISOvsmow_corr.std()
					val.append((i,"Control",ISO,sd,ISO_mem,sd_mem,ISO_drift,sd_drift,ISO_smow,sd_smow))

			else:
				if i == 'W22':
					dat = df.loc[i]
					ISOmean = dat.loc["Control W22"][col1]
					ISOmem_corr = dat.loc["Control W22"][col2]
					ISOdrift_corr = dat.loc["Control W22"][col3]
					ISOvsmow_corr = dat.loc["Control W22"][col4]
					ISO,ISO_mem,ISO_drift,ISO_smow = ISOmean.mean(),ISOmem_corr.mean(),ISOdrift_corr.mean(),ISOvsmow_corr.mean()
					sd,sd_mem,sd_drift,sd_smow = ISOmean.std(),ISOmem_corr.std(),ISOdrift_corr.std(),ISOvsmow_corr.std()
					val.append((i,"Control W22",ISO,sd,ISO_mem,sd_mem,ISO_drift,sd_drift,ISO_smow,sd_smow))
				
				else:
					dat = df.loc[i]
					ISOmean = dat[col1]
					ISOmem_corr = dat[col2]
					ISOdrift_corr = dat[col3]
					ISOvsmow_corr = dat[col4]
					ISO,ISO_mem,ISO_drift,ISO_smow = ISOmean.mean(),ISOmem_corr.mean(),ISOdrift_corr.mean(),ISOvsmow_corr.mean()
					sd,sd_mem,sd_drift,sd_smow = ISOmean.std(),ISOmem_corr.std(),ISOdrift_corr.std(),ISOvsmow_corr.std()
					val.append((i,dat.index.values[0][0],ISO,sd,ISO_mem,sd_mem,ISO_drift,sd_drift,ISO_smow,sd_smow))

		
	return pd.DataFrame(val,columns=["Identifier 1",
		"Identifier 2",
		"{}_raw".format(general_label),
		"stdev. raw",
		"{} memory".format(general_label),
		"stdev. memory",
		"{} drift".format(general_label),
		"stdev. drift",
		"{} vsmow".format(general_label),
		"stdev. vsmow"])

def OLSR(x,y,lims,ax=None):
	xi=np.arange(lims[0],lims[1],0.25)
	slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
	
	ax.plot(xi,slope*xi+intercept,color='black',linewidth=0.5)
	return (ax,(slope,intercept))

def MemoryCorrection(df,coeff_dictionary,isotope = "O"):

	if isotope == "O":
		col1 = "d(18_16)mem_corrected"
		col2 = "d(18_16)Mean"
	else: 
		col1 = "d(D_H)mem_corrected"
		col2 = "d(D_H)Mean"
	
	coefficients = coeff_dictionary["mem_coeffs"] # coefficients indexed from 1
	
	def makecoefflist(coeffs):
		initlist= [coeffs[i] for i in coeffs]
		newlist = []
		lengths = []
		previousline = []
		counter=0
		
		for j in range(3):
			lengths.append(10)
		for i in range(23):
			lengths.append(4)
		
		for l in lengths:
			counter +=1
			if counter <=3:
				previous = counter*l
			else:
				previous = 40-16 + counter*l
			for i in range(l):
				newlist.append(initlist[i])
				previousline.append(previous)
		newlist =[1,1,1,1] + newlist
		previousline = [7,8,9,10] + previousline
		

		return newlist,previousline
	
	newlist, previouslines = makecoefflist(coefficients)
	
	df["coeffs"] = newlist
	df["previous"] = previouslines
	
	def correctLine(line):

		line_temp = df.iloc[line]
		x_orig = line_temp[col2]
		x_previous = df.iloc[int(line_temp["previous"])-7][col2]
		
		result= x_orig + (1-line_temp["coeffs"])*(x_orig-x_previous)
		
		return result
	
	corrected_values = []
	for i in range(len(df)):
		corrected_values.append(correctLine(i))
		
	df[col1] = corrected_values
	return df[["Line",col2,col1]]
	

def memCorrPlot(df,isotope= "O"):

	if isotope == "O":
		col1 = "d(18_16)Mean"
		col2 = "d(18_16)mem_corrected"
		posxy=((7,-3.5,"HAUS1"),(7,-31.5,"HAUS2"),(7,-18,"TAP"),(3,-7,"W22"))
	else:
		col1 = "d(D_H)Mean"
		col2 = "d(D_H)mem_corrected"
		posxy=((7,0,"HAUS1"),(7,-200,"HAUS2"),(7,-110,"TAP"),(3,-25,"W22"))

	df.sort_index(axis= 0)

	fig, (ax1,ax2,ax3,ax4) = plt.subplots(4,1,figsize = (6,9))

	y1 = df.loc["HAUS1"][col1]
	x =[i[1] for i in df.loc["HAUS1"][col1].index.values]
	y2 = df.loc["HAUS1"][col2]
	ax1.plot(x,y1,'o',color = 'black',markerfacecolor='w')
	ax1.plot(x,y2,'o',color = 'black',markersize = 3, markerfacecolor='black')
	OLSR(x,y2,(0,10),ax1)[0]
	ax1.text(posxy[0][0],posxy[0][1],posxy[0][2])

	y1 = df.loc["HAUS2"][col1]
	y2 = df.loc["HAUS2"][col2]
	ax2.plot(x,y1,'o',color = 'black',markerfacecolor='w')
	ax2.plot(x,y2,'o',color = 'black',markersize = 3, markerfacecolor='black')
	OLSR(x,y2,(0,10),ax2)[0]
	ax2.text(posxy[1][0],posxy[1][1],posxy[1][2])

	ytap = df.loc["TAP"].loc["Standard"]
	ytap_corr = ytap[col2]
	ax3.plot(x,ytap[col1],'o',color = 'black',markerfacecolor='w')
	ax3.plot(x,ytap_corr,'o',color = 'black',markersize = 3, markerfacecolor='black')
	OLSR(x,ytap_corr,(0,10),ax3)[0]
	ax3.text(posxy[2][0],posxy[2][1],posxy[2][2])

	y = df.loc["W22"][col1]
	x =[i[1] for i in df.loc["W22"][col1].index.values]
	y_corr = df.loc["W22"][col2]
	ax4.plot(x,y,'o',color = 'black',markerfacecolor='w')
	ax4.plot(x,y_corr,'o',color = 'black',markersize = 3, markerfacecolor='black')
	OLSR(x,y_corr,(0,4),ax4)[0]
	ax4.text(posxy[3][0],posxy[3][1],posxy[3][2])
	ax4.set_xlim(0,5)

	for ax in (ax1,ax2,ax3,ax4):
		ax.grid()
		ax.tick_params(direction='in',top=True,right=True)

	plt.show()

def getDriftcorrection(df,isotope = "O"):

	if isotope == "O":
		isotope = "d(18_16)mem_corrected"
	else :
		isotope = "d(D_H)mem_corrected"
	x = df.loc['TAP']['Line']
	y = df.loc['TAP'][isotope]

	slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)

	return {'slope':slope,'intercept':intercept}


def driftCorrPlot(df,isotope = "O"):

	if isotope =="O":
		col = 'd(18_16)mem_corrected'
		col2 = 'd(18_16)drift_corrected'
		posxy = (20,-16.4)
	else:
		col = 'd(D_H)mem_corrected'
		col2 = 'd(D_H)drift_corrected'
		posxy = (20,-90)

	x = df.loc['TAP']['Line']
	y1 = df.loc['TAP'][col]
	y2 = df.loc['TAP'][col2]

	fig,ax = plt.subplots()

	OLSR(x,y1,(0,132),ax)[0]

	params = OLSR(x,y1,(0,132),ax)[1]

	ax.plot(x,
			y1,
			'o',
			color = 'black',
			markersize = 3, 
			markerfacecolor='white')

	ax.plot(x,
			y2,
			'o',
			color = 'black',
			markersize = 3, 
			markerfacecolor='black')

	ax.text(posxy[0],posxy[1],"y {:.4f} = x {:+.4f}".format(params[0],params[1]))
	ax.grid()
	ax.tick_params(direction='in',top=True,right=True)
	plt.tight_layout()
	plt.show()

def driftCorrect(df,isotope = "O"):
	
	if isotope == "O":
		col1 = "d(18_16)mem_corrected"
		col2 = "d(18_16)drift_corrected"
	else :
		col1 = "d(D_H)mem_corrected"
		col2 = "d(D_H)drift_corrected"

	params = getDriftcorrection(df,isotope)


	df[col2] = df[col1] - df["Line"]*params["slope"]
	
	return df


def normToVSMOW(df,isotope = "O"):
	if isotope == "O":
		col = "d(18_16)drift_corrected"
		defined = {"HAUS1":0.6,
		   "HAUS2":-29.88,
		   "TAP":-13.4}
	else :
		col = "d(D_H)drift_corrected"
		defined= {"HAUS1":3.7,
		   "HAUS2":-229.8,
		   "TAP":-95.2}


	xhaus1= df.loc["HAUS1"][col].values
	xhaus2= df.loc["HAUS2"][col].values
	xtap = df.loc["TAP"].loc["Standard"][col].values
	yhaus1 = np.ones(len(xhaus1))*defined["HAUS1"]
	yhaus2 = np.ones(len(xhaus2))*defined["HAUS2"]
	ytap = np.ones(len(xtap))*defined["TAP"]

	yvals = []
	xvals = []


	for i,j in zip(xhaus1,yhaus1):
		xvals.append(i)
		yvals.append(j)
	for i,j in zip(xhaus2,yhaus2):
		xvals.append(i)
		yvals.append(j)
	for i,j in zip(xtap,ytap):
		xvals.append(i)
		yvals.append(j)
	
	slope, intercept, r_value, p_value, std_err = stats.linregress(xvals,yvals)

	return {'slope':slope,'intercept':intercept,"xvals":xvals,"yvals":yvals}

def VSMOWCorrPlot(df,isotope = "O"):

	if isotope=="O":
		posxy = (-15,-25)
	else:
		posxy = (-60,-200)

	params = normToVSMOW(df,isotope)

	fig, ax = plt.subplots(figsize= (6,4))
	xi = np.arange(-250,50,1)

	ax.plot(params["xvals"],params["yvals"],'o',color = 'black',markerfacecolor='w')
	ax.plot(xi,xi*params["slope"]+params["intercept"], color = 'black')

	ax.text(posxy[0],posxy[1],"y {:.3f} = x {:+.3f}".format(params["slope"],params["intercept"]))
	ax.grid()
	ax.tick_params(direction='in',top=True,right=True)
	plt.tight_layout()
	plt.show()

def VSMOWcorrect(df,isotope = "O"):

	if isotope == "O":
		col1="d(18_16)vsmow_corrected"
		col2="d(18_16)drift_corrected"
	else:
		col1="d(D_H)vsmow_corrected"
		col2="d(D_H)drift_corrected"

	params = normToVSMOW(df,isotope)

	df[col1] = df[col2]*params["slope"]+params["intercept"]
	return df

def getSD(df):
	sd1 = df.iloc[0]["stdev. memory"]
	
	sd2 = df.loc[1]["stdev. memory"]
	sd3 = df.iloc[22]["stdev. memory"]
	return [sd1,sd2,sd3]

def updateCoeffs(dictionary,df,mem_coeffs):
	mydict = dictionary

	coefficients = [i for i in mem_coeffs["x"]]

	mydict["mem_coeffs"]=dict([(i,j) for i,j in zip(range(0,len(coefficients)),coefficients)])

	return mydict

def updateSD(dictionary,df,isotope = "O"):
	mydict = dictionary
	stdevs = getMeanSDs(df,isotope)
	sds = getSD(stdevs)
	combined = sqrt(sds[0]**2+sds[1]**2+sds[2]**2)

	mydict["combined_sd"]=combined

	return mydict

from scipy.optimize import differential_evolution
from scipy.optimize import basinhopping,minimize
from scipy.optimize import brute

def Optimize(df,isotope="O"):

	if isotope == "O":
		col = "d(18_16)Mean"
	else:
		col = "d(D_H)Mean"

	def getAllcoeffs(df):
		df = df
		stds=[]
		for loc1 in ["TAP","HAUS1","HAUS2"]:
			loc2 = "Standard"
			stds.append(np.array( df.loc[loc1].loc[loc2][col].values))
	
		std1_0 = df.loc["TAP"].loc["Conditioning"][col].values[-1:]
		std2_0 =stds[0][-1:]
		std3_0 = stds[1][-1:]
	
		return stds,(std1_0,std2_0,std3_0)

	def F(x):
		vals = getAllcoeffs(df)
		val_temp=[]
		for i,j in zip(vals[0],vals[1]):
			temp = []
			for k in range(0,len(x)):
				xi = i[k]+(1-x[k])*(i[k] - j)
				temp.append(xi)
			
		
			temp = np.array(temp)	
			val_temp.append(temp)
		val_temp = np.array(val_temp)
	
		A = np.std(val_temp[0])
		B = np.std(val_temp[1])
		C = np.std(val_temp[2])
		#print(vals[0][0],val_temp[0])
	
		SD = np.sqrt(A**2+B**2+C**2)
		return(SD)

	#def jacobian(f,x):
		#out = approx_fprime(x,f,0.001)

		#return out
	x0 = get_initial_coeffs(0.5) 

	#print(jacobian(F,x0))
	bounds = [(0.5,1),(0.5,1),(0.5,1),(0.5,1),(0.8,1),
	(0.8,1),(0.8,1),(0.9,1),(0.9,1),(0.999999,1)]

	cons = ({'type':"ineq",
		'fun': lambda x: np.array(x[1]-x[0])},
		{'type':"ineq",
		'fun': lambda x: np.array(x[2]-x[1])},
		{'type':"ineq",
		'fun': lambda x: np.array(x[3]-x[2])},
		{'type':"ineq",
		'fun': lambda x: np.array(x[4]-x[3])},
		{'type':"ineq",
		'fun': lambda x: np.array(x[5]-x[4])},
		{'type':"ineq",
		'fun': lambda x: np.array(x[6]-x[5])},
		{'type':"ineq",
		'fun': lambda x: np.array(x[7]-x[6])},
		{'type':"ineq",
		'fun': lambda x: np.array(x[8]-x[7])},
		{'type':"eq",
		'fun': lambda x: np.array(x[9]-1)},
		{'type':"ineq",
		'fun': lambda x: np.array(x[0]-0.75)})

	xnew = minimize(F,x0 = x0,bounds = bounds,constraints = cons)
	#xnew = minimize(F,x0,method='BFGS')

	#x0 = get_initial_coeffs(0.5) 

	#xnew = basinhopping(F,x0 = x0)

	print(xnew)
	return xnew

import csv

def initialTreatment(filename,iso = "O"):

	dat = readRaw(filename)
	checkEmpty(dat)
	checkVolume(dat)

	run = runSummary(dat)

	df = IsotopeSelect(dat,iso)
	df = ignoreFirstSeven(df)
	SD_coeffs = initMemCoeffs()
	memory_optimized = Optimize(df,isotope = iso)

	SD_coeffs = updateCoeffs(SD_coeffs,df,memory_optimized)

	lines = [ (i+1,SD_coeffs["mem_coeffs"][i]) for i in SD_coeffs["mem_coeffs"]]

	with open('ipynb/coeffs{}.csv'.format(iso), 'w') as coeffFile:
		writer = csv.writer(coeffFile)

	
		writer.writerows(lines)
	coeffFile.close()

	return run,df,SD_coeffs

def readCoeffs(filename):
    SD_coeffs = initMemCoeffs()
    updated = pd.read_csv(filename,header = None)

    SD_coeffs["mem_coeffs"]= dict([(i,j) for i,j in zip(updated[0],updated[1])])
        
    return SD_coeffs


def secondTreatment(run,df,iso = "O",option ="Plot"):

	SD_coeffs = readCoeffs("ipynb/coeffs{}.csv".format(iso))

	df = MemoryCorrection(df,SD_coeffs,isotope = iso)

	df = driftCorrect(df,isotope = iso)
	df = VSMOWcorrect(df,isotope = iso)

	SD_Coeffs = updateSD(SD_coeffs,df, isotope = iso)

	if option == "Plot":
		plotSummary(run)
		memCorrPlot(df,isotope = iso)
		driftCorrPlot(df,isotope = iso)
		VSMOWCorrPlot(df,isotope = iso)

	return run,df,SD_coeffs




