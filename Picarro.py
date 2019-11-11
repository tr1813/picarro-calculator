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

def getMeanSDs(df):
	val=[]
	
	for i in getSampleNames(df):
		if (i == "HAUS1") or (i=="HAUS2"):
			dat = df.loc[[i,'Standard']]
			O18mean = dat.iloc[0:10]["d(18_16)Mean"]
			O18mem_corr = dat.iloc[0:10]["d(18_16)mem_corrected"]
			O18drift_corr = dat.iloc[0:10]["d(18_16)drift_corrected"]
			O18vsmow_corr = dat.iloc[0:10]["d(18_16)vsmow_corrected"]
			O18,O18_mem,O18_drift,O18_smow = O18mean.mean(),O18mem_corr.mean(),O18drift_corr.mean(),O18vsmow_corr.mean()
			sd,sd_mem,sd_drift,sd_smow = O18mean.std(),O18mem_corr.std(),O18drift_corr.std(),O18vsmow_corr.std()
			val.append((i,"Standard",O18,sd,O18_mem,sd_mem,O18_drift,sd_drift,O18_smow,sd_smow))
			
		else:
			if i == "TAP":
				dat = df.loc[i]
				O18 = dat.loc["Conditioning"]["d(18_16)Mean"]
				O18mean = dat.loc["Conditioning"]["d(18_16)Mean"]
				O18mem_corr = dat.loc["Conditioning"]["d(18_16)mem_corrected"]
				O18drift_corr = dat.loc["Conditioning"]["d(18_16)drift_corrected"]
				O18vsmow_corr = dat.loc["Conditioning"]["d(18_16)vsmow_corrected"]
				O18,O18_mem,O18_drift,O18_smow = O18mean.mean(),O18mem_corr.mean(),O18drift_corr.mean(),O18vsmow_corr.mean()
				sd,sd_mem,sd_drift,sd_smow = O18mean.std(),O18mem_corr.std(),O18drift_corr.std(),O18vsmow_corr.std()
				val.append((i,"Conditioning",O18,sd,O18_mem,sd_mem,O18_drift,sd_drift,O18_smow,sd_smow))

				
				O18 = dat.loc["Standard"]["d(18_16)Mean"]
				O18mean = dat.loc["Standard"]["d(18_16)Mean"]
				O18mem_corr = dat.loc["Standard"]["d(18_16)mem_corrected"]
				O18drift_corr = dat.loc["Standard"]["d(18_16)drift_corrected"]
				O18vsmow_corr = dat.loc["Standard"]["d(18_16)vsmow_corrected"]
				O18,O18_mem,O18_drift,O18_smow = O18mean.mean(),O18mem_corr.mean(),O18drift_corr.mean(),O18vsmow_corr.mean()
				sd,sd_mem,sd_drift,sd_smow = O18mean.std(),O18mem_corr.std(),O18drift_corr.std(),O18vsmow_corr.std()
				val.append((i,"Standard",O18,sd,O18_mem,sd_mem,O18_drift,sd_drift,O18_smow,sd_smow))
				
				
				js = np.arange(0,3,1)
				for j in js:
					O18 = dat.loc['Control'].iloc[j*4:(j+1)]
					O18mean = dat.loc["Control"].iloc[j*4:(j+1)*4]["d(18_16)Mean"]
					O18mem_corr = dat.loc["Control"].iloc[j*4:(j+1)*4]["d(18_16)mem_corrected"]
					O18drift_corr = dat.loc["Control"].iloc[j*4:(j+1)*4]["d(18_16)drift_corrected"]
					O18vsmow_corr = dat.loc["Control"].iloc[j*4:(j+1)*4]["d(18_16)vsmow_corrected"]
					O18,O18_mem,O18_drift,O18_smow = O18mean.mean(),O18mem_corr.mean(),O18drift_corr.mean(),O18vsmow_corr.mean()
					sd,sd_mem,sd_drift,sd_smow = O18mean.std(),O18mem_corr.std(),O18drift_corr.std(),O18vsmow_corr.std()
					val.append((i,"Control {}".format(j+1),O18,sd,O18_mem,sd_mem,O18_drift,sd_drift,O18_smow,sd_smow))

			else:
				if i == 'W22':
					dat = df.loc[i]
					O18mean = dat.loc["Control W22"]["d(18_16)Mean"]
					O18mem_corr = dat.loc["Control W22"]["d(18_16)mem_corrected"]
					O18drift_corr = dat.loc["Control W22"]["d(18_16)drift_corrected"]
					O18vsmow_corr = dat.loc["Control W22"]["d(18_16)vsmow_corrected"]
					O18,O18_mem,O18_drift,O18_smow = O18mean.mean(),O18mem_corr.mean(),O18drift_corr.mean(),O18vsmow_corr.mean()
					sd,sd_mem,sd_drift,sd_smow = O18mean.std(),O18mem_corr.std(),O18drift_corr.std(),O18vsmow_corr.std()
					val.append((i,"Control W22",O18,sd,O18_mem,sd_mem,O18_drift,sd_drift,O18_smow,sd_smow))
				
				else:
					dat = df.loc[i]
					O18mean = dat["d(18_16)Mean"]
					O18mem_corr = dat["d(18_16)mem_corrected"]
					O18drift_corr = dat["d(18_16)drift_corrected"]
					O18vsmow_corr = dat["d(18_16)vsmow_corrected"]
					O18,O18_mem,O18_drift,O18_smow = O18mean.mean(),O18mem_corr.mean(),O18drift_corr.mean(),O18vsmow_corr.mean()
					sd,sd_mem,sd_drift,sd_smow = O18mean.std(),O18mem_corr.std(),O18drift_corr.std(),O18vsmow_corr.std()
					val.append((i,dat.index.values[0][0],O18,sd,O18_mem,sd_mem,O18_drift,sd_drift,O18_smow,sd_smow))

		
	return pd.DataFrame(val,columns=["Identifier 1","Identifier 2","O18_raw","stdev. raw","O18 memory","stdev. memory","O18 drift","stdev. drift","O18 vsmow","stdev. vsmow"])

def coeffCorrect(df,line,memcoeffs):

	mem_coeffs= memcoeffs
	x_after =df.iloc[line].values[1]
	x_before =df.iloc[line-1].values[1]
	inj = df.index.values[line][2]
	coeff = mem_coeffs[inj]
	
	corrected = x_after+(1-coeff)*(x_after-x_before)
	return corrected

def OLSR(x,y,lims,ax=None):
	xi=np.arange(lims[0],lims[1],0.25)
	slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
	
	ax.plot(xi,slope*xi+intercept,color='black',linewidth=0.5)
	return (ax,(slope,intercept))

def computeMemCorrected(df,memcoeffs): #where memcoeffs is a dictionary of coefficients
	memCorrected = []
	for i in (0,1,2):
		memCorrected.append(df.iloc[i].values[1])
	for i in np.arange(3,126,1):
		memCorrected.append(coeffCorrect(df,i,memcoeffs))
	
	return memCorrected

def memCorrPlot(df):

	df.sort_index(axis= 0)

	fig, (ax1,ax2,ax3,ax4) = plt.subplots(4,1,figsize = (6,9))

	y1 = df.loc["HAUS1"]["d(18_16)Mean"]
	x =[i[1] for i in df.loc["HAUS1"]["d(18_16)Mean"].index.values]
	y2 = df.loc["HAUS1"]["d(18_16)mem_corrected"]
	ax1.plot(x,y1,'o',color = 'black',markerfacecolor='w')
	ax1.plot(x,y2,'o',color = 'black',markersize = 3, markerfacecolor='black')
	OLSR(x,y2,(0,10),ax1)[0]
	ax1.text(7,-3.5,"HAUS1")

	y1 = df.loc["HAUS2"]["d(18_16)Mean"]
	y2 = df.loc["HAUS2"]["d(18_16)mem_corrected"]
	ax2.plot(x,y1,'o',color = 'black',markerfacecolor='w')
	ax2.plot(x,y2,'o',color = 'black',markersize = 3, markerfacecolor='black')
	OLSR(x,y2,(0,10),ax2)[0]
	ax2.text(7,-31.5,"HAUS2")

	ytap = df.loc["TAP"].loc["Standard"]
	ytap_corr = ytap["d(18_16)mem_corrected"]
	ax3.plot(x,ytap["d(18_16)Mean"],'o',color = 'black',markerfacecolor='w')
	ax3.plot(x,ytap_corr,'o',color = 'black',markersize = 3, markerfacecolor='black')
	OLSR(x,ytap_corr,(0,10),ax3)[0]
	ax3.text(7,-18,"TAP")

	y = df.loc["W22"]["d(18_16)Mean"]
	x =[i[1] for i in df.loc["W22"]["d(18_16)Mean"].index.values]
	y_corr = df.loc["W22"]["d(18_16)mem_corrected"]
	ax4.plot(x,y,'o',color = 'black',markerfacecolor='w')
	ax4.plot(x,y_corr,'o',color = 'black',markersize = 3, markerfacecolor='black')
	OLSR(x,y_corr,(0,4),ax4)[0]
	ax4.text(3,-7,"W22")
	ax4.set_xlim(0,4)

	for ax in (ax1,ax2,ax3,ax4):
		ax.grid()
		ax.tick_params(direction='in',top=True,right=True)

	plt.show()

def getDriftcorrection(df):
	x = df.loc['TAP']['Line']
	y = df.loc['TAP']['d(18_16)mem_corrected']

	slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)

	return {'slope':slope,'intercept':intercept}


def driftCorrPlot(df):
	x = df.loc['TAP']['Line']
	y = df.loc['TAP']['d(18_16)mem_corrected']

	fig,ax = plt.subplots()

	OLSR(x,y,(0,132),ax)[0]

	params = OLSR(x,y,(0,132),ax)[1]

	ax.plot(x,
			y,
			'o',
			color = 'black',
			markersize = 3, 
			markerfacecolor='black')

	ax.text(80,-16.8,"y {:.4f} = x {:+.4f}".format(params[0],params[1]))
	ax.grid()
	ax.tick_params(direction='in',top=True,right=True)
	plt.tight_layout()
	plt.show()

def driftCorrect(df):
	
	params = getDriftcorrection(df)
	df["d(18_16)drift_corrected"] = df["d(18_16)mem_corrected"] - df["Line"]*params["slope"]
	
	return df


def normToVSMOW(df):

	defined = {"HAUS1":0.6,
		   "HAUS2":-29.88,
		   "TAP":-13.4}

	xhaus1= df.loc["HAUS1"]["d(18_16)drift_corrected"].values
	xhaus2= df.loc["HAUS2"]["d(18_16)drift_corrected"].values
	xtap = df.loc["TAP"].loc["Standard"]["d(18_16)drift_corrected"].values
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

def VSMOWCorrPlot(df):

	params = normToVSMOW(df)

	fig, ax = plt.subplots(figsize= (6,4))
	xi = np.arange(-50,30,1)

	ax.plot(params["xvals"],params["yvals"],'o',color = 'black',markerfacecolor='w')
	ax.plot(xi,xi*params["slope"]+params["intercept"], color = 'black')
	ax.set_xlim(-40,0)
	ax.set_ylim(-40,20)


	ax.text(-10,-30,"y {:.3f} = x {:+.3f}".format(params["slope"],params["intercept"]))
	ax.grid()
	ax.tick_params(direction='in',top=True,right=True)
	plt.tight_layout()
	plt.show()

def VSMOWcorrect(df):

	params = normToVSMOW(df)

	df["d(18_16)vsmow_corrected"] = df["d(18_16)drift_corrected"]*params["slope"]+params["intercept"]
	return df

def getSD(df):
	sd1 = df.iloc[0]["stdev. memory"]
	
	sd2 = df.loc[1]["stdev. memory"]
	sd3 = df.iloc[22]["stdev. memory"]
	return [sd1,sd2,sd3]

def updateSD(dictionary,df):
	mydict = dictionary
	stdevs = getMeanSDs(df)
	sds = getSD(stdevs)
	combined = sqrt(sds[0]**2+sds[1]**2+sds[2]**2)

	mydict["combined_sd"]=combined

	return mydict


def initialTreatment(filename,option=None):

	dat = readRaw(filename)
	checkEmpty(dat)
	checkVolume(dat)

	run = runSummary(dat)

	SD_coeffs = initMemCoeffs()

	d18O = IsotopeSelect(dat,"O")
	d18O = ignoreFirstSeven(d18O)
	d18O["d(18_16)mem_corrected"] = computeMemCorrected(d18O,SD_coeffs["mem_coeffs"])
	d18O = driftCorrect(d18O)
	d18O = VSMOWcorrect(d18O)

	SD_coeffs = updateSD(SD_coeffs,d18O)
	if option == "Plot":
		plotSummary(run)
		memCorrPlot(d18O)
		driftCorrPlot(d18O)
		VSMOWCorrPlot(d18O)

	return run,d18O,SD_coeffs




