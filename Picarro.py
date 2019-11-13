import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from math import exp,log,sqrt
from scipy.optimize import minimize,differential_evolution
import os
import csv


class Isotope(object):
	"""docstring for Isotope
	A object of data"""
	def __init__(self,filename):
		
		self.filename = filename
		self.log = []
		self.raw = None
		self.summary = None
		self.corr = None
		self.coeffs = None
		self.SD_coeffs = None
		self.memory = None
		self.drift = None
		self.vsmow = None
		self.drift_params = None
		self.VSMOW_params = None
		self.run_overview = None
		self.run_overview = None
		self.defined_O = None
		self.defined_H = None
		self.combined_sd = None
		self.constraint = None
		self.HAUS1_raw_sd = None
		self.HAUS2_raw_sd = None
		self.TAP_raw_sd = None
		self.HAUS1_sd = None
		self.HAUS2_sd = None
		self.TAP_sd = None

	def Log(self):
		for i in self.log:
			print(i+"\n")
	def readRaw(self):
		df = pd.read_csv(self.filename)

		renamed = dict([(i,i.strip()) for i in df]) # strip white space from column names
		df.rename(columns=renamed, inplace=True)	# and rename all headers of the dataFrame with new 'clean names'
		df["Identifier 1"] = [i.strip() for i in df["Identifier 1"].values]
		df["Identifier 2"] = [i.strip() for i in df["Identifier 2"].values]
		df["Line"] = [int(i) for i in df["Line"]]
		df.set_index(["Identifier 1","Identifier 2","Inj Nr"], inplace=True)

		self.raw = df

	def checkEmpty(self):

		if np.any(np.isnan(self.raw["H2O_Mean"].values)) == 'True':  # checks a boolean array which evaluates to True if
															# if one cell is empty
		
			self.log.append("Warning: there is at least one empty line!")
		else: 
			self.log.append("No empty cells. Proceeding...")

	def checkVolume(self):

		conds = [self.raw.H2O_Mean.values < 17000 , self.raw.H2O_Mean.values > 23000]	# create the conditions for ignoring a certain analysis line
		choices = [1, 2]													# what the error code value will become

		out = np.select(conds, choices, default=0)							# mapping the error code on the conditions

		self.raw["Error Code"]=out												# updating the dataFrame

		above = self.raw.where(self.raw["Error Code"] == 2)
		below = self.raw.where(self.raw["Error Code"] == 1)
		warning = pd.concat([above,below]).dropna()

		lines =[i for i in warning["Line"].values]
		for i in lines:
			self.log.append("Warning: H20 value outside bounds on line {:.0f}...".format(i))

	def runSummary(self):

		self.summary = self.raw[["Line",
					   "H2O_Mean",
					   "Good",
					   "Error Code",
					   "d(18_16)Mean",
					   "d(D_H)Mean"]]
		
		self.summary = self.summary.sort_index(axis=0)

	def plotSummary(self):

		df = self.summary
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

	def IsotopeSelect(self,isotope="O"):
		if isotope == "O":
			self.corr =  self.raw[["Line","d(18_16)Mean","Ignore","Error Code"]]
			self.corr = self.corr.where(self.corr["Line"]>6).dropna()
		else:
			if isotope == "H":
				self.corr =  self.raw[["Line","d(D_H)Mean","Ignore","Error Code"]]
				self.corr = self.corr.where(self.corr["Line"]>6).dropna()
			else:
				self.log.append("You must choose a valid isotope to select")


	def initMemCoeffs(self):

		def get_k(i):
			k = exp(-log(i)/9)
			return k

		def get_initial_coeffs(i):
			coeffs=[]
			c=0
			k = get_k(i)
			while c <10:
					
				coeffs.append(round(i*k**c,7))			
				c+=1
			return coeffs

		injection = np.arange(0,10,1)
		coeffs =get_initial_coeffs(0.849029498880193) 

		self.coeffs =dict([(i+1,j) for i,j in zip(injection,coeffs)])


		self.SD_coeffs={
			"combined_sd":"None",
			"date":'None',
			"mem_coeffs":self.coeffs
		}


	def Optimize(self,isotope="O",method = 'default'):

		self.log.append("Checking a coeffs whether a coeffs file already exists")

		if os.path.isfile("ipynb/coeffs{}_{}.csv".format(isotope,self.filename[31:47]))==True:
			self.log.append("coeffs file already exists!")
		else:
			self.log.append("coeffs file does not exist yet. Switch to Optimization...")



		self.log.append("Now running the Optimization algorithm to minimize the combined standard deviation")

		if isotope == "O":
			col = "d(18_16)Mean"
		else:
			if isotope =="H":
				col = "d(D_H)Mean"
			else:
				self.log.append("You must choose a valid isotope")

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

		def F2(x):
			vals = getAllcoeffs(self.corr)
			val_temp=[]
			for i,j in zip(vals[0],vals[1]):
				temp = []
				for k in range(0,len(x)):
					xi = i[k]+(1-x[k])*(i[k] - j)
					temp.append(xi)
				
			
				temp = np.array(temp)	
				val_temp.append(temp)
			val_temp = np.array(val_temp)
		
			self.HAUS1_sd = np.std(val_temp[0])
			self.HAUS2_sd = np.std(val_temp[1])
			self.TAP_sd = np.std(val_temp[2])
		
			SD = (self.HAUS1_sd**2+self.HAUS2_sd**2+self.TAP_sd**2)**(0.5)

			return(SD)

		def F(x):
			vals = getAllcoeffs(self.corr)
			val_temp=[]
			for i,j in zip(vals[0],vals[1]):
				temp = []
				for k in range(0,len(x)):
					xi = i[k]+(1-x[k])*(i[k] - j)
					temp.append(xi)
				
			
				temp = np.array(temp)	
				val_temp.append(temp)
			val_temp = np.array(val_temp)
		
			self.HAUS1_raw_sd = np.std(val_temp[0])
			self.HAUS2_raw_sd = np.std(val_temp[1])
			self.TAP_raw_sd = np.std(val_temp[2])
		
			SD = (self.HAUS1_raw_sd+self.HAUS2_raw_sd+self.TAP_raw_sd)**(0.5)

			return(SD)

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

		x0= np.array([self.coeffs[i] for i in self.coeffs])

		self.log.append("Running SQSLP algorithm")
		xnew = minimize(F,x0 = x0,bounds = bounds,constraints = cons) # the main workhorse?
		self.log.append("Done")


		
		local = dict([(i,j) for i,j in zip(range(1,len(xnew["x"])+1),xnew["x"])])
		self.coeffs = local
		self.combined_raw_sd = F2(x0)
		self.combined_sd = F2(xnew["x"])

		def writeCoeffsfile(coeffs,iso = isotope):

			lines = [ (i+1,coeffs[i]) for i in coeffs]

			with open('ipynb/coeffs{}_{}.csv'.format(iso,self.filename[31:47]), 'w') as coeffFile:
				writer = csv.writer(coeffFile)

	
				writer.writerows(lines)
			coeffFile.close()

		writeCoeffsfile(self.coeffs)


	def getUpdatedSD(self):

		df= pd.DataFrame({"old":(self.combined_raw_sd),
			"updated":(self.combined_sd)})
		return df


	def MemoryCorrection(self,isotope = "O"):

		if isotope == "O":
			col1 = "d(18_16)mem_corrected"
			col2 = "d(18_16)Mean"
		else: 
			col1 = "d(D_H)mem_corrected"
			col2 = "d(D_H)Mean"
		
		coefficients = [self.coeffs[i] for i in self.coeffs] # coefficients indexed from 1
		
		def makecoefflist(coeffs):
			initlist= coeffs
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
		
		self.corr["coeffs"] = newlist
		self.corr["previous"] = previouslines
		
		def correctLine(line):

			line_temp = self.corr.iloc[line]
			x_orig = line_temp[col2]
			x_previous = self.corr.iloc[int(line_temp["previous"])-7][col2]
			
			result= x_orig + (1-line_temp["coeffs"])*(x_orig-x_previous)
			
			return result
		
		corrected_values = []
		for i in range(len(self.corr)):
			corrected_values.append(correctLine(i))
			
		self.corr[col1] = corrected_values
		self.log.append('Successfully corrected values for memory effects')
		self.memory = self.corr[["Line",col2,col1,"Error Code"]]

	def OLSR(self,x,y,lims,ax=None):
		xi=np.arange(lims[0],lims[1],0.25)
		slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
		
		ax.plot(xi,slope*xi+intercept,color='black',linewidth=0.5)
		return (ax,(slope,intercept))

	def memCorrPlot(self,isotope= "O"):
		self.log.append("Plotting the results of memory correction")

		if isotope == "O":
			col1 = "d(18_16)Mean"
			col2 = "d(18_16)mem_corrected"
			posxy=((7,-3.5,"HAUS1"),(7,-31.5,"HAUS2"),(7,-18,"TAP"),(3.2,-7,"W22"))
		else:
			col1 = "d(D_H)Mean"
			col2 = "d(D_H)mem_corrected"
			posxy=((7,0,"HAUS1"),(7,-200,"HAUS2"),(7,-110,"TAP"),(3.2,-25,"W22"))

		df = self.memory

		fig, (ax1,ax2,ax3,ax4) = plt.subplots(4,1,figsize = (6,9))

		y1 = df.loc["HAUS1"][col1]
		x =[i[1] for i in df.loc["HAUS1"][col1].index.values]
		y2 = df.loc["HAUS1"][col2]
		ax1.plot(x,y1,'o',color = 'black',markerfacecolor='w')
		ax1.plot(x,y2,'o',color = 'black',markersize = 3, markerfacecolor='black')
		self.OLSR(x,y2,(0,10),ax1)[0]
		ax1.text(posxy[0][0],posxy[0][1],posxy[0][2])

		y1 = df.loc["HAUS2"][col1]
		y2 = df.loc["HAUS2"][col2]
		ax2.plot(x,y1,'o',color = 'black',markerfacecolor='w')
		ax2.plot(x,y2,'o',color = 'black',markersize = 3, markerfacecolor='black')
		self.OLSR(x,y2,(0,10),ax2)[0]
		ax2.text(posxy[1][0],posxy[1][1],posxy[1][2])

		ytap = df.loc["TAP"].loc["Standard"]
		ytap_corr = ytap[col2]
		ax3.plot(x,ytap[col1],'o',color = 'black',markerfacecolor='w')
		ax3.plot(x,ytap_corr,'o',color = 'black',markersize = 3, markerfacecolor='black')
		self.OLSR(x,ytap_corr,(0,10),ax3)[0]
		ax3.text(posxy[2][0],posxy[2][1],posxy[2][2])

		y = df.loc["W22"][col1]
		x =[i[1] for i in df.loc["W22"][col1].index.values]
		y_corr = df.loc["W22"][col2]
		ax4.plot(x,y,'o',color = 'black',markerfacecolor='w')
		ax4.plot(x,y_corr,'o',color = 'black',markersize = 3, markerfacecolor='black')
		self.OLSR(x,y_corr,(0,4),ax4)[0]
		ax4.text(posxy[3][0],posxy[3][1],posxy[3][2])
		ax4.set_xlim(0,5)

		for ax in (ax1,ax2,ax3,ax4):
			ax.grid()
			ax.tick_params(direction='in',top=True,right=True)

		plt.show()

	def driftCorrect(self,isotope = "O"):
		self.drift = self.memory.copy()

		def getDriftcorrection(df,isotope = "O"):

			if isotope == "O":
				name = "d(18_16)mem_corrected"
			else :
				name = "d(D_H)mem_corrected"
			x = df.loc['TAP']['Line']
			y = df.loc['TAP'][name]

			slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)

			return  {'slope':slope,'intercept':intercept}
		
		self.drift_params = getDriftcorrection(self.drift,isotope)

		df = self.drift

		if isotope == "O":
			self.drift["d(18_16)drift_corrected"] = self.drift["d(18_16)mem_corrected"]-self.drift["Line"]*self.drift_params["slope"]
		else :
			self.drift["d(D_H)drift_corrected"] = self.drift["d(D_H)mem_corrected"]-self.drift["Line"]*self.drift_params["slope"]


		
		self.log.append("Successfully corrected for drift")

	def driftCorrPlot(self,isotope = "O"):

		if isotope =="O":
			col = 'd(18_16)mem_corrected'
			col2 = 'd(18_16)drift_corrected'
			posxy = (20,-16.4)
		else:
			col = 'd(D_H)mem_corrected'
			col2 = 'd(D_H)drift_corrected'
			posxy = (20,-90)

		x = self.drift.loc['TAP']['Line']
		y1 = self.drift.loc['TAP'][col]
		y2 = self.drift.loc['TAP'][col2]

		fig,ax = plt.subplots()

		self.OLSR(x,y1,(0,132),ax)[0]

		params = self.OLSR(x,y1,(0,132),ax)[1]

		ax.plot(x,
				y1,
				'o',
				color = 'black',
				markersize = 3, 
				markerfacecolor='white',label =" no drift correction")

		ax.plot(x,
				y2,
				'o',
				color = 'black',
				markersize = 3, 
				markerfacecolor='black',label="drift corrected")

		ax.legend()

		ax.text(posxy[0],posxy[1],"y {:.4f} = x {:+.4f}".format(params[0],params[1]))
		ax.grid()
		ax.tick_params(direction='in',top=True,right=True)
		plt.tight_layout()
		plt.show()

	def VSMOWcorrect(self,isotope = "O"):
		self.vsmow = self.drift.copy()

		if isotope == "O":
			col1="d(18_16)vsmow_corrected"
			col2="d(18_16)drift_corrected"

		else:
			col1="d(D_H)vsmow_corrected"
			col2="d(D_H)drift_corrected"


		def normToVSMOW(df,isotope = "O"):
			if isotope == "O":
				col = "d(18_16)drift_corrected"
				self.defined_O = {"HAUS1":0.6,
				   "HAUS2":-29.88,
				   "TAP":-13.4}
				xhaus1= df.loc["HAUS1"][col].values
				xhaus2= df.loc["HAUS2"][col].values
				xtap = df.loc["TAP"].loc["Standard"][col].values
				yhaus1 = np.ones(len(xhaus1))*self.defined_O["HAUS1"]
				yhaus2 = np.ones(len(xhaus2))*self.defined_O["HAUS2"]
				ytap = np.ones(len(xtap))*self.defined_O["TAP"]
			else :
				col = "d(D_H)drift_corrected"
				self.defined_H= {"HAUS1":3.7,
				   "HAUS2":-229.8,
				   "TAP":-95.2}
				xhaus1= df.loc["HAUS1"][col].values
				xhaus2= df.loc["HAUS2"][col].values
				xtap = df.loc["TAP"].loc["Standard"][col].values
				yhaus1 = np.ones(len(xhaus1))*self.defined_H["HAUS1"]
				yhaus2 = np.ones(len(xhaus2))*self.defined_H["HAUS2"]
				ytap = np.ones(len(xtap))*self.defined_H["TAP"]


			
			

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

			self.VSMOW_params = {'slope':slope,'intercept':intercept,"xvals":xvals,"yvals":yvals}


		normToVSMOW(self.vsmow,isotope)

		params = self.VSMOW_params

		results = []

		for i in self.vsmow[col2].values:
			results.append(i*params["slope"]+params["intercept"])

		self.vsmow[col1] = results
		self.log.append("Sucessfully calibrated to VSMOW scale")

	def VSMOWCorrPlot(self,isotope = "O"):

		if isotope=="O":
			posxy = (-15,-25)
		else:
			posxy = (-60,-200)

		params = self.VSMOW_params

		fig, ax = plt.subplots(figsize= (6,4))
		xi = np.arange(-250,50,1)

		ax.plot(params["xvals"],params["yvals"],'o',color = 'black',markerfacecolor='w')
		ax.plot(xi,xi*params["slope"]+params["intercept"], color = 'black')

		ax.text(posxy[0],posxy[1],"y {:.3f} = x {:+.3f}".format(params["slope"],params["intercept"]))
		ax.grid()
		ax.tick_params(direction='in',top=True,right=True)
		plt.tight_layout()
		plt.show()

	def getMeanSDs(self,isotope = "O"):

		df = self.vsmow.copy()

		if isotope == "O":
			col1="d(18_16)Mean"
			col2="d(18_16)mem_corrected"
			col3="d(18_16)drift_corrected"
			col4="d(18_16)vsmow_corrected"
			general_label= "d18O"
			limit = 0.1
		else:
			col1="d(D_H)Mean"
			col2="d(D_H)mem_corrected"
			col3="d(D_H)drift_corrected"
			col4="d(D_H)vsmow_corrected"
			general_label= "d2H"
			limit = 0.7

		val=[]


		def getSampleNames(df):
			return df.index.levels[0].values  
		
		for i in getSampleNames(self.vsmow):
			if (i == "HAUS1") or (i=="HAUS2"):
				dat = df.where(df["Error Code"] == 0)
				dat = dat.loc[[i,'Standard']]
				ISOmean = dat.iloc[0:10][col1]
				ISOmem_corr = dat.iloc[0:10][col2]
				ISOdrift_corr = dat.iloc[0:10][col3]
				ISOvsmow_corr = dat.iloc[0:10][col4]
				count,ISO,ISO_mem,ISO_drift,ISO_smow = ISOmean.count(),ISOmean.mean(),ISOmem_corr.mean(),ISOdrift_corr.mean(),ISOvsmow_corr.mean()
				sd,sd_mem,sd_drift,sd_smow = ISOmean.std(),ISOmem_corr.std(),ISOdrift_corr.std(),ISOvsmow_corr.std()
				val.append((i,"_Standard",ISO,sd,ISO_mem,sd_mem,ISO_drift,sd_drift,ISO_smow,sd_smow,count))
				if sd_smow >= limit:
					self.log.append("Warning: high standard deviation on sample {}".format(i))

				
			else:
				if i == "TAP":
					dat = df.where(df["Error Code"] == 0)
					dat = dat.loc[i]

					ISO = dat.loc["Conditioning"][col1]
					ISOmean = dat.loc["Conditioning"][col1]
					ISOmem_corr = dat.loc["Conditioning"][col2]
					ISOdrift_corr = dat.loc["Conditioning"][col3]
					ISOvsmow_corr = dat.loc["Conditioning"][col4]
					count,ISO,ISO_mem,ISO_drift,ISO_smow = ISOmean.count(),ISOmean.mean(),ISOmem_corr.mean(),ISOdrift_corr.mean(),ISOvsmow_corr.mean()
					sd,sd_mem,sd_drift,sd_smow = ISOmean.std(),ISOmem_corr.std(),ISOdrift_corr.std(),ISOvsmow_corr.std()
					val.append((i,"_Conditioning",ISO,sd,ISO_mem,sd_mem,ISO_drift,sd_drift,ISO_smow,sd_smow,count))
					if sd_smow >= limit:
						self.log.append("Warning: high standard deviation on sample {}".format(i))


					
					ISO = dat.loc["Standard"][col1]
					ISOmean = dat.loc["Standard"][col1]
					ISOmem_corr = dat.loc["Standard"][col2]
					ISOdrift_corr = dat.loc["Standard"][col3]
					ISOvsmow_corr = dat.loc["Standard"][col4]
					count,ISO,ISO_mem,ISO_drift,ISO_smow = ISOmean.count(),ISOmean.mean(),ISOmem_corr.mean(),ISOdrift_corr.mean(),ISOvsmow_corr.mean()
					sd,sd_mem,sd_drift,sd_smow = ISOmean.std(),ISOmem_corr.std(),ISOdrift_corr.std(),ISOvsmow_corr.std()
					val.append((i,"_Standard",ISO,sd,ISO_mem,sd_mem,ISO_drift,sd_drift,ISO_smow,sd_smow,count))
					if sd_smow >= limit:
						self.log.append("Warning: high standard deviation on sample {}".format(i))

					
					
					js = np.arange(0,3,1)
					for j in js:

						ISO = dat.loc['Control'].iloc[j*4:(j+1)]
						ISOmean = dat.loc["Control"].iloc[j*4:(j+1)*4][col1]
						ISOmem_corr = dat.loc["Control"].iloc[j*4:(j+1)*4][col2]
						ISOdrift_corr = dat.loc["Control"].iloc[j*4:(j+1)*4][col3]
						ISOvsmow_corr = dat.loc["Control"].iloc[j*4:(j+1)*4][col4]
						count,ISO,ISO_mem,ISO_drift,ISO_smow = ISOmean.count(),ISOmean.mean(),ISOmem_corr.mean(),ISOdrift_corr.mean(),ISOvsmow_corr.mean()
						sd,sd_mem,sd_drift,sd_smow = ISOmean.std(),ISOmem_corr.std(),ISOdrift_corr.std(),ISOvsmow_corr.std()
						val.append((i,"_Control {}".format(j+1),ISO,sd,ISO_mem,sd_mem,ISO_drift,sd_drift,ISO_smow,sd_smow,count))
					if sd_smow >= limit:
						self.log.append("Warning: high standard deviation on sample {}".format(i))

				else:
					if i == 'W22':

						dat = df.where(df["Error Code"] == 0)
						dat = dat.loc[i]
						ISOmean = dat.loc["Control W22"][col1]
						ISOmem_corr = dat.loc["Control W22"][col2]
						ISOdrift_corr = dat.loc["Control W22"][col3]
						ISOvsmow_corr = dat.loc["Control W22"][col4]
						count,ISO,ISO_mem,ISO_drift,ISO_smow = ISOmean.count(),ISOmean.mean(),ISOmem_corr.mean(),ISOdrift_corr.mean(),ISOvsmow_corr.mean()
						sd,sd_mem,sd_drift,sd_smow = ISOmean.std(),ISOmem_corr.std(),ISOdrift_corr.std(),ISOvsmow_corr.std()
						val.append((i,"_Control W22",ISO,sd,ISO_mem,sd_mem,ISO_drift,sd_drift,ISO_smow,sd_smow,count))
						if sd_smow >= limit:
							self.log.append("Warning: high standard deviation on sample {}".format(i))

					
					else:
						dat = df.where(df["Error Code"] == 0)
						dat = dat.loc[i]
						ISOmean = dat[col1]
						ISOmem_corr = dat[col2]
						ISOdrift_corr = dat[col3]
						ISOvsmow_corr = dat[col4]
						count,ISO,ISO_mem,ISO_drift,ISO_smow = ISOmean.count(),ISOmean.mean(),ISOmem_corr.mean(),ISOdrift_corr.mean(),ISOvsmow_corr.mean()
						sd,sd_mem,sd_drift,sd_smow = ISOmean.std(),ISOmem_corr.std(),ISOdrift_corr.std(),ISOvsmow_corr.std()
						val.append((i,dat.index.values[0][0],ISO,sd,ISO_mem,sd_mem,ISO_drift,sd_drift,ISO_smow,sd_smow,count))
						if sd_smow >= limit:
							self.log.append("Warning: high standard deviation on sample {}".format(i))

		def runOverview():

			self.run_overview = pd.DataFrame(val,columns=["Identifier 1",
					"Identifier 2",
					"{}_raw".format(general_label),
					"stdev. raw",
					"{} memory".format(general_label),
					"stdev. memory",
					"{} drift".format(general_label),
					"stdev. drift",
					"{} vsmow".format(general_label),
					"{} stdev. vsmow".format(general_label),
					"{} counts".format(general_label)]) 

		runOverview()


	def getFinalValues(self):

		return self.run_overview.iloc[:,[0,1,8,9,10]]

	def checkStandards(self,isotope="O"):

		df = self.getFinalValues()
		check = df.where(df["Identifier 2"]== "_Standard").dropna()

		if isotope == "O":
			self.constraint = 0.1
			standards = self.defined_O.copy()
		else:
			self.constraint = 0.7
			standards = self.defined_H.copy()
		diffs = []
		for i,j in zip(check.iloc[:,2],standards):
			diff = abs(i-standards[j])
			diffs.append(diff)

			string = "Measured:{:4.2f} --- Standard: {:4.2f} --- Difference: {:+4.2f} ".format(i,standards[j],diff)
			self.log.append(string)
			print(string)

		diffs = np.array(diffs)
		if np.any(diffs >self.constraint)==True:
			self.log.append("Warning! At least one of the standards has standard deviation superior to {}".format(self.constraint))
		else:
			self.log.append("Standard deviations of in-house standards seem to be fine")




def Run(iso,filename):

	RUN = Isotope(filename)
	RUN.readRaw()
	RUN.checkEmpty()
	RUN.checkVolume()
	RUN.runSummary()
	RUN.IsotopeSelect(iso)
	RUN.initMemCoeffs()
	RUN.Optimize(iso,method = 'default')
	RUN.MemoryCorrection(iso)
	RUN.driftCorrect(iso)
	RUN.VSMOWcorrect(iso)
	RUN.getMeanSDs(iso)
	RUN.checkStandards(iso)

	return RUN

def FullRun(filename):

	print("Running the corrections for Oxygen \n ... \n ...")
	X = Run("O",filename)
	print('Done! \n ... \n ...\n ... \n ...')

	print("Running the corrections for Deuterium \n ... \n ...")
	Y =Run("H",filename)
	print('Done!')

	return X,Y

def Merge(IsoO,IsoH):

	dfO = IsoO.getFinalValues()
	dfH = IsoH.getFinalValues()

	df = pd.merge(dfO,dfH, on = ["Identifier 1","Identifier 2"])

	return df

def OverviewPlot(IsoO,IsoH):

	merged = Merge(IsoO,IsoH)[0:19]

	fig,ax = plt.subplots()
	xi = merged["d18O vsmow"]
	yi = merged["d2H vsmow"]
	xi_err = merged["d18O stdev. vsmow"]
	yi_err = merged["d2H stdev. vsmow"]
	ax.plot(xi,yi,'o',label="results")
	ax.errorbar(xi,yi,yerr = yi_err  ,xerr= xi_err,
			marker = '.',
			markersize = 0,
			lw = 0,
			elinewidth = 1,
			ecolor = 'black',label="error")

	x = np.arange(ax.get_xlim()[0]-1,ax.get_xlim()[1]+1,0.1)
	ax.set_xlim(ax.get_xlim()[0],ax.get_xlim()[1])
	ax.set(xlabel="$\delta^{18}$O [‰] (VSMOW)",ylabel="$\delta^{2}$H [‰] (VSMOW)")
	ax.plot(x,x*8+10,label='GMWL')
	ax.fill_between(x,x*8+2,x*8+18,alpha = 0.5,color = "orange")

	ax.legend()
	ax.grid()

	ax.tick_params(direction='in',top=True,right=True)
	plt.title("Dual Isotope space plot of Results")
	plt.tight_layout()
	plt.show()

def DatatoCSV(IsoO,IsoH):
	merged = Merge(IsoO,IsoH)[0:19]
	merged.to_csv('ipynb/data.csv')



