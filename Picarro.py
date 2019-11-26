#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys; sys.path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from math import exp,log,sqrt
from scipy.optimize import minimize,differential_evolution
import os
import csv
import datetime
import time


class Isotope(object):
	"""docstring for Isotope
	A object of data"""
	def __init__(self):

		self.filename = None
		self.dir = None
		self.log = []
		self.raw = None
		self.noempty = None
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
		self.run_id = None

	def Log(self):
		for i in self.log:
			print(i+"\n")





	def readRaw(self,filename):
		self.filename = filename
		df = pd.read_csv(self.filename)

		renamed = dict([(i,i.strip()) for i in df]) # strip white space from column names
		df.rename(columns=renamed, inplace=True)	# and rename all headers of the dataFrame with new 'clean names'
		df["Identifier 1"] = [i.strip() for i in df["Identifier 1"].values]
		if len(df['Identifier 1'].values[0]) < 1:
			raise Exception('Sample description missing, load it manually')
		df["Identifier 2"] = [i.strip() for i in df["Identifier 2"].values]
		df["Line"] = [int(i) for i in df["Line"]]
		self.run_id = float(self.filename[-19:-11])

		df["RUN_ID"] = self.run_id * np.ones(len(df))
		df.set_index(["Identifier 1","Identifier 2","Inj Nr"], inplace=True)

		self.raw = df
	def DummyCheck(self):

		if len(self.raw)>132:
			raise Exception("Your file has 133 or more lines. I give up now...")


	def checkEmpty(self):

		df = self.raw.copy()
		df = df.where(df["H2O_Mean"]!="              ").dropna()
		df["H2O_Mean"] = pd.to_numeric(df["H2O_Mean"])
		for i in df:
			try:
				df[i] = pd.to_numeric(df[i])
			except:
				print("Cannot convert column {} to numeric type".format(i))

		self.noempty = df

		diff = len(self.raw) -len(self.noempty)
		if diff != 0:  # checks a boolean array which evaluates to True if
															# if one cell is empty

			self.log.append("Warning: there are some empty lines! {} to be precise".format(diff))
		else:
			self.log.append("No empty cells. Proceeding...")

	def checkVolume(self):

		conds = [self.noempty.H2O_Mean.values < 17000 , self.noempty.H2O_Mean.values > 23000]	# create the conditions for ignoring a certain analysis line
		choices = [1, 2]													# what the error code value will become

		out = np.select(conds, choices, default=0)							# mapping the error code on the conditions

		self.noempty["Error Code"]=out												# updating the dataFrame

		above = self.noempty.where(self.noempty["Error Code"] == 2)
		below = self.noempty.where(self.noempty["Error Code"] == 1)
		warning = pd.concat([above,below]).dropna()

		lines =[i for i in warning["Line"].values]
		for i in lines:
			self.log.append("Warning: H20 value outside bounds on line {:.0f}...".format(i))


	def setPrimaryKey(self):
		#def f(timecode):
			#t = datetime.datetime.strptime(timecode,'   %Y/%m/%d %H:%M:%S')
			#return int(time.mktime(t.timetuple()))

		#pk = [int(f(i)) for i in self.raw["Time Code"].values]
		pk = [int(str(i).strip(" ")[2:]) for i in self.noempty["Analysis"].values]

		self.noempty["key"] = pk



	def runSummary(self):

		self.summary = self.noempty[["Line",
						"H2O_Mean",
						"Good",
						"Error Code",
						"d(18_16)Mean",
						"d(D_H)Mean"]]

		self.summary = self.summary.sort_index(axis=0)

	def plotSummary(self):

		df = self.summary
		fig,(ax1,ax2,ax3) = plt.subplots(3,1,figsize= (6,8),sharex= True)

		ax1.plot(df["Line"],df["H2O_Mean"],'o',color = 'black',markerfacecolor='w', label="Injection Volumes")
		ax1.legend(loc='best')
		ax1.set_ylim(15000,30000)
		ax1.set_ylabel("H2O in ppmV")
		ax1.axhline(20000,color = 'red')
		ax1.fill_between(np.arange(-2,133,1),17000,23000,color='red',alpha=0.5)
		ax1.set_xlim(0,132)

		ax2.plot(df["Line"],df["d(18_16)Mean"],'o',color = 'black',markerfacecolor='w',label= "Raw values")
		ax2.legend(loc='best')
		ax2.set_ylim(-50,10)
		ax2.set_ylabel("$\delta^{18}O$ permil [VSMOW]")

		ax3.plot(df["Line"],df["d(D_H)Mean"],'o',color = 'black',markerfacecolor='w',label= "Raw values")
		ax3.legend(loc='best')
		ax3.set_ylim(-300,50)
		ax3.set_ylabel("$\delta^{18}O$ permil [VSMOW]")
		ax3.set_xlabel("Analysis line")


		for ax in (ax1,ax2,ax3):
			ax.grid()
			ax.tick_params(direction='in',top=True,right=True)

		plt.tight_layout()
		plt.show()

	def IsotopeSelect(self,isotope="O"):



		if isotope == "O":
			self.corr =  self.noempty[["key","Line","d(18_16)Mean","Ignore","Error Code","RUN_ID"]]
			self.corr = self.corr.where(self.corr["Line"]>6).dropna()
			self.constraint = 0.1
		else:
			if isotope == "H":
				self.corr =  self.noempty[["key","Line","d(D_H)Mean","Ignore","Error Code","RUN_ID"]]
				self.corr = self.corr.where(self.corr["Line"]>6).dropna()
				self.constraint = 0.8
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
		#self.dir = "../../Picarro_data/PICARRO_run_"+self.run_id

		self.log.append("Checking a coeffs whether a coeffs file already exists")

		#if os.path.isdir(self.dir)==False:
			#print('creating a directory to store the data')
			#os.makedirs(self.dir)
			#if os.path.isfile(self.dir+"/coeffs{}_{}.csv".format(isotope,self.filename[31:46]))==True:
				#self.log.append("coeffs file already exists!")
		#else:
			#self.log.append("coeffs file does not exist yet. Switch to Optimization...")



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
				stds.append(np.array(df.loc[(loc1,loc2),col].values))

			std1_0 = df.loc[("TAP","Conditioning"),col].values[-1:]
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

			SD = (self.HAUS1_raw_sd**(0.5)+self.HAUS2_raw_sd**(0.5)+self.TAP_raw_sd**(0.5))**(0.5)

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

			with open(self.dir+'/coeffs{}_{}.csv'.format(iso,self.filename[31:46]), 'w') as coeffFile:
				writer = csv.writer(coeffFile)


				writer.writerows(lines)
			coeffFile.close()

		#writeCoeffsfile(self.coeffs)


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

		def makecoefflist():
			df = self.corr["Line"]
			coeffs = [self.coeffs[i] for i in self.coeffs]

			conditioning = df[0:4].values
			coeff_conditioning = [1,1,1,1]

			standards = df[4:34].values-1-(df[4:34].values-1) %10

			samples = df[34:].values -1-(df[34:].values-1) %4
			coeff_standards = coeffs * 3
			previous = list(conditioning) + list(standards) + list(samples)
			coeff_samples = []
			index_diff = (df.values[34:]-1)%4
			for i in index_diff:
			    coeff_samples.append(coeffs[int(i)])
			new_list =  coeff_conditioning + coeff_standards + coeff_samples
			return previous,new_list

		previouslines,newlist  = makecoefflist()

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
		self.memory = self.corr[["key","Line",col2,col1,"Error Code","RUN_ID"]]

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
			posxy=((7,-3.5,"HAUS1"),(7,-31.5,"HAUS2"),(7,-18,"TAP"),(3.2,-7,"CONTROL"))
		else:
			col1 = "d(D_H)Mean"
			col2 = "d(D_H)mem_corrected"
			posxy=((7,0,"HAUS1"),(7,-200,"HAUS2"),(7,-110,"TAP"),(3.2,-25,"CONTROL"))

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

		y = df.iloc[34:38][col1] # for the control run (either W22 or HUSN)
		x =np.arange(1,5)
		y_corr = df.iloc[34:38][col2]
		ax4.plot(x,y,'o',color = 'black',markerfacecolor='w')
		ax4.plot(x,y_corr,'o',color = 'black',markersize = 3, markerfacecolor='black')
		self.OLSR(x,y_corr,(0,4),ax4)[0]
		posx = ax4.get_xlim()[1]+0.5
		posy = ax4.get_ylim()[0]
		ax4.text(posx,posy,df.iloc[34:38].index[0][0])
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
			xy1 = df.loc['TAP',['Line',name]][0:4]
			xy2 = df.loc['TAP',['Line',name]][10:]

			xy = xy1.append(xy2)

			slope, intercept, r_value, p_value, std_err = stats.linregress(xy["Line"],xy[name])

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
			posxy = (60,-16.5)
		else:
			col = 'd(D_H)mem_corrected'
			col2 = 'd(D_H)drift_corrected'
			posxy = (60,-93)

		xy1 = self.drift.loc['TAP',['Line',col,col2]][0:4]
		xy2 = self.drift.loc['TAP',['Line',col,col2]][10:]
		xy = xy1.append(xy2)


		fig,ax = plt.subplots()

		self.OLSR(xy["Line"],xy[col],(0,132),ax)[0]

		params = self.OLSR(xy["Line"],xy[col],(0,132),ax)[1]

		ax.plot(xy["Line"],xy[col],
				'o',
				color = 'black',
				markersize = 5,
				markerfacecolor='white',label =" no drift correction")

		ax.plot(xy["Line"],xy[col2],
				'o',
				color = 'black',
				markersize = 5,
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
			limx = np.arange(-40,10,1)
		else:
			posxy = (-60,-200)
			limx = np.arange(-250,50,1)

		params = self.VSMOW_params

		fig, ax = plt.subplots(figsize= (6,4))
		xi = limx

		ax.plot(params["xvals"],params["yvals"],'o',color = 'black',markerfacecolor='w')
		ax.plot(xi,xi*params["slope"]+params["intercept"], color = 'black')

		ax.text(posxy[0],posxy[1],"y {:.3f} = x {:+.3f}".format(params["slope"],params["intercept"]))
		ax.grid()
		ax.tick_params(direction='in',top=True,right=True)
		plt.tight_layout()
		plt.show()

	def getMeanSDs(self,isotope = "O"):

		df = self.vsmow.copy()
		self.run_id = self.vsmow["RUN_ID"].values[0]

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
			return set([i[0] for i in df.index.values])

		def runCheck2(df1,col4):
			df = df1.copy()

			length = len(df)

			df = df[col4]

			stds = []
			stds.append((df.std(),"no missing measurements"))

			for i in range(len(df)):
				statement1 = "missing measurement {}".format((i)%length)
				stds.append((df.iloc[[i,(i+1)%length,(i+2)%length]].std(),statement1))
				statement2 ="missing measurements {} and {}".format((i)%length,(i+1)%length)
				stds.append((df.iloc[[i,(i+1)%length]].std(),statement2))
				if i <2:
					statement3 ="missing measurements {} and {}".format((i-1)%length+1,(i+1)%length+1)
					stds.append((df.iloc[[i,(i+2)%length]].std(),statement3))

			def checks(mylist):

				if mylist[0][0]<self.constraint:
					self.log.append("Standard dev is good")
					print("Standard dev is good")
				else:
					conds = [mylist[1][0],mylist[4][0],mylist[7][0],mylist[9][0]]
					if min(conds) < self.constraint:
						val = ((conds.index(min(conds))-1)%4)+1
						self.log.append("Standard dev too high get rid of measurement {}".format(val))
						print("Standard dev too high get rid of measurement {}".format(val))
						df1.drop(index = val ,level = 1, inplace = True)

					else:
						conds2 = [(mylist[2][0],2),(mylist[3][0],3),(mylist[5][0],5),(mylist[6][0],6),(mylist[8][0],8),(mylist[10][0],10)]
						conds_temp = [i[0] for i in conds2]
						if min(conds_temp) < self.constraint:



							value= conds_temp.index(min(conds_temp))
							new_val = conds2[value][1]
							print("get rid of measurements {}".format(stds[new_val][1][-8:]))
							self.log.append("get rid of measurements {}".format(stds[new_val][1][-8:]))
							vals_to_drop = [stds[new_val][1][-7],stds[new_val][1][-1]]
							df1.drop(index = (int(vals_to_drop[0])-1)%4 ,level = 1, inplace = True)
							df1.drop(index = (int(vals_to_drop[1])-1)%4 ,level = 1, inplace = True)
						else:
							print("too high std. deviation")
							self.log.append("too high std. deviation")

			if len(df) == 4:
				checks(stds)
			return df1


		for i in getSampleNames(self.vsmow):
			if (i == "HAUS1") or (i=="HAUS2"):
				dat = df.where(df["Error Code"] == 0)
				dat = dat.loc[[i,'Standard']]
				key = max(dat.iloc[0:10]["key"])
				ISOmean = dat.iloc[0:10][col1]
				ISOmem_corr = dat.iloc[0:10][col2]
				ISOdrift_corr = dat.iloc[0:10][col3]
				ISOvsmow_corr = dat.iloc[0:10][col4]
				count,ISO,ISO_mem,ISO_drift,ISO_smow = ISOmean.count(),ISOmean.mean(),ISOmem_corr.mean(),ISOdrift_corr.mean(),ISOvsmow_corr.mean()
				sd,sd_mem,sd_drift,sd_smow = ISOmean.std(),ISOmem_corr.std(),ISOdrift_corr.std(),ISOvsmow_corr.std()
				val.append((key,i,"_Standard",ISO,sd,ISO_mem,sd_mem,ISO_drift,sd_drift,ISO_smow,sd_smow,count,self.run_id))
				print(key)
				if sd_smow >= limit:
					self.log.append("Warning: high standard deviation on sample {}".format(i))


			else:
				if i == "TAP":
					dat = df.where(df["Error Code"] == 0)
					dat = dat.loc[i]
					key = max(dat.loc["Conditioning"]["key"])
					ISO = dat.loc["Conditioning"][col1]
					ISOmean = dat.loc["Conditioning"][col1]
					ISOmem_corr = dat.loc["Conditioning"][col2]
					ISOdrift_corr = dat.loc["Conditioning"][col3]
					ISOvsmow_corr = dat.loc["Conditioning"][col4]
					count,ISO,ISO_mem,ISO_drift,ISO_smow = ISOmean.count(),ISOmean.mean(),ISOmem_corr.mean(),ISOdrift_corr.mean(),ISOvsmow_corr.mean()
					sd,sd_mem,sd_drift,sd_smow = ISOmean.std(),ISOmem_corr.std(),ISOdrift_corr.std(),ISOvsmow_corr.std()
					val.append((key,i,"_Conditioning",ISO,sd,ISO_mem,sd_mem,ISO_drift,sd_drift,ISO_smow,sd_smow,count,self.run_id))
					print(key)
					if sd_smow >= limit:
						self.log.append("Warning: high standard deviation on sample {}".format(i))


					key = max(dat.loc["Standard"]["key"])
					ISO = dat.loc["Standard"][col1]
					ISOmean = dat.loc["Standard"][col1]
					ISOmem_corr = dat.loc["Standard"][col2]
					ISOdrift_corr = dat.loc["Standard"][col3]
					ISOvsmow_corr = dat.loc["Standard"][col4]
					count,ISO,ISO_mem,ISO_drift,ISO_smow = ISOmean.count(),ISOmean.mean(),ISOmem_corr.mean(),ISOdrift_corr.mean(),ISOvsmow_corr.mean()
					sd,sd_mem,sd_drift,sd_smow = ISOmean.std(),ISOmem_corr.std(),ISOdrift_corr.std(),ISOvsmow_corr.std()
					val.append((key,i,"_Standard",ISO,sd,ISO_mem,sd_mem,ISO_drift,sd_drift,ISO_smow,sd_smow,count,self.run_id))
					print(key)
					if sd_smow >= limit:
						self.log.append("Warning: high standard deviation on sample {}".format(i))



					js = np.arange(0,3,1)
					for j in js:


						ISO = dat.loc['Control'].iloc[j*4:(j+1)]

						ISOmean = dat.loc["Control"].iloc[j*4:(j+1)*4-1][col1]
						key = min(dat.loc["Control"].iloc[j*4:(j+1)*4-1]["key"])
						ISOmem_corr = dat.loc["Control"].iloc[j*4:(j+1)*4-1][col2]
						ISOdrift_corr = dat.loc["Control"].iloc[j*4:(j+1)*4-1][col3]
						ISOvsmow_corr = dat.loc["Control"].iloc[j*4:(j+1)*4-1][col4]
						count,ISO,ISO_mem,ISO_drift,ISO_smow = ISOmean.count(),ISOmean.mean(),ISOmem_corr.mean(),ISOdrift_corr.mean(),ISOvsmow_corr.mean()
						sd,sd_mem,sd_drift,sd_smow = ISOmean.std(),ISOmem_corr.std(),ISOdrift_corr.std(),ISOvsmow_corr.std()
						val.append((key,i,"_Control {}".format(j+1),ISO,sd,ISO_mem,sd_mem,ISO_drift,sd_drift,ISO_smow,sd_smow,count,self.run_id))
					if sd_smow >= limit:
						self.log.append("Warning: high standard deviation on sample {}".format(i))

				else:
					if i == 'W22' or i=='HUSN':
						dat = df.where(df["Error Code"] == 0)
						dat = dat.loc[i]
						try:
							key = dat.loc[[i,"Control W22"],"key"][0]
							ISOmean = dat.loc["Control W22"][col1]
							ISOmem_corr = dat.loc["Control W22"][col2]
							ISOdrift_corr = dat.loc["Control W22"][col3]
							ISOvsmow_corr = dat.loc["Control W22"][col4]
							count,ISO,ISO_mem,ISO_drift,ISO_smow = ISOmean.count(),ISOmean.mean(),ISOmem_corr.mean(),ISOdrift_corr.mean(),ISOvsmow_corr.mean()
							sd,sd_mem,sd_drift,sd_smow = ISOmean.std(),ISOmem_corr.std(),ISOdrift_corr.std(),ISOvsmow_corr.std()
							val.append((key,i,"_Control W22",ISO,sd,ISO_mem,sd_mem,ISO_drift,sd_drift,ISO_smow,sd_smow,count,self.run_id))
						except:
							key = dat.loc[[i,"Control HUSN"],"key"][0]
							ISOmean = dat.loc["Control HUSN"][col1]
							ISOmem_corr = dat.loc["Control HUSN"][col2]
							ISOdrift_corr = dat.loc["Control HUSN"][col3]
							ISOvsmow_corr = dat.loc["Control HUSN"][col4]
							count,ISO,ISO_mem,ISO_drift,ISO_smow = ISOmean.count(),ISOmean.mean(),ISOmem_corr.mean(),ISOdrift_corr.mean(),ISOvsmow_corr.mean()
							sd,sd_mem,sd_drift,sd_smow = ISOmean.std(),ISOmem_corr.std(),ISOdrift_corr.std(),ISOvsmow_corr.std()
							val.append((key,i,"_Control HUSN",ISO,sd,ISO_mem,sd_mem,ISO_drift,sd_drift,ISO_smow,sd_smow,count,self.run_id))
						if sd_smow >= limit:
							self.log.append("Warning: high standard deviation on sample {}".format(i))


					else:
						dat = df.where(df["Error Code"] == 0)
						print("Checking: {} ...".format(i))
						self.log.append("Checking: {} ...".format(i))

						key = dat.loc[[i],"key"].dropna()[0]
						dat = dat.loc[i]
						dat = runCheck2(dat,col4)

						ISOmean = dat[col1]
						ISOmem_corr = dat[col2]
						ISOdrift_corr = dat[col3]
						ISOvsmow_corr = dat[col4]
						count,ISO,ISO_mem,ISO_drift,ISO_smow = ISOmean.count(),ISOmean.mean(),ISOmem_corr.mean(),ISOdrift_corr.mean(),ISOvsmow_corr.mean()
						sd,sd_mem,sd_drift,sd_smow = ISOmean.std(),ISOmem_corr.std(),ISOdrift_corr.std(),ISOvsmow_corr.std()
						val.append((key,i,dat.index.values[0][0],ISO,sd,ISO_mem,sd_mem,ISO_drift,sd_drift,ISO_smow,sd_smow,count,self.run_id))
						print(key)
						if sd_smow >= limit:
							self.log.append("Warning: high standard deviation on sample {}".format(i))

		def runOverview():

			df = pd.DataFrame(val,columns=["key","Identifier 1",
					"Identifier 2",
					"{}_raw".format(general_label),
					"stdev. raw",
					"{} memory".format(general_label),
					"stdev. memory",
					"{} drift".format(general_label),
					"stdev. drift",
					"{} vsmow".format(general_label),
					"{} stdev. vsmow".format(general_label),
					"{} counts".format(general_label),
					"RUN_ID"])
			df.set_index('key', inplace = True)
			df.sort_index(inplace = True)

			df["position"] = df.index.values-(df.index.values[0]-1)
			self.run_overview = df




		runOverview()

	def getFinalValues(self):

		return self.run_overview.iloc[:,[0,1,11,12,8,9,10]]

	def checkStandards(self,isotope="O"):

		df = self.getFinalValues()
		check = df.iloc[-3:]

		if isotope == "O":
			self.constraint = 0.1
			standards = self.defined_O.copy()
		else:
			self.constraint = 0.8
			standards = self.defined_H.copy()
		diffs = []
		for i,j in zip(check.iloc[:,0],standards):
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

	def getDir(self):

		return self.dir

class Merged(object):
	"""docstring for Isotope
	A object of data"""
	def __init__(self):

		self.merge = None
		self.trimmed = None
		self.O18 = None
		self.D = None
		self.coeffs = None
		self.non_replicate = None
		self.non_triplicate = None
		self.non_perfect = None
		self.high_std = None
		self.gmwl = None
		self.nickname = None
		self.fullrerun = None

	def setMerge(self,O18,D):
		self.O18 = O18
		self.D = D
		self.merge = Merge(self.O18,self.D)

	def setCoeffs(self):
		d = {"O":[i for i in self.O18.coeffs.values()],"H":[j for j in self.D.coeffs.values()]}

		df = pd.DataFrame(d,index = np.arange(1,11,1))
		self.coeffs = df

	def TrimStandards(self):

		ls = ["HAUS1", "HAUS2", "TAP", "W22","HUSN"]
		df = self.merge.copy()
		for i in ls:
			df = df[df["Identifier 1"] !=i]

		self.trimmed = df

	def setNickName(self):
			self.nickname = str(input("""Please set the run's nickname:  \n
			Suggested format: YYYY DD MM UserInitials RunName RunNumber \n
			Example: 2000 01 01 CDG Paris Run01\n \n """))


	def suggestedReruns(self):

		def CheckGMWL(df):

			cond1 = df["d2H vsmow"] - df["d18O vsmow"]*8 - 18 < 0
			cond2 = df["d2H vsmow"] - df["d18O vsmow"]*8 - 2 > 0

			result = []
			for i,j in zip(cond1.values,cond2.values):
				k = i and j == True
				if k == True:
					g = "inside"
					result.append(g)
				else:
					g = 'outside'
					result.append(g)

			df["inside GMWL"] = result

			return df.where(df["inside GMWL"] == 'outside').dropna()



		def checkSTDv(df):
			df1 = df.where(df["d18O stdev. vsmow"]>0.1)
			df2 = df.where(df["d2H stdev. vsmow"]>0.8)

			df3 = df1.dropna().append(df2.dropna()).drop_duplicates()

			return df3


		def checkCount(df,thresh):

			df1 = df.where(df["d18O counts"]<=2)
			df2 = df.where(df["d2H counts"]<=2)
			df3 = df1.dropna().append(df2.dropna()).drop_duplicates()

			return df3

		self.non_triplicate = checkCount(self.trimmed,2)
		self.non_perfect = checkCount(self.trimmed,3)

		self.high_std = checkSTDv(self.merge)
		self.gmwl = CheckGMWL(self.merge)

		print("Checking for triplicates...")

		if len(self.non_triplicate) > 0:

			print("Some samples were not triplicated")
			print(self.non_triplicate["Identifier 1"])
			print('\n')

		else:

			print("This was a good run")
			print('\n')

		print("Checking for high standard deviations ...")
		if len(self.high_std)  > 0:
				print("Suggested reruns for following samples, which had high standard deviations")
				print(self.high_std["Identifier 1"])
				print('\n')

		else:
			print("This was a good run")
			print('\n')

		print("Checking for samples lying outside of the GWML ...")
		if len(self.gmwl)  > 0:
				print("Suggested reruns for following samples, which were outside of the GMWL")
				print(self.gmwl["Identifier 1"])
				print('\n')
		else:
			print("Nothing to report")
			print('\n')

	def setFullRerun(self):
		reruns = self.high_std.append(self.non_triplicate)
		reruns = reruns.append(self.gmwl)

		self.fullrerun = reruns[~reruns.index.duplicated(keep='last')]
		self.fullrerun.set_index(["Identifier 1", "Identifier 2"], inplace = True)

def RunFromDB(by=None,date=None,name=None,run_num=None,conn=None,iso="O"):

	run = Isotope()
	run.raw = getRawData(by=by,date=date,name=name,run=run_num,conn=conn)

	def DummyCheck2(df):
		"""Checks the length of the dataframe,
		If longer than 132 lines, returns an exception"""

		if len(df) >132:
			raise Exception("The dataframe contains more than 132 lines. Please narrow down the search. Quitting...")

	try:
		DummyCheck2(run.raw)
	except:
		raise

	# Makes all the other checks and produces a run summary
	for f in (run.checkEmpty,run.checkVolume,run.setPrimaryKey,run.runSummary):
		f()


	run.IsotopeSelect(iso)
	run.initMemCoeffs()
	run.Optimize(iso,method = 'default')
	# Finish by computing all the corrections
	for f in (run.MemoryCorrection,run.driftCorrect,run.VSMOWcorrect,run.getMeanSDs):
		f(iso)

	return run

def Run(iso,filename):

	try:
		run = Isotope()
		run.readRaw(filename)
		for f in (run.DummyCheck,run.checkEmpty,run.checkVolume,run.setPrimaryKey,run.runSummary):
			f()
	except:
		raise

	run.IsotopeSelect(iso)
	run.initMemCoeffs()
	run.Optimize(iso,method = 'default')

	for f in (run.MemoryCorrection,run.driftCorrect,run.VSMOWcorrect,run.getMeanSDs):
		f(iso)

	return run

def FullRun(filename = None,name=None,run_num=None,conn=None, date = None, by= None):

	run = Merged()

	if filename is not None:

		run.setNickName()

		print("Running the corrections for Oxygen \n ... \n ...")
		O = Run(iso= "O",filename=filename)
		print('Done! \n ... \n ...\n ... \n ...')

		print("Running the corrections for Deuterium \n ... \n ...")
		D =Run(iso="H",filename=filename)
		print('Done!')

		#OverviewDatatoCSV(X,Y)
		print("now merging the O and H isotope data \n ... \n ...")

	else:
		print("Running the corrections for Oxygen \n ... \n ...")
		O = RunFromDB(date = date,by=by,name=name,run_num=run_num,conn=conn,iso = "O")
		run.nickname = O.raw["NickName"].values[0]
		print('Done! \n ... \n ...\n ... \n ...')

		print("Running the corrections for Deuterium \n ... \n ...")
		D = RunFromDB(date = date,by=by,name=name,run_num=run_num,conn=conn,iso = "H")
		print('Done!')

		#OverviewDatatoCSV(X,Y)
		print("now merging the O and H isotope data \n ... \n ...")

	run.setMerge(O,D)

	for f in (run.setCoeffs,run.TrimStandards,run.suggestedReruns,run.setFullRerun):
		f()

	return run

def Merge(IsoO,IsoH):

	dfO = IsoO.getFinalValues()
	dfH = IsoH.getFinalValues()

	df = pd.merge(dfO,dfH, on = ["key","Identifier 1", "Identifier 2","RUN_ID","position"])

	return df

def OverviewPlot(frun):

	fig,ax = plt.subplots()

	# Plotting the corrected values for each measured sample of the run
	xi,yi = frun.trimmed["d18O vsmow"],frun.trimmed["d2H vsmow"]
	xi_err,yi_err  = frun.trimmed["d18O stdev. vsmow"],frun.trimmed["d2H stdev. vsmow"]
	ax.plot(xi,yi,'o',label="results")
	ax.errorbar(xi,yi,yerr = yi_err  ,xerr= xi_err,
			marker = '.',
			markersize = 0,
			lw = 0,
			elinewidth = 1,
			ecolor = 'black',label="error")


	ax.set_xlim(ax.get_xlim()[0],ax.get_xlim()[1])

	# Plotting the Global Meteoric Water line
	x = np.arange(ax.get_xlim()[0]-1,ax.get_xlim()[1]+1,0.1)
	ax.plot(x,x*8+10,label='GMWL')
	ax.fill_between(x,x*8+2,x*8+18,alpha = 0.5,color = "orange")

	ax.set(xlabel="$\delta^{18}$O [‰] (VSMOW)",ylabel="$\delta^{2}$H [‰] (VSMOW)")
	ax.tick_params(direction='in',top=True,right=True)
	ax.legend()
	ax.grid()

	plt.title("Dual Isotope space plot of Results \n {}".format(frun.nickname))
	plt.tight_layout()
	plt.show()

def OverviewDatatoCSV(IsoO,IsoH):
	dir_path = IsoO.getDir()
	merged = Merge(IsoO,IsoH)
	summary = IsoO.summary
	log = IsoO.log

	def writelog():
		f = open(dir_path+'/log.txt','w')

		for line in log:
			f.write(line+'\n')

		f.close()

	if os.path.isdir(dir_path) == False:
		print("Directory does not yet exist.")

		os.makedirs(dir_path)

		print("Writing file")

		merged.to_csv(dir_path+'/data.csv')

		summary.to_csv(dir_path+'/run_summary.csv')
		writelog()



	else:
		print("Directory already exists")
		if os.path.isfile(dir_path+'/data.csv') == True:
			print("Summary data file already exists!")
			if os.path.isfile(dir_path+'/run_summary.csv'):
				print("Run summary already exists")
			else:
				print("writing run summary file")
				summary.to_csv(dir_path+'/run_summary.csv')
				writelog()

			writelog()
		else:
			print("writing data file")
			merged.to_csv(dir_path+'/data.csv')
			if os.path.isfile(dir_path+'/run_summary.csv'):
				print("Run summary already exists")
			else:
				print("writing run summary file")
				summary.to_csv(dir_path+'/run_summary.csv')
				writelog()

def getRawData(conn,by,date=None,name=None,run=None):
	if by == 'keyword':
		statement = """SELECT r.'Analysis',r.'Line',r.'Identifier 1',r.'Identifier 2',r.'Inj Nr',
	    	r.'d(18_16)Mean',r.'d(D_H)Mean',r.'Good',r.'Ignore',r.'Error Code',r.'RUN_ID',r.'H2O_Mean',rlk.NickName
	    	from rawrun r, runlookup rlk
	    	WHERE (rlk.NickName like '%{0}%{1}%' and r.RUN_ID = rlk.RUN_ID);""".format(name,run)
		df = pd.read_sql_query(statement,conn)
	else:
		if by == 'date':
			statement = """SELECT r.'Analysis',r.'Line',r.'Identifier 1',r.'Identifier 2',r.'Inj Nr',
	    	r.'d(18_16)Mean',r.'d(D_H)Mean',r.'Good',r.'Ignore',r.'Error Code',r.'RUN_ID',r.'H2O_Mean'
	    	from rawrun r, runlookup rlk
	    	WHERE (r.RUN_ID = {} and r.RUN_ID = rlk.RUN_ID);""".format(date)
			df = pd.read_sql_query(statement,conn)
		else:
			raise Exception("by= must be either date (format: yyyymmdd) or keyword")

	df["Identifier 1"] = [i.strip() for i in df["Identifier 1"].values]
	df["Identifier 2"] = [i.strip() for i in df["Identifier 2"].values]
	df["Line"] = [int(i) for i in df["Line"]]
	df.set_index(["Identifier 1","Identifier 2","Inj Nr"], inplace = True)

	return df

def getCorrectedFromDB(conn,by,date=None,name=None,run=None):
	if by == 'keyword':
		statement = """SELECT r.'key',r.'Identifier 1',r.'Identifier 2',r.'RUN_ID',r.'position',
    	r.'d18O vsmow',r.'d18O stdev. vsmow',r.'d18O counts',
    	r.'d2H vsmow',r.'d2H stdev. vsmow',r.'d2H counts',rlk.NickName
    	from runs r, runlookup rlk
    	WHERE (rlk.NickName like '%{0}%{1}%' and r.RUN_ID = rlk.RUN_ID);""".format(name,run)
		df = pd.read_sql_query(statement,conn)
	else:
		if by == 'date':
			statement =  """SELECT r.'key',r.'Identifier 1',r.'Identifier 2',r.'RUN_ID',r.'position',
    		r.'d18O vsmow',r.'d18O stdev. vsmow',r.'d18O counts',
    		r.'d2H vsmow',r.'d2H stdev. vsmow',r.'d2H counts',rlk.NickName
    		from runs r, runlookup rlk
    		WHERE (r.RUN_ID = {} and r.RUN_ID = rlk.RUN_ID);""".format(date)
			df = pd.read_sql_query(statement,conn)
		else:
			raise Exception("by= must be either date (format: yyyymmdd) or keyword")

	df["Identifier 1"] = [i.strip() for i in df["Identifier 1"].values]
	df["Identifier 2"] = [i.strip() for i in df["Identifier 2"].values]
	df.set_index(["key"], inplace = True)
	ls = ["HAUS1", "HAUS2", "TAP", "W22","HUSN"]
	for i in ls:
		df = df[df["Identifier 1"] !=i]

	return df

def PlotRunFromDB(conn,by,date=None,name=None,run=None):
	X = Merged()
	X.trimmed = getCorrectedFromDB(conn,by,date=date,name=name,run=run)
	OverviewPlot(X)
