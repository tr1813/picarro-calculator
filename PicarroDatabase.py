#!/usr/bin/python
# -*- coding: latin-1 -*-

import sys; sys.path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import Picarro as pica
import sqlite3
import glob
from sqlite3 import Error
import pickle




def listFiles(dir):
	return glob.glob("{}/*.csv".format(dir))

def CreateConnection(db_file):
	""" create a database connection to the SQLite database
		specified by db_file
	:param db_file: database file
	:return: Connection object or None
	https://www.sqlitetutorial.net/sqlite-python/create-tables/
	"""
	conn = None
	try:
		conn = sqlite3.connect(db_file)
		return conn
	except Error as e:
		print(e)

	return conn

def CreateTable(conn,statement):
	""" create a table from the create_table_sql statement
	:param conn: Connection object
	:param create_table_sql: a CREATE TABLE statement
	:return:

	https://www.sqlitetutorial.net/sqlite-python/create-tables/
	"""
	try:
		c = conn.cursor()

		c.execute(statement)
	except Error as e:
		print(e)

def AddSummaryRun(run,conn):

	statement = """CREATE TABLE runs
		('key' int,
		'Identifier 1' varchar(255),
		'Identifier 2' varchar(255),
		'RUN_ID' int,
		'position' int,
		'd18O vsmow' float(2),
		'd18O stdev. vsmow' float(2),
		'd18O counts' int,
		'd2H vsmow' float(2),
		'd2H stdev. vsmow' float(2),
		'd2H counts' int,
		'inside GMWL' varchar(255) ,
		PRIMARY KEY ('key')
		);
		"""
	CreateTable(conn,statement)


	c = conn.cursor()
	def getrows(df):
		rows = []

		for i,j in zip(df.index,df.values):
			rows.append(tuple([i]+list(j)))
		return rows

	rows = getrows(run.merge)
	sql ="""INSERT INTO runs
			('key',
			'Identifier 1',
			'Identifier 2',
			'RUN_ID',
			'position',
			'd18O vsmow',
			'd18O stdev. vsmow',
			'd18O counts',
			'd2H vsmow',
			'd2H stdev. vsmow',
			'd2H counts',
			'inside GMWL')
			VALUES (?,?,?,?,
			?,?,?,?,
			?,?,?,?);"""

	for row in rows:

		try:
			c.execute(sql,tuple(row))
		except Error as e:
			print(e)
	conn.commit()

def AddRun(filename,nickname,conn):

	statement1="""CREATE TABLE runlookup ('NickName' varchar(255),
		'RUN_ID' int,
		PRIMARY KEY ('RUN_ID')
		);"""

	CreateTable(conn,statement1)
	run = pd.read_csv(filename)

	run_id = float(filename[-19:-11])

	run["RUN_ID"] = run_id * np.ones(len(run))
	print(run_id)

	sql = """ INSERT INTO runlookup (
		'NickName',
		'RUN_ID')
		VALUES (?,?);"""

	c = conn.cursor()

	try:
		c.execute(sql,(nickname,run_id))
	except Error as e:
		print(e)

	conn.commit()



def AddRaw(filename,conn):

	statement=""" CREATE TABLE rawrun (
		'Line' int,
		'Analysis' varchar,
		'Time Code' timestamp,
		'Port' varchar,
		'Inj Nr' int,
		'd(18_16)Mean' float(5),
		'd(D_H)Mean' float(5),
		'H2O_Mean' float(5),
		'Ignore' int,
		'Good' int,
		'Identifier 1' varchar,
		'Identifier 2' varchar,
		'Gas Configuration',
		'Timestamp Mean' float(1),
		'd(18_16)_SD' float(5),
		'd(D_H)_SD' float(5),
		'H2O_SD' float(5),
		'd(18_16)_Sl' float(5),
		'd(D_H)_Sl' float(5),
		'H2O_Sl' float(5),
		'baseline_shift' float(5),
		'slope_shift' float(5),
		'residuals' float(5),
		'baseline_curvature' float(5),
		'interval' float(5),
		'ch4_ppm' float(5),
		'h16od_adjust' float(5),
		'h16od_shift' float(5),
		'n2_flag' int,
		'DAS Temp' float(5),
		'Tray' int,
		'Sample' int,
		'Job' int,
		'Method' varchar,
		'Error Code' int,
		'RUN_ID' int,
		PRIMARY KEY ('Timestamp Mean')
		);"""

	RAW = pd.read_csv(filename)

	run_id = float(filename[-19:-11])

	RAW["RUN_ID"] = run_id * np.ones(len(RAW))

	CreateTable(conn,statement)

	sql=""" INSERT INTO rawrun(
		'Line',
		'Analysis',
		'Time Code',
		'Port',
		'Inj Nr',
		'd(18_16)Mean',
		'd(D_H)Mean',
		'H2O_Mean',
		'Ignore',
		'Good',
		'Identifier 1',
		'Identifier 2',
		'Gas Configuration',
		'Timestamp Mean',
		'd(18_16)_SD',
		'd(D_H)_SD',
		'H2O_SD',
		'd(18_16)_Sl',
		'd(D_H)_Sl',
		'H2O_Sl',
		'baseline_shift',
		'slope_shift',
		'residuals',
		'baseline_curvature',
		'interval',
		'ch4_ppm',
		'h16od_adjust',
		'h16od_shift',
		'n2_flag',
		'DAS Temp',
		'Tray',
		'Sample',
		'Job',
		'Method',
		'Error Code','RUN_ID') VALUES(
		?,?,?,?,?,
		?,?,?,?,?,
		?,?,?,?,?,
		?,?,?,?,?,
		?,?,?,?,?,
		?,?,?,?,?,
		?,?,?,?,?,?
		);"""

	c = conn.cursor()

	for row in RAW.values:
		try:
			c.execute(sql,tuple(row))
		except Error as e:
			print(e)
	conn.commit()


def checkforrawdata(path_to_watch):
	newfile=list()
	try:
		with open(os.path.join(path_to_watch,'filelist.txt'),'rb') as fp:
			ls_old=pickle.load(fp)
			ls_new=list(f for f in glob.glob(os.path.join(path_to_watch,'*.csv')))
		if ls_old == ls_new:
			print('There is no new raw data in directory:',path_to_watch)

		else:
			added= [f for f in ls_new if not f in ls_old]
			removed= [f for f in ls_old if not f in ls_new]
			if removed:
				ls_new=ls_new+removed
				with open(os.path.join(path_to_watch,'filelist.txt'),'wb') as fp:
					pickle.dump(ls_new,fp)
				print('missing file:',removed)
			if added:
				newfile=list(i for i in added)
				with open(os.path.join(path_to_watch,'filelist.txt'),'wb') as fp:
					pickle.dump(ls_new,fp)
				print('new file',added)
				return newfile
	except EOFError:
		l=list()
		with open(os.path.join(path_to_watch,'filelist.txt'),'wb') as fp:
			pickle.dump(l,fp)
		print('Could not find FileList...created empty FileList in directory:',path_to_watch)

	except FileNotFoundError:
		l=list()
		with open(os.path.join(path_to_watch,'filelist.txt'),'wb') as fp:
			pickle.dump(l,fp)
		newfile=list(f for f in glob.glob(os.path.join(path_to_watch,'*.csv')))
		with open(os.path.join(path_to_watch,'filelist.txt'),'wb') as fp:
			pickle.dump(newfile,fp)
		print('Could not find FileList...created empty FileList in directory:',path_to_watch)
		print('Added all existing files in directory',path_to_watch)
		print('Existing files:',newfile)
		return newfile


def ReplaceName(conn,RUN_ID,newname):
    """ create a table from the create_table_sql statement
    :param conn: Connection object
    :param RUN_ID: the run id, an eight digit integer with format yyyymmdd
	:param newname: the new nick name
    """

    statement= """UPDATE runlookup
    SET RUN_ID = {0}, NickName = '{1}'
    WHERE RUN_ID = {0};""".format(RUN_ID,newname)

    try:
        c = conn.cursor()

        c.execute(statement)
    except Error as e:
            print(e)
    conn.commit()


def FullRunUpdate(filename):

	PATH = r"J:\c715\Picarro\Results\Database\data.db"

	conn_local = CreateConnection(PATH)

	# Add raw files to database
	AddRaw(filename,conn_local)

	# Add corrected values to database, automatically does the optimization procedure
	frun = pica.FullRun(filename)
	AddSummaryRun(frun,conn_local)

	# Add Nickname to Run look up table
	AddRun(filename,frun.nickname,conn_local)

	return frun




### In case the database is rebuilt, and the run look up table needs to be done again, uncomment the next bit of code.


#strips = ["~$","J:\\c715\\Picarro\\Results\\Results 2019\\Picarro-Results-",".xlsx","-","_"]

#nicknames = []

#for i in glob.glob(r"J:\c715\Picarro\Results\Results 2019\*.xlsx"):
#    full = i
#    for j in strips:

#        full = full.replace(j," ")
#    #print(full)
#    if full.startswith(" 2019") == True:
#        nicknames.append(full)

#RUN_IDS = glob.glob(r"J:\c715\Picarro\Results\Results 2019\Raw data\*.csv")

#for i,j in zip(nicknames,RUN_IDS):
#    pdb.AddRun(j,i,conn)
