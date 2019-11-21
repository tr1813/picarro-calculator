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

def AddSummaryRun(filename,conn):

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
		'inside GMWL' varchar(255),
		PRIMARY KEY ('key')
		);
		"""
	CreateTable(conn,statement)

	RUN = pica.FullRun(filename)

	c = conn.cursor()
	def getrows(df):
		rows = []
		
		for i,j in zip(df.index,df.values):
			rows.append(tuple([i]+list(j)))
		return rows
		
	rows = getrows(RUN.merge)
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

def AddRun(filename,conn):

	statement1="""CREATE TABLE runlookup ('Analysis start' timestamp,
		'RUN_ID' int,
		PRIMARY KEY ('RUN_ID')
		);"""

	CreateTable(conn,statement1)
	run = pd.read_csv(filename)

	run_id = float(filename[-19:-11]+filename[-10:-4])

	run["RUN_ID"] = run_id * np.ones(len(run))
	print(run_id)

	#sql = """ INSERT INTO runlookup (
	#	'Analysis start',
	#	'RUN_ID')
	#	VALUES (?,?);"""

	#c = conn.cursor()

	#try:
	#	c.execute(sql,tuple(row))
	#except Error as e:
	#	print(e)
	
	#conn.commit()



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