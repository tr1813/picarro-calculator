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

def CreateTable(conn):
	""" create a table from the create_table_sql statement
	:param conn: Connection object
	:param create_table_sql: a CREATE TABLE statement
	:return:

	https://www.sqlitetutorial.net/sqlite-python/create-tables/
	"""
	try:
		c = conn.cursor()
		statement = """CREATE TABLE runs 
			('key' int,
			'Identifier 1' varchar(255),
			'Identifier 2' varchar(255),
			'd18O vsmow' float(2),
			'd18O stdev. vsmow' float(2),
			'd18O counts' int,
			'd2H vsmow' float(2),
			'd2H stdev. vsmow' float(2),
			'd2H counts' int,
			PRIMARY KEY ('key')
			);
			"""
		c.execute(statement)
	except Error as e:
		print(e)

def AddSummaryRun(filename,conn):

	RUN = pica.FullRun(filename)

	c = conn.cursor()

	def ReplaceLine():
		statements = []
		list1 = [i for i in RUN.merge.index]
		list2 = [j for j in RUN.merge]
		list3 = [k for k in RUN.merge.values]
		for i,k in zip(list1,list3):
			statements.append("""INSERT INTO runs
				('key',
				'Identifier 1',
				'Identifier 2',
				'd18O vsmow',
				'd18O stdev. vsmow',
				'd18O counts',
				'd2H vsmow',
				'd2H stdev. vsmow',
				'd2H counts') 
				VALUES ({0},'{1}','{2}',{3},{4},{5},{6},{7},{8});""".format(i,
					k[0],k[1],k[2],k[3],k[4],k[5],k[6],k[7]))
			
		return statements

	statements = ReplaceLine()

	for i in statements[:-1]:
		try:
			c.execute(i)
		except Error as e:
			print(e)
	conn.commit()