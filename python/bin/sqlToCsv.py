#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Author: FCST Development Team
# License: TBD
# Copyright (c) 2014, MIT License
"""
Python script to extract data from the fcst database
"""
import sqlite3 as lite
import sys
import csv


def findEntry(cur, requirements, modelName, tolerance):
    """Function for finding model entry given a tolerance.
    Gets model meta data from the head, looks at entries and
    compares them againts your requirments.
    
    cur - database cursor 
    """
    sqlCmd = "SELECT * FROM HEAD WHERE MODEL_NAME = '" + modelName + "'"
    cur.execute(sqlCmd)
    rows = cur.fetchall()
    
    for row in rows:
        reqMet = False
        #we take a row, does this row meet requirements
        for req in requirements.keys():
            #For each requirement try find in row
            for i in range(len(row)-1):
                if row[i] == req:
                    diff = abs((float(row[i+1]) - float(requirements[req]))/float(requirements[req]))
                    if diff < tolerance:
                        reqMet = True 
                    else:
                        reqMet = False
                    break
                    
            #If we didn't find a rew in that row then go to next row
            if not reqMet:
                break
        
        if reqMet:
            #We found it in that row, get name
            return row[-1]
                
    return "NONE"

def outputTable(model, requirements, tolerance, dbAddr, dumpLocation):
    con = None
    
    try:
        con = lite.connect(dbAddr)
        cur = con.cursor()     
        table = findEntry(cur, requirements, model, tolerance)
        
        
        
        if table == "NONE":
            print "Could not find entry"
        else:
            print "Outputting table for " + table
            
            sqlCmd='select * from ' + table + ' order by x'
            cur.execute(sqlCmd)
            rows = cur.fetchall()
            
            dumpFileAddr = dumpLocation + table + ".csv"
            writer = csv.writer(open(dumpFileAddr, 'w'))
            for row in rows:
                writer.writerow(row)               
        
            print "Done!"
        
        
    except lite.Error, e:
        
        print "Error %s:" % e.args[0]
        sys.exit(1)
        
    finally:
        
        if con:
            con.close()


if __name__ == "__main__":
    
    dumpLocation="."
    dbAddr='./AFCCValidation04_100.00'
    reqs = {"x_O2": 0.1, "phi_m": -0.1, "phi_s": 0.625000}
    #modelName = "numericalionomer"
    modelName = "wateragglomeratenumerical"
    selectionTolerance = 0.11
    
    
    outputTable(modelName,reqs, selectionTolerance, dbAddr, dumpLocation)
    
    