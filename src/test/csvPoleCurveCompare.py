#! /usr/bin/env python2.7
"""
Script to compare two polarization curves.

Author: Philip Wardlaw, Dec 2014

Multiple functionals can be compared per static point.
Left column of potentiostatic/galvanostatic points must match exactly.

Assumed data format:

    # This is a comment
    # All comment lines are ignored
    # First uncommented line is the meta data line
    Cell Potential<delimiter>Functional1  
    # All uncommented lines proceeding the meta line are data lines
    
    # Lines of pure white space are also ignored
    
    1<delimiter>2
    3<delimiter>4
    
    # Number of colums for each line should be consistent
    # Delimiter can be a tab or ,
    


"""
import csv
import sys


def ignore(row):
  "Returns true if line should be ignored"
  answer = False
  
  if len(row) == 0:
    "Line is just white space or has no delimiters"
    answer = True
  else:
     firstChars = row[0].strip()
     if len(firstChars) != 0:
         if firstChars[0] == '#':
            "Line's first non whitespace character is a '#'"
            answer = True

      
    
  return answer


def readRow(row, cnvrt2Float = False):
  "Read currect csv row, stripping white space and converting to float if necessary"
  data = []
  
  for el in row:
    el = el.strip()
    
    if cnvrt2Float:
        try:
            data.append(float(el))
        except ValueError as e:
            print "Error converting '" + el + "' to a floating point number"
            raise ValueError
      
    else:
        data.append(el)
      
      
  return data

def readData(fileName, delim = '\t'):
    "Read data from file, returns metaData and data"
    metaData = []
    data     = []
   
    reader = csv.reader(open(fileName), delimiter = delim)
    
    readMetaData = False
    lineNo = 1
    for row in reader: 
        #Check to see if row should be ignored (if it is a comment)
        if not ignore(row):
            if readMetaData is False:
                rowData = readRow(row)
                metaData =  rowData
                readMetaData = True
            else:
                try:
                    rowData = readRow(row, True)
                    data.append(rowData)
                except ValueError:
                    print "Error reading data row ", row , "at line number", lineNo, "of file", fileName
                    raise ValueError
            
        lineNo +=1
	   
        
    return metaData,data


"""
 Main function for the Python script:
 
"""
if __name__ == "__main__":
    
    
    tolerance = 0.01 #  1%
        

    assert len(sys.argv) ==3, "Invalid number of arguments"
    
    #Read the data:
    meta0, data0 = readData(sys.argv[1])
    meta1, data1 = readData(sys.argv[2])
    
    #Start tests:
    # First compare the metadata list (compares each element in the meta data list)
    # So MetaData should be identical:    
    if meta0 != meta1:
        print "Metadata miss match"
        exit(1)
    
    # Then, make sure there are the same number of rows of data:
    if len(data0) != len(data1):
        print "Data length miss match"
        exit(1)
    
    size_ = len(data0)
    
    # Now, start comparing every element:
    for i in range(size_):       
        
        if data0[i][0] != data1[i][0]:
            print "Potentiostatic/galvanostatic values do not match - they must match exactly!"
            exit(1)
        
        if len(data0[i]) != len(meta0):
            print "Number of functionals/columns in data row " + str(i) + " inconsistent with meta data"
            exit(1)
        
        if len(data0[i]) != len(data1[i]):
            print "Number of functionals reported per potentiostatic/galvanostatic point do not match"
            exit(1)
         
        for j in range(1, len(data0[i])):
            
            diff = abs((data1[i][j]-data0[i][j])/data1[i][j])

            if diff > tolerance:
                print "Functional values do not match within tolerance"
                exit(1)
                
    print "Files match within tolerance"          
    exit(0)