'''
Created on Jun 10, 2013

@author: wardlawp

Script for removing "no_name" entries from the "Output Variables" and "System management"
sections of the default xml parameter file.


Use:
python ModifyXML.py <path to defaul.xml>



Expected XML format:

<ParameterHandler>
    ...
    <System_20management>
        <Solution_20variables>
            ...
        </Solution_20variables>
        <Equations>
            ...
        </Equations>
    <System_20management>
    ...
    <Output_20Variables>
        ...
    </Output_20Variables>
</ParameterHandler>

'''
from xml.etree.ElementTree import ElementTree
from optparse import OptionParser
import sys

if __name__ == "__main__":
    
    # Setup of the command line options:
    usage = "usage: %prog [options] filename.dat"
    parser = OptionParser(usage)
    dumby ,fileName = parser.parse_args(sys.argv[1:])
    fileName = str(fileName[0])
    #print fileName
    
    tree = ElementTree()
    tree.parse(fileName) 
    sysManTrunk = tree.findall('System_20management')[0] #managment trunk
    
    rootSections = [tree.findall('Output_20Variables')[0],sysManTrunk.findall('Solution_20variables')[0],sysManTrunk.findall('Equations')[0] ]
    for r in rootSections:
        for subSection in list(r):
            #print subSection.tag
            for subSubSection in list(subSection):
                #print "   " + subSubSection .tag
                if subSubSection.tag == "value":
                    if subSubSection.text == "no_name":
                        #print "Removing subsection: " + str(subSection)
                        r.remove(subSection)
    
    
    tree.write(fileName)