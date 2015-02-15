#! /usr/bin/env python


# @brief Script to run apropriate tests.

import os
import distutils.core
from shutil import rmtree, copyfile




"""Avaiable tests dictionary in the format no_of_test : name_of_test"""
tests = {0:"default Generator.dat with lot of comments and explanations", 
	1:"RHIC pt_pi, eta_pi; tecm = 200GeV; Lambda2=1", 
	2:"RHIC pt_pi, eta_pi, t1, t2; tecm = 200GeV; Lambda2=1", 
	3:"RHIC pt_pi, eta_pi, t1, t2; tecm = 200GeV; Lambda2=1.6", 
	4:"RHIC pt_pi, eta_pi; tecm = 500GeV; Lambda2=1", 
	5:"RHIC pt_pi, eta_pi, t1, t2; tecm = 500GeV; Lambda2=1", 
	6:"RHIC pt_pi, eta_pi, t1, t2; tecm = 500GeV; Lambda2=1.6",
	7:"LHC pt_pi, eta_pi; tecm = 7TeV, 1st; Lambda2=1.2", 
	8:"LHC pt_pi, eta_pi; tecm = 7TeV, 1st; Lambda2=1.6", 
	9:"LHC pt_pi, eta_pi; tecm = 7TeV, 2nd; Lambda2=1.2",
	10:"LHC pt_pi, eta_pi; tecm = 7TeV, 2nd; Lambda2=1.6", 
	11:"LHC pt_K, eta_K; tecm = 7TeV, 1st; Lambda2=1.2",
	12:"LHC pt_K, eta_K; tecm = 7TeV, 1st; Lambda2=1.6", 
	13:"LHC pt_K, eta_K; tecm = 7TeV, 2nd; Lambda2=1.2",
	14:"LHC pt_K, eta_K; tecm = 7TeV, 2nd; Lambda2=1.6",
	15:"2to5; y_pi, tecm = 200GeV", 
	16:"CPS, N=5, y_pi, tecm = 200GeV", 
	17:"2to5; y_pi, t, tecm = 200GeV",
	18:"CPS, N=5, y_pi, t, tecm = 200GeV",
	19:"CPS, N=5, Exploration Cuts, y_pi, t, tecm = 200GeV",
	20:"RHIC pt_pi, eta_pi, t1, t2; tecm = 200GeV; Lambda2=1.6; LS method of Phase Space generation",
	21:"RHIC pt_pi, eta_pi, t1, t2; tecm = 200GeV; Lambda2=1.6; 2toN (N=4) method of Phase Space generation",
	22:"RHIC pt_pi, eta_pi, t1, t2; tecm = 200GeV; Lambda2=1.6; nCells = 1000 = nSampl, y in [-8;8]",
	23:"RHIC pt_pi, eta_pi, t1, t2; tecm = 200GeV; Lambda2=1.6; nCells = 10000, nSampl = 1000, y in [-8;8]",
	24:"RHIC pt_pi, eta_pi, t1, t2; tecm = 200GeV; Lambda2=1.6; nCells = 10000, nSampl = 10000, y in [-8;8]",
	25:"RHIC pt_pi, eta_pi, t1, t2; tecm = 200GeV; Lambda2=1.6; nCells = 10000, nSampl = 1000, y in [-2;2]",
	26:"RHIC pt_pi, eta_pi, t1, t2; tecm = 200GeV; Lambda2=1.6; nCells = 10000 = nSampl, y in [-2;2]"
	}


def prepareTest( number, testDir = './Tests', testName = 'test', configFile = 'Generator.dat' ):
    """Prepare configuration file by picking one of the test file from testDir
    @param testDir  dir containing tests
    @param testName basename of test    
    @param configFile configuration file for generator
    """
    #save old config file
    copyfile(configFile, "OLD"+configFile)
    #copy configuration files
    fromDirectory = testDir + '/' + testName + str(number)
    copyfile(fromDirectory, configFile)
    return testDir    
    
    
def rmDir( directory = "./" ):
    """Remove directory and all its content
    @param directory base directory for project    
    """
    rmtree( directory )


def runMake( option, runDir = './'):
    """Run make with option in given directory
    @param option option for make
    @param runDir directory in which make will be executed    
    """
    os.system( 'make -C ' + runDir + ' ' + option)


def showTests( testDict ):
    """Show tests in dictionary
    @param testDict dictionary with tests in the format  no_of_test : name_of_test   
    """
    print("#########################") 
    print("AVAIBLE TESTS:")
    print("#########################") 
    for key in testDict:
        print( str(key) + ' -- '+ str(testDict[key]) )
    print("#########################")    

def pickTest( testDict ):
    """Allows user to pick option from the keys of dictionary and returns it
    @param testDict dictionary with tests in the format  no_of_test : name_of_test
    """
    finish = False    
    while not finish:        
        showTests(testDict)
        input_var = input("Enter option: ")
        print ("you entered " + str(input_var))
        if input_var in testDict.keys():
                finish = True        
    return input_var           
    




def main():
    """Simple test suit for GenEx. It copy Core GenEx files and selected 
    configuration files to one test directory and then run it and remove test directory"""
    testNo = pickTest( tests )
    print("Preparing generator...")
    prepareTest(testNo)
    print("...DONE")
    print("Start test...")
    runMake('run')
    print("...DONE")
    print("Cleaning dir...")
    runMake('clean')
    print("...DONE")
    



if __name__ == "__main__":
    main()
