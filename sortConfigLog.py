#!/usr/bin/env python

# @brief Sort configuration log removing duplicates 



import sys
import os


#Python procedures to sort config log and remove duplicates

def read_config(filename):
  '''
     Read config log and returns dictionary { vriable : value }
     @param filename configuration file in the format variable = value
  '''
  f = open(filename)
  config_dict = {}
  for lines in f:
    if not lines.startswith("#"):
      items = lines.strip().split('= ', 1)
      if items == ['']:
        continue
      config_dict[items[0]] = (items[1])
      
  return config_dict



def print_sorted_dict( d ):
  '''
     Print sorted dict
     @param d dictionary 
  '''
  for key in sorted(d.iterkeys()):
    print("%s: %s" % (key, d[key]))



def save_sorted_dict_in_file( d, filename ):
  '''
     Save sorted dict d into file
     @param d  dictionary to be sorted and saved
     @param filename file to save dictionary
  '''
  file_ = open(filename, 'w')
  tmp = ""
  for key in sorted(d.iterkeys()):
    tmp = key + " = " + d[key] + '\n'
    file_.write( tmp )
  file_.close()


def save_dict_in_file( d, filename ):
  '''
     Save dict d into file
     @param d  dictionary to be saved
     @param filename file to save dictionary
  '''
  file_ = open(filename, 'w')
  tmp = ""
  for key in d.iterkeys():
    tmp = key + " = " + d[key] + '\n'
    file_.write( tmp )
  file_.close()


def read_config_file_option(filename, option):
  '''
    Read config file of the type 
    option = value
    and extract value for option
    
    @param filename configuration file
    @param option option to be extracted
  '''  
  f = open(filename, 'r')
  config_dict = {}
  for lines in f:
    if not lines.startswith("#"):
      items = lines.strip().split('= ', 1)
      if items == ['']:
        continue
      config_dict[items[0].strip()] = (items[1].split('#', 1)[0].strip())
  value = config_dict[option]
  return value
  



def main(argv):
  '''
     Main method
  '''
  #read Generator.dat and extract name of config log
  logname =  read_config_file_option('Generator.dat', 'ConfigLogFile' )
  #sortedLogname = "sorted" + logname
  log = read_config(logname)
  #print_sorted_dict( log )
  #save_sorted_dict_in_file( log, sortedLogname )
  #save_dict_in_file( log, sortedLogname )
  #Remove oryginal config log and save the new one
  os.remove(logname)
  save_sorted_dict_in_file( log, logname )
  #save_dict_in_file( log, logname )
  #print( "Read: " + logname + " and write result in: " + sortedLogname)





if __name__ == "__main__":
  main(sys.argv)
