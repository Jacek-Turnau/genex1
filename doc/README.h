/*!
 @mainpage



Program name:  "GenEx" Generator



@authors                                                 

Radoslaw Kycia
rkycia@pk.edu.pl
Cracow University of Technology
Krakow, Poland

Janusz Chwastowski
jchwastowski@pk.edu.pl                       
Cracow University of Technology
Krakow, Poland

Rafal Staszewski
rafal.staszewski@gmail.com                       
Institute of Nuclear Physics, Polish Academy of Sciences 
Krakow, Poland                                                   

Jacek Turnau 
jacek.turnau@ifj.edu.pl                        
Institute of Nuclear Physics, Polish Academy of Sciences 
Krakow, Poland   

\brief

Generator provides universal scheme for HEP generation. It is configured for models of exclusive production of hadron pairs.
Currently it implements model of Lebiedowicz and Szczurek [1] (see also [2]) of reaction (see [3]) \f$2p->2p + \pi^{+} + \pi^{-}\f$ and \f$2p->2p + K^{+} + K^{-}\f$, but it is not restricetd to these models.


@section bibliography Bibliography:

[1] P. Lebiedowicz, A. Szczurek ,,Exclusive \f$2p->2p + \pi^{+} + \pi^{-}\f$ reaction: from the threshold to LHC'', arXiv:0912.0190

[2] L.A. Harland-Lang, V.A. Khoze , M.G. Ryskin, W.J. ,,Stirling The phenomenology of central exclusive production at hadron colliders'', arXiv:1204.4803 [hep-ph] 

[3] A. Szczurek, P. Lebiedowicz, arXiv:1110.4787[hep-ph]

[4] S. Jadach ,,Rapidity Generator for Monte-Carlo Calculations of Cylindrical Phase Space'' Comp. Phys. Comm. 9 (1975) 297

[5] S. Jadach ,,Foam: A General-Purpose Cellular Monte Carlo Event Generator'', arXiv:physics/0203033 [physics.comp-ph]

[6] ROOT framework, (access date: 04/01/2013): http://root.cern.ch/drupal/

[7] TFoam documentation (access date: 04/01/2013): http://root.cern.ch/root/html/TFoam.html

[8] Gnu Scientific Library, (access date: 04/01/2013):  http://www.gnu.org/software/gsl/

[9] Doxygen documentation, (access date: 04/01/2013): http://www.stack.nl/~dimitri/doxygen/


@section dependence Dependence

Program depends highly on ROOT framework [6]. Every properly installed ROOT with FOAM generator(default for new ROOTs)
should work. For FOAM consult [5].

For compiling (GNU) make system and c++ compiler is required - if ROOT was installed properly from 
the source then all necessary things to run generator should be installed.


There is a simple (GNU) Makefile available. Before run please check if locations to ROOT folder, ROOT includes and
ROOT libraries in your linking system were set properly.

\warning
1. Underlaying random number generator TRandom3 from ROOT package has the period = 2^19937-1. It is big but finite number. It should be large enough for any today's applications, however anyone who seriously consider generation of much larger set of events should change the generator, change the seed/use the presistency or use hardware random number generators.
2. It is assumed that first two particles are leading ones.

@section commands  Commands

- To compile and run type:  \b make \b all 
  
  Executable is in ./bin directory

- To compile and run type:  \b make \b run

- To run example test type:  \b make \b test
  
  Restoring oryginal Generator.dat (with comments) is by selection 0 option.

- To restore folder to original state type: \b make \b clean

- To count lines of the program type:  \b make \b line_count

- To generate Doxygen documentation [9] type: \b make \b Generate-doc


@section dirStruct  Directory structure

 - \b ./ - main directory contains Makefile, configuration file \b Generator.dat and scripts that perform various operations.
 - \b /bin - contains main.x generator program appears after compilation. It is wiped by \b make \b clean command.
 - \b /include - contains header files.
 - \b /src - contains source files and the main generator file \b main.cxx
 - \b /build - is temporary directory for building. It is created by main run command and deleted with all the content by \b make \b clean command.
 - \b /EVENTS - contains files with saved events form the run.
 - \b /MODEL_DATA - contains configuration (text) files for different models (e.g. matrix elements) used in generation.
 - \b /ACCEPTANCE_CUTS - contains files with cuts.
 - \b /Tests - contains sorted generator configuration files sortedConfigLog.log which are used in the second mode of generator running \b make \b test 

@section setup Setup

The main configuration file is \b Generator.dat. It is heavily commented.

*/
