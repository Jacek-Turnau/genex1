/// Simple ROOT script to browse/manipulate histograms created in generator
//
//There is  memory leak but it is the only way that browser "survive" the end of the script.


#include <iostream>

using namespace std;

browseHistograms()
{

 TFile *file = new TFile("histograms.root", "");
 
 cout << "Content of the histograms.root:" << endl;
 file->ls();


 TBrowser * browser = new TBrowser();

}
