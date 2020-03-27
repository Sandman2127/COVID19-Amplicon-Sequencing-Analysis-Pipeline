#!/usr/bin/env python

import os, sys, re, shutil
import argparse
parser = argparse.ArgumentParser(description='Build Individual Report for Each Patient')
#parser.add_argument('-inputDir', type=str, help='The absolute path to the input directory', required=True)
#parser.add_argument('-templateHtmlDir', type=str, help='The absolute path to the input directory', required=True)

#args = parser.parse_args()
#inputDir = args.inputDir
#templateHtmlDir = args.templateHtmlDir
cwd = os.getcwd()
match = re.search(r'(.*/COVID19-Amplicon-Sequencing-Analysis-Pipeline/lib).*',str(os.path.dirname(__file__)))
if match:
    lib = match.group(1)
else:
    print("lib directory not found, use the full path to the buildIndivReport.py, exiting now")
    exit
reportDir = lib + "/indiv_report"

def main():
    #TODO: establish position on the OS
    print("The lib for COVID19 amplicon analysis directory is:")
    print(lib)
    print("Your current working directory is:")
    print(cwd)
    print("The report directory to copy is:")
    print(reportDir)
    #TODO: copy the entire template HtmlDir into the current working directory and check that it worked
    Nreportdir = cwd + "/indiv_report"
    copystdOut = "Copying: " + reportDir + " to " + Nreportdir
    print(copystdOut)
    shutil.copytree(reportDir, Nreportdir)
    if os.path.isdir(Nreportdir):
        print("The new report directory exists, the copy was successful")
    else:
        print("The new report directory doesnt exist, the copy was unsuccessful, exiting now")
        exit
    #TODO: copy the appropriate png files into the correct location: i.e. CWD/templateHtmlDir(indiv_report)/images/
    origAmpliconBarPlot = cwd + "/Amplicon_Barplot.png"
    ampliconBarPlot = Nreportdir + "/images/Amplicon_Barplot.png"
    origAmpliconPositionPlot = cwd + "/Amplicon_Pos_Hist.png"
    ampliconPositionPlot = Nreportdir + "/images/Amplicon_Pos_Hist.png"
    shutil.copy(origAmpliconBarPlot,ampliconBarPlot)
    shutil.copy(origAmpliconPositionPlot,ampliconPositionPlot)
    if os.path.isfile(ampliconBarPlot) and os.path.isfile(ampliconPositionPlot) :
        print("The analysis images were successfull copied, proceeding to build report ")
    else:
        print("The analysis images were not copied succesfully, the copy was unsuccessful, exiting now")
        exit
    #TODO: fill in all individual information into the template result.html :
    #{_predicted_result_} {_sample_name_}

if __name__ == "__main__":
    main()



