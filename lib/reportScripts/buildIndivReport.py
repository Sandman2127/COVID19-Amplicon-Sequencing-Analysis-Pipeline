#!/usr/bin/env python

import os, sys, re, shutil
import argparse
# parser = argparse.ArgumentParser(description='Build Individual Report for Each Patient')
# parser.add_argument('-samOutSummary', type=str, help='The sample_summary.txt file from analyzeSam.py <sampleID>,<>', required=True)

# #args = parser.parse_args()
# #inputDir = args.inputDir
# #templateHtmlDir = args.templateHtmlDir
cwd = os.getcwd()
match = re.search(r'(.*/COVID19-Amplicon-Sequencing-Analysis-Pipeline/lib).*',str(os.path.dirname(__file__)))
if match:
    lib = match.group(1)
else:
    print("lib directory not found, use the full path to the buildIndivReport.py, exiting now")
    exit
reportDir = lib + "/indiv_report"


def replace_dict(template,dictA):  # pass it a dictionary with two prepped items it passes back replaced template
    for key,value in dictA.items():
        template=template.replace(key,str(value))
    return template


def parseSummaryFile():
    sampSumfile=cwd + "/sample_summary.txt"
    if os.path.isfile(sampSumfile):
        with open(sampSumfile,'r+') as fileIn:
            line = fileIn.read()
            line = line.strip()
            line = line.split(",")
            SampleID,Diagnosis=str(line[0]),str(line[1])
            listR = [SampleID,Diagnosis]
            return listR 
    else:
        print("Sample summary file from analyzeSam.py doesnt exist the report will not include a diagnosis")
        listE = [str(0),str(Unknown)]
        return listE

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
    #TODO: fill in all individual information into the template result.html:
    parsedResults = parseSummaryFile()
    sampID,diagnosis = parsedResults[0],parsedResults[1]
    predictedDiagnosis={"{_predicted_result_}":diagnosis}
    sampleID={"{_sample_name_}":sampID}

    replacement_pairs=[predictedDiagnosis,sampleID]
    
    #TODO: readin template and replace template sections:
    REPORT_TEMPLATE = Nreportdir + "/result.html" 
    with open(REPORT_TEMPLATE,'r') as f:
        template=f.read()
        f.close()

    for dict_pair in replacement_pairs:
        template=replace_dict(template,dict_pair)

    fout = open("result.html",'w+')
    fout.write(template)
    fout.close()

    fout_loc = cwd + "/result.html"
    if os.path.isfile(fout_loc):
        os.remove(REPORT_TEMPLATE)
        shutil.copy(fout_loc,REPORT_TEMPLATE)

if __name__ == "__main__":
    main()



