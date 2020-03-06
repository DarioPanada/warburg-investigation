'''
Script contains utilities to aggregate multiple reports

For consistency, reports should be run for the same number of epochs
'''
import csv
import os
import operator
from core.Toolkit import depickleFromLite
import matplotlib
import numpy as np
from models.alpha04c.scratches.SensitivityTestingTemplate import getSensitivityTestingLatexTemplate
from numpy.polynomial import Polynomial
matplotlib.use("Qt4Agg")
import matplotlib.pyplot as plt
from TabularResultsAnalyzer import getCorrelationWithVariable
import pandas as pd

'''
Given a model object, saves a scatter of the actual vs expected growth curve.
:param rep - The model object
'''
def displayExpectedVsActualCurve(rep):
    cancerCellNum = rep.output["agentNums"]["cancerCells"]
    pWarburgSwitch = rep.properties['agents']['cancerCells']['pWarburgSwitch']


    cancerCellNum = [x*120.*5000/10**9 for x in cancerCellNum]
    # plt.scatter([x/24.*2 for x in range(len(cancerCellNum))], cancerCellNum, label="In-Silico")

    # Empirical data (liu et al)
    xs = [0,13, 24]
    ys = [0, 50, 150]
    p = Polynomial.fit(xs,ys, 2)

    if len(cancerCellNum) < 300:
        cancerCellNum = cancerCellNum + [cancerCellNum[-1] for _ in range(300 - len(cancerCellNum))]

    modelValues = []

    for i,t in enumerate([x/24.*2 for x in range(len(cancerCellNum))]):

        modelValue = cancerCellNum[i]
        modelValues.append(modelValue)

    cancerCellNum = [x*120.*5000/10**9 for x in cancerCellNum]
    plt.scatter([x/24.*2 for x in range(len(modelValues))], modelValues, label="In-Silico")
    plt.plot(*p.linspace(), label="Empirical (Liu et al.)", color="orange")

    plt.xlabel("Time (days)")
    plt.ylabel("Volume (mm3)")
    plt.legend()
    plt.savefig("out/comparison.jpg")

'''
Given a model object, returns the mean error of the growth curve.
:param rep - The model object
:return The error (or none if the model halted before epoch 300)
'''
def getME(rep):
    cancerCellNum = rep.output["agentNums"]["cancerCells"]
    pWarburgSwitch = rep.properties['agents']['cancerCells']['pWarburgSwitch']


    cancerCellNum = [x*120.*5000/10**9 for x in cancerCellNum]
    # plt.scatter([x/24.*2 for x in range(len(cancerCellNum))], cancerCellNum, label="In-Silico")

    # Empirical data (liu et al)
    xs = [0,13, 24]
    ys = [0, 50, 150]
    p = Polynomial.fit(xs,ys, 2)

    diffs = []

    if len(cancerCellNum) < 300:
        print("%s has %s epochs recorded, padding" % (rep.properties["name"], str(len(cancerCellNum))))
        cancerCellNum = cancerCellNum + [cancerCellNum[-1] for _ in range(300 - len(cancerCellNum))]
        print(len(cancerCellNum))



    for i,t in enumerate([x/24.*2 for x in range(len(cancerCellNum))]):

        modelValue = cancerCellNum[i]
        empiricalValue = p(t)
        diffs.append(abs(modelValue-empiricalValue))


    score = np.mean(diffs)

    return score

'''
Given a list of reports, produces a csv where each row is a simulation and each
column the error at a cerain time-point.
:param rep - The reports
:param outputName - The name of the output report csv
'''
def getErrorSeries(rep, outputName):
    allErrs = []
    print("hi")
    padCount = 0
    for expName,k in rep.iteritems():
        cancerCellNum = k.output["agentNums"]["cancerCells"]
        pWarburgSwitch = k.properties['agents']['cancerCells']['pWarburgSwitch']
        minimumOxygenConcentration = k.properties['agents']['cancerCells']['minimumOxygenConcentration']

        cancerCellNum = [x*120.*5000/10**9 for x in cancerCellNum]
        # plt.scatter([x/24.*2 for x in range(len(cancerCellNum))], cancerCellNum, label="In-Silico")

        # Empirical data (liu et al)
        xs = [0,13, 24]
        ys = [0, 50, 150]
        p = Polynomial.fit(xs,ys, 2)


        if len(cancerCellNum) < 300:
            cancerCellNum = cancerCellNum + [cancerCellNum[-1] for _ in range(300 - len(cancerCellNum))]
            print(expName, pWarburgSwitch, minimumOxygenConcentration)
            padCount+=1


        errs = [expName, pWarburgSwitch, minimumOxygenConcentration]

        for i,t in enumerate([x/24.*2 for x in range(len(cancerCellNum))]):

            modelValue = cancerCellNum[i]
            empiricalValue = p(t)
            err = modelValue-empiricalValue
            errs.append(err)


        allErrs.append(errs)
    print("hi")
    try:
        print("HERE")
        df = pd.DataFrame(allErrs)
        df.to_csv("../out/errorSummaries/%s" % outputName, index=False, header=["name","pWarburgSwitch", "minimumOxygenConcentration"]+range(0,300))
        print("Padded: %s" % str(padCount))
    except Exception as e:
        print(str(e))


'''
Given an integer i > 0, always returns the same color. Useful to plot properties of different
experiments on plots so as to be able to correlate them. (Eg: Alive/Dead cell scatters for 5 experiments)
:param i : The index
:return : The color
'''
def getColor(i):
    colors = ["r","g","b","k","m","c", "olive", "deeppink", "chocolate", "peru",
              "teal", "tomato", "lawngreen", "sandybrown", "teal", "fuchsia",
              "maroon", "sienna", "palegreen", "greenyellow", "salmon", "darksalmon",
              "lightsalmon"]

    index = i % len(colors)
    return colors[index]

'''
Loads all pickles lite representing average model properties for a report. Ie as created by the function 
pickleAverageModel. 

:param reportsDir (optional, defaults to ../reports). Directory where the reports are contained.
:param reports (optional, defaults to all). If set to "all" then all reports all loaded. Can also be set to
a list of report names (or part of names). For example, setting to [130818_1535, 130818_1349] will load all reports
where the name includes such string.
:param keyword (optional, defaults to avg) - The keyword that is used to identify the appropriate .pickle file in the
reports directory. (Eg, setting this to avg will load any pickle file containing avg in the name). There should
only be one such file.
:return A dictionary mapping average report names to model lite objects.
'''
def loadAvgReports(reportsDir="../reports", reports="all", keyword="avg"):
    retrievedReports = os.listdir(reportsDir)

    if reports != "all":
        retrievedReports = [re for re in retrievedReports if reports in re]

    # Concatenating names of directories containing reports with main reports directory
    retrievedReports = ["%s/%s" % (reportsDir, r) for r in retrievedReports]
    # Exploring each such directory and only keeping the path to the avg pickle
    retrievedReports = [["%s/%s" % (r,x) for x in os.listdir(r) if x.endswith(".pickle") and keyword in x] for r in retrievedReports if "tmp" not in r]
    # Flattening the list
    retrievedReports = reduce(operator.add, retrievedReports)

    finalReports = dict()

    for p in retrievedReports:
        name = "_".join(p.split("/")[-2].split("_"))[:-1]
        name = "_".join(name.split("_")[:-1])
        pickle = depickleFromLite(p)
        finalReports[name] = pickle

    return finalReports

'''
Works the assumption the experiment has been run once and there is a single pickleLite object in /pickles, compares
multiple reports valid under such an assumption.

:param reportsDir (optional, defaults to ../reports). Directory where the reports are contained.
:param reports (optional, defaults to all). If set to "all" then all reports all loaded. Can also be set to
a list of report names (or part of names). For example, setting to [130818_1535, 130818_1349] will load all reports
where the name includes such string.
:param keyword (optional, defaults to avg) - The keyword that is used to identify the appropriate .pickle file in the
reports directory. (Eg, setting this to avg will load any pickle file containing avg in the name). There should
only be one such file.
:return A dictionary mapping average report names to model lite objects.
'''
def loadFinalReport(reportsDir="/home/dario/phdreports", reports="all", keyword="avg"):
    retrievedReports = os.listdir(reportsDir)

    if reports != "all":
        retrievedReports = [re for re in retrievedReports if reports in re]

    # Concatenating names of directories containing reports with main reports directory
    retrievedReports = ["%s/%s/pickles" % (reportsDir, r) for r in retrievedReports]
    # Exploring each such directory and only keeping the path to the avg pickle
    retrievedReports = [["%s/%s" % (r,x) for x in os.listdir(r) if x.endswith(".pickle")] for r in retrievedReports if "tmp" not in r and os.path.isdir(r)]
    # Flattening the list
    retrievedReports = reduce(operator.add, retrievedReports)

    finalReports = dict()

    for p in retrievedReports:
        name = "_".join(p.split("/")[-3].split("_")[:-1])
        pickle = depickleFromLite(p)
        finalReports[name] = pickle

    return finalReports

'''
Returns a scatter with time-series plots with of total number of cancer cells for each model.
:param models - A dictionary of models (as generated by loadAvgReports) where keys are report names and values
and model lite objects.
:param render (optional, defaults to False), set to true to display the scatter
:return - The scatter object'''
def compareTotalCells(models, render=False):
    totalCells = {
        k:v.output["agentNums"]["cancerCells"] for k,v in models.iteritems()
    }
    fig = plt.figure()
    mng = plt.get_current_fig_manager()

    xsize = ysize = 15

    fig.set_size_inches(xsize,ysize)

    for i,(k,v) in enumerate(totalCells.iteritems()):
        plt.scatter(range(len(v)), v, label=k, color=getColor(i))

    plt.title("Cancer Cells in Time")
    plt.xlabel("Epoch")
    plt.ylabel("Total Number of Cancer Cells")
    plt.legend(loc="best", bbox_to_anchor=(0.7, 1), prop={'size': 10}, ncol=2)

    if render:
        plt.show()

    return plt

'''
Compares cancer agent number for alive, dead and total cancer cell numbers across multiple reports.
:param models - A dictionary of models (as generated by loadAvgReports) where keys are report names and values
and model lite objects.
:param render - (optional, defaults to false) If set to true, scatters are rendered
:return aliveCells - A dictionary where keys are report names and values are lists of number of alive cancer cells at each epoch
:return deadCells - A dictionary where keys are report names and values are lists of number of dead cancer cells at each epoch
:return totalCells - A dictionary where keys are report names and values are lists of number of total cancer  cells at each epoch 
'''
def compareCancerAgentNums(models, render=False):

    totalCells = {
        k:v.output["agentNums"]["cancerCells"] for k,v in models.iteritems()
    }

    deadCells = {
        k:v.output["agentNums"]["deadCancerCells"] for k,v in models.iteritems()
    }

    aliveCells = {
        k:v.output["agentNums"]["aliveCancerCells"] for k,v in models.iteritems()
    }

    if render:

        f,axarr = plt.subplots(1,3, sharey=True)
        for i,(k,v) in enumerate(totalCells.iteritems()):
            axarr[0].scatter(range(len(v)), v, label=k, color=getColor(i))
        axarr[0].set_xlabel("Epoch")
        axarr[0].set_ylabel("Number of Total Cancer Cells")
        axarr[0].set_title("Number of Total Cancer Cells in Time")



        for i,(k,v) in enumerate(aliveCells.iteritems()):
            axarr[1].scatter(range(len(v)), v, label=k, color=getColor(i))
        axarr[1].set_xlabel("Epoch")
        axarr[1].set_ylabel("Number of Alive Cancer Cells")
        axarr[1].set_title("Number of Alive Cancer Cells in Time")


        for i,(k,v) in enumerate(deadCells.iteritems()):
            axarr[2].scatter(range(len(v)), v, label=k, color=getColor(i))
        axarr[2].set_xlabel("Epoch")
        axarr[2].set_ylabel("Number of Dead Cancer Cells")
        axarr[2].set_title("Number of Dead Cancer Cells in Time")
        axarr[2].legend(loc='upper right', ncol=2)
        plt.show()

    return totalCells, deadCells, aliveCells

'''
Compares number of endothelial cells across multiple reports.
:param models - A dictionary of model name => model lite
:param render - (optional, defaults to false) If set to true, displays the scatter
:return  A dictionary of model names to number of endothelial cells in time
'''
def compareNumEndothelial(models, render=False):
    numEndothelial = {
        k:v.output["agentNums"]["tipCells"] for k,v in models.iteritems()
    }

    if render:
        plt.figure()
        for k,v in numEndothelial.iteritems():
            plt.scatter(range(len(v)), v, label=k)

        plt.xlabel("Epoch")
        plt.ylabel("Number of Endothelial Cells")
        plt.title("Number of Endothelial Cells in Time")
        plt.legend(loc=2)
        plt.show()

    return numEndothelial

'''
Compares oxygen concentrations in time
Compares cancer cell properties in time.
:param models - A dictionary of model name => model lite
:param render - (optional, defaults to false) If set to true, displays the scatter
:return Avg oxygen - A dictionary of model names => avg oxygen concentrations
:return Max oxygen - A dictionary of model names => max oxygen concentrations
:return Min oxygen - A dictionary of model => min oxygen concentrations
'''
def compareOxygenConcentrations(models, render=False):
    maxOxygen = {
        k:v.output["cancerCellProperties"]["maxOxygen"] for k,v in models.iteritems()
    }

    minOxygen = {
        k:v.output["cancerCellProperties"]["minOxygen"] for k,v in models.iteritems()
    }

    avgOxygen = {
        k:v.output["cancerCellProperties"]["avgOxygen"] for k,v in models.iteritems()
    }

    if render:
        plt.figure()
        f, axarr = plt.subplots(1,3)

        for k,v in maxOxygen.iteritems():
            axarr[0][0].scatter(range(len(v)), v, label=k)
        axarr[0][0].set_xlabel("Epoch")
        axarr[0][0].set_ylabel("Oxygen Concentration")
        axarr[0][0].set_title("Maximum Oxygen Concentration in Time")
        axarr[0][0].legend(loc=2)

        for k,v in minOxygen.iteritems():
            axarr[0][1].scatter(range(len(v)), v, label=k)
        axarr[0][1].set_xlabel("Epoch")
        axarr[0][1].set_ylabel("Oxygen Concentrations")
        axarr[0][1].set_title("Minimum Oxygen Concentration in Time")
        axarr[0][1].legend(loc=2)

        for k,v in avgOxygen.iteritems():
            axarr[0][2].scatter(range(len(v)), v, label=k)
        axarr[0][2].set_xlabel("Epoch")
        axarr[0][2].set_ylabel("Oxygen Concentration")
        axarr[0][2].set_title("Average Oxygen Concentration in Time")
        axarr[0][2].legend(loc=2)

        plt.show()

    return maxOxygen, minOxygen, avgOxygen

'''
Compares cancer cell properties in time.
:param models - A dictionary of model name => model lite
:param render - (optional, defaults to false) If set to true, displays the scatter
:return  Avg hif - A dictionary of model names => hif expression rate
:return  Avg metabolic rate - A dictionary of model names => metabolic rate
:return  Avg proliferation rate - A dictionary of model names => proliferation rate
:return  Avg vegf secretion rate - A dictionary of model names => vegf secretion rate
'''
def compareCancerCellProperties(models, render=False):
    hifRates = {
        k:v.output["cancerCellProperties"]["avgHif"] for k,v in models.iteritems()
    }

    avgVegf = {
        k:v.output["cancerCellProperties"]["avgVegf"] for k,v in models.iteritems()
    }

    avgMetabolicRate = {
        k:v.output["cancerCellProperties"]["avgMetabolicRate"] for k,v in models.iteritems()
    }

    avgPSynthesis = {
        k:v.output["cancerCellProperties"]["avgPSynthesis"] for k,v in models.iteritems()
    }



    if render:
        plt.figure()
        f, axarr = plt.subplots(2,2)

        for k,v in hifRates.iteritems():
            axarr[0][0].scatter(range(len(v)), v, label=k)
        axarr[0][0].set_xlabel("Epoch")
        axarr[0][0].set_ylabel("Avg HIF Rate")
        axarr[0][0].set_title("Average HIF Expression Rates in Time")
        axarr[0][0].legend(loc=2)

        for k,v in avgVegf.iteritems():
            axarr[0][1].scatter(range(len(v)), v, label=k)
        axarr[0][1].set_xlabel("Epoch")
        axarr[0][1].set_ylabel("Avg VEGF Rate")
        axarr[0][1].set_title("Average VEGF Expression Rates in Time")
        axarr[0][1].legend(loc=2)

        for k,v in avgMetabolicRate.iteritems():
            axarr[1][0].scatter(range(len(v)), v, label=k)
        axarr[1][0].set_xlabel("Epoch")
        axarr[1][0].set_ylabel("Avg Metabolic Rate")
        axarr[1][0].set_title("Average Metabolic Expression Rates in Time")
        axarr[1][0].legend(loc=2)

        for k,v in avgPSynthesis.iteritems():
            axarr[1][1].scatter(range(len(v)), v, label=k)
        axarr[1][1].set_xlabel("Epoch")
        axarr[1][1].set_ylabel("Avg P Synthesis Rate")
        axarr[1][1].set_title("Average P Synthesis Expression Rates in Time")
        axarr[1][1].legend(loc=2)

        plt.show()
        print("HERE")

    return hifRates

'''
Compares tumour size across multiple reports.
:param models - A dictionary of model name => model lite
:param render - (optional, defaults to false) If set to true, displays the scatter
:return  A dictionary of model names to tumour volumes in time
'''
def compareTumourSizes(models, render=False):
    tumourVolumes = {
        k:v.output["maxDistances"] for k,v in models.iteritems()
    }

    if render:
        plt.figure()
        for k,v in tumourVolumes.iteritems():
            plt.scatter(range(len(v)), v, label=k)

        plt.xlabel("Epoch")
        plt.ylabel("Average Tumour Volume")
        plt.title("Average Tumour Volume in Time")
        plt.legend(loc=2)
        plt.show()

    return tumourVolumes

'''
Compares average vegf stimulus in time across multiple reports.
:param models - A dictionary of model name => model lite
:param render - (optional, defaults to false) If set to true, displays the scatter
:return  A dictionary of model names to average vegf stimulus in time
'''
def compareVegfStimulus(models, render=False):
    vegfStimuli = {
        k:v.output["endothelialCellProperties"]["avgVegf"] for k,v in models.iteritems()
    }

    if render:
        plt.figure()
        for k,v in vegfStimuli.iteritems():
            plt.scatter(range(len(v)), v, label=k)

        plt.title("Average VEGF Stimulus in Time")
        plt.xlabel("Epoch")
        plt.ylabel("Average VEGF Stimulus in Time")
        plt.legend(loc=2)
        plt.show()

    return vegfStimuli

'''
Match each experiment contained in a csv to a pickle lite of the average experiment end-state.
:param experimentPath - Relative path to the experiment csv
:param reports (optional, defaults to all). If set to "all" then all reports all loaded. Can also be set to
a list of report names (or part of names). For example, setting to [130818_1535, 130818_1349] will load all reports
where the name includes such string.
:return a dictionary where keys are experiment names, values are tuples where the first element is the list
of parameters (first element being the experiment name) and the second element is the de-pickled model end-state. 
(Derived from a  pickle lite) or None if no report was found.
An additional key "headers" is present that points to a list containing experiment header values.
'''
def matchExperimentReports(experimentPath, reports="all"):

    matchedExps = dict()

    experiments = getExperimentsFromFile_(experimentPath)
    reports = loadFinalReport(reports=reports)

    matchedExps["header"] = experiments[0]

    del experiments[0]

    reportNames = reports.keys()
    for e in experiments:
        expName = e[0]

        matchedReports = [r for r in reportNames if expName == r]

        if len(matchedReports) > 1:
            print("Error, experiment %s matches multiple reports: %s, SKIPPING" % (expName, ",".join(matchedReports)))
            continue

        if len(matchedReports) == 0:
            print("Experiment %s doesn't match any report!" % expName)
            report = None
        else:
            report = matchedReports[0]
            print("Matching %s to report %s" % (expName, report))
            report = reports[report]

        matchedExps[expName] = (e,report)

    return matchedExps

'''
Produces a csv file where each row represents an experiment with an added column that represents the number of epochs
 the simulation ran for and an added column that represents the output (as in
final number of cancer cells agents). Also ME is added as a column. Rows are sorted by ME in descending order
:param matchedReports a dictionary where keys are the experiment name and values are tuples where the first value
in a list of experiment parameters and the second value is a de-pickled pickle lite. 
There is an expectation that one of the keys is "header" and points to a list of column names.
(As produced by matchExperimentReports)
:param outPath The relative path of the output csv INCLUSIVE OF THE CSV NAME (eg: ../out/outA.csv)
:return - A list of where each element represents an experiment's properties PLUS the final tumour volume
'''
def outputExperimentsFinalNumAgentsAndME(matchedReports, outPath):

    expsWithOutput = []

    with open(outPath, 'a') as out:

        header = matchedReports["header"]
        header.append("numEpochs")
        header.append("finalNumCancerAgents")
        header.append("meanError")
        out.write(",".join(header))
        out.write("\n")
        del matchedReports["header"]

        for k,experiment in matchedReports.iteritems():
            model = experiment[1]

            if model is None:
                continue

            numAgents = model.output["agentNums"]["cancerCells"][-1]

            experimentProperties = experiment[0]
            experimentProperties.append(len(model.output["agentNums"]["cancerCells"]))
            experimentProperties.append(numAgents)
            experimentProperties.append(getME(model))

            expsWithOutput.append(experimentProperties)

        expsWithOutput.sort(key = lambda x : x[-1], reverse=True)

        for e in expsWithOutput:
            out.write(",".join(map(str,e)))
            out.write("\n")


    return expsWithOutput

'''
Loads all experiments from csv to a list of lists
:param experimentPath - Relative path to the experiment csv
:return A list of lists where each element represent an experiment's values. THE FIRST ROW IS THE HEADER
'''
def getExperimentsFromFile_(experimentPath):

    with open(experimentPath, 'r') as f:
        exps = csv.reader(f)
        return list(exps)

'''
Given a set of reports, generates a scatter comparing growth curves, creates
a csv adding final num agents to the experiments and calculates pearson correlation.

Given an out directory, within it a subdirectory with the experiment's name will be created. Inside it,
the following three files will be placed:

- growthCurves.png
- report.csv (A csv of experiments with a new column for finalAgentNum)
- correlation.txt (Containing the pearson correlation value)

:param reports - A list of depickled models, as generated by loadAvgReports or loadFinalReport
:param experimentsFile - Relative path to the original file detailing the experiments which generated the reports
:param targetVariable - The variable we are correlating
:param outDir - (optional, defaults to ../out) Relative path of the out dir where a subdirectory for this specific
experiment will be created
:return latex string
'''
def generateGraphAndCorrelationForExperiment(reports, experimentsFile, targetVariable, outDir="../out"):

    # Setup
    print("---")
    expName = "_".join(reports.keys()[0].split("_")[1:])
    print("I am starting work on " + expName)
    outPath = "%s/%s" % (outDir, expName)

    if os.path.isdir(outPath):
        print("Out path %s already exists!" % outPath)
        print("---")
        return

    os.mkdir(outPath)

    # Generating graph
    plt = compareTotalCells(reports)
    plt.savefig("%s/growthCurves_%s.png" % (outPath, targetVariable))

    # Generating output csv
    matchedReports = matchExperimentReports(experimentsFile)
    reportPath = "%s/report.csv" % outPath

    headerIndex = matchedReports["header"].index(targetVariable)

    outputExperimentsFinalNumAgents(matchedReports, reportPath)

    # Getting min and max values for experiment
    expVals = [e[0][headerIndex] for e in matchedReports.values()]
    minVal = min(expVals)
    maxVal = max(expVals)

    # Calculating correlation
    correlation = getCorrelationWithVariable(reportPath, targetVariable)
    with open("%s/correlation.txt" % outPath, 'w') as f:
        f.write(targetVariable + ":" + str(correlation))
        f.close()


    latex = getSensitivityTestingLatexTemplate(targetVariable, "growthCurves_%s.png" % targetVariable, correlation, minVal, maxVal)

    with open("%s/latex.txt" % outPath, 'w') as f:
        f.write(latex)
        f.close()
    print("---")

    return latex


'''
Generates correlation reports for multiple experiment.

:param experiments - A list of experiment objects containing keys
keyword - Used to identify reports
experimentFile - Relative path of the original experiment file
targetVariable - Name of target variable for correlation

Eg:

{
    "keyword":"boer",
    "experimentFile":"../experiments/experiments_170519_boer.csv", 
    "targetVariable":"baseOxygenEmissionRate"
}

Further creates a latex.txt file in outDir with the combined latex of all correlations.
'''
def bulkGenerateGraphAndCorrelationForExperiment(experiments, outDir="../out"):

    latexes = []

    for experiment in experiments:
        reports = loadAvgReports(reports=experiment["keyword"])
        latex = generateGraphAndCorrelationForExperiment(reports, experiment["experimentFile"], experiment["targetVariable"])

        latexes.append(latex)

    with open(outDir + "/latex.txt", 'w') as f:
        for latex in latexes:
            f.write(latex+"\n")

        f.close()

'''
Given multiple deplickled reports, produces a scatter
with dt value on the x-axis and number of epochs for completion
on the y-axis.
:param reports - Array of depickled reports
'''
def scatterNumEpochsFromReports(reports):
    points = []

    for report in reports.values():
        points.append((report.properties["diffusion"]["dt"], report.current_epoch))

    points.sort(key=lambda x : x[0])
    plt.scatter([x[0] for x in points], [x[1] for x in points])
    plt.xlabel("dt")
    plt.ylabel("Number of Epochs for Max Growth")
    plt.title("dt impact on number of epochs to achieve max tumour growth")
    plt.show()

if __name__ == "__main__":

    experiments = [
            {
                "keyword" : "_bomr_",
                "experimentFile" : "../experiments/experiments_200519_bomr.csv",
                "targetVariable" : "baseOxygenMetabolicRate"
            },
            {
                "keyword" : "hypoxic_threshold",
                "experimentFile" : "../experiments/experiments_200519_hypoxicThreshold.csv",
                "targetVariable" : "hypoxicThreshold"
            },
            {
                "keyword" : "minGlucoseNonWarburg",
                "experimentFile" : "../experiments/experiments_200519_minGlucoseNonWarburg.csv",
                "targetVariable" : "minGlucoseNonWarburg"
            },
            {
                "keyword" : "minGlucoseUptakeRate",
                "experimentFile" : "../experiments/experiments_200519_minGlucoseUptakeRate.csv",
                "targetVariable" : "minGlucoseUptakeRate"
            },
            {
                "keyword" : "minOxygenConcentration",
                "experimentFile" : "../experiments/experiments_200519_minOxygenConcentration.csv",
                "targetVariable" : "minimumOxygenConcentration"
            },
            {
                "keyword" : "numendo",
                "experimentFile" : "../experiments/experiments_200519_numEndo.csv",
                "targetVariable" : "numEndothelialCells"
            },
            {
                "keyword" : "ultraHypoxicThreshold",
                "experimentFile" : "../experiments/experiments_200519_ultraHypoxicThreshold.csv",
                "targetVariable" : "ultraHypoxicThreshold"
            },
        ]
    '''[
            {
                "keyword":"boer",
                "experimentFile":"../experiments/experiments_170519_boer.csv",
                "targetVariable":"baseOxygenEmissionRate"
            },
            {
                "keyword":"pwarburg",
                "experimentFile":"../experiments/experiments_170519_pwarburg.csv",
                "targetVariable":"pWarburgSwitch"
            },
            {
                "keyword":"iterations",
                "experimentFile":"../experiments/experiments_170519_diffusionIterations.csv",
                "targetVariable":"diffusionSolveIterations"
            },
            {
                "keyword":"edd",
                "experimentFile":"../experiments/experiments_170519_edd.csv",
                "targetVariable":"endothelialDivisionDelay"
            },
            {
                "keyword":"vegf",
                "experimentFile":"../experiments/experiments_170519_minvegf.csv",
                "targetVariable":"minimumVegfConcentration"
            },
            {
                "keyword":"psynthesis",
                "experimentFile":"../experiments/experiments_170519_minpsynthesis.csv",
                "targetVariable":"minPSynthesis"
            }
        ]'''

    #bulkGenerateGraphAndCorrelationForExperiment(experiments)
    # Useful command to copy all pngs: find . -name *.png -exec cp {} ~/phdwritings/parameter-summary/ \;

    reports_pws_80_moc_18 = loadAvgReports(reportsDir="/home/dario/phdreports", reports="0.8_moc_18")
    reports_pws_005_moc_18 = loadAvgReports(reportsDir="/home/dario/phdreports", reports="0.05_moc_18")
    #exploreClcCc(reports)

    #scatterNumEpochsFromReports(reports)

    #generateGraphAndCorrelationForExperiment(reports, "../experiments/experiments_170519_boer.csv", "baseOxygenEmissionRate")

    #totalCells, deadCells, aliveCells = compareCancerAgentNums(reports, render=True)
    # reports = {
    #     reports_pws_80_moc_18.keys()[0] : reports_pws_80_moc_18[reports_pws_80_moc_18.keys()[0]],
    #     reports_pws_20_moc_18.keys()[0] : reports_pws_20_moc_18[reports_pws_20_moc_18.keys()[0]]
    # }
    #compareTotalCells(reports, render=True)
    #endothelials = compareNumEndothelial(reports, render=True)
    #compareCancerCellProperties(reports, render=True)
    #compareTumourSizes(reports, render=True)
    #compareVegfStimulus(reports, render=True)
    #compareOxygenConcentrations(reports, render=True)
    #matchedReports = matchExperimentReports("../experiments/experiments_071019_hyp_testing_full.csv")
    #expsWithOutput = outputExperimentsFinalNumAgentsAndME(matchedReports, "../out/experiments_071019_hyp_testing_full.csv")
    displayExpectedVsActualCurve(reports_pws_80_moc_18[reports_pws_80_moc_18.keys()[1]])
    print(reports_pws_80_moc_18.keys()[1])
    print(reports_pws_005_moc_18.keys()[0])
    print(getME(reports_pws_80_moc_18[reports_pws_80_moc_18.keys()[1]]))
