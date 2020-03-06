from random import randint

from panaxea.core.Environment import ObjectGrid3D, NumericalGrid3D
from panaxea.core.Model import Model
from model.agents.CancerCell import CancerCell
from model.agents.EndothelialCell import TrunkCell, TipCell
from model.helpers.AgentCounter import AgentCounter
from model.helpers.CancerCellWatcher import CancerCellWatcher
from model.helpers.ExitConditionWatcher import ExitConditionWatcher
from model.helpers.GlucoseConcentrationWatcher import GlucoseConcentrationWatcher
from model.helpers.GlucoseDiffusionHelper import GlucoseDiffusionHelper
from model.helpers.OxygenConcentrationWatcher import OxygenConcentrationWatcher
from model.helpers.OxygenDiffusionHelper import OxygenDiffusionHelper
from model.helpers.TumourVolumeWatcher import TumourVolumeWatcher
from model.helpers.VegfDiffusionHelper import VegfDiffusionHelper
from model.helpers.VegfStimulusWatcher import VegfStimulusWatcher
from model.utils.OxygenHIFRelationsGenerator import OxygenHIFRelationsGenerator
from model.agents.HealthyCell import HealthyCell
from model.helpers.DeathCauseWatcher import DeathCauseWatcher


def generate_properties(properties_values):

    properties = dict()

    maxAgentDensity, oxygenDiffusivity, vegfDiffusivity, G1, S, G2, M, endothelialDivisionDelay, \
    minimumVegfConcentration, baseOxygenEmissionRate, minHIF, maxHIF, ultraHypoxicThreshold, \
    hypoxicThreshold, minPSynthesis, minimumOxygenConcentration, baseHifRate, numCancerCells, \
    numEndothelialCells, glucoseDiffusivity, glucoseSecretionRate, minGlucoseUptakeRate, maxGlucoseUptakeRate, pWarburgSwitch, \
    enhancedHypoxicThreshold, minGlucoseWarburg, baseOxygenMetabolicRate, minGlucoseNonWarburg, \
    healthyTissueOxygenUptakeRate, drugName, drugSecretionRate, drugDiffusivity, drugUptakeRate, \
    minDrugDosage, dt, diffusionSolveIterations, envSize, = properties_values

    initialAgentSetup = dict()

    initialAgentSetup["numCancerCells"] = numCancerCells
    initialAgentSetup["numEndothelialCells"] = numEndothelialCells

    properties["initialAgentSetup"] = initialAgentSetup

    # Environment properties, including names, etc.
    envNames = dict()
    envNames["agentEnvName"] = "agentEnv"
    envNames["oxygenEnvName"] = "oxygenEnv"
    envNames["vegfEnvName"] = "vegfEnv"
    envNames["glucoseEnvName"] = "glucoseEnv"
    envNames["drugEnvName"] = "drugEnv"

    properties["envNames"] = envNames

    properties["envSize"] = envSize

    properties["maxAgentDensity"] = maxAgentDensity

    # Properties of agents
    agents = dict()
    agents["baseCellCycleLength"] = {
        "G1": G1,
        "S": S,
        "G2": G2,
        "M": M
    }

    cancerCells = dict()

    cancerCells["HIFRange"] = [0.0, maxHIF]

    cancerCells["domains"] = {
        "ultraHypoxic": ultraHypoxicThreshold,
        "warburgHypoxic": enhancedHypoxicThreshold,
        "hypoxic": hypoxicThreshold
    }

    ohrg = OxygenHIFRelationsGenerator(minHIF=minHIF, maxHIF=maxHIF,
                                       ultraHypoxiaThreshold=ultraHypoxicThreshold, hypoxiaThreshold=hypoxicThreshold,
                                       enhancedHypoxicThreshold=enhancedHypoxicThreshold,
                                       baseOxygenMetabolicRate=baseOxygenMetabolicRate, minPSynthesis=minPSynthesis)

    ultraHypoxiaCoeffs, hypoxiaCoeffs = ohrg.getOxygenToHif()
    warburgHypoxicCoeffs = ohrg.getOxygenToHifWarburg()
    cancerCells["pWarburgSwitch"] = pWarburgSwitch
    cancerCells["baseHifRate"] = baseHifRate
    cancerCells["minGlucoseUptakeRate"] = minGlucoseUptakeRate
    cancerCells["maxGlucoseUptakeRate"] = maxGlucoseUptakeRate
    cancerCells["minGlucoseWarburg"] = minGlucoseWarburg
    cancerCells["minGlucoseNonWarburg"] = minGlucoseNonWarburg
    cancerCells["drugUptakeRate"] = drugUptakeRate
    cancerCells["minHIF"] = minHIF

    cancerCells["oxygenToHifCoeffs"] = {
        "hypoxic": hypoxiaCoeffs,
        "warburg": warburgHypoxicCoeffs,
        "ultraHypoxic": ultraHypoxiaCoeffs
    }

    def drugReactFunction(model, currentPos):
        drugAtPos = model.environments[model.properties["envNames"]["drugEnvName"]].grid[currentPos]
        pass

    cancerCells["drugReactFunction"] = drugReactFunction

    cancerCells["hifToMetabolicRateCoeffs"] = ohrg.getHifToMetabolicRate()

    # Minimum probability of progressing into synthesis
    cancerCells["minPSynthesis"] = minPSynthesis
    cancerCells["hifToProliferationRateCoeffs"] = ohrg.getHifToPSynthesis()

    cancerCells["hifToVegfSecretionRateCoeffs"] = ohrg.getHifToVegf()

    # Minimum oxygen concentration for survival
    cancerCells["minimumOxygenConcentration"] = minimumOxygenConcentration
    cancerCells["maxVegfSecretionRate"] = 10

    agents["cancerCells"] = cancerCells

    agents["healthyTissues"] = {
        "oxygenUptakeRate": healthyTissueOxygenUptakeRate
    }

    endothelialCells = dict()

    agents["endothelialCells"] = endothelialCells

    # Minimum vegf concentration to sprout angiogenesis
    endothelialCells["minimumVegfConcentration"] = minimumVegfConcentration
    # Minimum age of endothelial cells for sprouting
    endothelialCells["divisionDelay"] = endothelialDivisionDelay
    # Measured in mmHg
    endothelialCells["baseOxygenEmissionRate"] = baseOxygenEmissionRate

    endothelialCells["glucoseSecretionRate"] = glucoseSecretionRate
    endothelialCells["drugSecretionRate"] = drugSecretionRate

    properties["agents"] = agents

    diffusion = dict()

    diffusion["oxygenDissociationCurveCoeffs"] = {
        "a1": -8.5322289 * 10.0 ** 3.0,
        "a2": 2.1214010 * 10.0 ** 3.0,
        "a3": -6.7073989 * 10.0 ** 1.0,
        "a4": 9.3596087 * 10.0 ** 5.0,
        "a5": -3.1346258 * 10.0 ** 4.0,
        "a6": 2.3961674 * 10.0 ** 3.0,
        "a7": -6.7104406 * 10.0 ** 1.0

    }

    diffusion["oxygenDiffusivity"] = oxygenDiffusivity
    diffusion["vegfDiffusivity"] = vegfDiffusivity
    diffusion["glucoseDiffusivity"] = glucoseDiffusivity
    diffusion["dt"] = dt
    diffusion["diffusionSolveIterations"] = diffusionSolveIterations

    properties["diffusion"] = diffusion

    properties["drug"] = {
        "name": drugName,
        "secretionRate": drugSecretionRate,
        "diffusivity": drugDiffusivity,
        "uptakeRate": drugUptakeRate,
        "minDrugDosage": minDrugDosage
    }

    return properties

def generate_model(numEpochs=5, properties=dict(), verbose=False, defaultOxygenConc=1, numCancerCells=1,
                   defaultGlucoseConc=1, endothelialRadius=1):

    model = Model(numEpochs, verbose=verbose)

    model.properties = properties

    xsize = ysize = zsize = model.properties["envSize"]

    # Adding environments
    agentEnv = ObjectGrid3D(model.properties["envNames"]["agentEnvName"], xsize, ysize, zsize, model)
    oxygenEnv = NumericalGrid3D(model.properties["envNames"]["oxygenEnvName"], xsize, ysize, zsize, model)
    vegfEnv = NumericalGrid3D(model.properties["envNames"]["vegfEnvName"], xsize, ysize, zsize, model)
    glucoseEnv = NumericalGrid3D(model.properties["envNames"]["glucoseEnvName"], xsize, ysize, zsize, model)
    drugEnv = NumericalGrid3D(model.properties["envNames"]["drugEnvName"], xsize, ysize, zsize, model)

    hc = 0
    tc = 0
    tic = 0
    # Adding agents

    for x in range(xsize):
        for y in range(ysize):
            for z in range(zsize):
                if x % 2 == z % 2:
                    agent = HealthyCell(model)
                    hc += 1
                    agent.AddAgentToGrid(model.properties["envNames"]["agentEnvName"], (x,y,z), model)
                    model.schedule.agents.add(agent)
                else:
                    for _ in range(model.properties["initialAgentSetup"]["numEndothelialCells"]):
                        agent = TipCell(model)
                        tic += 1
                        agent.AddAgentToGrid(model.properties["envNames"]["agentEnvName"], (x,y,z), model)
                        model.schedule.agents.add(agent)



    print("hc %s tc %s tic %s" % (hc, tc, tic))

    for _ in range(model.properties["initialAgentSetup"]["numCancerCells"]):
        for x in range(8,12):
            for y in range(8,12):
                for z in range(8,12):
                    c = CancerCell(model)

                    state = ["G1", "S", "G2", "M"][randint(0,3)]
                    stateLength = model.properties["agents"]["baseCellCycleLength"][state]
                    progressInState = randint(0, stateLength-1)

                    c.currentState = state
                    c.progressInState = progressInState


                    c.AddAgentToGrid(model.properties["envNames"]["agentEnvName"], (x,y,z), model)
                    model.schedule.agents.add(c)

    # Adding diffusion helpers
    model.schedule.helpers.append(GlucoseDiffusionHelper(model))

    odh = OxygenDiffusionHelper(model)
    model.schedule.helpers.append(odh)
    model.schedule.helpers.append(VegfDiffusionHelper(model))
    #model.schedule.helpers.append(DrugDiffusionHelper(model))

    snapshotInterval = 10

    model.schedule.helpers.append(AgentCounter(model))
    model.schedule.helpers.append(CancerCellWatcher(model, distributionInterval=snapshotInterval))
    model.schedule.helpers.append(VegfStimulusWatcher(model))
    model.schedule.helpers.append(TumourVolumeWatcher(model))
    model.schedule.helpers.append(OxygenConcentrationWatcher(model, interval=snapshotInterval))
    model.schedule.helpers.append(GlucoseConcentrationWatcher(model, interval=snapshotInterval))
    model.schedule.helpers.append(DeathCauseWatcher(model, interval=snapshotInterval))
    #model.schedule.helpers.append(OxygenSourceSinkWatcher(model, odh))

    model.schedule.helpers.append(ModelPicklerLite(properties["pickles"]["targetDir"], pickleEvery=1000, prefix=properties["pickles"]["prefix"]))

    def numAgentsExitCondition(model):
        return len([a for a in model.schedule.agents if a.__class__.__name__ == "CancerCell"]) > 400000

    def noCancerCells(model):
        return len([c for c in model.schedule.agents if c.__class__.__name__ == "CancerCell" and not c.dead]) == 0

    model.schedule.helpers.append(ExitConditionWatcher([numAgentsExitCondition, noCancerCells]))

    return model
