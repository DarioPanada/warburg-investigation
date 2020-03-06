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


def generate_properties(p):
    """
    Given a dictionary of parameter values, converts these to a dictionary
    object with with the structure expected by our model implementation,
    including derived properties such as coefficients for functions
    describing HIF-mediated properties and HIF expression.
    
    Parameters
    ----------
    p : dict
        The properties dictionary where key/value paris should be as obtained
        from the csv file

    Returns
    -------
    dict
        The properties dictionary
    """

    properties = dict()

    initialAgentSetup = dict()

    initialAgentSetup["numCancerCells"] = p["numCancerCells"]
    initialAgentSetup["numEndothelialCells"] = p["numEndothelialCells"]

    properties["initialAgentSetup"] = initialAgentSetup

    # Environment properties, including names, etc.
    envNames = dict()
    envNames["agentEnvName"] = "agentEnv"
    envNames["oxygenEnvName"] = "oxygenEnv"
    envNames["vegfEnvName"] = "vegfEnv"
    envNames["glucoseEnvName"] = "glucoseEnv"
    envNames["drugEnvName"] = "drugEnv"

    properties["envNames"] = envNames

    properties["envSize"] = p["envSize"]

    properties["maxAgentDensity"] = p["maxAgentDensity"]

    # Properties of agents
    agents = dict()
    agents["baseCellCycleLength"] = {
        "G1": p["G1"],
        "S": p["S"],
        "G2": p["G2"],
        "M": p["M"]
    }

    cancerCells = dict()

    cancerCells["HIFRange"] = [0.0, p["maxHIF"]]

    cancerCells["domains"] = {
        "ultraHypoxic": p["ultraHypoxicThreshold"],
        "warburgHypoxic": p["enhancedHypoxicThreshold"],
        "hypoxic": p["hypoxicThreshold"]
    }

    ohrg = OxygenHIFRelationsGenerator(minHIF=p["minHIF"], maxHIF=p["maxHIF"],
                                       ultraHypoxiaThreshold=
                                       p["ultraHypoxicThreshold"],
                                       hypoxiaThreshold=p["hypoxicThreshold"],
                                       enhancedHypoxicThreshold=
                                       p["enhancedHypoxicThreshold"],
                                       baseOxygenMetabolicRate=
                                       p["baseOxygenMetabolicRate"],
                                       minPSynthesis=p["minPSynthesis"])

    ultraHypoxiaCoeffs, hypoxiaCoeffs = ohrg.getOxygenToHif()
    warburgHypoxicCoeffs = ohrg.getOxygenToHifWarburg()
    cancerCells["pWarburgSwitch"] = p["pWarburgSwitch"]
    cancerCells["baseHifRate"] = p["baseHifRate"]
    cancerCells["minGlucoseUptakeRate"] = p["minGlucoseUptakeRate"]
    cancerCells["maxGlucoseUptakeRate"] = p["maxGlucoseUptakeRate"]
    cancerCells["minGlucoseWarburg"] = p["minGlucoseWarburg"]
    cancerCells["minGlucoseNonWarburg"] = p["minGlucoseNonWarburg"]
    cancerCells["minHIF"] = p["minHIF"]

    cancerCells["oxygenToHifCoeffs"] = {
        "hypoxic": hypoxiaCoeffs,
        "warburg": warburgHypoxicCoeffs,
        "ultraHypoxic": ultraHypoxiaCoeffs
    }

    cancerCells["hifToMetabolicRateCoeffs"] = ohrg.getHifToMetabolicRate()

    # Minimum probability of progressing into synthesis
    cancerCells["minPSynthesis"] = p["minPSynthesis"]
    cancerCells["hifToProliferationRateCoeffs"] = ohrg.getHifToPSynthesis()

    cancerCells["hifToVegfSecretionRateCoeffs"] = ohrg.getHifToVegf()

    # Minimum oxygen concentration for survival
    cancerCells["minimumOxygenConcentration"] = p["minimumOxygenConcentration"]
    cancerCells["maxVegfSecretionRate"] = 10

    agents["cancerCells"] = cancerCells

    agents["healthyTissues"] = {
        "oxygenUptakeRate": p["healthyTissueOxygenUptakeRate"]
    }

    endothelialCells = dict()

    agents["endothelialCells"] = endothelialCells

    # Minimum vegf concentration to sprout angiogenesis
    endothelialCells["minimumVegfConcentration"] = \
        p["minimumVegfConcentration"]
    # Minimum age of endothelial cells for sprouting
    endothelialCells["divisionDelay"] = p["endothelialDivisionDelay"]
    # Measured in mmHg
    endothelialCells["baseOxygenEmissionRate"] = p["baseOxygenEmissionRate"]

    endothelialCells["glucoseSecretionRate"] = p["glucoseSecretionRate"]

    properties["agents"] = agents

    diffusion = dict()

    diffusion["oxygenDiffusivity"] = p["oxygenDiffusivity"]
    diffusion["vegfDiffusivity"] = p["vegfDiffusivity"]
    diffusion["glucoseDiffusivity"] = p["glucoseDiffusivity"]
    diffusion["dt"] = p["dt"]
    diffusion["diffusionSolveIterations"] = p["diffusionSolveIterations"]

    properties["diffusion"] = diffusion

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
