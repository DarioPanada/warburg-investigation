from panaxea.core.Environment import ObjectGrid3D, NumericalGrid3D
from panaxea.core.Model import Model
from panaxea.toolkit.Toolkit import ModelPicklerLite
from random import randint

from model.agents.CancerCell import CancerCell
from model.agents.EndothelialCell import TipCell
from model.agents.HealthyCell import HealthyCell
from model.helpers.AgentCounter import AgentCounter
from model.helpers.CancerCellWatcher import CancerCellWatcher
from model.helpers.DeathCauseWatcher import DeathCauseWatcher
from model.helpers.ExitConditionWatcher import ExitConditionWatcher
from model.helpers.GlucoseConcentrationWatcher import \
    GlucoseConcentrationWatcher
from model.helpers.GlucoseDiffusionHelper import GlucoseDiffusionHelper
from model.helpers.OxygenConcentrationWatcher import OxygenConcentrationWatcher
from model.helpers.OxygenDiffusionHelper import OxygenDiffusionHelper
from model.helpers.TumourVolumeWatcher import TumourVolumeWatcher
from model.helpers.VegfDiffusionHelper import VegfDiffusionHelper
from model.helpers.VegfStimulusWatcher import VegfStimulusWatcher
from model.utils.OxygenHIFRelationsGenerator import OxygenHIFRelationsGenerator


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

    initial_agent_setup = dict()

    initial_agent_setup["numCancerCells"] = p["numCancerCells"]
    initial_agent_setup["numEndothelialCells"] = p["numEndothelialCells"]

    properties["initialAgentSetup"] = initial_agent_setup

    # Environment properties, including names, etc.
    env_names = dict()
    env_names["agentEnvName"] = "agentEnv"
    env_names["oxygenEnvName"] = "oxygenEnv"
    env_names["vegfEnvName"] = "vegfEnv"
    env_names["glucoseEnvName"] = "glucoseEnv"
    env_names["drugEnvName"] = "drugEnv"

    properties["env_names"] = env_names

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

    cancer_cells = dict()

    cancer_cells["HIFRange"] = [0.0, p["maxHIF"]]

    cancer_cells["domains"] = {
        "ultraHypoxic": p["ultraHypoxicThreshold"],
        "warburgHypoxic": p["enhancedHypoxicThreshold"],
        "hypoxic": p["hypoxicThreshold"]
    }

    ohrg = OxygenHIFRelationsGenerator(
        minHIF=p["minHIF"],
        maxHIF=p["maxHIF"],
        ultraHypoxiaThreshold=p["ultraHypoxicThreshold"],
        hypoxiaThreshold=p["hypoxicThreshold"],
        enhancedHypoxicThreshold=p["enhancedHypoxicThreshold"],
        baseOxygenMetabolicRate=p["baseOxygenMetabolicRate"],
        minPSynthesis=p["minPSynthesis"])

    ultra_hypoxia_coeffs, hypoxia_coeffs = ohrg.getOxygenToHif()
    warburg_hypoxic_coeffs = ohrg.getOxygenToHifWarburg()
    cancer_cells["pWarburgSwitch"] = p["pWarburgSwitch"]
    cancer_cells["baseHifRate"] = p["baseHifRate"]
    cancer_cells["minGlucoseUptakeRate"] = p["minGlucoseUptakeRate"]
    cancer_cells["maxGlucoseUptakeRate"] = p["maxGlucoseUptakeRate"]
    cancer_cells["minGlucoseWarburg"] = p["minGlucoseWarburg"]
    cancer_cells["minGlucoseNonWarburg"] = p["minGlucoseNonWarburg"]
    cancer_cells["minHIF"] = p["minHIF"]

    cancer_cells["oxygenToHifCoeffs"] = {
        "hypoxic": hypoxia_coeffs,
        "warburg": warburg_hypoxic_coeffs,
        "ultraHypoxic": ultra_hypoxia_coeffs
    }

    cancer_cells["hifToMetabolicRateCoeffs"] = ohrg.getHifToMetabolicRate()

    # Minimum probability of progressing into synthesis
    cancer_cells["minPSynthesis"] = p["minPSynthesis"]
    cancer_cells["hifToProliferationRateCoeffs"] = ohrg.getHifToPSynthesis()

    cancer_cells["hifToVegfSecretionRateCoeffs"] = ohrg.getHifToVegf()

    # Minimum oxygen concentration for survival
    cancer_cells["minimumOxygenConcentration"] = p[
        "minimumOxygenConcentration"]
    cancer_cells["maxVegfSecretionRate"] = 10

    agents["cancer_cells"] = cancer_cells

    agents["healthyTissues"] = {
        "oxygenUptakeRate": p["healthyTissueOxygenUptakeRate"]
    }

    endothelial_cells = dict()

    agents["endothelial_cells"] = endothelial_cells

    # Minimum vegf concentration to sprout angiogenesis
    endothelial_cells["minimumVegfConcentration"] = \
        p["minimumVegfConcentration"]
    # Minimum age of endothelial cells for sprouting
    endothelial_cells["divisionDelay"] = p["endothelialDivisionDelay"]
    # Measured in mmHg
    endothelial_cells["baseOxygenEmissionRate"] = p["baseOxygenEmissionRate"]

    endothelial_cells["glucoseSecretionRate"] = p["glucoseSecretionRate"]

    properties["agents"] = agents

    diffusion = dict()

    diffusion["oxygenDiffusivity"] = p["oxygenDiffusivity"]
    diffusion["vegfDiffusivity"] = p["vegfDiffusivity"]
    diffusion["glucoseDiffusivity"] = p["glucoseDiffusivity"]
    diffusion["dt"] = p["dt"]
    diffusion["diffusionSolveIterations"] = p["diffusionSolveIterations"]

    properties["diffusion"] = diffusion

    return properties


def generate_model(properties, numEpochs):
    """
    Generates and sets up a model object for warburg investigation. Includes
    assigning properties, instantiating agents and environments, etc.

    Parameters
    ----------
    properties : dict
        A dictionary of model properties, as generated by the
        generate_properties function
    numEpochs : number
        Number of epochs the model should run for

    Returns
    -------
    Model
        A configured instance of the model object.
    """

    model = Model(numEpochs)

    model.properties = properties

    xsize = ysize = zsize = model.properties["envSize"]

    # Adding environments
    ObjectGrid3D(
        model.properties["envNames"]["agentEnvName"],
        xsize, ysize, zsize, model)
    NumericalGrid3D(
        model.properties["envNames"]["oxygenEnvName"],
        xsize,
        ysize,
        zsize,
        model)
    NumericalGrid3D(
        model.properties["envNames"]["vegfEnvName"],
        xsize,
        ysize,
        zsize,
        model)
    NumericalGrid3D(
        model.properties["envNames"]["glucoseEnvName"],
        xsize,
        ysize,
        zsize,
        model)
    NumericalGrid3D(
        model.properties["envNames"]["drugEnvName"],
        xsize,
        ysize,
        zsize,
        model)

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
                    agent.add_agent_to_grid(model.properties["envNames"][
                                                "agentEnvName"], (x, y, z),
                                            model)
                    model.schedule.agents.add(agent)
                else:
                    for _ in range(model.properties["initialAgentSetup"][
                                       "numEndothelialCells"]):
                        agent = TipCell(model)
                        tic += 1
                        agent.add_agent_to_grid(model.properties["envNames"][
                                                    "agentEnvName"], (x, y, z),
                                                model)
                        model.schedule.agents.add(agent)

    print("hc %s tc %s tic %s" % (hc, tc, tic))

    for _ in range(model.properties["initialAgentSetup"]["numCancerCells"]):
        for x in range(8, 12):
            for y in range(8, 12):
                for z in range(8, 12):
                    c = CancerCell(model)

                    state = ["G1", "S", "G2", "M"][randint(0, 3)]
                    state_length = \
                        model.properties["agents"]["baseCellCycleLength"][
                            state]
                    progress_in_state = randint(0, state_length - 1)

                    c.current_state = state
                    c.progress_in_state = progress_in_state

                    c.add_agent_to_grid(model.properties["envNames"][
                                            "agentEnvName"], (x, y, z), model)
                    model.schedule.agents.add(c)

    # Adding diffusion helpers
    model.schedule.helpers.append(GlucoseDiffusionHelper(model))

    odh = OxygenDiffusionHelper(model)
    model.schedule.helpers.append(odh)
    model.schedule.helpers.append(VegfDiffusionHelper(model))

    snapshot_interval = 10

    model.schedule.helpers.append(AgentCounter(model))
    model.schedule.helpers.append(
        CancerCellWatcher(model, distributionInterval=snapshot_interval))
    model.schedule.helpers.append(VegfStimulusWatcher(model))
    model.schedule.helpers.append(TumourVolumeWatcher(model))
    model.schedule.helpers.append(
        OxygenConcentrationWatcher(model, interval=snapshot_interval))
    model.schedule.helpers.append(
        GlucoseConcentrationWatcher(model, interval=snapshot_interval))
    model.schedule.helpers.append(
        DeathCauseWatcher(model, interval=snapshot_interval))
    model.schedule.helpers.append(ModelPicklerLite(model.properties["outDir"]))

    def num_agents_exit_condition(model):
        return len([a for a in model.schedule.agents if
                    a.__class__.__name__ == "CancerCell"]) > 400000

    def no_cancer_cells(model):
        return len([c for c in model.schedule.agents if
                    c.__class__.__name__ == "CancerCell" and not c.dead]) == 0

    model.schedule.helpers.append(
        ExitConditionWatcher([num_agents_exit_condition, no_cancer_cells]))

    return model
