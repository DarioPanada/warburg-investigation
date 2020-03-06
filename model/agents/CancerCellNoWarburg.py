from core.Steppables import Agent
import random
from numpy.polynomial import Polynomial


class CancerCellNoWarburg(Agent, object):
    def __init__(self, model,warburgSwitch=False):
        super(CancerCellNoWarburg, self).__init__()

        self.currentState = "G1"
        self.progressInState = 0
        self.cellCycleLength = model.properties["agents"]["baseCellCycleLength"]
        self.cellCycleOrder = ["G1", "S", "G2", "M"]
        self.agentEnvName = model.properties["envNames"]["agentEnvName"]

        self.oxygenEnvName = model.properties["envNames"]["oxygenEnvName"]

        cancerCellProps = model.properties["agents"]["cancerCells"]

        self.baseHifRate = cancerCellProps["baseHifRate"]
        self.oxygenToHifCoeffsHypoxic = cancerCellProps["oxygenToHifCoeffs"]["hypoxic"]
        self.oxygenToHifCoeffsUltraHypoxic = cancerCellProps["oxygenToHifCoeffs"]["ultraHypoxic"]
        self.oxygenToHifCoeffsHypoxicWarburg = cancerCellProps["oxygenToHifCoeffs"]["warburg"]
        self.oxygenHypoxicDomain = cancerCellProps["domains"]["hypoxic"]
        self.oxygenUltraHypoxicDomain = cancerCellProps["domains"]["ultraHypoxic"]
        self.oxygenWarburgHypoxicDomain = cancerCellProps["domains"]["warburgHypoxic"]

        self.currentHifRate = self.baseHifRate

        self.hifToMetabolicRateCoeffs = cancerCellProps["hifToMetabolicRateCoeffs"]
        self.currentMetabolicRate = 1

        self.hifToProliferationRateCoeffs = cancerCellProps["hifToProliferationRateCoeffs"]
        self.minPSynthesis = cancerCellProps["minPSynthesis"]
        self.currentPSynthesis = self.minPSynthesis

        self.hifToVegfSecretionRateCoeffs = cancerCellProps["hifToVegfSecretionRateCoeffs"]
        self.currentVegfSecretionRate = 1

        self.minimumOxygenConcentration = cancerCellProps["minimumOxygenConcentration"]
        self.minGlucoseWarburg = cancerCellProps["minGlucoseWarburg"]

        self.HIFRange = cancerCellProps["HIFRange"]

        self.minGlucoseUptakeRate = cancerCellProps["minGlucoseUptakeRate"]
        self.maxGlucoseUptakeRate = cancerCellProps["maxGlucoseUptakeRate"]
        self.glucoseUptakeRate = self.minGlucoseUptakeRate

        self.pWarburgSwitch = cancerCellProps["pWarburgSwitch"]
        self.warburgSwitch = warburgSwitch

        self.dead = False


    def progressCell_(self, model):
        # time to divide
        if self.currentState == self.cellCycleOrder[-1] and self.progressInState == self.cellCycleLength[
            self.currentState]:
            currentPos = self.environmentPositions[self.agentEnvName]
            targetPos = None

            if len(model.environments[self.agentEnvName].grid[(currentPos[0], currentPos[1], currentPos[2])]) < \
                    model.properties["maxAgentDensity"]:
                targetPos = currentPos
            else:
                mooreTarget = model.environments[self.agentEnvName].getLeastPopulatedMooreNeigh(currentPos)
                if len(model.environments[self.agentEnvName].grid[(mooreTarget[0], mooreTarget[1], mooreTarget[2])]) < \
                        model.properties["maxAgentDensity"]:
                    targetPos = mooreTarget

            if targetPos is not None:
                # Creating new cancer cell and adding it at current position in the designated environment
                c = CancerCellNoWarburg(model, warburgSwitch=self.warburgSwitch)
                c.AddAgentToGrid(self.agentEnvName, targetPos, model)
                model.schedule.agentsToSchedule.add(c)

                # Resetting current cell
                self.progressInState = 0
                self.currentState = self.cellCycleOrder[0]
        # progress in cell life-cyle
        elif self.progressInState == self.cellCycleLength[self.currentState]:
            if self.currentState == "G1" and random.random() > self.currentPSynthesis:
                pass
            else:
                self.progressInState = 0
                # Get index of current state and set next state
                currentIndex = self.cellCycleOrder.index(self.currentState)
                self.currentState = self.cellCycleOrder[currentIndex + 1]
        else:
            self.progressInState = self.progressInState + 1

    def stepMain(self, model):
        currentPos = self.environmentPositions[self.agentEnvName]
        self.oxygenAtPos = model.environments[self.oxygenEnvName].grid[currentPos]

        if self.warburgSwitch:
            exit("Cell has warburg phenotype on, this shouldn't be possible in this agent type!")

        if not self.dead:
            # Oxygen is a limiting factor is the warburg switch is turned on, otherwise glucose is, glucose
            # as limiting factor to be implemented
            if self.oxygenAtPos > self.minimumOxygenConcentration and not self.warburgSwitch:

                self.updateHifAndMediated_(model)
                self.progressCell_(model)
            else:
                self.dead = True

    def updateHifAndMediated_(self, model):
        self.updateHifExpressionRate_(model)
        self.updateMetabolicRate_()
        self.updatePSynthesis_()
        self.updateVegfSecretionRate_()

    def updatePSynthesis_(self):

        p = Polynomial(coef=self.hifToProliferationRateCoeffs, domain=self.HIFRange)

        pSynthesis = p(self.currentHifRate)

        self.currentPSynthesis = pSynthesis

    def updateVegfSecretionRate_(self):

        p = Polynomial(coef=self.hifToVegfSecretionRateCoeffs, domain=self.HIFRange)

        vegfRate = p(self.currentHifRate)

        self.currentVegfSecretionRate = 10*max(0, min(1, vegfRate))

    def updateMetabolicRate_(self):

        p = Polynomial(coef=self.hifToMetabolicRateCoeffs, domain=self.HIFRange)

        metabolicRate = p(self.currentHifRate)

        self.currentMetabolicRate = metabolicRate


    def updateHifExpressionRate_(self, model):
        newRate = self.calculateHifExpressionRateFromOxygen_(self.oxygenAtPos)
        self.currentHifRate = self.baseHifRate * newRate

    def calculateHifExpressionRateFromOxygen_(self, oxygenAtPos):

        if oxygenAtPos > self.oxygenHypoxicDomain:
            return 1

        if oxygenAtPos > self.oxygenUltraHypoxicDomain:
            coef = self.oxygenToHifCoeffsHypoxic
            domain = [self.oxygenUltraHypoxicDomain,self.oxygenHypoxicDomain]
        else:
            coef = self.oxygenToHifCoeffsUltraHypoxic
            domain = [0.0, self.oxygenUltraHypoxicDomain]

        p = Polynomial(coef=coef, domain=domain)
        return p(oxygenAtPos)

