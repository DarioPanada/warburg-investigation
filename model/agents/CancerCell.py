import random
from core.Steppables import Agent
from numpy.polynomial import Polynomial


class CancerCell(Agent, object):
    def __init__(self, model, warburgSwitch=False):
        super(CancerCell, self).__init__()

        self.currentState = "G1"
        self.progressInState = 0
        self.cellCycleLength = model.properties["agents"][
            "baseCellCycleLength"]
        self.cellCycleOrder = ["G1", "S", "G2", "M"]
        self.agentEnvName = model.properties["envNames"]["agentEnvName"]
        self.oxygenEnvName = model.properties["envNames"]["oxygenEnvName"]
        self.glucoseEnvName = model.properties["envNames"]["glucoseEnvName"]
        self.drugEnvName = model.properties["envNames"]["drugEnvName"]

        cancerCellProps = model.properties["agents"]["cancerCells"]

        self.baseHifRate = cancerCellProps["baseHifRate"]
        self.oxygenToHifCoeffsHypoxic = cancerCellProps["oxygenToHifCoeffs"][
            "hypoxic"]
        self.oxygenToHifCoeffsUltraHypoxic = \
        cancerCellProps["oxygenToHifCoeffs"]["ultraHypoxic"]
        self.oxygenToHifCoeffsHypoxicWarburg = \
        cancerCellProps["oxygenToHifCoeffs"]["warburg"]
        self.oxygenHypoxicDomain = cancerCellProps["domains"]["hypoxic"]
        self.oxygenUltraHypoxicDomain = cancerCellProps["domains"][
            "ultraHypoxic"]
        self.oxygenWarburgHypoxicDomain = cancerCellProps["domains"][
            "warburgHypoxic"]

        self.currentHifRate = cancerCellProps["minHIF"]

        self.hifToMetabolicRateCoeffs = cancerCellProps[
            "hifToMetabolicRateCoeffs"]
        self.currentMetabolicRate = 1

        self.hifToProliferationRateCoeffs = cancerCellProps[
            "hifToProliferationRateCoeffs"]
        self.minPSynthesis = cancerCellProps["minPSynthesis"]
        self.currentPSynthesis = self.minPSynthesis

        self.hifToVegfSecretionRateCoeffs = cancerCellProps[
            "hifToVegfSecretionRateCoeffs"]
        self.currentVegfSecretionRate = 1

        self.minimumOxygenConcentration = cancerCellProps[
            "minimumOxygenConcentration"]
        self.minGlucoseWarburg = cancerCellProps["minGlucoseWarburg"]
        self.minGlucoseNonWarburg = cancerCellProps["minGlucoseNonWarburg"]

        self.HIFRange = cancerCellProps["HIFRange"]

        self.minGlucoseUptakeRate = cancerCellProps["minGlucoseUptakeRate"]
        self.maxGlucoseUptakeRate = cancerCellProps["maxGlucoseUptakeRate"]
        self.glucoseUptakeRate = self.minGlucoseUptakeRate

        self.pWarburgSwitch = cancerCellProps["pWarburgSwitch"]
        self.warburgSwitch = warburgSwitch

        self.dead = False
        self.quiescent = False

        self.drugAtPos = 0

        # To allow for different drug effects to be tested, this will be
        # loaded as part of the config.
        self.maxVegf = cancerCellProps["maxVegfSecretionRate"]
        self.age = 0
        self.minHIF = cancerCellProps["minHIF"]

    def progressCell_(self, model):
        # time to divide
        if self.currentState == self.cellCycleOrder[
            -1] and self.progressInState == self.cellCycleLength[
            self.currentState]:
            currentPos = self.environmentPositions[self.agentEnvName]
            targetPos = None

            if len(model.environments[self.agentEnvName].grid[
                       (currentPos[0], currentPos[1], currentPos[2])]) < \
                    model.properties["maxAgentDensity"]:
                targetPos = currentPos
            else:
                mooreTarget = model.environments[
                    self.agentEnvName].getLeastPopulatedMooreNeigh(currentPos)
                if len(model.environments[self.agentEnvName].grid[
                           (mooreTarget[0], mooreTarget[1], mooreTarget[
                               2])]) < \
                        model.properties["maxAgentDensity"]:
                    targetPos = mooreTarget

            if targetPos is not None:
                # Creating new cancer cell and adding it at current position
                # in the designated environment
                c = CancerCell(model)
                c.AddAgentToGrid(self.agentEnvName, targetPos, model)
                model.schedule.agentsToSchedule.add(c)

                # Resetting current cell
                self.progressInState = 0
                self.currentState = self.cellCycleOrder[0]
                self.warburgSwitch = False
                self.age = 0

        # progress in cell life-cyle
        elif self.progressInState == self.cellCycleLength[self.currentState]:
            if self.currentState == "G1" and random.random() > \
                    self.currentPSynthesis:
                pass
            else:
                self.progressInState = 0
                # Get index of current state and set next state
                currentIndex = self.cellCycleOrder.index(self.currentState)
                self.currentState = self.cellCycleOrder[currentIndex + 1]
        else:
            self.progressInState = self.progressInState + 1

    def step_main(self, model):
        currentPos = self.environmentPositions[self.agentEnvName]
        self.oxygenAtPos = model.environments[self.oxygenEnvName].grid[
            currentPos]
        self.glucoseAtPos = model.environments[self.glucoseEnvName].grid[
            currentPos]
        self.drugAtPos = model.environments[self.drugEnvName].grid[currentPos]

        '''# If a cell is quiescent and the drug dropped, we have it return 
        active
        if self.quiescent and self.drugAtPos < self.minDrugDosage:
            self.quiescent = False'''

        # No logic executed on dead or quiescent cells
        if not self.dead and not self.quiescent:
            self.age = self.age + 1

            if random.random() < self.pWarburgSwitch and not \
                    self.warburgSwitch:
                self.warburgSwitch = True
                self.glucoseUptakeRate = self.maxGlucoseUptakeRate
            # If the cell is not warburg, then in order for it to survive it
            # must:
            # Have access to minimum oxygen concentration
            # Have access to minimum glucose for non warburg cells
            # If the cell is a warburg cell then, in order for it to
            # survive, it must:
            # Have access to minimum glucose for warburg cells
            if not self.decideDie_():

                # If a cell is a warburg cell and the dosage at position is
                # greater than the minimum threshold,
                # it becomes quiescent. (We return so that no further cell
                # logic is executed)

                # TEMPORARILY DISABLING DRUGS, THIS SHOULD BE TURNED INTO A
                # FLAG THAT CAN BE SET FROM PARAMS

                '''if self.drugAtPos >= self.minDrugDosage and 
                self.warburgSwitch:
                    self.quiescent = True
                    return'''

                self.reactToDrug_(model, currentPos)
                self.updateHifAndMediated_(model)
                self.progressCell_(model)
            else:
                # A dead cell is not quiescent
                self.dead = True
                self.quiescent = False

    def decideDie_(self):
        if self.warburgSwitch and self.glucoseAtPos < self.minGlucoseWarburg:
            self.causeOfDeath = {
                "cause": "Lack of glucose",
                "glucoseAtPos": self.glucoseAtPos,
                "warburg": self.warburgSwitch,
                "age": self.age
            }
            return True
        elif not self.warburgSwitch:
            if self.oxygenAtPos < self.minimumOxygenConcentration:
                self.causeOfDeath = {
                    "cause": "Lack of oxygen",
                    "oxygenAtPos": self.oxygenAtPos,
                    "warburg": self.warburgSwitch,
                    "age": self.age
                }
                return True
            elif self.glucoseAtPos < self.minGlucoseNonWarburg:
                self.causeOfDeath = {
                    "cause": "Lack of glucose",
                    "glucoseAtPos": self.glucoseAtPos,
                    "warburg": self.warburgSwitch,
                    "age": self.age
                }
                return True

        return False

    def updateHifAndMediated_(self, model):
        self.updateHifExpressionRate_(model)
        self.updateMetabolicRate_()
        self.updatePSynthesis_()
        self.updateVegfSecretionRate_()

    def updatePSynthesis_(self):

        p = Polynomial(coef=self.hifToProliferationRateCoeffs,
                       domain=self.HIFRange)

        pSynthesis = p(self.currentHifRate)

        self.currentPSynthesis = pSynthesis

    def updateVegfSecretionRate_(self):

        p = Polynomial(coef=self.hifToVegfSecretionRateCoeffs,
                       domain=self.HIFRange)

        vegfRate = p(self.currentHifRate)

        self.currentVegfSecretionRate = self.maxVegf * max(0, min(1, vegfRate))

    def updateMetabolicRate_(self):

        p = Polynomial(coef=self.hifToMetabolicRateCoeffs,
                       domain=self.HIFRange)

        metabolicRate = p(self.currentHifRate)

        self.currentMetabolicRate = metabolicRate

    def updateHifExpressionRate_(self, model):

        if not self.warburgSwitch:
            newRate = self.calculateHifExpressionRateFromOxygen_(
                self.oxygenAtPos)
        else:
            newRate = self.calculateHifExpressionRateFromOxygenWarburg_(
                self.oxygenAtPos)

        maxShiftPerEpoch = 0.1

        currentRate = self.currentHifRate

        if currentRate > newRate:
            newRate = max(0, currentRate - min(maxShiftPerEpoch,
                                               currentRate - newRate))
        else:
            newRate = min(
                currentRate + min(maxShiftPerEpoch, newRate - currentRate), 16)

        self.currentHifRate = self.baseHifRate * newRate

    def calculateHifExpressionRateFromOxygen_(self, oxygenAtPos):

        if oxygenAtPos > self.oxygenHypoxicDomain:
            return 1

        if oxygenAtPos > self.oxygenUltraHypoxicDomain:
            coef = self.oxygenToHifCoeffsHypoxic
            domain = [self.oxygenUltraHypoxicDomain, self.oxygenHypoxicDomain]
        else:
            coef = self.oxygenToHifCoeffsUltraHypoxic
            domain = [0.0, self.oxygenUltraHypoxicDomain]

        p = Polynomial(coef=coef, domain=domain)
        return p(oxygenAtPos)

    def calculateHifExpressionRateFromOxygenWarburg_(self, oxygenAtPos):
        if oxygenAtPos > self.oxygenWarburgHypoxicDomain:
            return self.minHIF

        if oxygenAtPos > self.oxygenUltraHypoxicDomain:
            coef = self.oxygenToHifCoeffsHypoxicWarburg
            domain = [self.oxygenUltraHypoxicDomain,
                      self.oxygenWarburgHypoxicDomain]
        else:
            coef = self.oxygenToHifCoeffsUltraHypoxic
            domain = [0.0, self.oxygenUltraHypoxicDomain]

        p = Polynomial(coef=coef, domain=domain)
        return p(oxygenAtPos)
