import unittest
import numpy as np
from core.Environment import NumericalGrid3D, ObjectGrid3D
from core.Model import Model
from models.alpha04c.agents.CancerCell import CancerCell
from models.alpha04c.models.model_070818 import generateProperties


class testCancerCell(unittest.TestCase):

    def test_cancer_cell(self):
        model = Model(5, verbose=True)
        properties = generateProperties()
        model.properties = properties

        a = CancerCell(model)

        # VEGF secretion rate should increase as HIF concentration increases
        delta = 0.1
        for h in np.arange(2,15, delta):
            a.currentHifRate = h
            a.__update_vegf_secretion_rate()
            vegfRate = a.currentVegfSecretionRate

            a.currentHifRate = h+delta
            a.__update_vegf_secretion_rate()
            self.assertGreater(a.currentVegfSecretionRate, vegfRate)
            self.assertLessEqual(0, a.currentVegfSecretionRate)
            self.assertGreaterEqual(a.currentVegfSecretionRate, vegfRate)

        # Probability of progressing into synthesis should increase as hif concentration increases
        delta = 0.1
        for h in np.arange(2,15, delta):
            a.currentHifRate = h
            a.__update_p_synthesis()
            pSyntheis = a.currentPSynthesis

            a.currentHifRate = h+delta
            a.__update_p_synthesis()
            self.assertGreaterEqual(a.currentPSynthesis, pSyntheis)
            self.assertLessEqual(0, a.currentPSynthesis)
            self.assertGreaterEqual(1, a.currentPSynthesis)

        # Metabolic rate should decreases as hif concentration increase
        delta = 0.1
        for h in np.arange(2,15, delta):
            a.currentHifRate = h
            a.__update_metabolic_rate()
            mr = a.currentMetabolicRate
            a.currentHifRate = h+delta
            a.__update_metabolic_rate()
            self.assertGreaterEqual(mr, a.currentMetabolicRate)
            self.assertLessEqual(0, a.currentMetabolicRate)
            self.assertGreaterEqual(1, a.currentMetabolicRate)

        # Expectation is that HIF expression rate increases as oxygen decreases, until again falling sharply
        for o in np.arange(0.2,1, 0.1):
            self.assertEquals(a.__calculate_hif_expression_rate_from_oxygen(o), 1)

        delta = 0.1
        # HIF expression increases between 0.2 and ~0.02
        for o in np.arange(0.02, 0.2, delta):
            self.assertGreaterEqual(a.__calculate_hif_expression_rate_from_oxygen(o), a.__calculate_hif_expression_rate_from_oxygen(o + delta))

        # HIF expression decreases from ~0.02
        for o in np.arange(0,0.02, delta):
            self.assertGreaterEqual(a.__calculate_hif_expression_rate_from_oxygen(o + delta), a.__calculate_hif_expression_rate_from_oxygen(o))

        oxygenEnvName = "oxygenEnv"
        agentEnvName = "agentEnv"
        model.properties["envNames"]["oxygenEnvName"] = oxygenEnvName
        model.properties["envNames"]["agentEnvName"] = agentEnvName
        envOxygen = NumericalGrid3D(oxygenEnvName, 10,10,10, model)
        envAgents = ObjectGrid3D(agentEnvName, 10,10,10, model)
        envOxygen.grid[(5,5,5)] = 1
        a.add_agent_to_grid(agentEnvName, (5,5,5), model)

        # Checking cancer cells die in low oxygen concentration
        a.step_main(model)
        self.assertEqual(0, [a for a in envAgents.grid[(5,5,5)] if a.dead].__len__())

        self.assertEqual((5,5,5), a.environment_positions[agentEnvName])

        envOxygen.grid[(5,5,5)] = 0.00001
        a.step_main(model)


        self.assertEqual(1, [a for a in envAgents.grid[(5,5,5)] if a.dead].__len__())
