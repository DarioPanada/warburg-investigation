import numpy as np
import unittest
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
        for h in np.arange(2, 15, delta):
            a.current_hif_rate = h
            a._update_vegf_secretion_rate()
            vegf_rate = a.current_vegf_secretion_rate

            a.current_hif_rate = h + delta
            a._update_vegf_secretion_rate()
            self.assertGreater(a.current_vegf_secretion_rate, vegf_rate)
            self.assertLessEqual(0, a.current_vegf_secretion_rate)
            self.assertGreaterEqual(a.current_vegf_secretion_rate, vegf_rate)

        # Probability of progressing into synthesis should increase as hif
        # concentration increases
        delta = 0.1
        for h in np.arange(2, 15, delta):
            a.current_hif_rate = h
            a._update_p_synthesis()
            pSyntheis = a.currentPSynthesis

            a.current_hif_rate = h + delta
            a._update_p_synthesis()
            self.assertGreaterEqual(a.currentPSynthesis, pSyntheis)
            self.assertLessEqual(0, a.currentPSynthesis)
            self.assertGreaterEqual(1, a.currentPSynthesis)

        # Metabolic rate should decreases as hif concentration increase
        delta = 0.1
        for h in np.arange(2, 15, delta):
            a.current_hif_rate = h
            a._update_metabolic_rate()
            mr = a.current_metabolic_rate
            a.current_hif_rate = h + delta
            a._update_metabolic_rate()
            self.assertGreaterEqual(mr, a.current_metabolic_rate)
            self.assertLessEqual(0, a.current_metabolic_rate)
            self.assertGreaterEqual(1, a.current_metabolic_rate)

        # Expectation is that HIF expression rate increases as oxygen
        # decreases, until again falling sharply
        for o in np.arange(0.2, 1, 0.1):
            self.assertEquals(a._calculate_hif_expression_rate_from_oxygen(o),
                              1)

        delta = 0.1
        # HIF expression increases between 0.2 and ~0.02
        for o in np.arange(0.02, 0.2, delta):
            self.assertGreaterEqual(
                a._calculate_hif_expression_rate_from_oxygen(o),
                a._calculate_hif_expression_rate_from_oxygen(o + delta))

        # HIF expression decreases from ~0.02
        for o in np.arange(0, 0.02, delta):
            self.assertGreaterEqual(
                a._calculate_hif_expression_rate_from_oxygen(o + delta),
                a._calculate_hif_expression_rate_from_oxygen(o))

        oxygen_env_name = "oxygenEnv"
        agent_env_name = "agentEnv"
        model.properties["envNames"]["oxygenEnvName"] = oxygen_env_name
        model.properties["envNames"]["agentEnvName"] = agent_env_name
        env_oxygen = NumericalGrid3D(oxygen_env_name, 10, 10, 10, model)
        env_agents = ObjectGrid3D(agent_env_name, 10, 10, 10, model)
        env_oxygen.grid[(5, 5, 5)] = 1
        a.add_agent_to_grid(agent_env_name, (5, 5, 5), model)

        # Checking cancer cells die in low oxygen concentration
        a.step_main(model)
        self.assertEqual(0, [a for a in env_agents.grid[(5, 5, 5)] if
                             a.dead].__len__())

        self.assertEqual((5, 5, 5), a.environment_positions[agent_env_name])

        env_oxygen.grid[(5, 5, 5)] = 0.00001
        a.step_main(model)

        self.assertEqual(1, [a for a in env_agents.grid[(5, 5, 5)] if
                             a.dead].__len__())
