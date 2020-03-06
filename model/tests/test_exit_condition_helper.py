import unittest

from core.Model import Model
from core.Steppables import Agent
from models.alpha04c.helpers.ExitConditionWatcher import ExitConditionWatcher

class testExitCondition(unittest.TestCase):

    def test_exit_conditions(self):

        class TickerAgent(Agent, object):

            def __init__(self):
                self.counter = 0

            def stepPrologue(self, model):
                self.counter = self.counter + 1

        epochs = 20
        model = Model(epochs, verbose=False)
        model.schedule.agents.add(TickerAgent())

        def tickerCondition(model):
            expiredTickers = [a for a in model.schedule.agents if a.counter >=4]
            return len(expiredTickers) > 0

        model.schedule.helpers.append(ExitConditionWatcher([tickerCondition]))

        model.run()
        self.assertEqual(model.currentEpoch, 4)

if __name__ == '__main__':
    unittest.main()