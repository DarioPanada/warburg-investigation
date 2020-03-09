import unittest
from core.Model import Model
from models.alpha04c.helpers.ExitConditionWatcher import ExitConditionWatcher
from panaxea.core.Steppables import Agent


class TestExitCondition(unittest.TestCase):

    def test_exit_conditions(self):
        class TickerAgent(Agent, object):

            def __init__(self):
                self.counter = 0

            def stepPrologue(self, model):
                self.counter = self.counter + 1

        epochs = 20
        model = Model(epochs, verbose=False)
        model.schedule.agents.add(TickerAgent())

        def ticker_condition(model):
            expired_tickers = [a for a in model.schedule.agents if
                               a.counter >= 4]
            return len(expired_tickers) > 0

        model.schedule.helpers.append(ExitConditionWatcher([ticker_condition]))

        model.run()
        self.assertEqual(model.current_epoch, 4)


if __name__ == '__main__':
    unittest.main()
