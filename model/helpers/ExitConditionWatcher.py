from panaxea.core.Steppables import Helper

'''
Helper checks if any custom exit condition for the simulation is true and, 
if so, sets the 'exit' flag of the model
to true
'''


class ExitConditionWatcher(Helper, object):
    '''
    :param exitConditions - List of functions that, taken an instance of the
    model, return true or false.
    These should be functions, *NOT* lambdas, as pickle has issues with
    lambdas!
    '''

    def __init__(self, exitConditions):
        self.exitConditions = exitConditions

    def step_prologue(self, model):
        # It is important this logic is in the prologue and picklers run in
        # the epilogue. Pickler will always run
        # if the exit flag is set to true, so as to capture the end-state of
        # the model.

        for condition in self.exitConditions:
            if condition(model):
                model.exit = True
