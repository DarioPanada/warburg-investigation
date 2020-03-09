from panaxea.core.Steppables import Helper


class ExitConditionWatcher(Helper, object):
    """
    At each epoch, checks whether any exit condition is met and if so sets
    the "exit" flag in the model to true.

    Attributes
    ----------
    exit_conditions : list
        A list of functions, each of which should take the current model
        instance and return a boolean.

        These should be functions, *NOT* lambdas, as pickle has issues with
        lambdas!
    """

    def __init__(self, exit_conditions):
        self.exit_conditions = exit_conditions

    def step_prologue(self, model):
        # It is important this logic is in the prologue and picklers run in
        # the epilogue. Pickler will always run
        # if the exit flag is set to true, so as to capture the end-state of
        # the model.

        for condition in self.exit_conditions:
            if condition(model):
                model.exit = True
