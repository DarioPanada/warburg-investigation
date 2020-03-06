from numpy.polynomial import Polynomial
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

'''
This class allows to generate coefficients for function relating oxygen concentration to hif expression rates
and hif expression rates to vegf secretion rate, proliferation rate and metabolic rate.
'''


class OxygenHIFRelationsGenerator():
    '''
    :param minHIF - Minimum HIF expression rates
    :param maxHIF - Maximum HIF expression rates
    :param ultraHypoxiaThreshold - A value x (0 < x < 1) such that oxygen concentrations 0 to x are considered ultra
    hypoxia. (HIF expression rate at 0 will be 0 and at x will be maxHIF)
    :param hypoxiaThreshold - A value x (ultraHypoxiaThreshold < hypoxiaThreshold < 1) such that HIF expression rates
    at hypoxiaThreshold will 1. HIF Expression rate will decrease from a value maxHIF at ultraHypoxia threshold
    until reaching minHIF at hypoxiaThreshold
    :param enhancedHypoxicThreshold - Hypoxic threshold to be used during warburg metabolism
    :param baseOxygenMetabolicRate - Oxygen uptake rate in normoxia
    :param minPSynthesis - (Optional, defaults to 0), Lower value for probability of progressing into synthesis
    '''

    def __init__(self, minHIF, maxHIF, ultraHypoxiaThreshold, hypoxiaThreshold, enhancedHypoxicThreshold, baseOxygenMetabolicRate, minPSynthesis=0):
        self.minHIF = minHIF
        self.maxHIF = maxHIF
        self.ultraHypoxiaThreshold = ultraHypoxiaThreshold
        self.hypoxiaThreshold = hypoxiaThreshold
        self.minPSynthesis = minPSynthesis
        self.enhancedHypoxicThreshold = enhancedHypoxicThreshold
        self.baseOxygenMetabolicRate = baseOxygenMetabolicRate

    '''
    Returns coefficients for hif to metabolic rate. The domain of the function is [0, maxHIF].
    :param degree - (Optional, defaults to 2) If degree > 2 function has to be changed to include more points)
    :param render - (Optional, defaults to False) To be used for testing, if set to true plots the polynomials
    :return Coefficients for hif to metabolic rate
    '''
    def getHifToMetabolicRate(self, degree=2, render=False):
        points = [
            (0,self.baseOxygenMetabolicRate),
            (self.maxHIF,0)
        ]

        xs = [p[0] for p in points]
        ys = [p[1] for p in points]


        p = Polynomial.fit(xs,ys, 2)
        # Getting data to reproduce polynomials, use of sort:
        # p = Polynomial(coef=[list of coeffs], domain=[xstart, xend])

        if render:
            plt.plot(*p.linspace(), label="Fit", color="orange")
            plt.scatter(xs, ys, label="Raw Points")
            plt.legend()
            plt.xlabel("HIF Expression Rates")
            plt.ylabel("Metabolic Rate")
            plt.title("Metabolic Rates across HIF Expression Rates")
            plt.show()

        return [round(c,2) for c in list(p.coef)]

    '''
    Returns coefficients for hif to probability of synthesis. The domain of the function is [0, maxHIF].
    :param degree - (Optional, defaults to 2) If degree > 2 function has to be changed to include more points)
    :param render - (Optional, defaults to False) To be used for testing, if set to true plots the polynomials
    :return Coefficients for hif to p synthesis
    '''
    def getHifToPSynthesis(self, degree=1, render=False):
        points = [
            (0,self.minPSynthesis),
            (self.maxHIF,1)
        ]

        xs = [p[0] for p in points]
        ys = [p[1] for p in points]


        p = Polynomial.fit(xs,ys, degree)
        # Getting data to reproduce polynomials, use of sort:
        # p = Polynomial(coef=[list of coeffs], domain=[xstart, xend]

        if render:
            plt.plot(*p.linspace(), label="Fit", color="orange")
            plt.scatter(xs, ys, label="Raw Points")
            plt.legend()
            plt.xlabel("HIF Expression Rates")
            plt.ylabel("Probability Synthesis")
            plt.title("Probability of Synthesis across HIF Expression Rates")
            plt.show()

        return [round(c,2) for c in list(p.coef)]

    '''
    Returns coefficients for hif to vegf relation. The domain of the function is [0, maxHIF].
    :param degree - (Optional, defaults to 2) If degree > 2 function has to be changed to include more points)
    :param render - (Optional, defaults to False) To be used for testing, if set to true plots the polynomials
    :return Coefficients for hif to vegf
    '''
    def getHifToVegf(self, degree=2, render=False):
        points = [
            (0, 0),
            (self.maxHIF, 1)
        ]

        xs = [p[0] for p in points]
        ys = [p[1] for p in points]

        p = Polynomial.fit(xs, ys, degree)
        # Getting data to reproduce polynomials, use of sort:
        # p = Polynomial(coef=[list of coeffs], domain=[xstart, xend])

        if render:
            plt.plot(*p.linspace(), label="Fit", color="orange")
            plt.scatter(xs, ys, label="Raw Points")
            plt.legend()
            plt.xlabel("HIF Expression Rates")
            plt.ylabel("VEGF Secretion Rate")
            plt.title("VEGF Secretion Rates across HIF Expression Rates")
            plt.show()

        return [round(c, 2) for c in list(p.coef)]

    '''
    Returns coefficients for ultra-hypoxia and hypoxia polynomials to be used by cancer cells if the warburg
    switch has been activated. At present, only hypoxic threshold is changed whereas ultra-hypoxic is left as-is.
    :param degree - (Optional, defaults to 2) If degree > 2 function has to be changed to include more points)
    :param render - (Optional, defaults to False) To be used for testing, if set to true plots the polynomials
    :return ultraHypoxiaCoeffs - Coefficients of polynomial for ultra-hypoxia (as a list)
    :return hypoxiaCoeffs - Coefficients of polynomial for hypoxia (as a list)
    '''
    def getOxygenToHifWarburg(self, degree=1, render=False):
        points = [
            # Points to 0.02
            (self.ultraHypoxiaThreshold, self.maxHIF),
            (self.enhancedHypoxicThreshold, self.minHIF),
        ]

        xs = [p[0] for p in points]
        ys = [p[1] for p in points]

        # Fitting points between 0.02 and 0.2
        xs= [xs[0], xs[1]]
        ys= [ys[0], ys[1]]

        p = Polynomial.fit(xs, ys, degree)


        # Getting data to reproduce polynomials, use of sort:
        # p = Polynomial(coef=[list of coeffs], domain=[xstart, xend])

        if render:
            plt.plot(*p.linspace(), label="Ascending", color="orange")
            plt.scatter(xs, ys, label="Raw Points")
            plt.legend()
            plt.xlabel("Oxygen Saturation (%)")
            plt.ylabel("Relative HIF Expression Rates")
            plt.title("Oxygen Saturation to Relative HIF Expression Rates")
            plt.show()

        return [round(c, 2) for c in list(p.coef)]

    '''
    Returns coefficients for ultra-hypoxia and hypoxia polynomials.
    :param degree - (Optional, defaults to 2) If degree > 2 function has to be changed to include more points)
    :param render - (Optional, defaults to False) To be used for testing, if set to true plots the polynomials
    :return ultraHypoxiaCoeffs - Coefficients of polynomial for ultra-hypoxia (as a list)
    :return hypoxiaCoeffs - Coefficients of polynomial for hypoxia (as a list)
    '''
    def getOxygenToHif(self, degree=1, render=False):
        points = [
            # Points below 0.02
            (0, 0),
            (self.ultraHypoxiaThreshold, self.maxHIF),
            # Points to 0.02
            (self.ultraHypoxiaThreshold, self.maxHIF),
            (self.hypoxiaThreshold, 1),
        ]

        xs = [p[0] for p in points]
        ys = [p[1] for p in points]

        # Fitting points between 0.02 and 0.2
        xsAscending = [xs[0], xs[1]]
        ysAscending = [ys[0], ys[1]]

        pA = Polynomial.fit(xsAscending, ysAscending, degree)

        # Fitting points between 0.02 and 0
        xsDescending = [xs[2], xs[3]]
        ysDescending = [ys[2], ys[3]]

        pD = Polynomial.fit(xsDescending, ysDescending, degree)

        # Getting data to reproduce polynomials, use of sort:
        # p = Polynomial(coef=[list of coeffs], domain=[xstart, xend])

        if render:
            plt.plot(*pA.linspace(), label="Ascending", color="orange")
            plt.plot(*pD.linspace(), label="Descending", color="green")
            plt.scatter(xs, ys, label="Raw Points")
            plt.legend()
            plt.xlabel("Oxygen Saturation (%)")
            plt.ylabel("Relative HIF Expression Rates")
            plt.title("Oxygen Saturation to Relative HIF Expression Rates")
            plt.show()

        return [round(c, 2) for c in list(pA.coef)], [round(c, 2) for c in list(pD.coef)]


if __name__ == "__main__":
    ohrg = OxygenHIFRelationsGenerator(1, 16, 0.02, 0.2, 0.4, 90)
    #ohrg.getOxygenToHif(render=True)
    #ohrg.getHifToVegf(render=True)
    #ohrg.getHifToPSynthesis(render=True)
    #ohrg.getHifToMetabolicRate(render=True)
    #ohrg.getOxygenToHifWarburg(render=True)
