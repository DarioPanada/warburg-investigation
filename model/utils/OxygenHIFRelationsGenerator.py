import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.switch_backend("agg")
from numpy.polynomial import Polynomial


class OxygenHIFRelationsGenerator():
    """
    This class allows to generate coefficients for function relating oxygen
    concentration to hif expression rates
    and hif expression rates to vegf secretion rate, proliferation rate and
    metabolic rate.

    Parameters
    ----------
    min_hif : int
        Minimum HIF expression rates
    max_hif : int
        Maximum HIF expression rates
    ultra_hypoxia_threshold : float
        A value x (0 < x < 1) such that oxygen
        concentrations 0 to x are considered ultra
        hypoxia. (HIF expression rate at 0 will be 0 and at x will be maxHIF)
    hypoxia_threshold : float
        A value x (ultraHypoxiaThreshold <
        hypoxiaThreshold < 1) such that HIF expression rates
        at hypoxiaThreshold will 1. HIF Expression rate will decrease from a
        value maxHIF at ultraHypoxia threshold
        until reaching minHIF at hypoxiaThreshold
    enhanced_hypoxic_threshold : float
        Hypoxic threshold to be used during
        warburg metabolism
    base_oxygen_metabolic_rate : float
        oxygen uptake rate in normoxia
    min_p_synthesis : float, Optional
        Lower value for probability of progressing into synthesis. Defaults
        to 0
    """

    def __init__(self, min_hif, max_hif, ultra_hypoxia_threshold,
                 hypoxia_threshold,
                 enhanced_hypoxic_threshold, base_oxygen_metabolic_rate,
                 min_p_synthesis=0):
        self.min_hif = min_hif
        self.max_hif = max_hif
        self.ultra_hypoxia_threshold = ultra_hypoxia_threshold
        self.hypoxia_threshold = hypoxia_threshold
        self.min_p_synthesis = min_p_synthesis
        self.enhanced_hypoxic_threshold = enhanced_hypoxic_threshold
        self.base_oxygen_metabolic_rate = base_oxygen_metabolic_rate

    def get_hif_to_metabolic_rate(self, render=False):
        """
        Returns coefficients for hif to metabolic rate. The domain of the
        function is [0, maxHIF].

        Parameters
        ----------
        render : bool, optional
            If set to true, will display the generated function

        Returns
        -------
        list
            Function coefficients
        """
        points = [
            (0, self.base_oxygen_metabolic_rate),
            (self.max_hif, 0)
        ]

        xs = [p[0] for p in points]
        ys = [p[1] for p in points]

        p = Polynomial.fit(xs, ys, 2)
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

        return [round(c, 2) for c in list(p.coef)]

    def get_hif_to_p_synthesis(self, degree=1, render=False):
        """
        Returns coefficients for hif to probability of synthesis. The domain
        of
        the function is [0, maxHIF].

        Parameters
        ----------
        degree : int, optional
            Degree of the generated polynomial. Optional, defaults to 1.
        render : bool, optional
            If set to true, will display the generated function

        Returns
        -------
        list
            Function coefficients
        """

        points = [
            (0, self.min_p_synthesis),
            (self.max_hif, 1)
        ]

        xs = [p[0] for p in points]
        ys = [p[1] for p in points]

        p = Polynomial.fit(xs, ys, degree)
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

        return [round(c, 2) for c in list(p.coef)]

    def get_hif_to_vegf(self, degree=2, render=False):
        """
        Returns coefficients for hif to probability of vegf secretion rate.
        The domain of the function is [0, maxHIF].

        Parameters
        ----------
        degree : int, optional
            Degree of the generated polynomial. Optional, defaults to 2.
        render : bool, optional
            If set to true, will display the generated function

        Returns
        -------
        list
            Function coefficients
        """
        points = [
            (0, 0),
            (self.max_hif, 1)
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

    def get_oxygen_to_hif_warburg(self, degree=1, render=False):
        """
        Returns coefficients for ultra-hypoxia and hypoxia polynomials to be
        used by cancer cells if the warburg
        switch has been activated. At present, only hypoxic threshold is
        changed whereas ultra-hypoxic is left as-is.

        Parameters
        ----------
        degree : int, optional
            Degree of the generated polynomial. Optional, defaults to 1.
        render : bool, optional
            If set to true, will display the generated function

        Returns
        -------
        list
            Coefficients of polynomial for ultra-hypoxia
        list
            Coefficients of polynomial for hypoxia
        """
        points = [
            # Points to 0.02
            (self.ultra_hypoxia_threshold, self.max_hif),
            (self.enhanced_hypoxic_threshold, self.min_hif),
        ]

        xs = [p[0] for p in points]
        ys = [p[1] for p in points]

        # Fitting points between 0.02 and 0.2
        xs = [xs[0], xs[1]]
        ys = [ys[0], ys[1]]

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

    def get_oxygen_to_hif(self, degree=1, render=False):
        """
        Returns coefficients for ultra-hypoxia and hypoxia polynomials to be
        used by cancer cells if the warburg
        switch has *NOT* been activated. At present, only hypoxic threshold is
        changed whereas ultra-hypoxic is left as-is.

        Parameters
        ----------
        degree : int, optional
            Degree of the generated polynomial. Optional, defaults to 1.
        render : bool, optional
            If set to true, will display the generated function

        Returns
        -------
        list
            Coefficients of polynomial for ultra-hypoxia
        list
            Coefficients of polynomial for hypoxia
        """
        points = [
            # Points below 0.02
            (0, 0),
            (self.ultra_hypoxia_threshold, self.max_hif),
            # Points to 0.02
            (self.ultra_hypoxia_threshold, self.max_hif),
            (self.hypoxia_threshold, 1),
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

        return [round(c, 2) for c in list(pA.coef)], [round(c, 2) for c in
                                                      list(pD.coef)]
