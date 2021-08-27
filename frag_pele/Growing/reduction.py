
class FragmentReduction():
    """

    """
    def __init__(self, iterations, start_growing_from=None):
        """

        Returns
        -------

        """
        self._lam_initial = self._define_lam(iterations, start_growing_from)
        self._inv_lam = 1-self._lam_initial

    def _define_lam(self, iterations, start_growing_from):
        if not start_growing_from:
            lam_initial = 1 / (iterations + 1)
        else:
            teoric_initial = 1 / (iterations + 1)
            i = 1
            while i <= iterations:
                lam_initial = (i) * teoric_initial
                i = i + 1
                if lam_initial > start_growing_from:
                    break
        return lam_initial


