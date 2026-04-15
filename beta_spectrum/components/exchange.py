from pathlib import Path
import numpy as np
from beta_spectrum.base import SpectrumComponent


class ExchangeCorrection(SpectrumComponent):
    """
    Atomic xchange correction using Hayen 2018 Table X coefficients (App. G)
    """

    def __init__(self, Z: int, filename=None, eps: float = 1e-6):
        self.Z = Z
        self.eps = eps

        if filename is None:
            filename = self._default_coeff_path()

        self.coeffs = self._load_coeffs(filename)

    def _load_coeffs(self, filename):
        data = np.genfromtxt(filename, delimiter=",", names=True)

        for row in data:
            if int(row["Z"]) == self.Z:
                return {
                    "a": row["a"],
                    "b": row["b"],
                    "c": row["c"],
                    "d": row["d"],
                    "e": row["e"],
                    "f": row["f"],
                    "g": row["g"],
                    "h": row["h"],
                    "i": row["i"],
                }
        raise ValueError(f"No exchange coefficients available for Z={self.Z}")

    def _default_coeff_path(self):
        return Path(__file__).resolve().parent.parent / "data" / "exchange_coeff.csv"

    def __call__(self, W: np.ndarray) -> np.ndarray:
        W = np.asarray(W)

        # Define physical cutoff
        W_cut = 1.005

        # Compute full expression ONLY above cutoff
        W_safe = np.maximum(W, W_cut)

        Wp = W_safe - 1.0

        a = self.coeffs["a"]
        b = self.coeffs["b"]
        c = self.coeffs["c"]
        d = self.coeffs["d"]
        e = self.coeffs["e"]
        f_coef = self.coeffs["f"]
        g = self.coeffs["g"]
        h = self.coeffs["h"]
        i = self.coeffs["i"]

        term1 = a / Wp
        term2 = b / (Wp**2)
        term3 = c * np.exp(-d * Wp)

        base = np.maximum(W_safe - f_coef, self.eps)
        term4 = e * np.sin((base**g + h)) / (W_safe**i)

        X_safe = 1.0 + term1 + term2 + term3 + term4

        return X_safe
