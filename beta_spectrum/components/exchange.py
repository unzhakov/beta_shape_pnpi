from pathlib import Path
from typing import Dict, Optional

import logging

import numpy as np

from beta_spectrum.base import SpectrumComponent


class ExchangeCorrection(SpectrumComponent):
    """
    Atomic xchange correction using Hayen 2018 Table X coefficients (App. G)
    """

    def __init__(
        self,
        Z: int,
        filename: Optional[str] = None,
        eps: float = 1e-6,
        logger: Optional[logging.Logger] = None,
    ):
        super().__init__(logger=logger)
        self.Z = Z
        self.eps = eps

        if filename is None:
            filename = str(self._default_coeff_path())

        self.coeffs = self._load_coeffs(filename)

        if self._logger:
            self._logger.debug(
                "ExchangeCorrection: Z=%d, coefficients loaded from %s",
                Z,
                filename,
            )

    def _load_coeffs(self, filename: str) -> Dict[str, float]:
        data = np.genfromtxt(filename, delimiter=",", names=True)

        for row in data:
            if int(row["Z"]) == self.Z:
                return {
                    "a": float(row["a"]),
                    "b": float(row["b"]),
                    "c": float(row["c"]),
                    "d": float(row["d"]),
                    "e": float(row["e"]),
                    "f": float(row["f"]),
                    "g": float(row["g"]),
                    "h": float(row["h"]),
                    "i": float(row["i"]),
                }
        raise ValueError(f"No exchange coefficients available for Z={self.Z}")

    def _default_coeff_path(self) -> Path:
        return Path(__file__).resolve().parent.parent / "data" / "exchange_coeff.csv"

    def __call__(self, W: np.ndarray) -> np.ndarray:
        if self._logger:
            self._logger.debug(
                "ExchangeCorrection: evaluating at %d points, W range=[%.4f, %.4f]",
                len(W),
                W.min(),
                W.max(),
            )

        W = np.asarray(W)

        # Define physical cutoff
        W_cut = 1.003

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

        result = np.asarray(X_safe, dtype=np.float64)

        if self._logger:
            self._logger.debug(
                "ExchangeCorrection: output range=[%.6e, %.6e]",
                result.min(),
                result.max(),
            )

        return result
