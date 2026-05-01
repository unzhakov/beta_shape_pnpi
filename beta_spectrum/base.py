# base.py
from __future__ import annotations
from abc import ABC, abstractmethod
import numpy as np

ArrayLike = np.ndarray


class SpectrumComponent(ABC):
    """
    Base class for multiplicative components of beta-spectrum
    """

    @abstractmethod
    def __call__(self, W: ArrayLike) -> ArrayLike:
        """
        Evaluate component  at total energy W (in m_e units)
        """
        pass
