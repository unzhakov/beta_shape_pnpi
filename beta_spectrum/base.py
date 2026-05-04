# base.py
from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from typing import Optional

import numpy as np

ArrayLike = np.ndarray


class SpectrumComponent(ABC):
    """
    Base class for multiplicative components of beta-spectrum.

    Parameters
    ----------
    logger : logging.Logger, optional
        Logger instance for debug/info output. If None, logging is disabled.
        Opt-in design: components are silent by default for library use.
    """

    def __init__(self, logger: Optional["logging.Logger"] = None) -> None:
        self._logger = logger

    @abstractmethod
    def __call__(self, W: ArrayLike) -> ArrayLike:
        """
        Evaluate component at total energy W (in m_e units).
        """
        pass
