import copy
import warnings
import numpy as np


__title__ = "HOPS"
__author__ = "D. I. G. Bennett, L. Varvelo"
__version__ = "1.2"

"""
Test1. description
"""


###########################################################
class HopsTrajectory:
    """
     HopsTrajectory is the class that a user should interface with to run a single 
     trajectory calculation. TEST 2 - class
    """

    def __init__(self):
        """
        TEST 3.
    
        :param kind: Optional "kind" of ingredients.
        :type kind: list[str] or None
        :raise lumache.InvalidKindError: If the kind is invalid.
        :return: The ingredients list.
        :rtype: list[str]
        """
        # Variables for the current state of the system


    def initialize(self, psi_0):
        """
        TEST 3.
    
        :param kind: Optional "kind" of ingredients.
        :type kind: list[str] or None
        :raise lumache.InvalidKindError: If the kind is invalid.
        :return: The ingredients list.
        :rtype: list[str]
        """
        return psi_0
