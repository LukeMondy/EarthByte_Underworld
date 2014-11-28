
# gLucifer visualisation module

'''
    We try to provide a more flexible interface to the database at least

'''

import setup
import live

# Initialisation
_gLucifer = False


def gLucifer():
    """
    Flag to determine if gLucifer is active
    """
    return _gLucifer


def setGLucifer( Flag ):
    """
    Sets the flag to determine if gLucifer is active
    """

    global _gLucifer

    _gLucifer = Flag
    return
