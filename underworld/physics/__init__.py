
# Physics module - replaces the various "module XML files" in Underworld:
#    LidDriven.xml
#    LidDrivenPIC.xml
#    RayleighTaylor.xml
#    ...
#
import setup

_problem_type = "None"


def problemType():
    """
    Returns a string representing the physics of the problem that has been set up

    """
    return (_problem_type)


def setProblemType( typeString ):
    """
    Sets the string representing the physics of the problem that has been set up

    """

    global _problem_type
    _problem_type = typeString

    return
