import underworld as _underworld
import os as _os


def options( optionsFilename="", optionsString="" ):
    """
    Injects options into the PETSc options database. Note that this is called by the
    Stokes block ksp solver as it runs and is therefore not guaranteed to affect anything which PETSc has
    already set up. But it does work for all the options associated with these solvers.

    This is implemented via PetscOptionsInsertFile() and PetscOptionsInsertString().
    It seems that the latter does not over-ride the former so use -log_summary and be careful !!

    Args:
        globalDict         (dict)   : The global dictionary
        optionsFilename (string) : A petsc options file
        optionsString   (string) : A petsc options string

    Returns:
        Nothing yet
    """

    _solverType = _underworld.solvers.solverType()

    if _solverType == "None":
        print "solvers.options(): You should first call a solver activation function "
        return

    if "Uzawa" in _solverType:
        print " !! Petsc Options cannot be set in code for the Uzawa solver"
        return

    globalDict = _underworld.dictionary.GetDictionary()

    # If it exists, pass it to PETSc (can't really check if it is sane)

    globalDict["components"]["stokesblockkspinterface"]["OptionsFile"] = optionsFilename
    if optionsFilename != "":
        print " *  Adding options file - {}.".format(optionsFilename)

    # The default of "" is already assumed by the ksp component
    globalDict["components"]["stokesblockkspinterface"]["OptionsString"] = optionsString
    if optionsString != "":
        print " *  Adding options string -> {}.".format(optionsString)


# Here are some "hard coded" options files as strings:


def petsc_options_StokesBlockKSP_accelerating_mg():
    """
    returns petsc options filename corresponding to cached copy of options-scr-mg-accelerating.opt

    Note: the petsc options in the dictionary is a one-shot affair and so
    we have to manage any merging of options strings outselves (at the moment)

    """

    petsc_opt_filename = _os.path.join(_os.path.dirname(__file__), "_options-scr-mg-accelerating.opt")

    return (petsc_opt_filename)
