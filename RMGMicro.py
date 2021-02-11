# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 16:41:32 2021

Github: https://github.com/tjczec01

@author: Travis J Czechorski

E-mail: tjczec01@gmail.com

"""

# Data sources
database(
    thermoLibraries = ['primaryThermoLibrary'],
    reactionLibraries = [],
    seedMechanisms = [],
    kineticsDepositories = ['training'],
    kineticsFamilies = 'default',
    kineticsEstimator = 'rate rules',
)

# List of species
species(
    label='ethane',
    reactive=True,
    structure=SMILES("CC"),
)

# Reaction systems
simpleReactor(
    temperature=(1350,'K'),
    pressure=(1.0,'bar'),
    initialMoleFractions={
        "ethane": 1.0,
    },
    terminationConversion={
        'ethane': 0.9,
    },
    terminationTime=(1e6,'s'),
)

simulator(
    atol=1e-16,
    rtol=1e-8,
)

model(
    toleranceKeepInEdge=0.0,
    toleranceMoveToCore=0.1,
    toleranceInterruptSimulation=0.1,
    maximumEdgeSpecies=100000,
    filterReactions=True,
)

options(
    units='si',
    generateOutputHTML=True,
    generatePlots=False,
    saveEdgeSpecies=True,
    saveSimulationProfiles=True,
)
# 5.2. 1,3-hexadiene pyrolysis
# This example models the pyrolysis of 1,3-hexadiene and demonstrates the effect of turning on the pressure-dependence module within RMG.

# Data sources
database(
    thermoLibraries = ['primaryThermoLibrary', 'GRI-Mech3.0'],
    reactionLibraries = [],
    seedMechanisms = [],
    kineticsDepositories = ['training'], 
    kineticsFamilies = 'default',
    kineticsEstimator = 'rate rules',
)

# Constraints on generated species
generatedSpeciesConstraints(
    maximumRadicalElectrons = 2,
)

# List of species
species(
    label='HXD13',
    reactive=True,
    structure=SMILES("C=CC=CCC"),
)
species(
    label='CH4',
    reactive=True,
    structure=SMILES("C"),
)
species(
    label='H2',
    reactive=True,
    structure=adjacencyList(
        """
        1 H u0 p0 {2,S}
        2 H u0 p0 {1,S}
        """),
)
species(
    label='N2',
    reactive=False,
    structure=InChI("InChI=1/N2/c1-2"),
)

# Reaction systems
simpleReactor(
    temperature=(1350,'K'),
    pressure=(1.0,'bar'),
    initialMoleFractions={
        "HXD13": 6.829e-4,
        "CH4": 0.104,
        "H2": 0.0156,
        "N2": 0.8797,
    },
    terminationConversion={
        'HXD13': 0.9,
    },
    terminationTime=(1e0,'s'),
)

simulator(
    atol=1e-16,
    rtol=1e-8,
)

model(
    toleranceKeepInEdge=0.0,
    toleranceMoveToCore=0.5,
    toleranceInterruptSimulation=0.5,
    maximumEdgeSpecies=100000
)

quantumMechanics(
    software='mopac',
    method='pm3',
    # fileStore='QMfiles', # relative to where it is run. Defaults within the output folder.
    scratchDirectory = None, # not currently used
    onlyCyclics = True,
    maxRadicalNumber = 0,
    )

pressureDependence(
    method='modified strong collision',
    maximumGrainSize=(0.5,'kcal/mol'),
    minimumNumberOfGrains=250,
    temperatures=(300,2000,'K',8),
    pressures=(0.01,100,'bar',5),
    interpolation=('Chebyshev', 6, 4),
)

options(
    units='si',
    generateOutputHTML=False,
    generatePlots=False,
)
# 5.3. Commented input file
# This is a fully commented input file with all optional blocks for new users to better understand the options of rmg input files

# Uncomment either of the blocks below to restart from a seed mechanism
#
# # Option 1: Specify the path to an RMG (version > 2.4.1) generated seed mechanism folder, which contains all of the
# # required files (core and edge seed, filters and mappings) in their default locations and names in the seed folder.
# restartFromSeed(path='seed')  # Location of the seed mechanism (with `filters` subfolder) to load for restarting
#
# # Option 2: Specify the paths of each of the required files individually.
# restartFromSeed(coreSeed='seed/seed'  # Path to core seed folder. Must contain `reactions.py` and `dictionary.txt`
#                 edgeSeed='seed/seed_edge'  # Path to edge seed folder containing `reactions.py` and `dictionary.txt`
#                 filters='seed/filters/filters.h5',
#                 speciesMap='seed/filters/species_map.yml')

# Data sources
database(
    # overrides RMG thermo calculation of RMG with these values.
    # libraries found at http://rmg.mit.edu/database/thermo/libraries/
    # if species exist in multiple libraries, the earlier libraries overwrite the
    # previous values
    thermoLibraries=['BurkeH2O2', 'primaryThermoLibrary', 'DFT_QCI_thermo', 'CBS_QB3_1dHR'],
    # overrides RMG transport calculations with these values.
    # if species exist in multiple libraries, the earlier libraries overwrite the previous values
    transportLibraries=['PrimaryTransportLibrary.py'],
    # overrides RMG kinetics estimation if needed in the core of RMG.
    # list of libraries found at http://rmg.mit.edu/database/kinetics/libraries/
    # libraries can be input as either a string or tuple of form ('library_name',True/False)
    # where a `True` indicates that all unused reactions will be automatically added
    # to the chemkin file at the end of the simulation. Placing just string values
    # defaults the tuple to `False`. The string input is sufficient in almost
    # all situations
    reactionLibraries=[('C3', False)],
    # seed mechanisms are reactionLibraries that are forced into the initial mechanism
    # in addition to species listed in this input file.
    # This is helpful for reducing run time for species you know will appear in
    # the mechanism.
    seedMechanisms=['BurkeH2O2inN2', 'ERC-FoundationFuelv0.9'],
    # lists specific families used to generate the model. 'default' uses a list of
    # families from RMG-Database/input/kinetics/families/recommended.py
    # a visual list of families is available in PDF form at RMG-database/families
    kineticsFamilies='default',
    # this is normally not changed in general RMG runs.  Usually used for testing with
    # outside kinetics databases
    kineticsDepositories='default',
    # specifies how RMG calculates rates.  currently, the only option is 'rate rules'
    kineticsEstimator='rate rules',
)

# List of species
# list initial and expected species below to automatically put them into the core mechanism.
# 'structure' can utilize method of SMILES("put_SMILES_here"),
# adjacencyList("""put_adj_list_here"""), or InChI("put_InChI_here")
# for molecular oxygen, use the smiles string [O][O] so the triplet form is used
species(
    label='butane',
    reactive=True,  # this parameter is optional if true
    structure=SMILES("CCCC"),
)
species(
    label='O2',
    structure=SMILES("[O][O]"),
)
species(
    label='N2',
    reactive=False,
    structure=adjacencyList("""
    1 N u0 p1 c0 {2,T}
    2 N u0 p1 c0 {1,T}
    """),
)
# You can list species not initially in reactor to make sure RMG includes them in the mechanism
species(
    label='QOOH',
    reactive=True,
    structure=SMILES("OOCC[CH]C")
)
species(
    label='CO2',
    reactive=True,
    structure=SMILES("O=C=O")
)

# Reaction systems
# currently RMG models only constant temperature and pressure as homogeneous batch reactors.
# two options are: simpleReactor for gas phase or liquidReactor for liquid phase
# use can use multiple reactors in an input file for each condition you want to test.
simpleReactor(
    # specifies reaction temperature with units
    temperature=(700, 'K'),
    # specifies reaction pressure with units
    pressure=(10.0, 'bar'),
    # list initial mole fractions of compounds using the label from the 'species' label.
    # RMG will normalize if sum/=1
    initialMoleFractions={
        "N2": 4,
        "O2": 1,
        "butane": 1. / 6.5,
    },
    # number of simulations used to explore variable temperature and pressure reactors
    nSims=6,
    # the following two values specify when to determine the final output model
    # only one must be specified
    # the first condition to be satisfied will terminate the process
    terminationConversion={
        'butane': .99,
    },
    terminationTime=(40, 's'),
    # the next two optional values specify how RMG computes sensitivities of
    # rate coefficients with respect to species concentrations.
    # sensitivity contains a list of species' labels to conduct sensitivity analysis on.
    # sensitvityThreshold is the required sensitiviy to be recorded in the csv output file
    # sensitivity=['CH4'],
    # sensitivityThreshold=0.0001,
)

# liquidReactor(
#     temperature=(500, 'K'),
#     initialConcentrations={
#         "N2": 4,
#         "O2": 1,
#         "CO": 1,
#     },
#     terminationConversion=None,
#     terminationTime=(3600, 's'),
#     sensitivity=None,
#     sensitivityThreshold=1e-3
# )
# liquid reactors also have solvents, you can specify one solvent
# list of solvents available at : http://rmg.mit.edu/database/solvation/libraries/solvent/
# solvation('water')

# determines absolute and relative tolerances for ODE solver and sensitivities.
# normally this doesn't cause many issues and is modified after other issues are
# ruled out
simulator(
    atol=1e-16,
    rtol=1e-8,
    # sens_atol=1e-6,
    # sens_rtol=1e-4,
)

# used to add species to the model and to reduce memory usage by removing unimportant additional species.
# all relative values are normalized by a characteristic flux at that time point
model(
    # determines the relative flux to put a species into the core.
    # A smaller value will result in a larger, more complex model
    # when running a new model, it is recommended to start with higher values and then decrease to converge on the model
    toleranceMoveToCore=0.1,
    # comment out the next three terms to disable pruning
    # determines the relative flux needed to not remove species from the model.
    # Lower values will keep more species and utilize more memory
    toleranceKeepInEdge=0.01,
    # determines when to stop a ODE run to add a species.
    # Lower values will improve speed.
    # if it is too low, may never get to the end simulation to prune species.
    toleranceInterruptSimulation=1,
    # number of edge species needed to accumulate before pruning occurs
    # larger values require more memory and will prune less often
    maximumEdgeSpecies=100000,
    # minimum number of core species needed before pruning occurs.
    # this prevents pruning when kinetic model is far away from completeness
    minCoreSizeForPrune=50,
    # make sure that the pruned edge species have existed for a set number of RMG iterations.
    # the user can specify to increase it from the default value of 2
    minSpeciesExistIterationsForPrune=2,
    # filter the reactions during the enlarge step to omit species from reacting if their
    # concentration are deemed to be too low
    filterReactions=False,
    # for bimolecular reactions, will only allow them to react if
    # filterThreshold*C_A*C_B > toleranceMoveToCore*characteristic_rate
    # and if filterReactions=True
    filterThreshold=1e8,
)

options(
    # provides a name for the seed mechanism produced at the end of an rmg run default is 'Seed'
    name='SeedName',
    # if True (default) every iteration it saves the current model as libraries/seeds
    # (and deletes the old one)
    # Unlike HTML this is inexpensive time-wise
    # note a seed mechanism will be generated at the end of a completed run and some incomplete
    # runs even if this is set as False
    generateSeedEachIteration=True,
    # If True the mechanism will also be saved directly as kinetics and thermo libraries in the database
    saveSeedToDatabase=False,
    # only option is 'si'
    units='si',
    # Draws images of species and reactions and saves the model output to HTML.
    # May consume extra memory when running large models.
    generateOutputHTML=True,
    # generates plots of the RMG's performance statistics. Not helpful if you just want a model.
    generatePlots=False,
    # saves mole fraction of species in 'solver/' to help you create plots
    saveSimulationProfiles=False,
    # gets RMG to output comments on where kinetics were obtained in the chemkin file.
    # useful for debugging kinetics but increases memory usage of the chemkin output file
    verboseComments=False,
    # gets RMG to generate edge species chemkin files. Uses lots of memory in output.
    # Helpful for seeing why some reaction are not appearing in core model.
    saveEdgeSpecies=False,
    # Sets a time limit in the form DD:HH:MM:SS after which the RMG job will stop. Useful for profiling on jobs that
    # do not converge.
    # wallTime = '00:00:00',
    # Forces RMG to import library reactions as reversible (default). Otherwise, if set to True, RMG will import library
    # reactions while keeping the reversibility as as.
    keepIrreversible=False,
    # Allows families with three products to react in the diverse direction (default).
    trimolecularProductReversible=True,
    # Allows a seed to be saved every n iterations.
    # The default of -1 causes the iteration to only be saved at the end of the RMG job
    saveSeedModulus=-1
)

# optional module allows for correction to unimolecular reaction rates at low pressures and/or temperatures.
pressureDependence(
    # two methods available: 'modified strong collision' is faster and less accurate than 'reservoir state'
    method='modified strong collision',
    # these two categories determine how fine energy is descretized.
    # more grains increases accuracy but takes longer
    maximumGrainSize=(0.5, 'kcal/mol'),
    minimumNumberOfGrains=250,
    # the conditions for the rate to be output over
    # parameter order is: low_value, high_value, units, internal points
    temperatures=(300, 2200, 'K', 2),
    pressures=(0.01, 100, 'bar', 3),
    # The two options for interpolation are 'PDepArrhenius' (no extra arguments) and
    # 'Chebyshev' which is followed by the number of basis sets in
    # Temperature and Pressure. These values must be less than the number of
    # internal points specified above
    interpolation=('Chebyshev', 6, 4),
    # turns off pressure dependence for molecules with number of atoms greater than the number specified below
    # this is due to faster internal rate of energy transfer for larger molecules
    maximumAtoms=15,
)

# optional block adds constraints on what RMG can output.
# This is helpful for improving the efficiency of RMG, but wrong inputs can lead to many errors.
generatedSpeciesConstraints(
    # allows exceptions to the following restrictions
    allowed=['input species', 'seed mechanisms', 'reaction libraries'],
    # maximum number of each atom in a molecule
    maximumCarbonAtoms=4,
    maximumOxygenAtoms=7,
    maximumNitrogenAtoms=0,
    maximumSiliconAtoms=0,
    maximumSulfurAtoms=0,
    # max number of non-hydrogen atoms
    # maximumHeavyAtoms=20,
    # maximum radicals on a molecule
    maximumRadicalElectrons=1,
    # maximum number of singlet carbenes (lone pair on a carbon atom) in a molecule
    maximumSingletCarbenes=1,
    # maximum number of radicals on a molecule with a singlet carbene
    # should be lower than maximumRadicalElectrons in order to have an effect
    maximumCarbeneRadicals=0,
    # If this is false or missing, RMG will throw an error if the more less-stable form of O2 is entered
    # which doesn't react in the RMG system. normally input O2 as triplet with SMILES [O][O]
    # allowSingletO2=False,
    # maximum allowed number of non-normal isotope atoms:
    # maximumIsotopicAtoms=2,
)

# optional block allows thermo to be estimated through quantum calculations
# quantumMechanics(
#     # the software package for calculations...can use 'mopac' or 'gaussian' if installed
#     software='mopac',
#     # methods available for calculations. 'pm2' 'pm3' or 'pm7' (last for mopac only)
#     method='pm3',
#     # where to store calculations
#     fileStore='QMfiles',
#     # where to store temporary run files
#     scratchDirectory=None,
#     # onlyCyclics allows linear molecules to be calculated using bensen group addivity....need to verify
#     onlyCyclics=True,
#     # how many radicals should be utilized in the calculation.
#     # If the amount of radicals is more than this, RMG will use hydrogen bond incrementation method
#     maxRadicalNumber=0,
# )

# optional block allows thermo to be estimated through ML estimator
# mlEstimator(
#     thermo=True,
#     # Name of folder containing ML architecture and parameters in database
#     name='main',
#     # Limits on atom numbers
#     minHeavyAtoms=1,
#     maxHeavyAtoms=None,
#     minCarbonAtoms=0,
#     maxCarbonAtoms=None,
#     minOxygenAtoms=0,
#     maxOxygenAtoms=None,
#     minNitrogenAtoms=0,
#     maxNitrogenAtoms=None,
#     # Limits on cycles
#     onlyCyclics=False,
#     onlyHeterocyclics=False, # If onlyHeterocyclics is True, the machine learning estimator is restricted to only
#                                heterocyclics species regardless of onlyCyclics setting.
#                                But onlyCyclics should also be True if onlyHeterocyclics is True.
#     minCycleOverlap=0,  # specifies the minimum number of atoms that must be shared between any two cycles.
#                           If minCycleOverlap is greater than zero, the machine learning estimator is restricted to
#                           only cyclic species with the specified minimum cyclic overlap regardless of onlyCyclics
#                           setting.
#     # If the estimated uncertainty of the thermo prediction is greater than
#     # any of these values, then don't use the ML estimate
#     H298UncertaintyCutoff=(3.0, 'kcal/mol'),
#     S298UncertaintyCutoff=(2.0, 'cal/(mol*K)'),
#     CpUncertaintyCutoff=(2.0, 'cal/(mol*K)')
# )