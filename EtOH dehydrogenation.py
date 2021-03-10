# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 01:52:14 2021

Github: https://github.com/tjczec01

@author: Travis J Czechorski

E-mail: tjczec01@gmail.com

"""

# Data sources

database(
    thermoLibraries=['surfaceThermoPt111', 'primaryThermoLibrary', 'thermo_DFT_CCSDTF12_BAC','DFT_QCI_thermo', "Klippenstein_Glarborg2016"], # 'surfaceThermoPt' is the default. Thermo data is derived using bindingEnergies for other metals 
    reactionLibraries = [('Surface/CPOX_Pt/Deutschmann2006', False), ("TEOS", False)], # when Ni is used change the library to Surface/Deutschmann_Ni 
    seedMechanisms = [],
    kineticsDepositories = ['training'],
    kineticsFamilies = ['surface','default'],
    kineticsEstimator = 'rate rules',

)


# mgosio2 = 100.3887 # g/mol
# sa = 80.4 # m^2/g
# cacm = sa*10000.0  # cm^2/g
# sacm = 1.0 / cacm  # g/cm^2
# saf = sacm/mgosio2  #  mol/cm^2
# print(saf)



#  https://doi.org/10.1002/cctc.202001007

# MgO = [50.25, 50.80, 51.00] # Mg-O
# Si = [20.3, 17.4, 20.3]
# O = [60.8, 56.2, 59.6]
# C = [2.2, 9.6, 6.4]

# Mgavg = sum(MgO)/len(MgO)
# Siavg = sum(Si)/len(Si)
# Oavg = sum(O)/len(O)
# Cavg = sum(C)/len(C)

catalystProperties(
bindingEnergies = {
'H': (50.30, 'eV/molecule'),
'O': (50.68333333333334, 'eV/molecule'),
'Si': (101.3, 'eV/molecule'),},
surfaceSiteDensity=(1.2389652366524949e-08, 'mol/cm^2'),
)

# List of species
species(
	label='C2H6O',
	reactive=True,
	structure=InChI("InChI=1/C2H6O/c1-2-3/h3H,2H2,1H3"),
)


species(
	label='C2H4',
	reactive=True,
	structure=SMILES("C=C"),  # structure=SMILES("CC")
)


species(
	label='C2OH4',
	reactive=True,
	structure=SMILES("CC=O"),  # structure=SMILES("CC")
)


species(
	label='1-Methoxypropan-2-ol',  # C4H10O2
	reactive=True,
	structure=SMILES("CC(O)COC"),  # structure=SMILES("CC")
)

species(
	label='Ethoxyethane',
	reactive=True,
	structure=SMILES("CCOCC"),  # structure=SMILES("CC")
)

species(
    label='H2O',
    reactive=True,
    structure=SMILES("O"),
)

species(
    label='H2',
    reactive=True,
    structure=SMILES("[H][H]"),
)

species(
    label='site',
    reactive=True,
    structure=adjacencyList("1 X u0"),
)


surfaceReactor(
    temperature=(723, 'K'),
    initialPressure=(1.0, 'atm'),
    initialGasMoleFractions={
        "C2H6O": 1.0,
        "C2H4": 0.0,
        "C2OH4": 0.0,
        "1-Methoxypropan-2-ol": 0.0,
        "Ethoxyethane": 0.0,
        "H20": 0.0,
        "H2": 0.0,
    },
    initialSurfaceCoverages={
        "site": 1.0,
    },
    surfaceVolumeRatio=(1.0e5, 'm^-1'),
    terminationConversion = {"C2H6O": 0.1,},
    terminationTime=(28800.0, 's'),
)


simulator(
    atol=1e-18,
    rtol=1e-12,
)

model(
    toleranceKeepInEdge=0.0,
    toleranceMoveToCore=1e-6,
    toleranceInterruptSimulation=0.1,
    maximumEdgeSpecies=100000,
)

options(
    units='si',
    generateOutputHTML=True,
    generatePlots=True, # Enable to make plots of core and edge size etc.. But takes a lot of the total runtime!
    saveEdgeSpecies=True,
    saveSimulationProfiles=True,
    verboseComments=True,
)


# Reaction systems
# simpleReactor(
#     temperature=(723,'K'),
#     pressure=(1.0,'atm'),
#     initialMoleFractions={
#         "EtOH": 1.0,
#     },
#     terminationConversion={
#         'EtOH': 0.1,
#     },
#     terminationTime=(28800,'s'),
# )


# simulator(
#     atol=1e-16,
#     rtol=1e-8,
# )

# model(
#     toleranceKeepInEdge=0.0,
#     toleranceMoveToCore=0.1,
#     toleranceInterruptSimulation=0.1,
#     maximumEdgeSpecies=100000,
#     filterReactions=True,
# )

# options(
#     units='si',
#     generateOutputHTML=True,
#     generatePlots=False,
#     saveEdgeSpecies=True,
#     saveSimulationProfiles=True,
# )

# database(
#     thermoLibraries = ['primaryThermoLibrary'],
#     reactionLibraries = [],
#     seedMechanisms = [],
#     kineticsDepositories = ['training'],
#     kineticsFamilies = 'default',
#     kineticsEstimator = 'rate rules',
# )

# # List of species
# species(
#     label='ethane',
#     reactive=True,
#     structure=SMILES("CC"),
# )

# Reaction systems
# simpleReactor(
#     temperature=(1350,'K'),
#     pressure=(1.0,'bar'),
#     initialMoleFractions={
#         "ethane": 1.0,
#     },
#     terminationConversion={
#         'ethane': 0.9,
#     },
#     terminationTime=(1e6,'s'),
# )


# surfaceReactor(
#     temperature=(723,'K'),
#     initialPressure=(1.0, 'atm'),
#     initialGasMoleFractions={
#         "CH4": 0.0500,
#         "O2": 0.1995,
#         "N2": 0.7505,
#     },
#     initialSurfaceCoverages={
#         "X": 1.0,
#     },
#     surfaceVolumeRatio=(1.0e4, 'm^-1'),
#     terminationConversion = { "CH4":0.9 },
#     terminationRateRatio=0.01
# )


# 5.2. 1,3-hexadiene pyrolysis
# This example models the pyrolysis of 1,3-hexadiene and demonstrates the effect of turning on the pressure-dependence module within RMG.

# Data sources
# database(
#     thermoLibraries = ['primaryThermoLibrary', 'GRI-Mech3.0'],
#     reactionLibraries = [],
#     seedMechanisms = [],
#     kineticsDepositories = ['training'], 
#     kineticsFamilies = 'default',
#     kineticsEstimator = 'rate rules',
# )
