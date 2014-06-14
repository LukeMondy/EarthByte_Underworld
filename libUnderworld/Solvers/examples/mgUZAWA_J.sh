#! /bin/bash

export UWPATH=`./getUWD.sh`
export UWEXEC=$UWPATH/build/bin/Underworld
OUT="outputMGUZAWA_J"
mkdir $OUT >& /dev/null

$UWEXEC $UWPATH/Experimental/InputFiles/+FlatEarth/MantleConvectionCrustPIC.xml \
       $UWPATH/Experimental/InputFiles/StokesPIC+AugmentedLagrangianUZ.xml \
       $UWPATH/Experimental/InputFiles/Experimental_Components/MultigridForRegularUzawa.xml \
       $UWPATH/Experimental/InputFiles/+FlatEarth/BaseApps/PeriodicLinearMesh.xml \
       -options_file $UWPATH/Experimental/InputFiles/options/options-uzawa-mg.opt \
       --components.constitutiveMatrix.incompressibility_Penalty=10.0 \
       --components.constitutiveMatrix.viscosity_weighting=false \
            --outputPath=$OUT                                       \
            --mgLevels=2                                                   \
            --dim=3                                                        \
            --elementResI=16                                               \
            --elementResK=16                                               \
            --elementResJ=16                  \
            --maxTimeSteps=1                                              \
            --gravity=1.0e6                                                \
            --components.temperatureDependence.eta0=1.0                    \
            --components.temperatureDependence.activationEnergy=0.0        \
            --components.temperatureDependence2.eta0=10000.0               \
            --components.temperatureDependence2.activationEnergy=0.0       \
            --components.lightMaterial.density=0.0                         \
            --nonLinearTolerance=1.0e-3                                    \
            --nonLinearMaxIterations=25                                    \
            --cBlockX1=-0.355                                                \
            --cBlockX2=0.355                                                \
            --cBlockZ1=-0.355                                                \
            --cBlockZ2=0.355                                                \
            --cBlockY1=0.90                                                \
            --cBlockY2=1.0                                                 \

            --minX=-0.5                                                    \
            --maxX=0.5                                                     \
            --minZ=-0.5                                                    \
            --maxZ=0.5                                                     \
            --minY=0.0                                                     \
            --maxY=1.0                                                     \

