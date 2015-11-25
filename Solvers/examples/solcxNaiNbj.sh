#! /bin/bash

#
#   Solve using BSSCR with Viscous Penalty Method matrix assembled from Assembly/src/Matrix_NaiNbj.c in Solvers Toolbox
#   This is set up via $UWPATH/Solvers/InputFiles/AugLagStokesSLE-NaiNbj.xml
#
# Options:
# -ksp_type bsscr -pc_type none   // these tell the KSPinterface Solver to use the BSSCR KSP. This is the current default so not necessary to specify explicitly.
# -ksp_k2_type <option>           // one of (NULL GG GMG DGMGD SLE)
# -ksp_k2_type SLE                // use a user defined matrix for the Augmented Lagrangian. e.g.  $UWPATH/Solvers/InputFiles/AugLagStokesSLE-NaiNbj.xml
# -ksp_k2_type NULL               // deactivate augmented method
# -ksp_k2_type GG                 // use G^T*G as augmented matrix.
# -ksp_k2_type GMG                // use G^T*M*G as augmented matrix. Where M is inverse of pressure mass matrix: 
#                                 // is set up in $UWPATH/Solvers/InputFiles/AugLagStokesSLE-GtMG.xml \
# -ksp_k2_type DGMGD              // use D*G^T*M*G*D as augmented matrix. Where D is scaling from K stiffness matrix.
# -build_const_pressure_nullspace // tell BSSCR to check for and remove the constant nullspace in the pressure space.
# -build_cb_pressure_nullspace    // tell BSSCR to check for and remove the checker board nullspace in the pressure space.
# -Q22_pc_type $UW                // specify the preconditioner for the Schur solve: one of (uw, uwscale, gtkg)
# -uzawastyle 0                   // set to 1 to imitate uzawa style solve: 0 is default
# -forcecorrection                // subtract force vector terms with penalty of a user-defined Penalty matrix when using SLE option with -ksp_k2_type from system RHS.
#
#
# Prefixes for sub KSPs in BSSCR
# -scr_          // for the schur outer solve (Pressure)
# -A11_          // for the schur inner solve
# -backsolveA11_ // for the final velocity solve
# -scr_pc_gtkg_  // for the two Laplacian solves when using the gtkg preconditioner
# -scrPCKSP_     // for preconditioner solve when imitating Uzawa
# --penaltyNumber=$PEN // using XML to pass the penalty in via the SLE currently.

count=0
PROCS=1

export UWPATH=`./getUWD.sh`
export UWEXEC="cgdb --args $UWPATH/build/bin/Underworld"
export UWEXEC="$UWPATH/build/bin/Underworld"

echo "| p its | v its | p solve time | constraint | gperror | NL its | avg P its | minp | maxp | minv | maxv | penalty | -Q22_pc_type | scale | scr | scr tol | scr norm type | A11 | A11 tol | res | MG | DIR |"

for UW in uwscale
do
for SCR in minres fgmres cg
do
for A11 in gmres
do
for SCRTOL in 1e-5
do
for A11TOL in 1e-6
do
echo "|-------+-------+------------+----------+------+------+------+------+---------+----------------+-------+-----+---------+---------------+-----+---------+-----+----|"
for PEN in 0.0 0.02 0.1 1.0 2.0 10.0 20.0 100.0 200.0 1000.0
do


#SCRP="unpreconditioned"
SCRP="default"

    if [ "$SCRP" = "unpreconditioned" ]
        then
	SCRNORMTYPE="-scr_ksp_norm_type unpreconditioned"
    else
	SCRNORMTYPE=" "
    fi

#MG=boomeramg
#MG="ml"
MG=lu
MGOP=" "
    if [ "$MG" = "gmg" ]
        then
	MGOP="$UWPATH/Solvers/InputFiles/MultigridForRegularSCR.xml -options_file ./options-scr-mg-accelerating.opt "
    fi
    if [ "$MG" = "boomeramg" ]                                            
        then                                                                                                                                                                                          
        MGOP=" -A11_pc_type hypre -XA11_pc_hypre_type boomerang -help "
    fi
    if [ "$MG" = "boomeramgfs" ]
        then
        MGOP="-A11_ksp_type gmres -A11_ksp_view -A11_pc_type fieldsplit -A11_pc_fieldsplit_block_size 2  -A11_fieldsplit_0_ksp_type richardson -A11_fieldsplit_0_pc_type hypre -A11_fieldsplit_0_pc_hypre_type boomeramg -A11_fieldsplit_1_ksp_type richardson -A11_fieldsplit_1_pc_type hypre -A11_fieldsplit_1_pc_hypre_type boomeramg "
    fi
    if [ "$MG" = "ml" ]
        then
        MGOP="-A11_pc_type ml -help -A11_pc_ml_PrintLevel 0 -A11_pc_ml_maxNlevels 2 "
    fi
    if [ "$MG" = "lu" ]
        then
        MGOP="-A11_pc_type lu -A11_ksp_type preonly "
    fi

ID=$SCR$A11
RES=32
RESX=$RES
RESY=$RES
PP=30

SCALE=-no_scale
#SCALE=-scaled
let "count+=1"
VC=7
VV=`echo "10^($VC)" | bc -l`

NAME="solcx_halfcell_vc${VC}_"
DIR="${NAME}_${RESX}x${RESY}_${BVISC}_${PP}"
OUT="$DIR/${PEN}_$count"
mkdir $DIR >& /dev/null
mkdir $OUT >& /dev/null

#                --components.MatrixTerm.incompressibility_Penalty=$PEN \

# --solCx_xc=0.50390625 \ gives end of material on half cell for res 128. Is 0.5+1/256
$UWEXEC $UWPATH/Underworld/SysTest/PerformanceTests/testVelicSolCx.xml \
    $UWPATH/Solvers/InputFiles/AugLagStokesSLE-NaiNbj.xml \
    $UWPATH/Solvers/InputFiles/kspinterface.xml \
    $UWPATH/Solvers/InputFiles/analyticVis.xml \
    $UWPATH/Solvers/InputFiles/quiet.xml  \
    $MGOP \
                --particlesPerCell=$PP \
  		--outputPath="./$OUT" \
  		--components.stokesEqn.isNonLinear=False \
  		--saveDataEvery=1 --checkpointEvery=2 --checkpointWritePath="./$OUT/Checkpoints" --checkpointAppendStep=1 \
                --components.MatrixTerm.viscosity_weighting=false \
                --components.FieldTest.normaliseByAnalyticSolution=False \
                --solCx_etaA=$VV \
                --solCx_xc=0.50390625 \
                --solCx_etaB=0.001 \
                --solCx_n=2.0 \
                --wavenumberY=2.0 \
                -Xuzawastyle 0 \
                -XscrPCKSP_ksp_type fgmres \
                -XscrPCKSP_ksp_converged_reason \
                -XscrPCKSP_ksp_view \
  		-ksp_k2_type SLE --penaltyNumber=$PEN \
  		-Q22_pc_type $UW \
  		-XQ22_pc_type gtkg -Xrestore_K $SCALE \
                -Xscr_pc_gtkg_ksp_view -Xscr_pc_gtkg_ksp_monitor -scr_pc_gtkg_ksp_rtol 1e-6 -scr_pc_gtkg_ksp_type cg \
  		-Xbuild_cb_pressure_nullspace \
  		-build_const_pressure_nullspace \
  		--mgLevels=5 \
  		-Xscr_ksp_max_it 1000 \
                -scr_ksp_type $SCR \
                -scr_ksp_view \
                $SCRNORMTYPE \
                -Xscr_ksp_left_pc \
  		-scr_ksp_rtol $SCRTOL \
                -Xscr_ksp_monitor_true_residual \
              -XA11_pc_type hypre -XA11_pc_hypre_type boomeramg -XA11_pc_hypre_boomeramg_print_statistics \
  		-A11_ksp_rtol $A11TOL \
                -A11_ksp_type $A11 \
              -XA11_pc_hypre_boomeramg_grid_sweeps_all 5 \
              -XA11_pc_hypre_boomeramg_tol 1e-3 \
                -A11_ksp_converged_reason \
               -XA11_ksp_norm_inf_monitor \
               -XA11_use_norm_inf_stopping_condition \
               -XA11_ksp_monitor_true_residual \
                -XA11_ksp_view \
                -backsolveA11_ksp_type fgmres \
                -backsolveA11_ksp_rtol 1.0e-6 \
  		--elementResI=$RES --elementResJ=$RES \
  		--maxTimeSteps=0 -Xdump_matvec -forcecorrection \
    > "./$OUT/output.txt" 2>&1

./getconv2.pl < "$OUT/output.txt"

echo " $PEN | $UW | $SCALE | $SCR | $SCRTOL | $SCRP  | $A11 | $A11TOL | $RESX | $MG |$DIR | $ID |"

done
done
done
done
done
done

# Mutually incompatible at present: 
# -build_cb_pressure_nullspace  
# -build_const_pressure_nullspace


# Choices
# -ksp_k2_type NULL GG GMG DGMGD SLE

function makepng(){
FF=solcx_halfcell_vc${VC}_${RESX}x${RESY}.txt
F=var.txt
T='SolCx Half Cell Jump\nViscosity Contrast = '$VC'\nmesh '${RESX}'x'${RESY}' eta0=1.0\n2x1 rectangle\nProcs = '$PROCS' SCR PC= '$UW' scaling '$SCALE
M='SolCx_halfcell_'${RESX}'x'${RESY}'_'$VC

cp $F $DIR/$FF
mv $F $DIR/$F
cd $DIR

for cg in cg$A11 fgmres$A11 minres$A11
do

CG=$( echo "$cg" | tr -s  '[:lower:]'  '[:upper:]' )
../getvars.pl $cg < $F > ./${cg}plot.txt
VAR=$(cat <<EOF
set term png size 800,700 font arial 9
set output '${M}_$CG.png'
set format y "%g";
set format x "%g";
set xlabel "Penalty Number"
set style line 6 lt rgb "blue" lw 1
set multiplot layout 2,2 title "$CG $T"
set log xy
unset key
unset logscale y
set style data linespoints
set ylabel "Pressure Iterations"
plot './${cg}plot.txt' u 12:1 ls 6
set ylabel "Total Velocity Iterations"
plot './${cg}plot.txt' u 12:2 ls 6
set ylabel "||G'u+Cp-h||"
set log xy
unset key
plot './${cg}plot.txt' u 12:4 ls 6
set ylabel "Time (secs)"
set log xy
unset key
unset logscale y
plot './${cg}plot.txt' u 12:3 ls 6
unset multiplot

EOF
)

echo "$VAR" > ${cg}.gp
gnuplot < ${cg}.gp

done

cd ..

}

makepng
