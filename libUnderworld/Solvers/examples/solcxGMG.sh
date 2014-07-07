#! /bin/bash

#
#   Solve using BSSCR GG/GMG augmented lagrangian method (optionally).
#   In the solver K <- K+penalty*K2 where K2 is G'*M^(-1)*G.
#   K2 is calculated in the BSSCR solver.
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
# -remove_constant_pressure_null_space  0      // tell BSSCR to check for and remove the constant nullspace in the pressure space.
# -remove_checkerboard_pressure_null_space 0   // tell BSSCR to check for and remove the checker board nullspace in the pressure space.
# -Q22_pc_type $UW                // specify the preconditioner for the Schur solve: one of (uw, uwscale, gtkg)
# -uzawa_style 0                  // set to 1 to imitate uzawa style solve: 0 is default
# -force_correction 1             // subtract force vector terms with penalty of a user-defined Penalty matrix when using SLE option with -ksp_k2_type from system RHS.
# -rescale_equations 0            // scale equations
# -k_scale_only 0                 // only scale the velocity terms in the system of equations
# -augmented_lagrangian 1         // use augmented lagrangian method
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
export UWEXEC="mpirun -n 8 $UWPATH/build/bin/Underworld"

echo "| p its | v its | p solve time | constraint | gperror | NL its | avg P its | minp | maxp | minv | maxv | penalty | -Q22_pc_type | scale | scr | scr tol | scr norm type | A11 | A11 tol |res | MG | DIR | ID | VC |" | tee var.txt
for VC in 3
do
for SC in 0
do
for UW in gkgdiag
do
for SCR in fgmres
do
for A11 in preonly
do
for SCRTOL in 1e-5
do
for A11TOL in 1e-7
do
echo "|-------+-------+------------+----------+------+------+------+------+---------+----------------+-------+-----+---------+---------------+-----+---------+-----+----+----|" | tee -a var.txt
#for PEN in 0.0 0.0001 0.05 0.1 1.0 5.0 10.0 20.0 50.0 100.0 200.0 500.0 1000.0 2000.0
#for PEN in 0.0 0.0001 0.05 0.1 1.0 5.0 10.0
#for PEN in 0.0 0.02 0.1 1.0 2.0 10.0 20.0 100.0 200.0 1000.0
#for PEN in 0.0 10.0 100.0 1000.0 10000.0
#10.0 100.0 1000.0
#for PENEXP in -4 -1 0 1 2 3 4 5
for PENEXP in 0
do
#dividing penalty by 4 to make equivalent to NaiNbj examples
#PEN=`echo "0.25*$PEN" | bc -l`
PEN=`echo "10^($PENEXP)" | bc -l`
#SCRP="unpreconditioned"
SCRP="default"
#SCRP="unpreconditioned"
    if [ "$SCRP" = "unpreconditioned" ]
        then
	SCRNORMTYPE="-scr_ksp_norm_type unpreconditioned"
    else
	SCRNORMTYPE=" "
    fi

#MG=boomeramg
#MG="ml"
MG=gmg
MGOP=" "
    if [ "$MG" = "gmg" ]
        then
	    #MGOP="$UWPATH/Solvers/InputFiles/MultigridForRegularSCR.xml -options_file ./options-scr-mg-accelerating.opt "
	    MGOP="$UWPATH/Solvers/InputFiles/MultigridForRegularSCR.xml "
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
        A11=preonly
        MGOP="-A11_pc_type lu -A11_ksp_type preonly "
    fi
    if [ "$MG" = "mumps" ]
        then
	A11=preonly
        MGOP="-options_file ./options-scr-mumps-petsc3.opt "
    fi

ID=$SCR$A11

RES=$1

RESX=$RES
RESY=$RES
PP=40

SCALETEXT=no_scale

if [ "$SC" = 1 ]                                            
then
    SCALETEXT=scaled
fi

#SC=0
SCALE="-rescale_equations $SC "
#SCALE=' '
#SCALE=-k_scale
#OUT="outputMGBSSCR_GMG_SolCx"
#OUT="TEST"
let "count+=1"
#VC=3
VV=`echo "10^($VC)" | bc -l`
#VCC=0
#VB=`echo "10^(-$VCC)" | bc -l`
VB=1.0
PCRES=15

#NAME="solcxGMG_vc${VC}_${A11TOL}_${SCRTOL}_${SCALE}_${UW}_ppc=${PP}_procs_${PROCS}_${MG}"
#NAME="solcxGMG"
NAME="crunsSolCx_conditionNumberMatrices"
#NAME="penTest"
DIR="${NAME}_${RESX}x${RESY}"
OUT="$DIR/cx_10e${VC}_${SCALETEXT}"
mkdir $DIR >& /dev/null
mkdir $OUT >& /dev/null

$UWEXEC $UWPATH/Underworld/SysTest/PerformanceTests/testVelicSolCx.xml \
    $UWPATH/Solvers/InputFiles/AugLagStokesSLE-GtMG.xml \
    $UWPATH/Solvers/InputFiles/VelocityMassMatrixSLE.xml \
    $UWPATH/Solvers/InputFiles/kspinterface.xml \
    $UWPATH/Solvers/InputFiles/analyticVis.xml \
    $UWPATH/Solvers/InputFiles/quiet.xml  \
    $MGOP \
                --particlesPerCell=$PP \
                --components.weights.resolutionX=$PCRES --components.weights.resolutionY=$PCRES --components.weights.resolutionZ=$PCRES \
  		--outputPath="./$OUT" \
  		--components.stokesEqn.isNonLinear=False \
       		--saveDataEvery=1 --checkpointEvery=2 --checkpointWritePath="./$OUT/Checkpoints" --checkpointAppendStep=1 \
                --components.FieldTest.normaliseByAnalyticSolution=False \
                --solCx_etaA=$VV \
                --solCx_xc=0.5 \
                --solCx_etaB=$VB \
                --solCx_n=2.0 \
                --wavenumberY=2.0 \
                -scr_ksp_set_min_it_converge 1 \
                -force_correction 1 -k_scale_only 1 \
                -uzawastyle 0 \
                -scrPCKSP_ksp_type fgmres \
                -XscrPCKSP_ksp_converged_reason \
                -XscrPCKSP_ksp_view \
  		-ksp_type bsscr -pc_type none -ksp_k2_type GMG -augmented_lagrangian 1 --penaltyNumber=$PEN \
  		-Q22_pc_type $UW \
  		-XQ22_pc_type gtkg -Xrestore_K $SCALE \
                -Xscr_pc_gtkg_ksp_view -Xscr_pc_gtkg_ksp_monitor -scr_pc_gtkg_ksp_rtol 1e-6 -scr_pc_gtkg_ksp_type cg \
  		-remove_checkerboard_pressure_null_space 0 \
  		-remove_constant_pressure_null_space 1 \
  		--mgLevels=5 \
  		-Xscr_ksp_max_it 1000 \
                -scr_ksp_type $SCR \
                -Xscr_ksp_view \
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
  		--maxTimeSteps=0 -XA11_ksp_view -XA11_mg_levels_ksp_view \
    -Xdump_matvec -matsuffix "_${RES}x${RES}_${SCALETEXT}_10e${VC}_cx_" -matdumpdir $OUT -solutiondumpdir $OUT \
   --components.stokesblockkspinterface.OptionsString="-A11_ksp_type cg -A11_ksp_rtol 1e-3 -A11_ksp_monitor -backsolveA11_ksp_type cg -log_summary" \

#    > "./$OUT/output.txt" 2>&1

#./getconv2.pl < "$OUT/output.txt"  | tee -a var.txt

cp "$OUT/window.00000.png" "$PNGNAME"
#echo " $PEN | $UW | $SC | $SCR | $SCRTOL | $SCRP  | $A11 | $A11TOL | $RES | $MG | $DIR | $ID |" | tee -a  var.txt

done
done
done
done
done
done
done
done

#               < /dev/null 2> "./$OUT/errors.txt" /dev/null 1> "./$OUT/output.txt"
#                -Xscr_ksp_norm_type natural \ <-- not supported when using KSPRelativeRhsConvergenced()

#    
#
#    $UWPATH/Solvers/InputFiles/MultigridForRegularSCR.xml \
#    -options_file ./options-scr-mg-accelerating.opt \

# Mutually incompatible at present: 
# -build_cb_pressure_nullspace  
# -build_const_pressure_nullspace


# Choices
# -ksp_k2_type NULL GG GMG DGMGD SLE
function makepng(){
FF=solcx_halfcell_vc${VC}_${RESX}x${RESY}.txt
F=var.txt
T='SolCx Half Cell Jump (GMG method)\nViscosity Contrast = '$VC'\nmesh '${RESX}'x'${RESY}' eta0=1.0\n2x1 rectangle\nProcs = '$PROCS' SCR PC= '$UW' scaling '$SCALE' Velocity Solve '$MG'\nA11TOl='$A11TOL' SCRTOL='$SCRTOL
M='SolCx_halfcell_'${RESX}'x'${RESY}'_'$VC

cp $F $DIR/$FF
mv $F $DIR/$F
cd $DIR

#for cg in cg$A11 fgmres$A11 minres$A11
for cg in fgmres$A11
do

#CG=${cg^^}
CG=$( echo "$cg" | tr -s  '[:lower:]'  '[:upper:]' )
../getvars.pl $cg < $F > ./${cg}plot.txt

VAR=$(cat <<EOF
set term png size 800,700 font arial 9
set output '${M}_$CG.png'
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

#makepng

