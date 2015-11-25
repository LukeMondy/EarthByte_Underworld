#! /bin/bash

#
#   Solve using UZAWA.

count=0
PROCS=1

export UWPATH=`./getUWD.sh`
export UWEXEC="cgdb --args $UWPATH/build/bin/Underworld"
#export UWEXEC="$UWPATH/build/bin/Underworld"

#echo "| p its | v its | p solve time | constraint | gperror | NL its | avg P its | minp | maxp | minv | maxv | penalty | -Q22_pc_type | scale | scr | scr tol | scr norm type | A11 | A11 tol |res | MG | DIR | ID |" | tee var.txt

for UW in uwscale
do
for SCR in cg
do
for A11 in fgmres
do
for SCRTOL in 1e-10
do
for A11TOL in 1e-9
do
#echo "|-------+-------+------------+----------+------+------+------+------+---------+----------------+-------+-----+---------+---------------+-----+---------+-----+----|" | tee -a var.txt

       $UWPATH/Experimental/InputFiles/StokesPIC+AugmentedLagrangianUZ.xml \
       $UWPATH/Experimental/InputFiles/Experimental_Components/MultigridForRegularUzawa.xml \
       -options_file $UWPATH/Experimental/InputFiles/options/options-uzawa-mg.opt \
       --components.constitutiveMatrix.incompressibility_Penalty=10.0 \
       --components.constitutiveMatrix.viscosity_weighting=false \

#MG=boomeramg
#MG="ml"
MG=lu
MGOP=" "
    if [ "$MG" = "lu" ]
        then
        A11=preonly
        MGOP="-Uzawa_velSolver_pc_type lu -Uzawa_velSolver_ksp_type preonly "
    fi
    if [ "$MG" = "mumps" ]
        then
	A11=preonly
        MGOP="-options_file ./options-scr-mumps-petsc3.opt "
    fi

ID=$SCR$A11
RES=128
RESX=$RES
RESY=$RES
PP=40

let "count+=1"
VC=9
VV=`echo "10^($VC)" | bc -l`
VC=0
VB=`echo "10^(-$VC)" | bc -l`
PCRES=15


NAME="XsolkxUZ"
DIR="${NAME}_${RESX}x${RESY}_${BVISC}_${PP}"
OUT="$DIR/x_$count"
mkdir $DIR >& /dev/null
mkdir $OUT >& /dev/null

$UWEXEC $UWPATH/Solvers/InputFiles/testVelicSolKx.xml \
    $UWPATH/Solvers/InputFiles/analyticVis.xml \
    $UWPATH/Solvers/InputFiles/quiet.xml  \
    $MGOP \
                -NN $OUT \
                --particlesPerCell=$PP \
  		--outputPath="./$OUT" \
  		--components.stokesEqn.isNonLinear=False \
  		--saveDataEvery=1 --checkpointEvery=1 --checkpointWritePath="./$OUT/Checkpoints" --checkpointAppendStep=1 \
                --components.FieldTest.normaliseByAnalyticSolution=False \
                --solKx_n=3.0   --wavenumberX=3.0 \
                --solKx_m=3.0   --wavenumberY=3.0 \
                --solKx_sigma=1.0 \
                --solKx_twiceB=1.0 \
  		--mgLevels=3 \
  		-A11_ksp_rtol $A11TOL \
                -A11_ksp_type $A11 \
                -A11_ksp_converged_reason \
  		--elementResI=$RES --elementResJ=$RES \
  		--maxTimeSteps=0 -Xdump_matvec -matsuffix "_${RES}x${RES}_${VV}_" \

#    > "./$OUT/output.txt" 2>&1

#./getconv2.pl < "$OUT/output.txt"  | tee -a var.txt

#echo " $PEN | $UW | $SCALE | $SCR | $SCRTOL | $SCRP  | $A11 | $A11TOL | $RES | $MG | $DIR | $ID |" | tee -a  var.txt

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

for cg in cg$A11 fgmres$A11 minres$A11
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

##makepng
