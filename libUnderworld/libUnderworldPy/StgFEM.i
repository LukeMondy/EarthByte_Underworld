/* -*- C -*-  (not really, but good for syntax highlighting) */

%module StgFEM


%{
/* Includes the header in the wrapper code */
#define SWIG_FILE_WITH_INIT
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
%}

%include "numpy.i"

%init %{
import_array();
%}

%import "StGermain.i"
%import "StgDomain.i"

%include "Discretisation/FeVariable.h"
%include "Discretisation/Element.h"
%include "Discretisation/ElementType.h"
%include "Discretisation/ElementType_Register.h"
%include "Discretisation/AnalyticFeVariable.h"
%include "Discretisation/AnalyticSolution.h"
%include "Discretisation/BilinearElementType.h"
%include "Discretisation/BilinearInnerElType.h"
%include "Discretisation/Biquadratic.h"
%include "Discretisation/C0Generator.h"
%include "Discretisation/C2Generator.h"
%include "Discretisation/ConstantElementType.h"
%include "Discretisation/ErrorFeVariable.h"
%include "Discretisation/FeEquationNumber.h"
%include "Discretisation/FeInterface.h"
%include "Discretisation/FeMesh.h"
%include "Discretisation/FeMesh_Algorithms.h"
%include "Discretisation/FeMesh_ElementType.h"
%include "Discretisation/FeSwarmVariable.h"
%include "Discretisation/FieldTest.h"
%include "Discretisation/FunctionSuite.h"
%include "Discretisation/Inner2DGenerator.h"
%include "Discretisation/LinearElementType.h"
%include "Discretisation/LinearTriangleElementType.h"
%include "Discretisation/LinkedDofInfo.h"
%include "Discretisation/OperatorFeVariable.h"
/*%include "Discretisation/P1Generator.h"*/
%include "Discretisation/PETScErrorChecking.h"
%include "Discretisation/ShapeFeVariable.h"
%include "Discretisation/shortcuts.h"
%include "Discretisation/TrilinearElementType.h"
%include "Discretisation/TrilinearInnerElType.h"
%include "Discretisation/Triquadratic.h"
%include "Discretisation/types.h"
%include "Discretisation/units.h"
%include "SLE/SystemSetup/Assembler.h"
%include "SLE/SystemSetup/EntryPoint.h"
%include "SLE/SystemSetup/FiniteElementContext.h"
%include "SLE/SystemSetup/ForceTerm.h"
%include "SLE/SystemSetup/MGOpGenerator.h"
%include "SLE/SystemSetup/MultigridSolver.h"
%include "SLE/SystemSetup/PETScMGSolver.h"
%include "SLE/SystemSetup/shortcuts.h"
%include "SLE/SystemSetup/SLE_Solver.h"
%include "SLE/SystemSetup/SolutionVector.h"
%include "SLE/SystemSetup/SROpGenerator.h"
%include "SLE/SystemSetup/StiffnessMatrix.h"
%include "SLE/SystemSetup/StiffnessMatrixTerm.h"
%include "SLE/SystemSetup/SystemLinearEquations.h"
%include "SLE/SystemSetup/ForceVector.h"
%include "SLE/SystemSetup/types.h"
%include "SLE/SystemSetup/units.h"
%include "SLE/ProvidedSystems/AdvectionDiffusion/AdvectionDiffusionSLE.h"
%include "SLE/ProvidedSystems/AdvectionDiffusion/LumpedMassMatrixForceTerm.h"
%include "SLE/ProvidedSystems/AdvectionDiffusion/Multicorrector.h"
%include "SLE/ProvidedSystems/AdvectionDiffusion/Residual.h"
%include "SLE/ProvidedSystems/AdvectionDiffusion/Timestep.h"
%include "SLE/ProvidedSystems/AdvectionDiffusion/types.h"
%include "SLE/ProvidedSystems/AdvectionDiffusion/UpwindParameter.h"
%include "SLE/ProvidedSystems/Energy/Energy_SLE.h"
%include "SLE/ProvidedSystems/Energy/Energy_SLE_Solver.h"
%include "SLE/ProvidedSystems/Energy/types.h"
%include "SLE/ProvidedSystems/StokesFlow/Stokes_SLE.h"
%include "SLE/ProvidedSystems/StokesFlow/Stokes_SLE_PenaltySolver.h"
%include "SLE/ProvidedSystems/StokesFlow/Stokes_SLE_UzawaSolver.h"
%include "SLE/ProvidedSystems/StokesFlow/types.h"
%include "SLE/ProvidedSystems/StokesFlow/UpdateDt.h"
%include "SLE/ProvidedSystems/StokesFlow/UzawaPreconditionerTerm.h"
%include "SLE/types.h"
%include "Assembly/DivergenceMatrixTerm.h"
%include "Assembly/GradientStiffnessMatrixTerm.h"
%include "Assembly/IsoviscousStressTensorTerm.h"
%include "Assembly/LaplacianStiffnessMatrixTerm.h"
%include "Assembly/MassMatrixTerm.h"
%include "Assembly/shortcuts.h"
%include "Assembly/ThermalBuoyancyForceTerm.h"
%include "Assembly/types.h"
%include "Assembly/units.h"
%include "Utils/IrregularMeshParticleLayout.h"
%include "Utils/SemiLagrangianIntegrator.h"
%include "Utils/types.h"


/* #  the following extends the FeVariable class such that we can extract a numpy style array via numpy.i magic */
%extend FeVariable
{

  void getAsNumpyArray( double** ARGOUTVIEW_ARRAY2, int* DIM1, int* DIM2 ){

     Variable* variable = $self->dofLayout->baseVariables[0];
     if(variable->parent != NULL)
        variable = variable->parent;
     *ARGOUTVIEW_ARRAY2 = variable->arrayPtr;
     *DIM1 = variable->arraySize;
     *DIM2 = variable->subVariablesCount;

  }

}

%extend FeMesh
{

  void getAsNumpyArray( double** ARGOUTVIEW_ARRAY2, int* DIM1, int* DIM2 ){

	Sync* sync = (Sync*)IGraph_GetDomain( (IGraph*)$self->topo, 0 );
	*DIM1 = Sync_GetNumDomains( sync );
	*DIM2 = $self->topo->nDims;

	/* need to offset the value because in stg it is allocated as a 2d array. note that it is contiguous. */
    *ARGOUTVIEW_ARRAY2 = (Pointer)( $self->vertices ) ;

  }

}


%pythoncode %{
	def SquareBracketsFinal(self, index):
		return self.numpyArray[index]

	def SquareBracketsPre(self, index):
	    # this routines setups array first
		self.numpyArray = self.getAsNumpyArray()
		# now that its setup, switch this overload to direct function.  this acts on the instance, not the class, aka monkeypatching
		import new
		self.GetIndexFunc = new.instancemethod(SquareBracketsFinal, self, None)
		return SquareBracketsFinal(self, index)

	def SquareBracketsFixed(self,index):
		return self.GetIndexFunc(index)

	FeVariable.GetIndexFunc = SquareBracketsPre
	FeVariable.__getitem__ = SquareBracketsFixed

	FeMesh.GetIndexFunc = SquareBracketsPre
	FeMesh.__getitem__ = SquareBracketsFixed
%}

/*

# add this to the FeVariable constructor
# %feature("pythonappend") FeVariable %{
#  #do something after C++ call
#  self.numpyArray = self.getAsNumpyArray()
# %}
*/