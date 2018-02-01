#if !defined(KRATOS_OVERSET_INTERPOLATION_METHOD_H_INCLUDED )
#define  KRATOS_OVERSET_INTERPOLATION_METHOD_H_INCLUDED

// External includes 
#include "boost/smart_ptr.hpp"

#include <mpi.h>
#include <iostream>
#include <unistd.h>
#include <vector>

#include "custom_utilities/DistributedAssignment.h"
#include "custom_utilities/Interpolator.h"
#include "custom_utilities/InterpolationInput.h"
#include "custom_utilities/InterpolationOutput.h"

namespace Kratos
{
namespace OversetAssembly
{

class InterpolationMethod
{
public:
    using OversetCommunicator = DistributedAssignment::Communication::MpiCommunicator;
    using Location = OversetCommunicator::Location;
    using KeyIssuer = DistributedAssignment::DistributedAssignment::DistributedKeyIssuer<Location>;
    using Key = KeyIssuer::Key;

    using DummyModelPartHolder = DistributedAssignment::DistributedAssignment::DummyContractor<Key>;
    using DummyModelPartHolderManager = DistributedAssignment::DistributedAssignment::DistributedContractorManager<DummyModelPartHolder,OversetCommunicator,DistributedAssignment::DistributedAssignment::DistributedKeyIssuer>;
    
    using InterpolatorManager = DistributedAssignment::DistributedAssignment::DistributedContractorManager<Interpolator,OversetCommunicator,DistributedAssignment::DistributedAssignment::DistributedKeyIssuer>;

    using InterpolationAssignmentManager = DistributedAssignment::DistributedAssignment::DistributedAssignmentManager<DummyModelPartHolder,Interpolator,InterpolationInput,InterpolationOutput,OversetCommunicator,DistributedAssignment::DistributedAssignment::DistributedKeyIssuer,DistributedAssignment::DistributedAssignment::DistributedKeyIssuer>;

    using InterpolationAssignmentInputData  = typename InterpolationAssignmentManager::template AssignmentDataType<InterpolationInput>;
    using InterpolationAssignmentOutputData = typename InterpolationAssignmentManager::template AssignmentDataType<InterpolationOutput>;

    using DoubleVariable = Variable<double>;
    using Array1dComponentVariable = VariableComponent<VectorComponentAdaptor<array_1d<double,3>>>;

public:
    InterpolationMethod() = delete;

    InterpolationMethod( OversetCommunicator & r_communicator, DummyModelPartHolderManager & r_dummy_model_part_holder_manager, const ModelPart & r_model_part )
        :   mrOversetCommunicator{r_communicator},
            mrDummyModelPartHolderManager{r_dummy_model_part_holder_manager},
            mpInterpolatorManager{nullptr},
            mpInterpolationAssignmentManager{nullptr}
    {
        //interpolator and manager
        mpInterpolatorManager = new InterpolatorManager{mrOversetCommunicator};
        mpInterpolatorManager->RegisterLocalContractor( * (new Interpolator{r_model_part}), "Interpolator" );
        mpInterpolatorManager->GenerateGlobalContractorsRegistry();

        //search assignment manager
        mpInterpolationAssignmentManager = new InterpolationAssignmentManager{r_communicator, mrDummyModelPartHolderManager, * mpInterpolatorManager};
    }

    virtual ~InterpolationMethod()
    {
        //interpolators and manager
        for( typename InterpolatorManager::ContractorPointerPairByContractorKey & r_interpolator_pair : mpInterpolatorManager->LocalContractorsPointer() )
            delete r_interpolator_pair.second;

        delete mpInterpolatorManager;

        //assignment mananger
        delete mpInterpolationAssignmentManager;
    }

    void AddVariableNeedEquationId( const DoubleVariable & r_variable )
    {
        for( const auto & r_local_interpolator_pointer_paired_by_key : mpInterpolatorManager->LocalContractorsPointer() )
            r_local_interpolator_pointer_paired_by_key.second -> AddVariableNeedEquationId(r_variable);
    }

    void AddVariableNeedValue( const DoubleVariable & r_variable )
    {
        for( const auto & r_local_interpolator_pointer_paired_by_key : mpInterpolatorManager->LocalContractorsPointer() )
            r_local_interpolator_pointer_paired_by_key.second -> AddVariableNeedValue(r_variable);
    }

    void AddVariableNeedDX( const DoubleVariable & r_variable )
    {   
        for( const auto & r_local_interpolator_pointer_paired_by_key : mpInterpolatorManager->LocalContractorsPointer() )
            r_local_interpolator_pointer_paired_by_key.second -> AddVariableNeedDX(r_variable);
    }

    void AddVariableNeedEquationId( const Array1dComponentVariable & r_variable )
    {
        for( const auto & r_local_interpolator_pointer_paired_by_key : mpInterpolatorManager->LocalContractorsPointer() )
            r_local_interpolator_pointer_paired_by_key.second -> AddVariableNeedEquationId(r_variable);
    }

    void AddVariableNeedValue( const Array1dComponentVariable & r_variable )
    {
        for( const auto & r_local_interpolator_pointer_paired_by_key : mpInterpolatorManager->LocalContractorsPointer() )
            r_local_interpolator_pointer_paired_by_key.second -> AddVariableNeedValue(r_variable);
    }

    void AddVariableNeedDX( const Array1dComponentVariable & r_variable )
    {   
        for( const auto & r_local_interpolator_pointer_paired_by_key : mpInterpolatorManager->LocalContractorsPointer() )
            r_local_interpolator_pointer_paired_by_key.second -> AddVariableNeedDX(r_variable);
    }

    void ClearAllInterpolationAssignments()
    {
        mpInterpolationAssignmentManager->ClearAllAssignments();
    }

    Key AddInterpolationAssignment( const Key & r_interpolator_key, const InterpolationInput & r_input )
    {
        const Key & r_dummy_model_part_holder_key = *(mrDummyModelPartHolderManager.LocalContractorsKey().begin());

        return mpInterpolationAssignmentManager->AddAssignment( r_dummy_model_part_holder_key, r_interpolator_key, r_input );
    }

    void ExecuteAllInterpolationAssignments()
    {
        mpInterpolationAssignmentManager->ExecuteAllDistributedAssignments();
    }

    void GetInterpolationResults( std::vector<InterpolationAssignmentOutputData> & r_outputs )
    {
        mpInterpolationAssignmentManager->GetResultsAtAssignor( r_outputs );
    }

private:
    OversetCommunicator & mrOversetCommunicator;
    DummyModelPartHolderManager & mrDummyModelPartHolderManager;
    InterpolatorManager * mpInterpolatorManager;
    InterpolationAssignmentManager * mpInterpolationAssignmentManager;
};

}//namespace OversetAssemlby
}//namespace Kratos
#endif