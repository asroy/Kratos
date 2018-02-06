#if !defined(KRATOS_OVERSET_TRILINOS_BLOCK_BUILDER_AND_SOLVER_PERIODIC_H )
#define  KRATOS_OVERSET_TRILINOS_BLOCK_BUILDER_AND_SOLVER_PERIODIC_H

#include "trilinos_application/custom_strategies/builder_and_solvers/trilinos_block_builder_and_solver_periodic.h"
#include "custom_utilities/OversetAssembler.h"

namespace Kratos
{

template<class TSparseSpace,
         class TDenseSpace,
         class TLinearSolver>
class OversetTrilinosBlockBuilderAndSolverPeriodic
    :   public TrilinosBlockBuilderAndSolverPeriodic< TSparseSpace,TDenseSpace, TLinearSolver >
{
    using BaseType = TrilinosBlockBuilderAndSolverPeriodic<TSparseSpace,TDenseSpace, TLinearSolver >;

public:
    KRATOS_CLASS_POINTER_DEFINITION( OversetTrilinosBlockBuilderAndSolverPeriodic );

    OversetTrilinosBlockBuilderAndSolverPeriodic
    (
        Epetra_MpiComm& Comm,
        int guess_row_size,
        typename TLinearSolver::Pointer pNewLinearSystemSolver,
        OversetAssembly::OversetAssembler & r_overset_assembler
    )
        :   BaseType(Comm, guess_row_size,pNewLinearSystemSolver),
            mrOversetAssembler{r_overset_assembler}
    {}

    ~OversetTrilinosBlockBuilderAndSolverPeriodic() override
    {}

    void SetUpSystem(ModelPart& r_model_part) override
    {

        std::cout<<"inside OversetTrilinosBlockBuilderAndSolverPeriodic::"<<__func__<<std::endl;
        
        if( & r_model_part != & (mrOversetAssembler.rModelPart()) )
        {
            std::cout<<__func__<<"wrong! not the same model part!"<<std::endl;
            exit(EXIT_FAILURE);
        }

        //call base type method
        BaseType::SetUpSystem(r_model_part);

        mrOversetAssembler.GetOversetConditionsDonorEquationsId();
    }
    
protected:
    OversetAssembly::OversetAssembler & mrOversetAssembler;

};

}//namespace Kratos
#endif