#if !defined(KRATOS_OVERSET_TRILINOS_RESIDUAL_BASED_INCREMENTAL_UPDATE_STATIC_SCHEME )
#define  KRATOS_OVERSET_TRILINOS_RESIDUAL_BASED_INCREMENTAL_UPDATE_STATIC_SCHEME

#include "trilinos_application/custom_strategies/schemes/trilinos_residualbased_incrementalupdate_static_scheme.h"
#include "custom_utilities/OversetAssembler.h"

namespace Kratos
{

template<class TSparseSpace,
         class TDenseSpace >
class OversetTrilinosResidualBasedIncrementalUpdateStaticScheme
    :   public TrilinosResidualBasedIncrementalUpdateStaticScheme <TSparseSpace,TDenseSpace>
{
    using BaseType = TrilinosResidualBasedIncrementalUpdateStaticScheme <TSparseSpace,TDenseSpace>;

public:
    KRATOS_CLASS_POINTER_DEFINITION( OversetTrilinosResidualBasedIncrementalUpdateStaticScheme );

    OversetTrilinosResidualBasedIncrementalUpdateStaticScheme (OversetAssembly::OversetAssembler & r_overset_assembler)
        :   BaseType(),
            mrOversetAssembler{r_overset_assembler}
    {}

    ~OversetTrilinosResidualBasedIncrementalUpdateStaticScheme() override
    {}

    void Update(ModelPart & r_model_part,
                typename BaseType::DofsArrayType & rDofSet,
                typename BaseType::TSystemMatrixType & A,
                typename BaseType::TSystemVectorType & Dx,
                typename BaseType::TSystemVectorType & b) override
    {
        std::cout<<"inside OversetTrilinosResidualBasedIncrementalUpdateStaticScheme::"<<__func__<<std::endl;
        
        if( & r_model_part != & (mrOversetAssembler.rModelPart()) )
        {
            std::cout<<__func__<<"wrong! not the same model part!"<<std::endl;
            exit(EXIT_FAILURE);
        }

        //call base type method
        BaseType::Update(r_model_part, rDofSet, A, Dx, b);

        //update hinge donor data
        mrOversetAssembler.InterpolateHingesDonorData();
    }
    
protected:
    OversetAssembly::OversetAssembler & mrOversetAssembler;
};

}//namespace Kratos
#endif