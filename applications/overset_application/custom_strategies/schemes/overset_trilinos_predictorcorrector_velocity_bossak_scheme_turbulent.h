#if !defined(KRATOS_OVERSET_TRILINOS_PREDICTOR_CORRECTOR_VELOCITY_BOSSAK_SCHEME_TURBULENT )
#define  KRATOS_OVERSET_TRILINOS_PREDICTOR_CORRECTOR_VELOCITY_BOSSAK_SCHEME_TURBULENT

#include "trilinos_application/custom_strategies/schemes/trilinos_predictorcorrector_velocity_bossak_scheme_turbulent.h"
#include "custom_utilities/OversetAssembler.h"

namespace Kratos
{

template<class TSparseSpace,
         class TDenseSpace >
class OversetTrilinosPredictorCorrectorVelocityBossakSchemeTurbulent
    :   public TrilinosPredictorCorrectorVelocityBossakSchemeTurbulent <TSparseSpace,TDenseSpace>
{
    using BaseType = TrilinosPredictorCorrectorVelocityBossakSchemeTurbulent <TSparseSpace,TDenseSpace>;

public:
    KRATOS_CLASS_POINTER_DEFINITION( OversetTrilinosPredictorCorrectorVelocityBossakSchemeTurbulent );

    OversetTrilinosPredictorCorrectorVelocityBossakSchemeTurbulent (double NewAlphaBossak, double MoveMeshStrategy, unsigned int DomainSize, OversetAssembly::OversetAssembler & r_overset_assembler)
        :   BaseType( NewAlphaBossak, MoveMeshStrategy, DomainSize ),
            mrOversetAssembler{r_overset_assembler}
    {}

    OversetTrilinosPredictorCorrectorVelocityBossakSchemeTurbulent (double NewAlphaBossak, double MoveMeshStrategy, unsigned int DomainSize, Process::Pointer pTurbulenceModel, OversetAssembly::OversetAssembler & r_overset_assembler)
        :   BaseType( NewAlphaBossak, MoveMeshStrategy, DomainSize, pTurbulenceModel ),
            mrOversetAssembler{r_overset_assembler}
    {}

    OversetTrilinosPredictorCorrectorVelocityBossakSchemeTurbulent (double NewAlphaBossak, double MoveMeshStrategy, unsigned int DomainSize, const Variable<int>& rPeriodicIdVar, OversetAssembly::OversetAssembler & r_overset_assembler)
        :   BaseType( NewAlphaBossak, MoveMeshStrategy, DomainSize, rPeriodicIdVar),
            mrOversetAssembler{r_overset_assembler}
    {}

    ~OversetTrilinosPredictorCorrectorVelocityBossakSchemeTurbulent() override
    {}

    void Update( ModelPart & r_model_part, typename BaseType::DofsArrayType & rDofSet, typename BaseType::TSystemMatrixType & A, typename BaseType::TSystemVectorType & Dx, typename BaseType::TSystemVectorType & b) override
    {
        std::cout<<"inside OversetTrilinosPredictorCorrectorVelocityBossakSchemeTurbulent::"<<__func__<<std::endl;
        
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