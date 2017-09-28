#if !defined(KRATOS_OVERSET_INTERPOLATOR_H_INCLUDED )
#define  KRATOS_OVERSET_INTERPOLATOR_H_INCLUDED

#include <stdio.h>
#include <string>
#include <vector>
#include <cmath>

#include "overset_application.h"

#include "custom_utilities/DistributedAssignment.h"
#include "custom_utilities/InterpolationInput.h"
#include "custom_utilities/InterpolationOutput.h"

namespace Kratos
{
namespace OversetAssembly
{

class Interpolator
{
public:
    using Location = DistributedAssignment::Communication::MpiLocation;
    using InterpolatorKey = DistributedAssignment::DistributedAssignment::DistributedKey<Location>;

    using PointType = Point<3>;

public:
    Interpolator() = delete;

    Interpolator(const ModelPart & r_model_part)
        :   mrModelPart{r_model_part}
    {}

    virtual ~Interpolator()
    {}

    void SetName( const std::string & r_name )
    { mName = r_name; }

    std::string Name() const
    { return mName; }
    
    void SetKey(const InterpolatorKey key)
    { mKey = key; }
    
    InterpolatorKey Key() const
    { return mKey; }

    void Execute( const InterpolationInput & r_input, InterpolationOutput & r_output )
    {
        const ModelPart::ElementType::IndexType element_id = r_input.mElementId;

        const ModelPart::ElementType & r_element = const_cast<ModelPart &> (mrModelPart).GetElement(element_id);
        
        const ModelPart::ElementType::GeometryType & r_geometry = r_element.GetGeometry();

        //check if node_id match
        {
            std::cout<<__func__<<": donor_node_id: "<<std::endl;
            
            std::size_t i = 0;
            for( const auto & r_node : r_geometry )
            {
                std::cout<<r_node.GetId()<<", "<<r_input.mNodesId[i]<<std::endl;

                if( r_node.GetId() != r_input.mNodesId[i] )
                {
                    std::cout<<__func__<<"wrong! element_id and Nodes_id not match!"<<std::endl;
                    exit(EXIT_FAILURE);
                }

                i++;
            }
        }

        //
        const PointType point{ r_input.mBarycentricCoordinate[0],
                               r_input.mBarycentricCoordinate[1],
                               r_input.mBarycentricCoordinate[2] };

        r_output.InterpolateDataFromElement( point, r_element );
    }

private:
    void Print( const DistributedAssignment::DataUtility::DataPrinter & r_printer ) const
    {
        std::cout << "{Interpolator: ";
        std::cout << "{Name: "<< mName <<"},",
        std::cout << "{ContractorKey: ";
        r_printer.Print(mKey);
        std::cout << "},";
        std::cout << "},";
    }

private:
    std::string mName;
    InterpolatorKey mKey;

    const ModelPart & mrModelPart;

    friend class DistributedAssignment::DataUtility::DataPrinter;
};

}//namespace OversetAssembly
}//namespace Kratos
#endif
