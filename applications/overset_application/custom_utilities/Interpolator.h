#if !defined(KRATOS_OVERSET_INTERPOLATOR_H_INCLUDED )
#define  KRATOS_OVERSET_INTERPOLATOR_H_INCLUDED

#include <stdio.h>
#include <string>
#include <vector>
#include <cmath>

#include "overset_application.h"

#include "custom_utilities/DistributedAssignment.h"
#include "custom_utilities/InterpolationInput.h"
#include "custom_conditions/HingeData.h"

namespace Kratos
{
namespace OversetAssembly
{

class Interpolator
{
public:
    using Location = DistributedAssignment::Communication::MpiLocation;
    using InterpolatorKey = DistributedAssignment::DistributedAssignment::DistributedKey<Location>;

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

    void Execute( const InterpolationInput & r_interpolation_input, HingeData & r_hinge_data )
    {
        const ModelPart::ElementType::IndexType element_id = r_interpolation_input.mElementId;

        const ModelPart::ElementType::GeometryType & r_geometry = const_cast<ModelPart &> (mrModelPart).GetElement(element_id).GetGeometry();

        //check
        {
            std::cout<<__func__<<": donor_node_id: "<<std::endl;
            
            std::size_t i = 0;
            for( const auto & r_node : r_geometry )
            {
                std::cout<<r_node.GetId()<<", "<<r_interpolation_input.mNodesId[i]<<std::endl;

                if( r_node.GetId() != r_interpolation_input.mNodesId[i] )
                {
                    std::cout<<__func__<<"wrong! element_id and Nodes_id not match!"<<std::endl;
                    exit(EXIT_FAILURE);
                }

                i++;
            }
        }

        Vector Ns(r_geometry.size());

        Point<3> point{ r_interpolation_input.mBarycentricCoordinate[0],
                        r_interpolation_input.mBarycentricCoordinate[1],
                        r_interpolation_input.mBarycentricCoordinate[2] };

        r_geometry.ShapeFunctionsValues(Ns, point);

        Vector r_coordinate(3);
        noalias(r_coordinate) = ZeroVector(3);
        for(std::size_t i = 0; i < r_geometry.size(); i++ )
            noalias(r_coordinate) += Ns[i]*r_geometry[i];

        r_hinge_data.mCoordinate[0] = r_coordinate[0];
        r_hinge_data.mCoordinate[1] = r_coordinate[1];
        r_hinge_data.mCoordinate[2] = r_coordinate[2];
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
