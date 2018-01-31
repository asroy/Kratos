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

    using PointType = Point;

    using DoubleVariable = Variable<double>;
    using Array1dComponentVariable = VariableComponent<VectorComponentAdaptor<array_1d<double,3>>>;
    using VariableKeyType = typename VariableData::KeyType;

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

    void AddVariableNeedEquationId( const DoubleVariable & r_variable )
    {
        mDoubleVariablesNeedEquationId.insert(r_variable);
    }

    void AddVariableNeedValue( const DoubleVariable & r_variable )
    {
        mDoubleVariablesNeedValue.insert(r_variable);
    }

    void AddVariableNeedDX( const DoubleVariable & r_variable )
    {   
        mDoubleVariablesNeedDX.insert(r_variable);
    }

    void Execute( const InterpolationInput & r_input, InterpolationOutput & r_output )
    {
        const ModelPart::ElementType::IndexType element_id = r_input.mElementId;

        const ModelPart::ElementType & r_element = const_cast<ModelPart &> (mrModelPart).GetElement(element_id);
        
        const ModelPart::ElementType::GeometryType & r_geometry = r_element.GetGeometry();

        //check if node_id match
        {
            std::size_t i = 0;
            for( const auto & r_node : r_geometry )
            {
                if( r_node.GetId() != r_input.mNodesId[i] )
                {
                    std::cout<<__func__<<"wrong! element_id and Nodes_id not match!"<<std::endl;
                    exit(EXIT_FAILURE);
                }
                i++;
            }
        }

        const PointType point{ r_input.mBarycentricCoordinate[0],
                               r_input.mBarycentricCoordinate[1],
                               r_input.mBarycentricCoordinate[2] };

        InterpolateVariablesFromElement( r_output, point, r_element );
    }

    void InterpolateVariablesFromElement( InterpolationOutput & r_output, const Point & r_point, const ModelPart::ElementType & r_element )
    {
        //
        r_output.Clear();

        //
        const ModelPart::ElementType::GeometryType & r_geometry = r_element.GetGeometry();
        const std::size_t num_node = r_geometry.size();

        //shape functions value and derivative
        Vector Ns(num_node);
        r_geometry.ShapeFunctionsValues(Ns, r_point);

        Matrix DNs_DEs(num_node,3);
        r_geometry.ShapeFunctionsLocalGradients(DNs_DEs, r_point);

        //Jacbian matrix
        Matrix jinv(3,3);
        r_geometry.InverseOfJacobian( jinv, r_point );

        //shape functions Gradient
        Matrix DNs_DXs = prod(DNs_DEs, jinv);

        //
        r_output.mNumberOfDonorNodes = num_node;

        r_output.mNs.clear();
        r_output.mDNsDXs.clear();

        r_output.mNs.resize(num_node);
        r_output.mDNsDXs.resize(num_node);

        for(std::size_t i = 0; i < num_node; i++)
        {
            r_output.mNs[i] = Ns[i];

            r_output.mDNsDXs[i].resize(3);

            r_output.mDNsDXs[i][0] = DNs_DXs(i,0);
            r_output.mDNsDXs[i][1] = DNs_DXs(i,1);
            r_output.mDNsDXs[i][2] = DNs_DXs(i,2);
        }

        // double variables equation nodes id 
        for( auto it = mDoubleVariablesNeedEquationId.begin(); it != mDoubleVariablesNeedEquationId.end(); it = std::next(it) )
        {
            const DoubleVariable & r_variable = * it;
            const VariableKeyType variable_key = r_variable.Key();

            std::vector<std::size_t> & r_nodes_equation_id = r_output.mDoubleVariablesNodesEquationId[variable_key];

            r_nodes_equation_id.clear();
            r_nodes_equation_id.resize(num_node);

            for(std::size_t i = 0; i < num_node; i++)
            {
                r_nodes_equation_id[i] = const_cast<ModelPart::NodeType &> (r_geometry[i]).GetDof(r_variable).EquationId();
            }
        }

        // double variables value at quadrature point
        for( auto it = mDoubleVariablesNeedValue.begin(); it != mDoubleVariablesNeedValue.end(); it = std::next(it) )
        {
            const DoubleVariable & r_variable =  * it;
            const VariableKeyType variable_key = r_variable.Key();

            Vector nodes_value = ZeroVector(num_node);
            for(std::size_t i = 0; i < num_node; i++)
                nodes_value[i] = r_geometry[i].FastGetSolutionStepValue(r_variable);

            r_output.mDoubleVariables[variable_key] = boost::numeric::ublas::inner_prod( Ns, nodes_value );
        }

        // double variables 1st derivatives
        for( auto it = mDoubleVariablesNeedDX.begin(); it != mDoubleVariablesNeedDX.end(); it = std::next(it) )
        {
            const DoubleVariable & r_variable = * it;
            const VariableKeyType variable_key = r_variable.Key();

            // value at quadrature point
            Vector nodes_value = ZeroVector(num_node);
            for(std::size_t i = 0; i < num_node; i++)
                nodes_value[i] = r_geometry[i].FastGetSolutionStepValue(r_variable);

            Vector Dvalue_DXs = prod( trans(DNs_DXs), nodes_value );

            std::vector<double> & r_gradient = r_output.mDoubleVariablesDXs[variable_key];

            r_gradient.clear();
            r_gradient.resize(3);
            r_gradient[0] = Dvalue_DXs[0];
            r_gradient[1] = Dvalue_DXs[1];
            r_gradient[2] = Dvalue_DXs[2];
        }

        //coordinate (for debugging)
        Vector r_coordinate(3);
        noalias(r_coordinate) = ZeroVector(3);
        for(std::size_t i = 0; i < num_node; i++ )
            noalias(r_coordinate) += Ns[i]*r_geometry[i];

        r_output.mCoordinate.clear();
        r_output.mCoordinate.resize(3);
        r_output.mCoordinate[0] = r_coordinate[0];
        r_output.mCoordinate[1] = r_coordinate[1];
        r_output.mCoordinate[2] = r_coordinate[2];
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

    std::set<DoubleVariable> mDoubleVariablesNeedEquationId; //variables that need to get equation id
    std::set<DoubleVariable> mDoubleVariablesNeedValue;      //variables whose value need to be interpolated
    std::set<DoubleVariable> mDoubleVariablesNeedDX;         //variables whose gradient need to be interpolated

    // std::set<Array1dComponentVariable> mArray1dComponentVariablesNeedEquationId;
    // std::set<Array1dComponentVariable> mArray1dComponentVariablesNeedValue;
    // std::set<Array1dComponentVariable> mArray1dComponentVariablesNeedDX;

    friend class DistributedAssignment::DataUtility::DataPrinter;
};

}//namespace OversetAssembly
}//namespace Kratos
#endif
