#if !defined(KRATOS_OVERSET_INTERPOLATIONOUTPUT_H_INCLUDED )
#define  KRATOS_OVERSET_INTERPOLATIONOUTPUT_H_INCLUDED

namespace Kratos
{
namespace OversetAssembly
{

class InterpolationOutput
{
public://make them public for debugging
    std::vector<std::size_t> mEquationsId;
    std::vector<double> mNs;
    std::vector<std::vector<double>> mDNsDXs;
    double mTemperature;
    std::vector<double> mTempGradient;
    std::vector<double> mCoordinate;//for debugging

public:
    InterpolationOutput()
    {}

    virtual ~InterpolationOutput()
    {}

    virtual void InterpolateDataFromElement( const Point & r_point, const ModelPart::ElementType & r_element )
    {
        const ModelPart::ElementType::GeometryType & r_geometry = r_element.GetGeometry();
        const std::size_t num_node = r_geometry.size();

        //mEquationsId
        mEquationsId.clear();
        mEquationsId.resize(num_node);
        for(std::size_t i = 0; i < num_node; i++)
		        mEquationsId[i] = const_cast<ModelPart::NodeType &> (r_geometry[i]).GetDof(TEMPERATURE).EquationId();			


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
        
        //field variable
        //TEMPERATURE
        Vector nodes_temp = ZeroVector(num_node);
        for(std::size_t i = 0; i < num_node; i++)
            nodes_temp[i] = r_geometry[i].FastGetSolutionStepValue(TEMPERATURE);

        double temp = boost::numeric::ublas::inner_prod( Ns, nodes_temp );

        Vector Dtemp_DXs = prod( trans(DNs_DXs), nodes_temp );

        //coordinate (for debugging)
        Vector r_coordinate(3);
        noalias(r_coordinate) = ZeroVector(3);
        for(std::size_t i = 0; i < num_node; i++ )
            noalias(r_coordinate) += Ns[i]*r_geometry[i];

        //
        mNs.clear();
        mDNsDXs.clear();

        mNs.resize(num_node);
        mDNsDXs.resize(num_node);
        
        for(std::size_t i = 0; i < num_node; i++)
        {
            mNs[i] = Ns[i];
            
            mDNsDXs[i].resize(3);

            mDNsDXs[i][0] = DNs_DXs(i,0);
            mDNsDXs[i][1] = DNs_DXs(i,1);
            mDNsDXs[i][2] = DNs_DXs(i,2);
        }

        mTemperature = temp;

        mTempGradient.clear();
        mTempGradient.resize(3);
        mTempGradient[0] = Dtemp_DXs[0];
        mTempGradient[1] = Dtemp_DXs[1];
        mTempGradient[2] = Dtemp_DXs[2];

        mCoordinate.clear();
        mCoordinate.resize(3);
        mCoordinate[0] = r_coordinate[0];
        mCoordinate[1] = r_coordinate[1];
        mCoordinate[2] = r_coordinate[2];

        {
            DistributedAssignment::DataUtility::DataPrinter printer;

            // std::cout<<__func__<<std::endl;
            // printf("%lg %lg %lg %lg\n",nodes_temp[0],nodes_temp[1],nodes_temp[2],nodes_temp[3]);
            // printf("mTempGradient\n");
            // printer.Print(mTempGradient);
        }
    }

private:
    virtual void Save( DistributedAssignment::DataUtility::Serializer & r_serializer ) const
    {
        r_serializer.Save(mEquationsId);
        r_serializer.Save(mNs);
        r_serializer.Save(mDNsDXs);
        r_serializer.Save(mTemperature);
        r_serializer.Save(mTempGradient);
        r_serializer.Save(mCoordinate);
    }

    virtual void Load( DistributedAssignment::DataUtility::Serializer & r_serializer )
    {
        r_serializer.Load(mEquationsId);
        r_serializer.Load(mNs);
        r_serializer.Load(mDNsDXs);
        r_serializer.Load(mTemperature);
        r_serializer.Load(mTempGradient);
        r_serializer.Load(mCoordinate);
    }

    virtual void Profile( DistributedAssignment::DataUtility::DataProfile & r_profile ) const
    {
        r_profile.SetIsTrivial(false);
    }

    virtual void Print( const DistributedAssignment::DataUtility::DataPrinter & r_printer ) const
    {
        std::cout << "{InterpolationOutput: ";
        r_printer.Print(mEquationsId);
        r_printer.Print(mNs);
        r_printer.Print(mDNsDXs);
        r_printer.Print(mTemperature);
        r_printer.Print(mTempGradient);
        r_printer.Print(mCoordinate);
        std::cout << "},";
    }

    friend class DistributedAssignment::DataUtility::Serializer;
    friend class DistributedAssignment::DataUtility::DataProfile;
    friend class DistributedAssignment::DataUtility::DataPrinter;
};

}//namespace OversetAssembly
}//namespace Kratos
#endif
