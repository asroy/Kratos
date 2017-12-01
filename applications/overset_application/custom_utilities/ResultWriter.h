//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//  Kratos default license: kratos/license.txt
//
//  Main authors:    YOUR_NAME_HERE
//

#if !defined(KRATOS_OVERSET_RESULT_WRITER_UTILITY_INCLUDED )
#define  KRATOS_OVERSET_RESULT_WRITER_UTILITY_INCLUDED

// System includes
#include <string>
#include <iostream> 

// Project includes 
#include "includes/define.h"
#include "includes/variables.h" 
#include "includes/model_part.h"
#include "includes/node.h"
#include "includes/element.h"

namespace Kratos
{
namespace OversetAssembly
{

class ResultWriter
{
public:
	
    KRATOS_CLASS_POINTER_DEFINITION(ResultWriter);

    ResultWriter(const ModelPart & r_model_part)
        :   mrModelPart(r_model_part)
    {}

    ~ResultWriter()
    {}

    void WriteVtkScalar(const std::string file_name, Variable<double> variable ) const
    {
        using GlobalToLocalNode = std::map<std::size_t,std::size_t>;

        GlobalToLocalNode node_global_to_local;

        //node global to local id
        {
            const ModelPart::NodesContainerType & r_nodes_pointer = const_cast<ModelPart &>(mrModelPart).Nodes();
            
            //loop over nodes
            std::size_t local_node_id = 0;
            for( ModelPart::NodesContainerType::ptr_const_iterator it_p_node = r_nodes_pointer.ptr_begin(); it_p_node != r_nodes_pointer.ptr_end(); it_p_node = std::next(it_p_node) )
            {
                std::size_t global_node_id = (* it_p_node)->GetId();
                node_global_to_local[global_node_id] = local_node_id;
                local_node_id++;
            }
        }

        //write
        {
            const ModelPart::NodesContainerType & r_nodes_pointer = const_cast<ModelPart &>(mrModelPart).Nodes();
            const ModelPart::ElementsContainerType & r_elements_pointer = const_cast<ModelPart &>(mrModelPart).Elements();
            
            FILE* fp;

            if ( ! (fp = fopen(file_name.c_str(), "w")) )
            {
                printf("wrong!: File cannot be opened\n");
                exit(EXIT_FAILURE);
            }

            fprintf(fp, "# vtk DataFile Version 2.0\nmesh\nASCII\nDATASET UNSTRUCTURED_GRID\n");

            //write nodes coordiante
            fprintf(fp, "POINTS %lu float\n", mrModelPart.NumberOfNodes());

            for( ModelPart::NodesContainerType::ptr_const_iterator it_p_node = r_nodes_pointer.ptr_begin(); it_p_node != r_nodes_pointer.ptr_end(); it_p_node = std::next(it_p_node) )
                fprintf(fp, "%lf %lf %lf\n", (* it_p_node)->X(), (* it_p_node)->Y(), (* it_p_node)->Z() );

            //write elements connectivity
            fprintf(fp, "CELLS %lu %lu\n",mrModelPart.NumberOfElements(), 5*mrModelPart.NumberOfElements() );

            for( ModelPart::ElementsContainerType::ptr_const_iterator it_p_element = r_elements_pointer.ptr_begin(); it_p_element != r_elements_pointer.ptr_end(); it_p_element = std::next(it_p_element) )
            {
                Element::GeometryType & r_geometry = (* it_p_element)->GetGeometry();
                
                fprintf(fp, "%lu ", r_geometry.size() );

                for( std::size_t i = 0; i < r_geometry.size(); i++ )
                    fprintf(fp, "%lu ", node_global_to_local[r_geometry[i].GetId()] );

                fprintf(fp, "\n");
            }

            //write element type
            fprintf(fp, "CELL_TYPES %lu\n", mrModelPart.NumberOfElements() );

            for( std::size_t i = 0; i < mrModelPart.NumberOfElements(); i++)
                fprintf(fp, "%d\n", 10 );

            //write nodes data
            fprintf(fp, "POINT_DATA %lu\n", mrModelPart.NumberOfNodes());
            fprintf(fp, "SCALARS scalars float 1\n");
            fprintf(fp, "LOOKUP_TABLE default\n");

            for( ModelPart::NodesContainerType::ptr_const_iterator it_p_node = r_nodes_pointer.ptr_begin(); it_p_node != r_nodes_pointer.ptr_end(); it_p_node = std::next(it_p_node) )
                fprintf(fp, "%lf\n", (* it_p_node)->GetSolutionStepValue(variable) );

            fclose(fp);
        }
    }

    void WriteVtkVector3(const std::string file_name, Variable<array_1d<double,3>> variable ) const
    {
        using GlobalToLocalNode = std::map<std::size_t,std::size_t>;

        GlobalToLocalNode node_global_to_local;

        //node global to local id
        {
            const ModelPart::NodesContainerType & r_nodes_pointer = const_cast<ModelPart &>(mrModelPart).Nodes();
            
            //loop over nodes
            std::size_t local_node_id = 0;
            for( ModelPart::NodesContainerType::ptr_const_iterator it_p_node = r_nodes_pointer.ptr_begin(); it_p_node != r_nodes_pointer.ptr_end(); it_p_node = std::next(it_p_node) )
            {
                std::size_t global_node_id = (* it_p_node)->GetId();
                node_global_to_local[global_node_id] = local_node_id;
                local_node_id++;
            }
        }

        //write
        {
            const ModelPart::NodesContainerType & r_nodes_pointer = const_cast<ModelPart &>(mrModelPart).Nodes();
            const ModelPart::ElementsContainerType & r_elements_pointer = const_cast<ModelPart &>(mrModelPart).Elements();
            
            FILE* fp;

            if ( ! (fp = fopen(file_name.c_str(), "w")) )
            {
                printf("wrong!: File cannot be opened\n");
                exit(EXIT_FAILURE);
            }

            fprintf(fp, "# vtk DataFile Version 2.0\nmesh\nASCII\nDATASET UNSTRUCTURED_GRID\n");

            //write nodes coordiante
            fprintf(fp, "POINTS %lu float\n", mrModelPart.NumberOfNodes());

            for( ModelPart::NodesContainerType::ptr_const_iterator it_p_node = r_nodes_pointer.ptr_begin(); it_p_node != r_nodes_pointer.ptr_end(); it_p_node = std::next(it_p_node) )
                fprintf(fp, "%lf %lf %lf\n", (* it_p_node)->X(), (* it_p_node)->Y(), (* it_p_node)->Z() );

            //write elements connectivity
            fprintf(fp, "CELLS %lu %lu\n",mrModelPart.NumberOfElements(), 5*mrModelPart.NumberOfElements() );

            for( ModelPart::ElementsContainerType::ptr_const_iterator it_p_element = r_elements_pointer.ptr_begin(); it_p_element != r_elements_pointer.ptr_end(); it_p_element = std::next(it_p_element) )
            {
                Element::GeometryType & r_geometry = (* it_p_element)->GetGeometry();
                
                fprintf(fp, "%lu ", r_geometry.size() );

                for( std::size_t i = 0; i < r_geometry.size(); i++ )
                    fprintf(fp, "%lu ", node_global_to_local[r_geometry[i].GetId()] );

                fprintf(fp, "\n");
            }

            //write element type
            fprintf(fp, "CELL_TYPES %lu\n", mrModelPart.NumberOfElements() );

            for( std::size_t i = 0; i < mrModelPart.NumberOfElements(); i++)
                fprintf(fp, "%d\n", 10 );

            //write nodes data
            fprintf(fp, "POINT_DATA %lu\n", mrModelPart.NumberOfNodes());
            fprintf(fp, "VECTORS vectors float\n");

            for( ModelPart::NodesContainerType::ptr_const_iterator it_p_node = r_nodes_pointer.ptr_begin(); it_p_node != r_nodes_pointer.ptr_end(); it_p_node = std::next(it_p_node) )
            {
                array_1d<double,3> array = (* it_p_node)->GetSolutionStepValue(variable);
                fprintf(fp, "%lf %lf %lf\n", array[0], array[1], array[2] );
            }

            fclose(fp);
        }
    }

private:
    const ModelPart & mrModelPart;
};

}  // namespace OversetAssembly
}  // namespace Kratos
#endif