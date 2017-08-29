#include <iostream>
#include <cmath>

// Project includes
#include "includes/define.h"
#include "includes/kratos_flags.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/cfd_variables.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "includes/serializer.h"
#include "utilities/geometry_utilities.h"

namespace Kratos {
    template< unsigned int TDim >
    class Eigen {
    public:
        double static Calculate_MatrixTrace3(
            const boost::numeric::ublas::bounded_matrix< double, TDim, TDim> matrix
        )
        {
            return matrix(0,0) + matrix(1,1) + matrix(2,2);
        }

        double static Calculate_MatrixTrace2(
            const boost::numeric::ublas::bounded_matrix< double, TDim, TDim> matrix
        )
        {
            return matrix(0,0) + matrix(1,1);
        }
        
        double static Calculate_DotProduct(
            const array_1d<double, TDim>& vec1,
            const array_1d<double, TDim>& vec2
        )
        {
            KRATOS_TRY
            
            if (TDim == 2)
                return Calculate_DotProduct2(vec1, vec2);
            else if (TDim == 3)
                return Calculate_DotProduct3(vec1, vec2);
            else
                KRATOS_THROW_ERROR(std::runtime_error,
                    "Eigen decomposition only supports 2D and 3D.","")
            
            KRATOS_CATCH("")

            return 0.0;
        }

        double static Calculate_DotProduct3(
            const array_1d<double, TDim>& vec1,
            const array_1d<double, TDim>& vec2
        )
        {
            return vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2];
        }

        double static Calculate_DotProduct2(
            const array_1d<double, TDim>& vec1,
            const array_1d<double, TDim>& vec2
        )
        {
            return vec1[0]*vec2[0] + vec1[1]*vec2[1];
        }

        void static Calculate_CrossProduct3(
            array_1d<double, TDim>& result,
            const array_1d<double, TDim>& vec1,
            const array_1d<double, TDim>& vec2
        )
        {
            result[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
            result[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2];
            result[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];
        }

        void static Calculate_CrossProduct2(
            array_1d<double, TDim>& result,
            const array_1d<double, TDim>& v1,
            const array_1d<double, TDim>& v2
        )
        {
            array_1d<double, 3> vec1;
            array_1d<double, 3> vec2;

            vec1[0] = v1[0];
            vec1[1] = v1[1];
            vec1[2] = 0.0;

            vec2[0] = v2[0];
            vec2[1] = v2[1];
            vec2[2] = 0.0;

            Calculate_CrossProduct3(result, vec1, vec2);
        }

        double static Calculate_MatrixDet3(
            const boost::numeric::ublas::bounded_matrix< double, TDim, TDim> a
        )
        {
            return a(0,0)*(a(1,1)*a(2,2) - a(1,2)*a(2,1)) - 
                a(0,1)*(a(1,0)*a(2,2) - a(2,0)*a(1,2)) +
                a(0,2)*(a(1,0)*a(2,1) - a(1,1)*a(2,0));
        }

        double static Calculate_MatrixDet2(
            const boost::numeric::ublas::bounded_matrix< double, TDim, TDim> a
        )
        {
            return a(0,0)*a(1,1) - a(0,1)*a(1,0);
        }

        void static Calculate_MatrixProduct3(
            boost::numeric::ublas::bounded_matrix< double, TDim, TDim>& result,
            const boost::numeric::ublas::bounded_matrix< double, TDim, TDim> a,
            const boost::numeric::ublas::bounded_matrix< double, TDim, TDim> b
        )
        {
            result(0,0) = a(0,0)*b(0,0) + a(0,1)*b(1,0) + a(0,2)*b(2,0);
            result(0,1) = a(0,0)*b(0,1) + a(0,1)*b(1,1) + a(0,2)*b(2,1);
            result(0,2) = a(0,0)*b(0,2) + a(0,1)*b(1,2) + a(0,2)*b(2,2);

            result(1,0) = a(1,0)*b(0,0) + a(1,1)*b(1,0) + a(1,2)*b(2,0);
            result(1,1) = a(1,0)*b(0,1) + a(1,1)*b(1,1) + a(1,2)*b(2,1);
            result(1,2) = a(1,0)*b(0,2) + a(1,1)*b(1,2) + a(1,2)*b(2,2);

            result(2,0) = a(2,0)*b(0,0) + a(2,1)*b(1,0) + a(2,2)*b(2,0);
            result(2,1) = a(2,0)*b(0,1) + a(2,1)*b(1,1) + a(2,2)*b(2,1);
            result(2,2) = a(2,0)*b(0,2) + a(2,1)*b(1,2) + a(2,2)*b(2,2);
        }

        void static Calculate_MatrixProduct2(
            boost::numeric::ublas::bounded_matrix< double, TDim, TDim>& result,
            const boost::numeric::ublas::bounded_matrix< double, TDim, TDim> a,
            const boost::numeric::ublas::bounded_matrix< double, TDim, TDim> b
        )
        {
            result(0,0) = a(0,0)*b(0,0) + a(0,1)*b(1,0);
            result(0,1) = a(0,0)*b(0,1) + a(0,1)*b(1,1);

            result(1,0) = a(1,0)*b(0,0) + a(1,1)*b(1,0);
            result(1,1) = a(1,0)*b(0,1) + a(1,1)*b(1,1);
        }        

        void static Calculate_EigenValues2(
            const boost::numeric::ublas::bounded_matrix< double, TDim, TDim> symm_matrix,
            array_1d<double, TDim>& eigen_values
        )
        {
            double temp_1 = 0.5*(symm_matrix(0,0)+symm_matrix(1,1));
            double temp_2 = 0.5*std::sqrt(
                                    std::pow(
                                        symm_matrix(0,0) - symm_matrix(1,1),
                                        2
                                    ) + 
                                    4*std::pow(
                                        symm_matrix(0,1),
                                        2
                                    )
                              );
            eigen_values[0] = temp_1 + temp_2;
            eigen_values[1] = temp_1 - temp_2;
        }

        void static Calculate_EigenValues3(
            const boost::numeric::ublas::bounded_matrix< double, TDim, TDim> symm_matrix,
            array_1d<double, TDim>& eigen_values
        )
        {
            //calculate  the trace of symm_matrix
            double q = Calculate_MatrixTrace3(symm_matrix)/3;
            
            boost::numeric::ublas::bounded_matrix< double, TDim, TDim> identity_matrix;
            identity_matrix(0,0) = 1;
            identity_matrix(0,1) = 0;
            identity_matrix(0,2) = 0;
            identity_matrix(1,0) = 0;
            identity_matrix(1,1) = 1;
            identity_matrix(1,2) = 0;
            identity_matrix(2,0) = 0;
            identity_matrix(2,1) = 0;
            identity_matrix(2,2) = 1;


            boost::numeric::ublas::bounded_matrix< double, TDim, TDim> temp;
            boost::numeric::ublas::bounded_matrix< double, TDim, TDim> result;
            temp = symm_matrix - q*identity_matrix;
            Calculate_MatrixProduct3(result,temp,temp);

            double p = std::sqrt(Calculate_MatrixTrace3(result)/6);
            noalias(temp) = (symm_matrix-q*identity_matrix)/p;
            
            double theta = std::acos(Calculate_MatrixDet3(temp)/2)/3;

            double t1 = p*2*std::cos(theta)+q;
            double t2 = p*2*std::cos(theta + 2.0*M_PI/3.0)+q;
            double t3 = p*2*std::cos(theta + 4.0*M_PI/3.0)+q;
            double swap;

            if (t2>t1)
            {
                swap = t1;
                t1 = t2;
                t2 = swap;
            }

            if (t3>t2){
                swap = t2;
                t2 = t3;
                t3 = swap;                
            }

            if (t2>t1)
            {
                swap = t1;
                t1 = t2;
                t2 = swap;                
            }

            eigen_values[0] = t1;
            eigen_values[1] = t2;
            eigen_values[2] = t3;
        }

        void static Calculate_EigenVector2(
            const boost::numeric::ublas::bounded_matrix< double, TDim, TDim> symm_matrix,
            const double eigen_value,
            array_1d<double, TDim>& eigen_vector
        )
        {
            double temp = std::sqrt(
                std::pow(symm_matrix(0,1), 2) + 
                std::pow(symm_matrix(0,0) - eigen_value, 2)
            );
            eigen_vector[0] = symm_matrix(0,1)/temp;
            eigen_vector[1] = -(symm_matrix(0,0)-eigen_value)/temp;
        }

        void static Calculate_EigenVector3(
            const boost::numeric::ublas::bounded_matrix< double, TDim, TDim> symm_matrix,
            const double eigen_value,
            array_1d<double, TDim>& eigen_vector
        )
        {
            array_1d<double, TDim> row0;
            array_1d<double, TDim> row1;
            array_1d<double, TDim> row2;

            row0[0] = symm_matrix(0,0)-eigen_value;
            row0[1] = symm_matrix(0,1);
            row0[2] = symm_matrix(0,2);

            row1[0] = symm_matrix(1,0);
            row1[1] = symm_matrix(1,1)-eigen_value;
            row1[2] = symm_matrix(1,2);

            row2[0] = symm_matrix(2,0);
            row2[1] = symm_matrix(2,1);
            row2[2] = symm_matrix(2,2)-eigen_value;

            
            array_1d<double, TDim> r0xr1;
            array_1d<double, TDim> r0xr2;
            array_1d<double, TDim> r1xr2;

            Calculate_CrossProduct3(r0xr1, row0, row1);
            Calculate_CrossProduct3(r0xr2, row0, row2);
            Calculate_CrossProduct3(r1xr2, row1, row2);

            double d0 = Calculate_DotProduct3(r0xr1, r0xr1);
            double d1 = Calculate_DotProduct3(r0xr2, r0xr2);
            double d2 = Calculate_DotProduct3(r0xr1, r0xr2);

            double dmax = d0;
            int imax = 0;
            if (d1>dmax) {
                dmax = d1;
                imax = 1;
            }
            if (d2>dmax){
                imax = 2;
            }

            if (imax==0)
                eigen_vector = r0xr1/std::sqrt(d0);
            if (imax==1)
                eigen_vector = r0xr2/std::sqrt(d1);
            if (imax==2)
                eigen_vector = r1xr2/std::sqrt(d2);
        }

        void static Calculate_EigenValues(
            const boost::numeric::ublas::bounded_matrix< double, TDim, TDim > symm_matrix,
            array_1d<double, TDim>& eigen_values
        )
        {
            KRATOS_TRY

            if (TDim == 2)
                Calculate_EigenValues2(symm_matrix, eigen_values);
            else if (TDim == 3)
                Calculate_EigenValues3(symm_matrix, eigen_values);
            else
                KRATOS_THROW_ERROR(std::runtime_error,
                    "Eigen decomposition only supports 2D and 3D.","")

            KRATOS_CATCH("")
        }

        void static Calculate_EigenVector(
            const boost::numeric::ublas::bounded_matrix< double, TDim, TDim > symm_matrix,
            const double eigen_value,
            array_1d<double, TDim>& eigen_vector
        )
        {
            KRATOS_TRY

            if (TDim == 2)
                Calculate_EigenVector2(symm_matrix, eigen_value, eigen_vector);
            else if (TDim == 3)
                Calculate_EigenVector3(symm_matrix, eigen_value, eigen_vector);
            else
                KRATOS_THROW_ERROR(std::runtime_error,
                    "Eigen decomposition only supports 2D and 3D.","")

            KRATOS_CATCH("")
        }        

        void static test_eigen22()
        {
            boost::numeric::ublas::bounded_matrix< double, 2, 2 > symm_matrix;
            symm_matrix(0,0) = -8;
            symm_matrix(0,1) = 2;
            symm_matrix(1,0) = 2;
            symm_matrix(1,1) = -5;

            array_1d<double, 2> eigen_values;
            array_1d<array_1d<double, 2>, 2> eigen_vectors;

            Calculate_EigenValues2(symm_matrix, eigen_values);
            std::cout<<eigen_values<<std::endl;

            Calculate_EigenVector2(symm_matrix, eigen_values[0], eigen_vectors[0]);
            Calculate_EigenVector2(symm_matrix, eigen_values[1], eigen_vectors[1]);

            std::cout<<eigen_vectors<<std::endl;            
        }

        void static test_eigen33()
        {
            boost::numeric::ublas::bounded_matrix< double, 3, 3 > symm_matrix;
            symm_matrix(0,0) = -8;
            symm_matrix(0,1) = 2;
            symm_matrix(0,2) = 3;
            symm_matrix(1,0) = 2;
            symm_matrix(1,1) = -5;
            symm_matrix(1,2) = 6;
            symm_matrix(2,0) = 3;
            symm_matrix(2,1) = 6;
            symm_matrix(2,2) = 9;
            
            array_1d<double, 3> eigen_values;
            array_1d<array_1d<double, 3>, 3> eigen_vectors;

            Calculate_EigenValues3(symm_matrix, eigen_values);
            std::cout<<eigen_values<<std::endl;

            Calculate_EigenVector3(symm_matrix, eigen_values[0], eigen_vectors[0]);
            Calculate_EigenVector3(symm_matrix, eigen_values[1], eigen_vectors[1]);
            Calculate_EigenVector3(symm_matrix, eigen_values[2], eigen_vectors[2]);

            std::cout<<eigen_vectors<<std::endl;

        }
    };
}