//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//  Kratos default license: kratos/license.txt
//
//  Main authors:    YOUR_NAME_HERE
//

#if !defined(KRATOS_PUREDIFFUSION_APPLICATION_H_INCLUDED )
#define  KRATOS_PUREDIFFUSION_APPLICATION_H_INCLUDED

// System includes
#include <string>
#include <iostream> 

// External includes:none in this case 

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"
#include "custom_elements/poisson_2d.h" //including the file for the element
#include "custom_elements/poisson_3d.h"
#include "includes/condition.h"         //we'll also need conditions for the point heat loads
#include "custom_conditions/pointsource_2d.h"
#include "custom_conditions/pointsource_3d.h"
#include "custom_conditions/poisson_overset_2d.h"
#include "custom_conditions/poisson_overset_3d.h"
#include "includes/ublas_interface.h"


namespace Kratos
{

 
	///@name Kratos Globals

	///@{ 


	// Variables definition 
	KRATOS_DEFINE_VARIABLE(double, POINT_HEAT_SOURCE)


	class KratosPureDiffusionApplication : public KratosApplication
	{
	public:

		/// Pointer definition of KratosPureDiffusionApplication
		KRATOS_CLASS_POINTER_DEFINITION(KratosPureDiffusionApplication);


		/// Default constructor.
		KratosPureDiffusionApplication();


		/// Destructor.
		virtual ~KratosPureDiffusionApplication(){} 


		virtual void Register();


		/// Turn back information as a string.
		virtual std::string Info() const
		{
			return "KratosPureDiffusionApplication";
		}


		/// Print information about this object.
		virtual void PrintInfo(std::ostream& rOStream) const
		{
			rOStream << Info();
			PrintData(rOStream);
		}


		///// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
      KRATOS_WATCH("in my application");
      KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size() );
			rOStream << "Variables:" << std::endl;
			KratosComponents<VariableData>().PrintData(rOStream);
			rOStream << std::endl;
			rOStream << "Elements:" << std::endl;
			KratosComponents<Element>().PrintData(rOStream);
			rOStream << std::endl;
			rOStream << "Conditions:" << std::endl;
			KratosComponents<Condition>().PrintData(rOStream);
    }	



	protected:


	private:

 		const Poisson2D  mPoisson2D; //and here is our element.
		const PointSource2D  mPointSource2D; //and our condition
		const PoissonOverset2D  mPoissonOverset2D; //and our condition

 		const Poisson3D  mPoisson3D; //and here is our element.
		const PointSource3D  mPointSource3D; //and our condition
		const PoissonOverset3D  mPoissonOverset3D; //and our condition
		
		/// Assignment operator.
		KratosPureDiffusionApplication& operator=(KratosPureDiffusionApplication  const& rOther);


		/// Copy constructor.
		KratosPureDiffusionApplication(KratosPureDiffusionApplication const& rOther);

	}; // Class KratosPureDiffusionApplication 


}  // namespace Kratos.


#endif // KRATOS_PUREDIFFUSION_APPLICATION_H_INCLUDED  defined
