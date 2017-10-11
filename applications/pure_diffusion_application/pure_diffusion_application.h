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
#include "includes/condition.h"
#include "includes/ublas_interface.h"

#include "custom_elements/PoissonElement3D.h"
#include "custom_conditions/PointSourceCondition3D.h"
#include "custom_conditions/PoissonHeatFluxCondition3D.h"
#include "custom_conditions/PoissonOversetCondition3D.h"


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

 		const PoissonElement3D  mPoissonElement3D;
		const PointSourceCondition3D  mPointSourceCondition3D;
		const PoissonHeatFluxCondition3D  mPoissonHeatFluxCondition3D;
		const PoissonOversetCondition3D  mPoissonOversetCondition3D;
		
		/// Assignment operator.
		KratosPureDiffusionApplication& operator=(KratosPureDiffusionApplication  const& rOther);

		/// Copy constructor.
		KratosPureDiffusionApplication(KratosPureDiffusionApplication const& rOther);

	}; // Class KratosPureDiffusionApplication 


}  // namespace Kratos.


#endif // KRATOS_PUREDIFFUSION_APPLICATION_H_INCLUDED  defined
