//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_OVERSET_APPLICATION_H_INCLUDED )
#define  KRATOS_OVERSET_APPLICATION_H_INCLUDED


// System includes
#include <string>
#include <iostream> 

// External includes 

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"

#include "custom_conditions/OversetCondition.h"


namespace Kratos
{

	// Variables definition 
  	KRATOS_DEFINE_VARIABLE(int, BLOCK_ID )

	class KratosOversetApplication : public KratosApplication
	{
	public:
		KRATOS_CLASS_POINTER_DEFINITION(KratosOversetApplication);

		KratosOversetApplication();

		virtual ~KratosOversetApplication()
		{}

		virtual void Register();

		virtual std::string Info() const
		{
			return "KratosOversetApplication";
		}

		virtual void PrintInfo(std::ostream& rOStream) const
		{
			rOStream << Info();
			PrintData(rOStream);
		}

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

	private:
		const OversetCondition mOversetCondition;

		KratosOversetApplication& operator=(KratosOversetApplication const& rOther);

		KratosOversetApplication(KratosOversetApplication const& rOther);

	}; // Class KratosOversetApplication 

}  // namespace Kratos.

#endif // KRATOS_OVERSET_APPLICATION_H_INCLUDED  defined 


