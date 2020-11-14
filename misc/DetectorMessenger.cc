#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

DetectorMessenger::DetectorMessenger(DetectorConstruction* det)
: Det(det)
{
	thisDir = new G4UIdirectory("/silicon/det/");
	thisDir->SetGuidance("silicon control");

	UpdateCmd = new G4UIcmdWithoutParameter("/silicon/det/update",this);
	UpdateCmd->SetGuidance("update geometry if any changes");
	UpdateCmd->AvailableForStates(G4State_Idle);
}

DetectorMessenger::~DetectorMessenger()
{
	delete thisDir;
	delete UpdateCmd;
}
void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
	if(command == UpdateCmd) Det->updateGeometry();
}
