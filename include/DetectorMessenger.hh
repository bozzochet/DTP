#include"globals.hh"
#include"G4UImessenger.hh"

class DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

class DetectorMessenger : public G4UImessenger
{
	public:
		DetectorMessenger(DetectorConstruction *);
		~DetectorMessenger();

	public:
		void SetNewValue(G4UIcommand * ,G4String);

	private:
		DetectorConstruction* Det;

		G4UIdirectory* thisDir;

		G4UIcmdWithoutParameter* UpdateCmd;
};
