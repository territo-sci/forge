/* 
 * Author: Victor Sand (victor.sand@gmail.com)
 *
 */

#include <Forge.h>
#include <string>
#include <cstdlib>

using namespace osp;

int main() {

  std::string in = "/home/vsand/OpenSpace/enlilTestData_256_256_256_32.vdf";
  std::string out = "/home/vsand/OpenSpace/bricks_256_256_256_32.bdf";

  Forge* forge = Forge::New();
  
  forge->SetInFilename(in);
  forge->SetOutFilename(out);
  forge->SetStructure(0); // flat structure
  forge->SetBrickDimensions(32, 32, 32);
  
  //if (!forge->Read()) exit(1);
  //if (!forge->Write()) exit(1);
  
  // TODO testing testing
  if (!forge->Construct()) exit(1);


  delete forge;

  exit(0);
}
