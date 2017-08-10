/* 
 * Author: Victor Sand (victor.sand@gmail.com)
 *
 */

#include <Forge.h>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iostream>
#include <boost/program_options.hpp>

using namespace osp;

namespace ForgeApp {
namespace options {
    boost::program_options::variables_map opts;
    const char *CONFIG = "config";
    const char *HELP = "help";
} // namespace options
} // namespace FurnaceApp

namespace options = boost::program_options;


void parseOptions(int ac, char *av[], options::variables_map &vm) {
    options::options_description desc("Options");
    desc.add_options()
            (ForgeApp::options::CONFIG, options::value<std::string>(), "configuration file")
            (ForgeApp::options::HELP, "print help");

    options::store(options::parse_command_line(ac, av, desc), vm);
}

int main(int ac, char *av[]) {

    // Very simple and short config file
    // TODO nicer implementation
    std::string config = "config/forgeConfig.txt";

    parseOptions(ac, av, ForgeApp::options::opts);

    if (ForgeApp::options::opts.count(ForgeApp::options::CONFIG)) {
        config = ForgeApp::options::opts[ForgeApp::options::CONFIG].as<std::string>();
    }

    std::ifstream in;
    in.open(config.c_str(), std::ifstream::in);
    if (!in.is_open()) {
        std::cout << "Could not open " << config << std::endl;
        exit(1);
    }

    std::string line;
    std::string inFilename = "notSet";
    std::string outFilename = "notSet";
    unsigned int xBrickDim = 0;
    unsigned int yBrickDim = 0;
    unsigned int zBrickDim = 0;
    while (std::getline(in, line)) {
        // Ignore empty lines and comments
        if (!line.empty() && line.at(0) != '#') {
            // Read variable name
            std::stringstream ss;
            ss.str(line);
            std::string var;
            ss >> var;
            // Read value
            if (var == "in_filename") {
                ss >> inFilename;
            } else if (var == "out_filename") {
                ss >> outFilename;
            } else if (var == "brick_dimensions") {
                ss >> xBrickDim;
                ss >> yBrickDim;
                ss >> zBrickDim;
            } else {
                std::cout << "Variable " << var << " not recognized" << std::endl;
                exit(1);
            }
        }
    }

    in.close();

    std::cout << "In filename: " << inFilename << std::endl;
    std::cout << "Out filename: " << outFilename << std::endl;
    std::cout << "Brick dimensions: " << xBrickDim << " " << yBrickDim
              << " " << zBrickDim << std::endl;

    Forge *forge = Forge::New();

    forge->SetInFilename(inFilename);
    forge->SetOutFilename(outFilename);
    forge->SetBrickDimensions(xBrickDim, yBrickDim, zBrickDim);

    // Construct TSP tree and calculate errors
    if (!forge->Construct()) {
        std::cerr << "Failed to construct TSP tree" << std::endl;
        exit(1);
    }

    delete forge;
    exit(0);
}
