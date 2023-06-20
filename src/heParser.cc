#include "heParser.h"
#include <fstream>

using namespace std;
//using namespace xfem;


xParseData getInfos(std::string parsename){
    xParseData parsedinfo;
    cout<<"Parse info in "<<parsename<<endl;

    parsedinfo.registerInt("order_approx_min",1);
    parsedinfo.registerInt("debug",0);
    parsedinfo.registerInt("XFEM");
    parsedinfo.registerInt("Discontinuity",1);
    

    // Octree parameters
   parsedinfo.registerInt("levelC",2);
   parsedinfo.registerInt("levelG",1);
   parsedinfo.registerInt("dimension",2);
//    parsedinfo.registerVector("octreeBoundingBoxMax", xtensor::xVector<>(10.,10.,0.));
//    parsedinfo.registerVector("octreeBoundingBoxMin", xtensor::xVector<>(-10.,-10.,0.));

    //Options
//    parsedinfo.registerInt("doErrorComputation",0);
//    parsedinfo.registerInt("doPostPro",0);

//  Generic interface parameters
    parsedinfo.registerInt("Generic_Interface", 1);
    parsedinfo.registerDouble("A",1.);
    parsedinfo.registerDouble("B",1.);
    parsedinfo.registerDouble("C",0.);
    parsedinfo.registerDouble("D",1.);

    parsedinfo.registerDouble("Resistivity");
    parsedinfo.registerDouble("Mat_Thickness");
    parsedinfo.registerDouble("Porosity");
    parsedinfo.registerDouble("Tortuosity");
    parsedinfo.registerDouble("Thermal_length");
    parsedinfo.registerDouble("Viscous_length");


    //Helmholtz related
    parsedinfo.registerString("Coupling_type", "A");
    parsedinfo.registerDouble("H_angular_frequency",6800);
    parsedinfo.registerDouble("H_wave_angle",0.);
    parsedinfo.registerDouble("Frequency_Hz",100.);
    parsedinfo.registerDouble("Cylinder_radius",0.1);

    parsedinfo.registerDouble("Nitsche_fluid");
    parsedinfo.registerDouble("Nitsche_solid");
    parsedinfo.registerDouble("Average_gamma");
    parsedinfo.registerDouble("Element_size");
    


    if (parsename=="") std::cout << "WARNING : info file NOT DEFINED, USING DEFAULT VALUES" << std::endl;
    else{
        std:: cout << "Parsing data info file" << std::endl;
        std::ifstream parsedinfofile (parsename.c_str());
        if(parsedinfofile)
            parsedinfo.parse(parsedinfofile);
        else
            std::cout << "WARNING : can\'t open file info file " << parsename.c_str() << " using default values " <<std::endl;
    }
    return parsedinfo;
}
