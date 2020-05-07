#ifndef FUNCTION_PARAMETERS_HPP
#define FUNCTION_PARAMETERS_HPP

#include <vector>
#include <string>
#include <fstream>

class FunctionParameters {
private:
    void readNextParams(std::ifstream &input);

public:
    struct ROI {
        ROI(unsigned int x, unsigned int y, unsigned int sizeX, unsigned int sizeY) :
        x(x), y(y), sizeX(sizeX), sizeY(sizeY) { }

        std::string toString();

        unsigned int x;
        unsigned int y;
        unsigned int sizeX;
        unsigned int sizeY;

        std::vector<double> params;
    };

    FunctionParameters() : function(""), inputDir(""), outputDir("") { };
    FunctionParameters(std::ifstream &input);

    void setParams(std::ifstream &input);
    std::string toString();

    std::string function;
    std::string inputDir;
    std::string outputDir;
    std::vector<ROI> rois;
};

#endif // FUNCTION_PARAMETERS_HPP