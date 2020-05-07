#include "FunctionParameters.hpp"

#include <sstream>

std::string FunctionParameters::ROI::toString()
{
    std::stringstream ss;

    ss << x << " " << y << " " << sizeX << " " << sizeY << " ";
    for (auto param : params) {
        ss << param << " ";
    }
    ss << "\n";

    return ss.str();
}

FunctionParameters::FunctionParameters(std::ifstream &input)
{
    readNextParams(input);
}

/**
 * - Tokenize the passed line, adding each token to the parameters vector.
 */
void FunctionParameters::readNextParams(std::ifstream &input)
{
    input >> function >> inputDir >> outputDir;

    std::string line;
    int x, y, sizeX, sizeY, param;
    std::getline(input, line);
    
    while (!isalpha(input.peek()) && !input.eof()) {
        std::getline(input, line);
        std::stringstream ss(line);

        ss >> x >> y >> sizeX >> sizeY;
        ROI roi(x, y, sizeX, sizeY);

        while (ss >> param) {
            roi.params.push_back(param);
        }
        
        rois.push_back(roi);
    }
}

void FunctionParameters::setParams(std::ifstream &input)
{
    rois.clear();
    readNextParams(input);
}

std::string FunctionParameters::toString()
{
    std::stringstream ss;

    ss << function << " " << inputDir << " " << outputDir << "\n";
    for (auto roi : rois) {
        ss << roi.toString();
    }

    return ss.str();
}