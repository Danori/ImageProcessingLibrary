#ifndef IMAGE_FUNCTIONS_HPP
#define IMAGE_FUNCTIONS_HPP

#include "RGBImage.hpp"
#include "FunctionParameters.hpp"

#include <map>

class ImageFunctions
{
private:
    // Consts for use in functions.
    static const int MAXRGB = 255;
    static const int MINRGB = 0;
    static constexpr double PI = 3.14159265358979323846;

    // Private helper functions.
    static unsigned int checkIntensity(const int value);
    static double euclideanDistance(int r, int g, int b, int cR, int cG, int cB);
    static std::map<int, int> initHistogramStretchMap(int lowerIntensity, int upperIntensity);
    static unsigned int findOptimalThreshold(RGBImage source, FunctionParameters::ROI roi);

public:
    // Public helper / utility functions.
    static RGBImage initHistogram(RGBImage source, FunctionParameters::ROI roi);

    // Image processing functions.
    static RGBImage addIntensity(RGBImage source, FunctionParameters params);
    static RGBImage binarize(RGBImage source, FunctionParameters params);
    static RGBImage scale(RGBImage source, FunctionParameters params);
    static RGBImage thresholdAdd(RGBImage source, FunctionParameters params);
    static RGBImage thresholdBinarize(RGBImage source, FunctionParameters params);
    static RGBImage smooth2D(RGBImage source, FunctionParameters params);
    static RGBImage smooth1D(RGBImage source, FunctionParameters params);
    static RGBImage smooth1DInc(RGBImage source, FunctionParameters params);
    static RGBImage colorBinarize(RGBImage source, FunctionParameters params);
    static RGBImage histogramStretch(RGBImage source, FunctionParameters params);
    static RGBImage optimalThreshold(RGBImage source, FunctionParameters params);
    static RGBImage stretchThreshold(RGBImage source, FunctionParameters params);
    static RGBImage colorHistogramStretch(RGBImage source, FunctionParameters params);
    static RGBImage greyEdgeDetection(RGBImage source, FunctionParameters params);
    static RGBImage colorEdgeDetection(RGBImage source, FunctionParameters params);
};

#endif // IMAGE_FUNCTIONS_HPP