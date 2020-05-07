#include "RGBImage.hpp"
#include "ImageFunctions.hpp"
#include "FunctionParameters.hpp"

#include <iostream>
#include <fstream>
#include <chrono>

int main(int argc, char *argv[])
{
    // Ensure parameters file was passed at command line.
    if (argc > 1) {
        std::ifstream input(argv[1]);

        if (input) {
            FunctionParameters params;

            while (!input.eof()) {
                params.setParams(input);
                RGBImage img(params.inputDir);

                try {
                    if (params.function == "addIntensity") {
                        ImageFunctions::addIntensity(img, params).outputImage(params.outputDir);
                    }
                    else if (params.function == "binarize") {
                        ImageFunctions::binarize(img, params).outputImage(params.outputDir);
                    }
                    else if (params.function == "scale") {
                        ImageFunctions::scale(img, params).outputImage(params.outputDir);
                    }
                    else if (params.function == "thresholdAdd") {
                        ImageFunctions::thresholdAdd(img, params).outputImage(params.outputDir);
                    }
                    else if (params.function == "thresholdBinarize") {
                        ImageFunctions::thresholdBinarize(img, params).outputImage(params.outputDir);
                    }
                    else if (params.function == "smooth2D") {
                        ImageFunctions::smooth2D(img, params).outputImage(params.outputDir);
                    }
                    else if (params.function == "smooth1D") {
                        ImageFunctions::smooth1D(img, params).outputImage(params.outputDir);
                    }
                    else if (params.function == "smooth1DInc") {
                        ImageFunctions::smooth1DInc(img, params).outputImage(params.outputDir);
                    }
                    else if (params.function == "colorBinarize") {
                        ImageFunctions::colorBinarize(img, params).outputImage(params.outputDir);
                    }
                    else if (params.function == "histogramStretch") {
                        ImageFunctions::histogramStretch(img, params).outputImage(params.outputDir);
                    }
                    else if (params.function == "optimalThreshold") {
                        ImageFunctions::optimalThreshold(img, params).outputImage(params.outputDir);
                    }
                    else if (params.function == "stretchThreshold") {
                        ImageFunctions::stretchThreshold(img, params).outputImage(params.outputDir);
                    }
                    else if (params.function == "colorHistogramStretch") {
                        ImageFunctions::colorHistogramStretch(img, params).outputImage(params.outputDir);
                    }
                    else if (params.function == "greyEdgeDetection") {
                        ImageFunctions::greyEdgeDetection(img, params);
                    }
                    else if (params.function == "colorEdgeDetection") {
                        ImageFunctions::colorEdgeDetection(img, params);
                    }
                }
                catch (std::exception e) {
                    std::cerr << e.what();
                }
            }
        }

        input.close();
    }
    else {
        std::cout << "Usage: ./ImageProcessing <Text document containing parameters as specified in the README.txt>\n";
    }
}