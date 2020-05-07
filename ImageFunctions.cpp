#include "ImageFunctions.hpp"
#include "HSIImage.hpp"

#include <iostream>
#include <cmath>
#include <algorithm>
#include <exception>

/**
 * Checks if a pixel value is outside the range [0,255]
 * If so, sets it to the maxima or minima.
 */
unsigned int ImageFunctions::checkIntensity(const int value)
{
    if (value > MAXRGB) {
        return MAXRGB;
    }

    if (value < MINRGB) {
        return MINRGB;
    }

    return value;
}

/**
 * Calculates the Euclidean distance between the points (r, g, b) and (cR, cG, cB) in the RGB space.
 */
double ImageFunctions::euclideanDistance(int r, int g, int b, int cR, int cG, int cB)
{
    return sqrt(static_cast<double>(((r-cR)*(r-cR)) + ((g-cG)*(g-cG)) + ((b-cB)*(b-cB))));
}

/**
 * Create a histogram of the intensities of the source image within the passed ROI.
 */
RGBImage ImageFunctions::initHistogram(RGBImage source, FunctionParameters::ROI roi)
{
    RGBImage histogram(256, 256, RGBImage::Format::PGM);
    std::vector<int> intensityCount(256);
    std::vector<double> pixelRatios(256);

    // Count up the intensities of every pixel within the ROI.
    for (int i = roi.y; i < roi.y + roi.sizeY; i++) {
        for (int j = roi.x; j < roi.x + roi.sizeX; j++) {
            intensityCount[source.getPixel(i, j)]++;
        }
    }

    // Find the maximum count, the histogram heights are a ratio of the maximum.
    int maxCount = *std::max_element(intensityCount.begin(), intensityCount.end());
    
    // Draw the histogram lines for every intensity.
    for (int i = 0; i < 256; i++) {
        pixelRatios[i] = static_cast<double>(intensityCount[i]) / maxCount;
        pixelRatios[i] *= 255;
        for (int j = 0; j < static_cast<int>(pixelRatios[i]); j++) {
            histogram.setPixel(255 - j, i, MAXRGB);
        }
    }

    return histogram;
}

/**
 * Creates a map to stretch intensities within the range [lowerIntensity, upperIntensity] to [0, 255]
 * Intensities outside these ranges map to 0 if below lowerIntensity, or 255 if above upperIntensity.
 */
std::map<int, int> ImageFunctions::initHistogramStretchMap(int lowerIntensity, int upperIntensity)
{
    std::map<int, int> intensityMap;

    // Map intensites below lowerIntensity to 0.
    for (int i = 0; i < lowerIntensity; i++) {
        intensityMap.insert(std::make_pair(i, 0));
    }

    // Map intensities within range [lowerIntensity, upperIntensity] to [0, 255] linearly.
    int rangeOfStretching = upperIntensity - lowerIntensity;
    double stretchVal = static_cast<double>(255) / rangeOfStretching;
    double newIntensity = 0;
    for (int i = lowerIntensity; i < upperIntensity; i++, newIntensity += stretchVal) {
        intensityMap.insert(std::make_pair(i, static_cast<int>(newIntensity)));
    }

    // Map intensies above upperIntensity to 255.
    for (int i = upperIntensity; i < 256; i++) {
        intensityMap.insert(std::make_pair(i, 255));
    }

    return intensityMap;
}

unsigned int ImageFunctions::findOptimalThreshold(RGBImage source, FunctionParameters::ROI roi)
{
    int oldThreshold = 0, newThreshold = 0, pixelSum = 0;

    // Sum the pixels.
    for (int i = roi.y; i < roi.y + roi.sizeY; i++) {
        for (int j = roi.x; j < roi.x + roi.sizeX; j++) {
            pixelSum += source.getPixel(i, j);
        }
    }

    // Use the mean of the ROI as a starting point.
    oldThreshold = pixelSum / (roi.sizeX * roi.sizeY);
    newThreshold = oldThreshold;

    do {
        oldThreshold = newThreshold;
        int numBackgroundPixels = 0, numForegroundPixels = 0;
        long backgroundSum = 0, foregroundSum = 0;

        // Partition pixels as background or foreground, sum their intensities.
        for (int i = roi.y; i < roi.y + roi.sizeY; i++) {
            for (int j = roi.x; j < roi.x + roi.sizeX; j++) {
                int pixel = source.getPixel(i, j);
                if (pixel < oldThreshold) {
                    numBackgroundPixels++;
                    backgroundSum += pixel;
                }
                else {
                    numForegroundPixels++;
                    foregroundSum += pixel;
                }        
            }
        }

        // Find the the mean intensity of both the background and foreground.
        double backgroundMean = static_cast<double>(backgroundSum) / numBackgroundPixels;
        double foregroundMean = static_cast<double>(foregroundSum) / numForegroundPixels;

        newThreshold = (backgroundMean + foregroundMean) / 2;
    } while (std::abs(newThreshold - oldThreshold) >= 1); // Continue until little to no change is made to the threshold.

    return newThreshold;
}

/**
 * Increases the brightness of the entire source image.
 */
RGBImage ImageFunctions::addIntensity(RGBImage source, FunctionParameters params)
{
    RGBImage output = source;

    for (auto roi : params.rois) {

        int value = static_cast<int>(roi.params[0]);

        switch (output.getFormat()) {
            case RGBImage::Format::PGM:
                for (int i = 0; i < output.getHeight(); i++) {
                    for (int j = 0; j < output.getWidth(); j++) {
                        output.setPixel(i, j, checkIntensity(source.getPixel(i, j) + value));
                    }
                }
            break;

            case RGBImage::Format::PPM:
                for (int i = 0; i < output.getHeight(); i++) {
                    for (int j = 0; j < output.getWidth(); j++) {
                        for (int k = 0; k < 3; k++) {
                            output.setPixel(i, j, k, checkIntensity(source.getPixel(i, j, k) + value));
                        }
                    }
                }
            break;
        }
    }

    return output;
}

/**
 * For every pixel, if the intensity is above the threshold value, set it to white.
 * Otherwise, set it to black.
 */
RGBImage ImageFunctions::binarize(RGBImage source, FunctionParameters params)
{
    if (source.getFormat() != RGBImage::Format::PGM) {
        throw std::runtime_error("Error: binarize function only compatible with .pgm images.\n");
    }

    RGBImage output = source;

    for (auto roi : params.rois) {
        int threshold = roi.params[0];

        for (int i = 0; i < output.getHeight(); i++) {
            for (int j = 0; j < output.getWidth(); j++) {
                if (source.getPixel(i, j) < threshold) {
                    output.setPixel(i, j, MINRGB);
                }
                else {
                    output.setPixel(i, j, MAXRGB);
                }
            }
        }
    }

    return output;
}

/**
 * Scale the image by a factor of 0.5 or 2.
 */
RGBImage ImageFunctions::scale(RGBImage source, FunctionParameters params)
{
    if (source.getFormat() != RGBImage::Format::PGM) {
        throw std::runtime_error("Error: scale function only compatible with .pgm images.\n");
    }

    RGBImage output = source;

    for (auto roi : params.rois) {
        double ratio = roi.params[0];

        if (ratio != 0.5 || ratio != 2.0) {
            throw std::runtime_error("Error: invalid scale ratio passed to scale function.\n");
        }

        int newHeight = static_cast<int>(static_cast<double>(source.getHeight()) * ratio);
        int newWidth = static_cast<int>(static_cast<double>(source.getWidth()) * ratio);

        output.resize(newWidth, newHeight);

        for (int i = 0; i < newHeight; i++) {
            for (int j = 0; j < newWidth; j++) {

                int i2 = static_cast<int>(floor(static_cast<double>(i/ratio)));
                int j2 = static_cast<int>(floor(static_cast<double>(j/ratio)));

                if (ratio == 2) {
                    output.setPixel(i, j, checkIntensity(source.getPixel(i2, j2)));
                }
                else if (ratio == 0.5) {
                    int value = source.getPixel(i2, j2) + source.getPixel(i2, j2+1) + source.getPixel(i2+1, j2) + source.getPixel(i2+1, j2+1);
                    output.setPixel(i, j, checkIntensity(value/4));
                }
            }
        }
    }

    return output;
}

/**
 * For every pixel, if the source intensity is above the threshold, leave it unchanged.
 * Otherwise, increase its intensity.
 */
RGBImage ImageFunctions::thresholdAdd(RGBImage source, FunctionParameters params)
{
    if (source.getFormat() != RGBImage::Format::PGM) {
        throw std::runtime_error("Error: thresholdAdd function only compatible with .pgm images.\n");
    }

    RGBImage output = source;

    for (auto roi : params.rois) {
        int threshold = static_cast<int>(roi.params[0]);
        int intensity = static_cast<int>(roi.params[1]);

        for (int i = 0; i < output.getHeight(); i++) {
            for (int j = 0; j < output.getWidth(); j++) {
                if (source.getPixel(i, j) < threshold) {
                    output.setPixel(i, j, checkIntensity(source.getPixel(i, j) + intensity));
                }
            }
        }
    }

    return output;
}

/**
 * For every pixel in each ROI, if the pixel is between thresholds T1 / T2, set the pixel to white.
 * Otherwise, set the pixel to black.
 */
RGBImage ImageFunctions::thresholdBinarize(RGBImage source, FunctionParameters params)
{
    if (source.getFormat() != RGBImage::Format::PGM) {
        throw std::runtime_error("Error: thresholdBinarize function only compatible with .pgm images.\n");
    }

    RGBImage output = source;

    for (auto roi : params.rois) {
        if (roi.params.size() != 2) {
            throw std::runtime_error("Error: Incorrect number of parameters passed to thresholdBinarize (requires 2).\n");
        }

        int lowerThreshold = roi.params[0];
        int upperThreshold = roi.params[1];

        if (upperThreshold < lowerThreshold) {
            int temp = lowerThreshold;
            lowerThreshold = upperThreshold;
            upperThreshold = temp;
        }

        for (int i = roi.y; i < roi.y + roi.sizeY; i++) {
            for (int j = roi.x; j < roi.x + roi.sizeX; j++) {
                if (lowerThreshold < source.getPixel(i, j) && source.getPixel(i, j) < upperThreshold) {
                    output.setPixel(i, j, MAXRGB);
                }
                else {
                    output.setPixel(i, j, MINRGB);
                }
            }
        }
    }

    return output;
}

/**
 * For every pixel in each ROI, calculate the average of the pixels in the neighborhood around the pixel
 * Set the pixel to that average. This will result in a smoothing effect.
 */
RGBImage ImageFunctions::smooth2D(RGBImage source, FunctionParameters params)
{
    if (source.getFormat() != RGBImage::Format::PGM) {
        throw std::runtime_error("Error: smooth2D function only compatible with .pgm images.\n");
    }

    RGBImage output = source;

    for (auto roi : params.rois) {
        if (roi.params.size() != 1) {
            std::cerr << "Error: Incorrect number of parameters passed to smooth2D (1 required).\n"
                << "Erraneous ROI: " << roi.toString();
            continue;
        }

        int windowSize = roi.params[0];
        int halfWindowSize = windowSize / 2;

        if (windowSize % 2 != 1) {
            std::cerr << "Error: Window size for smooth2D function not odd.\n"
                << "Erraneous ROI: " << roi.toString();
            continue;
        }

        for (int i = roi.y; i < roi.y + roi.sizeY; i++) {
            for (int j = roi.x; j < roi.x + roi.sizeX; j++) {
                int numValidPixels = 0;
                int sum = 0;
                for (int l = i - halfWindowSize; l < i + halfWindowSize + 1; l++) {
                    for (int k = j - halfWindowSize; k < j + halfWindowSize + 1; k++) {
                        if (source.validPixel(l, k)) {
                            numValidPixels++;
                            sum += source.getPixel(l, k);
                        }
                    }
                }

                output.setPixel(i, j, sum / numValidPixels);
            }
        }
    }


    return output;
}

/**
 * Similar to smooth2D, however decomposses the 2D pass into two seperate 1D passes.
 * This change results in better time complexity in terms of the windowsize.
 */
RGBImage ImageFunctions::smooth1D(RGBImage source, FunctionParameters params)
{
    if (source.getFormat() != RGBImage::Format::PGM) {
        throw std::runtime_error("Error: smooth1D function only compatible with .pgm images.\n");
    }

    RGBImage output = source;

    for (auto roi : params.rois) {
        if (roi.params.size() != 1) {
            std::cerr << "Error: Incorrect number of parameters passed to smooth1D (1 required).\n"
                << "Erraneous ROI: " << roi.toString();
            continue;
        }

        int windowSize = roi.params[0];
        int halfWindowSize = windowSize / 2;

        if (windowSize % 2 != 1) {
            std::cerr << "Error: Window size for smooth1D function not odd.\n"
                << "Erraneous ROI: " << roi.toString();
            continue;
        }

        for (int i = roi.y; i < roi.y + roi.sizeY; i++) {
            for (int j = roi.x; j < roi.x + roi.sizeX; j++) {
                int sum = 0;
                int numValidPixels = 0;
                for (int k = j - halfWindowSize; k < j + halfWindowSize + 1; k++) {
                    if (source.validPixel(i, k)) {
                        numValidPixels++;
                        sum += source.getPixel(i, k);
                    }
                }

                output.setPixel(i, j, sum / numValidPixels);
            }
        }

        RGBImage source = output;

        for (int i = roi.y; i < roi.y + roi.sizeY; i++) {
            for (int j = roi.x; j < roi.x + roi.sizeX; j++) {
                int numValidPixels = 0;
                int sum = 0;
                for (int k = i - halfWindowSize; k < i + halfWindowSize + 1; k++) {
                    if (source.validPixel(k, j)) {
                        numValidPixels++;
                        sum += source.getPixel(k, j);
                    }
                }

                output.setPixel(i, j, sum / numValidPixels);
            }
        }
    }

    return output;
}

/**
 * Similar to smooth1D, however now defines the the sums in terms of the prior
 * calculated sum for every pixel. This makes the time complexity completely
 * independent of the windowsize.
 */
RGBImage ImageFunctions::smooth1DInc(RGBImage source, FunctionParameters params)
{
    if (source.getFormat() != RGBImage::Format::PGM) {
        throw std::runtime_error("Error: smooth1DInc function only compatible with .pgm images.\n");
    }

    RGBImage output = source;

    for (auto roi : params.rois) {
        if (roi.params.size() != 1) {
            std::cerr << "Error: Incorrect number of parameters passed to smooth1DInc (1 required).\n"
                << "Erraneous ROI: " << roi.toString();
            continue;
        }

        int windowSize = roi.params[0];
        int halfWindowSize = windowSize / 2;

        if (windowSize % 2 != 1) {
            std::cerr << "Error: Window size for smooth1DInc function not odd.\n"
                << "Erraneous ROI: " << roi.toString();
            continue;
        }

        for (int i = roi.y + halfWindowSize; i < roi.y + roi.sizeY - halfWindowSize; i++) {
            int prevSum = 0;
            for (int j = roi.x + halfWindowSize; j < roi.x + roi.sizeX - halfWindowSize; j++) {
                if (j == roi.x + halfWindowSize) {
                    for (int k = j - halfWindowSize; k < j + halfWindowSize + 1; k++) {
                        prevSum += source.getPixel(i, k);
                    }
                    output.setPixel(i, j, prevSum / windowSize);
                }
                else {
                    prevSum = prevSum + source.getPixel(i, j + halfWindowSize) - source.getPixel(i, j - halfWindowSize - 1);
                    output.setPixel(i, j, prevSum / windowSize);
                }
            }
        }

        RGBImage source = output;

        for (int j = roi.x + halfWindowSize; j < roi.x + roi.sizeX - halfWindowSize; j++) {
            int prevSum = 0;
            for (int i = roi.y + halfWindowSize; i < roi.y + roi.sizeY - halfWindowSize; i++) {
                if (i == roi.y + halfWindowSize) {
                    for (int k = i - halfWindowSize; k < i + halfWindowSize + 1; k++) {
                        prevSum += source.getPixel(k, j);
                    }
                    output.setPixel(i, j, prevSum / windowSize);
                }
                else {
                    prevSum = prevSum + source.getPixel(i + halfWindowSize, j) - source.getPixel(i - halfWindowSize - 1, j);
                    output.setPixel(i, j, prevSum / windowSize);
                }
            }
        }
    }

    return output;
}

/**
 * For every pixel in each ROI, calculates the euclidean distance of each pixel to a
 * color C (defined by cR cG cB), and if the pixel is within a threshold distance T,
 * increase the intensity of each RGB channel by intensity I.
 */
RGBImage ImageFunctions::colorBinarize(RGBImage source, FunctionParameters params)
{
    if (source.getFormat() != RGBImage::Format::PPM) {
        throw std::runtime_error("Error: colorBinarize function only compatible with .ppm images.\n");
    }

    RGBImage output = source;

    for (auto roi: params.rois) {
        if (roi.params.size() != 5) {
            std::cerr << "Error: Incorrect number of parameters passed to colorBinarize (5 required).\n"
                << "Erraneous ROI: " << roi.toString();
            continue;
        }

        int cR = roi.params[0], cG = roi.params[1], cB = roi.params[2];
        double threshold = roi.params[3];
        int intensity = roi.params[4];

        for (int i = roi.y; i < roi.y + roi.sizeY; i++) {
            for (int j = roi.x; j < roi.x + roi.sizeX; j++) {
                int r = source.getPixel(i, j, RGBImage::Channel::RED);
                int g = source.getPixel(i, j, RGBImage::Channel::GREEN);
                int b = source.getPixel(i, j, RGBImage::Channel::BLUE);

                // If the distance between the two points in the RGB space is less than the
                if (euclideanDistance(r, g, b, cR, cG, cB) < threshold) {
                    output.setPixel(i, j, RGBImage::Channel::RED,
                        checkIntensity(source.getPixel(i, j, RGBImage::Channel::RED) + intensity));
                    output.setPixel(i, j, RGBImage::Channel::GREEN,
                        checkIntensity(source.getPixel(i, j, RGBImage::Channel::GREEN) + intensity));
                    output.setPixel(i, j, RGBImage::Channel::BLUE, 
                        checkIntensity(source.getPixel(i, j, RGBImage::Channel::BLUE) + intensity));
                }
                else {
                    output.setPixel(i, j, RGBImage::Channel::RED,
                        MINRGB);
                    output.setPixel(i, j, RGBImage::Channel::GREEN,
                        MINRGB);
                    output.setPixel(i, j, RGBImage::Channel::BLUE,
                        MINRGB);
                }
            }
        }
    }

    return output;
}

/**
 * This function takes two parameters for each ROI: lowerT and upperT, and
 * stretches the intensity values within that range to 0-255 linearly. Any
 * intensities outside the lowerT-upperT range are set to 0 or 255 accordingly.
 * It also has the functionality to output histogram images before and after
 * the stretching to visualize the stretching programmatically.
 */
RGBImage ImageFunctions::histogramStretch(RGBImage source, FunctionParameters params)
{
    if (source.getFormat() != RGBImage::Format::PGM) {
        throw std::runtime_error("Error: histogramStretch function only compatible with .pgm images.\n");
    }

    RGBImage output = source;

    int roiCount = 0;
    for (auto roi : params.rois) {
        if (roi.params.size() != 2) {
            std::cerr << "Error: Incorrect number of parameters passed to histogramStretch (2 required).\n"
                << "Erraneous ROI: " << roi.toString();
            continue;
        }

        int lowerIntensity = roi.params[0];
        int upperIntensity = roi.params[1];

        // Set intensities accordingly if they are out of order.
        if (upperIntensity < lowerIntensity) {
            int temp = lowerIntensity;
            lowerIntensity = upperIntensity;
            upperIntensity = temp;
        }

        // Create output directory from the user's entered output directory.
        int dotIndex = params.outputDir.find_last_of('.');
        std::string beforeHistogramDir = params.outputDir.substr(0, dotIndex) + 
            "BeforeROI" + std::to_string(roiCount) + ".pgm";

        // Create the before histogram, output the image.
        initHistogram(source, roi).outputImage(beforeHistogramDir);

        // Create the map from unstretched intensity values to stretched intensity values.
        std::map<int, int> intensityMap = initHistogramStretchMap(lowerIntensity, upperIntensity);

        // Adjust the intensities based on the map.
        for (int i = roi.y; i < roi.y + roi.sizeY; i++) {
            for (int j = roi.x; j < roi.x + roi.sizeX; j++) {
                output.setPixel(i, j, intensityMap[source.getPixel(i, j)]);
            }
        }

        // Create the after histogram, output the image.
        std::string afterHistogramDir = params.outputDir.substr(0, dotIndex) +
            "AfterROI" + std::to_string(roiCount) + ".pgm";
        initHistogram(output, roi).outputImage(afterHistogramDir);

        roiCount++;
    }

    return output;
}

/**
 * This function uses the optimal thresholding algorithm discussed in class to
 * find a threshold which will binarize the image and retain as much detail as
 * possible in each ROI. The algorithm uses the mean intensity within the ROI
 * as a starting point, then partitions the pixels in the ROI as foreground or
 * background pixels. It then finds the mean of foreground and background pixels,
 * takes the average of the two, and iterates until little to no change is made
 * to the the threshold, at which point the ROI is binarized based on the final
 * threshold.
 */
RGBImage ImageFunctions::optimalThreshold(RGBImage source, FunctionParameters params)
{
    if (source.getFormat() != RGBImage::Format::PGM) {
        throw std::runtime_error("Error: optimalThreshold function compatible with only .pgm images.\n");
    }

    RGBImage output = source;

    for (auto roi : params.rois) {
        int threshold = findOptimalThreshold(source, roi);

        // Binarize the image.
        for (int i = roi.y; i < roi.y + roi.sizeY; i++) {
            for (int j = roi.x; j < roi.x + roi.sizeX; j++) {
                if (source.getPixel(i, j) < threshold) {
                        output.setPixel(i, j, MINRGB);
                    }
                    else {
                        output.setPixel(i, j, MAXRGB);
                    }     
            }
        }
    }

    return output;
}

/**
 * This function combines the two above functions by first performing optimal
 * thresholding to partition the pixels into background and foreground pixels,
 * and then stretching the foreground and background pixels seperately. The lowerT
 * and upperT used for stretching are the minimum and maximum intensity values for
 * the foreground, and similarly the minimum and maximum intensity values for the 
 * background.
 */
RGBImage ImageFunctions::stretchThreshold(RGBImage source, FunctionParameters params)
{
    if (source.getFormat() != RGBImage::Format::PGM) {
        throw std::runtime_error("Error: stretchThreshold function compatible with only .pgm images.\n");
    }

    RGBImage output = source;

    for (auto roi : params.rois) {
        int threshold = findOptimalThreshold(source, roi);

        // Map the (i, j) positions of each pixel to whether that pixel is background or foreground.
        // Find the minimum and maximum for both the background and foreground.
        int backgroundMin = 256, backgroundMax = -1;
        int foregroundMin = 256, foregroundMax = -1;
        std::map<std::pair<int, int>, bool> isBackground;
        for (int i = roi.y; i < roi.y + roi.sizeY; i++) {
            for (int j = roi.x; j < roi.x + roi.sizeX; j++) {
                int pixel = source.getPixel(i, j);
                if (pixel < threshold) {
                    // Iteratively find the background minimum.
                    if (pixel < backgroundMin) {
                        backgroundMin = pixel;
                    }

                    // Iteratively find the background maximum.
                    if (pixel > backgroundMax) {
                        backgroundMax = pixel;
                    }

                    isBackground.insert(std::make_pair(std::make_pair(i, j), true));
                }
                else {
                    // Iteratively find the foreground minimum.
                    if (pixel < foregroundMin) {
                        foregroundMin = pixel;
                    }

                    // Iteratively find the foreground maximum.
                    if (pixel > foregroundMax) {
                        foregroundMax = pixel;
                    }

                    isBackground.insert(std::make_pair(std::make_pair(i, j), false));
                }     
            }
        }

        // Create intensity maps to stretch the histograms.
        std::map<int, int> backgroundIntensityMap = initHistogramStretchMap(backgroundMin, backgroundMax);
        std::map<int, int> foregroundIntensityMap = initHistogramStretchMap(foregroundMin, foregroundMax);

        // Adjust the intensities based on whether the pixel is foreground or background, based on the respective map.
        for (int i = roi.y; i < roi.y + roi.sizeY; i++) {
            for (int j = roi.x; j < roi.x + roi.sizeX; j++) {
                if (isBackground[std::pair<int, int>(i, j)]) {
                    output.setPixel(i, j, backgroundIntensityMap[output.getPixel(i, j)]);
                }  
                else {
                    output.setPixel(i, j, foregroundIntensityMap[output.getPixel(i, j)]);
                }
            }
        }
    }

    return output;
}

/**
 * This function is similar to the histogramStretch function above, however this
 * function is compatible with .ppm formats. In addition to the lowerT and upperT
 * parameters, there is an optional third parameter: RGB. If provided, the
 * function will perform the stretching only on the respective RGB channel in that
 * ROI. If not provided it will perform the stretching on all three channels with
 * the same lowerT and upperT. The RGB parameter may be 0, 1, or 2, which
 * corresponds to R, G, or B, respectively.
 */
RGBImage ImageFunctions::colorHistogramStretch(RGBImage source, FunctionParameters params)
{
    if (source.getFormat() != RGBImage::Format::PPM) {
        throw std::runtime_error("Error: colorHistogramStretch function only compatible with .ppm images.\n");
    }

    RGBImage output = source;

    int roiCount = 0;
    for (auto roi : params.rois) {
        if (roi.params.size() != 2 && roi.params.size() != 3) {
            std::cerr << "Error: Incorrect number of parameters passed to colorHistogramStretch (2-3 required).\n"
                << "Erraneous ROI: " << roi.toString();
            continue;
        }

        int lowerIntensity = roi.params[0];
        int upperIntensity = roi.params[1];
        RGBImage::Channel channel = RGBImage::Channel::INVALID;

        // If a specific channel is provided by the user, assign channel to it
        if (roi.params.size() == 3) {
            int rgb = roi.params[2];

            switch (rgb) {
                case RGBImage::Channel::RED:
                    channel = RGBImage::Channel::RED;
                break;

                case RGBImage::Channel::GREEN:
                    channel = RGBImage::Channel::GREEN;
                break;

                case RGBImage::Channel::BLUE:
                    channel = RGBImage::Channel::BLUE;
                break;

                default:
                    std::cerr << "Error: Unrecognized rgb parameter passed to colorHistogramStretch function.\n"
                        << "Valid values include: 0 (Red), 1 (Green), 2 (Blue)\n"
                        << "Applying histogram stretch to all channels by default.\n";
                break;
            }
        }

        // Set intensities accordingly if they are out of order.
        if (upperIntensity < lowerIntensity) {
            int temp = lowerIntensity;
            lowerIntensity = upperIntensity;
            upperIntensity = temp;
        }

        // Create the map from unstretched intensity values to stretched intensity values.
        std::map<int, int> intensityMap = initHistogramStretchMap(lowerIntensity, upperIntensity);

        // Adjust the intensities based on the map.
        for (int i = roi.y; i < roi.y + roi.sizeY; i++) {
            for (int j = roi.x; j < roi.x + roi.sizeX; j++) {
                // If channel option was passed, adjust only individual channel.
                if (channel != RGBImage::Channel::INVALID) {
                    output.setPixel(i, j, channel, 
                        intensityMap[output.getPixel(i, j, channel)]);
                }
                else {
                    output.setPixel(i, j, RGBImage::Channel::RED, 
                        intensityMap[output.getPixel(i, j, RGBImage::Channel::RED)]);
                    output.setPixel(i, j, RGBImage::Channel::GREEN,
                        intensityMap[output.getPixel(i, j, RGBImage::Channel::GREEN)]);
                    output.setPixel(i, j, RGBImage::Channel::BLUE,
                        intensityMap[output.getPixel(i, j, RGBImage::Channel::BLUE)]);
                }
            }
        }

        roiCount++;
    }

    return output;
}

/**
 * This function uses the sobel filter to perform edge detection on .pgm images.
 * It allows the user to select either 3x3 or 5x5 sobel filters. Using the sobel
 * filter it calculates the gradient for every pixel as well as the direction of
 * the edge. The function outputs three images: an image to visualize the gradient
 * magnitudes, an image to binarize this prior image using a threshold value entered
 * by the user, and an image which further binarizes this image by edge direction.
 */
RGBImage ImageFunctions::greyEdgeDetection(RGBImage source, FunctionParameters params)
{
    if (source.getFormat() != RGBImage::Format::PGM) {
        throw std::runtime_error("Error: greyEdgeDetection function only compatible with .pgm images.\n");
    }

    // Create output images.
    RGBImage edgeIntensity = source;
    RGBImage edgeThreshold = source;
    RGBImage edgeDirectionThreshold = source;

    const double radiansToDegrees = 180 / PI;

    std::vector<int> xSobel;
    std::vector<int> ySobel;
    double normalizeVal;

    for (auto roi : params.rois) {
        int sobelSize = roi.params[0];
        int halfSobelSize = sobelSize / 2;
        double threshold = roi.params[1];
        double userDirection = roi.params[2];

        // Assign sobel appropriately for the ROIs.
        if (sobelSize == 5) {
            xSobel = {-4,  -5,   0,   5,   4,
                      -8,  -10,  0,  10,   8,
                      -10, -20,  0,  20,  10,
                      -8,  -10,  0,  10,   8,
                      -4,  -5,   0,   5,   4};

            ySobel = { 4,   8,  10,   8,   4,
                       5,  10,  20,  10,   5,
                       0,   0,   0,   0,   0,
                      -5, -10, -20, -10,  -5,
                      -4,  -8, -10,  -8,  -4,};

            normalizeVal = 30293;
        }
        else {
            if (sobelSize != 3) {
                std::cerr << "Error: Invalid sobel size in greyEdgeDetection. Options: 3, 5. Using sobel size 3 by default.\n";
                sobelSize = 3;
                halfSobelSize = sobelSize / 2;
            }

            xSobel = {-1,  0,  1,
                      -2,  0,  2,
                      -1,  0,  1};

            ySobel = { 1,  2,  1, 
                       0,  0,  0,
                      -1, -2, -1};

            normalizeVal = 1443;
        }

        for (int i = roi.y; i < roi.y + roi.sizeY; i++) {
            for (int j = roi.x; j < roi.x + roi.sizeX; j++) {
                int Gx = 0, Gy = 0, pos = 0;
                double length = 0, direction = 0;

                // Sum up the gradient values.
                for (int si = -halfSobelSize; si <= halfSobelSize; si++) {
                    for (int sj = -halfSobelSize; sj <= halfSobelSize; sj++, pos++) {
                        if (source.validPixel(i + si, j + sj)) {
                            Gx += source.getPixel(i + si, j + sj) * xSobel[pos];
                            Gy += source.getPixel(i + si, j + sj) * ySobel[pos]; 
                        }
                    }
                }

                // Calculate the gradient magnitude.
                length = (sqrt((Gx * Gx) + (Gy * Gy)) / normalizeVal) * 255;

                // Calculate the gradient direction.
                if (Gx != 0) {
                    direction = static_cast<int>(atan(Gy / Gx) * radiansToDegrees) + 90;
                }
                else {
                    direction = static_cast<int>(atan(Gy / 0.0001) * radiansToDegrees) + 90;
                }

                // Assign intensity image pixel accordingly.
                edgeIntensity.setPixel(i, j, static_cast<int>(length));

                // Binarize the threshold image.
                if (length < threshold) {
                    edgeThreshold.setPixel(i, j, MINRGB);
                    edgeDirectionThreshold.setPixel(i, j, MINRGB);
                }
                else {
                    edgeThreshold.setPixel(i, j, MAXRGB);

                    // Further binarize the direction output image.
                    if (userDirection - 10 < direction && direction < userDirection + 10) {
                        edgeDirectionThreshold.setPixel(i, j, MAXRGB);
                    }
                    else {
                        edgeDirectionThreshold.setPixel(i, j, MINRGB);
                    }
                }
            }
        }
    }

    edgeIntensity.outputImage(params.outputDir + "Intensity.pgm");
    edgeThreshold.outputImage(params.outputDir + "Threshold.pgm");
    edgeDirectionThreshold.outputImage(params.outputDir + "Direction.pgm");

    return edgeIntensity;
}

/**
 * This function is similar to the greyEdgeDetection function above, however this function
 * is compatible with .ppm formats. The function applies the edge detection and thresholding
 * to each RGB channel individually, and then performs a bitwise OR to combine these edge
 * detection operations into one image. The function also performs an RGB-to-HSI conversion,
 * and then performs the same edge detection to the intensity channel. The function outputs
 * grey images to visualize the edge detection for each RGB channel. The function also converts
 * the HSI image back to RGB and outputs that as well.
 */
RGBImage ImageFunctions::colorEdgeDetection(RGBImage source, FunctionParameters params)
{
    if (source.getFormat() != RGBImage::Format::PPM) {
        throw std::runtime_error("Error: colorEdgeDetection function only compatible with .ppm images.\n");
    }

    // Create output images, and HSI image
    RGBImage edgeThresholdR(source, RGBImage::Channel::RED);
    RGBImage edgeThresholdG(source, RGBImage::Channel::GREEN);
    RGBImage edgeThresholdB(source, RGBImage::Channel::BLUE);
    RGBImage edgeThresholdCombined = source;
    HSIImage hsi = source;
    HSIImage outputHSI = source;

    std::vector<int> xSobel = {-1,  0,  1,
                               -2,  0,  2,
                               -1,  0,  1};

    std::vector<int> ySobel = {-1, -2, -1, 
                                0,  0,  0,
                                1,  2,  1};

    const double radiansToDegrees = 180 / PI;
    const double normalizeVal = 1443;

    for (auto roi : params.rois) {
        int threshold = roi.params[0];

        for (int i = roi.y; i < roi.y + roi.sizeY; i++) {
            for (int j = roi.x; j < roi.x + roi.sizeX; j++) {
                int Gxr = 0, Gxg = 0, Gxb = 0, pos = 0;
                int Gyr = 0, Gyg = 0, Gyb = 0;
                int Gxi = 0, Gyi = 0;
                double lenR = 0, lenG = 0, lenB = 0, lenHSI = 0;

                // Sum up gradient values.
                for (int si = -1; si <= 1; si++) {
                    for (int sj = -1; sj <= 1; sj++, pos++) {
                        if (source.validPixel(i + si, j + sj)) {
                            Gxr += source.getPixel(i + si, j + sj, RGBImage::Channel::RED) * xSobel[pos];
                            Gyr += source.getPixel(i + si, j + sj, RGBImage::Channel::RED) * ySobel[pos];

                            Gxg += source.getPixel(i + si, j + sj, RGBImage::Channel::GREEN) * xSobel[pos];
                            Gyg += source.getPixel(i + si, j + sj, RGBImage::Channel::GREEN) * ySobel[pos];

                            Gxb += source.getPixel(i + si, j + sj, RGBImage::Channel::BLUE) * xSobel[pos];
                            Gyb += source.getPixel(i + si, j + sj, RGBImage::Channel::BLUE) * ySobel[pos];

                            Gxi += hsi.getIntensity(i + si, j + sj) * xSobel[pos];
                            Gyi += hsi.getIntensity(i + si, j + sj) * ySobel[pos];
                        }
                    }
                }

                // Calculate the gradient magnitude.
                lenR = (sqrt((Gxr * Gxr) + (Gyr * Gyr)) / normalizeVal) * 255;
                lenG = (sqrt((Gxg * Gxg) + (Gyg * Gyg)) / normalizeVal) * 255;
                lenB = (sqrt((Gxb * Gxb) + (Gyb * Gyb)) / normalizeVal) * 255;
                lenHSI = (sqrt((Gxi * Gxi) + (Gyi * Gyi)) / normalizeVal) * 255;

                // Binarize each output image.
                if (lenR < threshold) {
                    edgeThresholdR.setPixel(i, j, MINRGB);
                }
                else {
                    edgeThresholdR.setPixel(i, j, MAXRGB);
                }

                if (lenG < threshold) {
                    edgeThresholdG.setPixel(i, j, MINRGB);
                }
                else {
                    edgeThresholdG.setPixel(i, j, MAXRGB);
                }

                if (lenB < threshold) {
                    edgeThresholdB.setPixel(i, j, MINRGB);
                }
                else {
                    edgeThresholdB.setPixel(i, j, MAXRGB);
                }

                if (lenHSI < threshold) {
                    outputHSI.setIntensity(i, j, MINRGB);
                }
                else {
                    outputHSI.setIntensity(i, j, MAXRGB);
                }
                outputHSI.setSaturation(i, j, 0);

                // Perform logical OR for combination image.
                int r = edgeThresholdR.getPixel(i, j);
                int g = edgeThresholdG.getPixel(i, j);
                int b = edgeThresholdB.getPixel(i, j);
                if (r | g | b) {
                    edgeThresholdCombined.setPixel(i, j, RGBImage::Channel::RED, MAXRGB);
                    edgeThresholdCombined.setPixel(i, j, RGBImage::Channel::GREEN, MAXRGB);
                    edgeThresholdCombined.setPixel(i, j, RGBImage::Channel::BLUE, MAXRGB);
                }
                else {
                    edgeThresholdCombined.setPixel(i, j, RGBImage::Channel::RED, MINRGB);
                    edgeThresholdCombined.setPixel(i, j, RGBImage::Channel::GREEN, MINRGB);
                    edgeThresholdCombined.setPixel(i, j, RGBImage::Channel::BLUE, MINRGB);
                }
            }
        }
    }

    edgeThresholdR.outputImage(params.outputDir + "R.pgm");
    edgeThresholdG.outputImage(params.outputDir + "G.pgm");
    edgeThresholdB.outputImage(params.outputDir + "B.pgm");
    edgeThresholdCombined.outputImage(params.outputDir + "Combined.ppm");
    outputHSI.convertToRGB().outputImage(params.outputDir + "HSI.ppm");

    return source;
}