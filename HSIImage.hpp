#ifndef HSI_IMAGE_HPP
#define HSI_IMAGE_HPP

#include "RGBImage.hpp"

#include <vector>

class HSIImage
{
private:
    const unsigned int pos(const unsigned int i, const unsigned int j);
    void clearChannels();
    void copyChannels(const HSIImage &image);
    void resizeChannels(const unsigned int length);

    std::vector<double> hueChannel;
    std::vector<double> satChannel;
    std::vector<double> intChannel;
    unsigned int height;
    unsigned int width;

public:
    HSIImage() : height(0), width(0)  { }
    HSIImage(const unsigned int height, const unsigned int width);
    HSIImage(RGBImage &image);
    HSIImage(const HSIImage &image);

    RGBImage convertToRGB();
    void convertFromRGB(RGBImage &image);

    void deleteImage();
    void copyImage(const HSIImage &image);
    void resize(const unsigned int height, const unsigned int width);

    bool validPixel(int i, int j);

    void setHue(const unsigned int i, const unsigned int j, const unsigned int value);
    void setSaturation(const unsigned int i, const unsigned int j, const unsigned int value);
    void setIntensity(const unsigned int i, const unsigned int j, const unsigned int value);

    double getHue(int i, int j);
    double getSaturation(int i, int j);
    double getIntensity(int i, int j);
};

#endif // HSI_IMAGE_HPP