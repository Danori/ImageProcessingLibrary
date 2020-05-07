#include "HSIImage.hpp"

#include <cmath>
#include <algorithm>
#include <iostream>

HSIImage::HSIImage(const unsigned int height, const unsigned int width)
{
    resize(height, width);
}

HSIImage::HSIImage(RGBImage &image)
{
    convertFromRGB(image);
}

HSIImage::HSIImage(const HSIImage &image)
{
    copyImage(image);
}

RGBImage HSIImage::convertToRGB()
{
    const double pi = 3.14159265358979323846;

    RGBImage output(height, width, RGBImage::Format::PPM);

    double x = 0, y = 0, z = 0;
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            double hue = (getHue(i, j) * pi) / 180;
            double saturation = getSaturation(i, j) / 100;
            double intensity = getIntensity(i, j) / 255;

            if (hue < ((2 * pi) / 3)) {
                x = intensity * (1 - saturation);
                y = intensity * (1 + ((saturation * cos(hue)) / (cos((pi / 3) - hue))));
                z = (3 * intensity) - (x + y);
                output.setPixel(i, j, RGBImage::Channel::BLUE, static_cast<unsigned int>(x * 255));
                output.setPixel(i, j, RGBImage::Channel::RED, static_cast<unsigned int>(y * 255));
                output.setPixel(i, j, RGBImage::Channel::GREEN, static_cast<unsigned int>(z * 255));
            }
            else if (((2 * pi) / 3) <= hue && hue < ((4 * pi) / 3)) {
                hue = hue - ((2 * pi) / 3);
                x = intensity * (1 - saturation);
                y = intensity * (1 + ((saturation * cos(hue)) / (cos((pi / 3) - hue))));
                z = (3 * intensity) - (x + y);
                output.setPixel(i, j, RGBImage::Channel::RED, static_cast<unsigned int>(x * 255));
                output.setPixel(i, j, RGBImage::Channel::GREEN, static_cast<unsigned int>(y * 255));
                output.setPixel(i, j, RGBImage::Channel::BLUE, static_cast<unsigned int>(z * 255));
            }
            else if (((4 * pi) / 3) <= hue && hue < (2 * pi)) {
                hue = hue - ((4 * pi) / 3);
                x = intensity * (1 - saturation);
                y = intensity * (1 + ((saturation * cos(hue)) / (cos((pi / 3) - hue))));
                z = (3 * intensity) - (x + y);
                output.setPixel(i, j, RGBImage::Channel::GREEN, static_cast<unsigned int>(x * 255));
                output.setPixel(i, j, RGBImage::Channel::BLUE, static_cast<unsigned int>(y * 255));
                output.setPixel(i, j, RGBImage::Channel::RED, static_cast<unsigned int>(z * 255));
            }
        }
    }

    return output;
}

void HSIImage::convertFromRGB(RGBImage &image)
{
    const double pi = 3.14159265358979323846;

    resize(image.getHeight(), image.getWidth());

    double hue = 0, saturation = 0, intensity = 0;
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            double r = image.getPixel(i, j, RGBImage::Channel::RED);
            double g = image.getPixel(i, j, RGBImage::Channel::GREEN);
            double b = image.getPixel(i, j, RGBImage::Channel::BLUE);

            double normR = r / (r + g + b);
            double normG = g / (r + g + b);
            double normB = b / (r + g + b);

            if (normB <= normG) {
                hue = acos((0.5 * ((r - g) + (r - b))) / (sqrt(((r - g) * (r - g)) + ((r- b) * (g - b)))));
            }
            else {
                hue = (2 * pi) - acos((0.5 * ((r - g) + (r - b))) / (sqrt(((r - g) * (r - g)) + ((r- b) * (g - b)))));
            }

            saturation = 1 - (3 * std::min({normR, normG, normB}));

            intensity = (r + g + b) / (3 * 255);

            setHue(i, j, (hue * 180) / pi);
            setSaturation(i, j, saturation * 100);
            setIntensity(i, j, intensity * 255);
        }
    }
}

void HSIImage::deleteImage()
{
    height = width = 0;

    clearChannels();
}

void HSIImage::copyImage(const HSIImage &image)
{
    height = image.height;
    width = image.width;

    copyChannels(image);
}

void HSIImage::resize(const unsigned int height, const unsigned int width)
{
    this->height = height;
    this->width = width;

    clearChannels();
    resizeChannels(height * width);
}

bool HSIImage::validPixel(int i, int j)
{
    return ((0 <= i) && (i < height) && (0 <= j) && (j < width));
}

void HSIImage::setHue(const unsigned int i, const unsigned int j, const unsigned int value)
{
    hueChannel[pos(i,j)] = value;
}

void HSIImage::setSaturation(const unsigned int i, const unsigned int j, const unsigned int value)
{
    satChannel[pos(i,j)] = value;
}

void HSIImage::setIntensity(const unsigned int i, const unsigned int j, const unsigned int value)
{
    intChannel[pos(i,j)] = value;
}

double HSIImage::getHue(int i, int j)
{
    return hueChannel[pos(i,j)];
}

double HSIImage::getSaturation(int i, int j)
{
    return satChannel[pos(i,j)];
}

double HSIImage::getIntensity(int i, int j)
{
    return intChannel[pos(i,j)];
}

void HSIImage::clearChannels()
{
    hueChannel.clear();
    satChannel.clear();
    intChannel.clear();
}

void HSIImage::copyChannels(const HSIImage &image)
{
    hueChannel = image.hueChannel;
    satChannel = image.satChannel;
    intChannel = image.intChannel;
}

void HSIImage::resizeChannels(const unsigned int length)
{
    hueChannel.resize(length);
    satChannel.resize(length);
    intChannel.resize(length);
}

const unsigned int HSIImage::pos(const unsigned int i, const unsigned int j)
{
    return (i * width) + j;
}