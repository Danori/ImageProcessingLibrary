#include "RGBImage.hpp"

#include <fstream>
#include <iostream>
#include <cstring>

RGBImage::RGBImage(const std::string &fileName)
{
    initImage(fileName);
}

RGBImage::RGBImage(const unsigned int height, const unsigned int width, const Format format)
{
    this->format = format;
    resize(height, width);
}

RGBImage::RGBImage(const RGBImage &image, Channel rgb)
{
    format = PGM;
    resize(image.height, image.width);

    switch (rgb) {
        case RED:
            greyChannel = image.redChannel;
        break;

        case GREEN:
            greyChannel = image.greenChannel;
        break;

        case BLUE:
            greyChannel = image.blueChannel;
        break;
    }
}

RGBImage::RGBImage(const RGBImage &image)
{
    copyImage(image);
}

void RGBImage::deleteImage()
{
    height = width = 0;

    clearChannels();
}

void RGBImage::copyImage(const RGBImage &image)
{
    height = image.height;
    width = image.width;
    format = image.format;

    copyChannels(image);
}

void RGBImage::resize(const unsigned int height, const unsigned int width)
{
    this->height = height;
    this->width = width;

    clearChannels();
    resizeChannels(height * width);
}

bool RGBImage::validPixel(int i, int j)
{
    return ((0 <= i) && (i < height) && (0 <= j) && (j < width));
}

void RGBImage::setPixel(const unsigned int i, const unsigned int j, const unsigned int value)
{
    greyChannel[pos(i,j)] = value;
}

void RGBImage::setPixel(const unsigned int i, const unsigned int j, const unsigned int rgb, const unsigned int value)
{
    switch (rgb) {
        case RED:
            redChannel[pos(i,j)] = value;
        break;

        case GREEN:
            greenChannel[pos(i,j)] = value;
        break;

        case BLUE:
            blueChannel[pos(i,j)] = value;
        break;
    }
}

const unsigned int RGBImage::getHeight()
{
    return height;
}

const unsigned int RGBImage::getWidth()
{
    return width;
}

const RGBImage::Format RGBImage::getFormat()
{
    return format;
}

unsigned int RGBImage::getPixel(int i, int j)
{
    return greyChannel[pos(i,j)];
}

unsigned int RGBImage::getPixel(int i, int j, const int rgb)
{
    switch (rgb) {
        case RED:
            return redChannel[pos(i,j)];
        break;

        case GREEN:
            return greenChannel[pos(i,j)];
        break;

        case BLUE:
            return blueChannel[pos(i,j)];
        break;
    }
}

void RGBImage::outputImage(const std::string &fileName)
{
    std::ofstream output(fileName);

    if (output) {
        switch (format) {
            case PGM:
                output << "P5\n";
                output << width << " " << height << "\n";
                output << "255\n";

                for (int i = 0; i < height; i++) {
                    for (int j = 0; j < width; j++) {
                        output << static_cast<unsigned char>(greyChannel[pos(i,j)]);
                    }
                }
            break;

            case PPM:
                output << "P6\n";
                output << width << " " << height << "\n";
                output << "255\n";

                for (int i = 0; i < height; i++) {
                    for (int j = 0; j < width; j++) {
                        output << static_cast<unsigned char>(redChannel[pos(i,j)]);
                        output << static_cast<unsigned char>(greenChannel[pos(i,j)]);
                        output << static_cast<unsigned char>(blueChannel[pos(i,j)]);
                    }
                }
            break;

            default:

            break;
        }
    }
}

void RGBImage::initImage(const std::string &fileName)
{
    deleteImage();

    std::ifstream input(fileName);
    std::string strBuffer;

    if (fileName.find(".pgm") != std::string::npos) {
        format = PGM;
    }
    else if (fileName.find(".ppm") != std::string::npos) {
        format = PPM;
    }
    else {
        std::cerr << "Error: ImageProcessing only compatible with .pgm or .ppm images.\n";
        return;
    }

    if (input) {
        std::getline(input, strBuffer);

        if (isdigit(input.peek())) {
            input >> width >> height;
        }
        else {
            std::getline(input, strBuffer);
            input >> width >> height;
        }

        resizeChannels(width * height);

        std::getline(input, strBuffer);
	    std::getline(input, strBuffer);

        switch (format) {
            case PGM: {
                char *buffer = new char[width * height];
                input.read(buffer, width * height);

                for (int i = 0; i < height; i++) {
                    for (int j = 0; j < width; j++) {
                        greyChannel[pos(i,j)] = static_cast<unsigned char>(buffer[pos(i,j)]);
                    }
                }

                delete[] buffer;
            }
            break;

            case PPM: {
                char *buffer = new char[width * height * 3];
                input.read(buffer, width * height * 3);

                for (int i = 0; i < height; i++) {
                    for (int j = 0; j < width; j++) {
                        redChannel[pos(i,j)] = static_cast<unsigned char>(buffer[pos(i*3,j*3)]);
                        greenChannel[pos(i,j)] = static_cast<unsigned char>(buffer[pos(i*3,j*3)+1]);
                        blueChannel[pos(i,j)] = static_cast<unsigned char>(buffer[pos(i*3,j*3)+2]);
                    }
                }

                delete[] buffer;
            }
            break;

            default:

            break;
        }
    }
    else {
        std::cerr << "Error: Failed to open " << fileName << " in initImage.\n";
    }
}

// Helper Functions

void RGBImage::clearChannels()
{
    redChannel.clear();
    greenChannel.clear();
    blueChannel.clear();
    greyChannel.clear();
}

void RGBImage::copyChannels(const RGBImage &image)
{
    redChannel = image.redChannel;
    greenChannel = image.greenChannel;
    blueChannel = image.blueChannel;
    greyChannel = image.greyChannel;
}

void RGBImage::resizeChannels(const unsigned int length)
{
    switch (format) {
        case PGM:
            greyChannel.resize(length);
        break;

        case PPM:
            redChannel.resize(length);
            greenChannel.resize(length);
            blueChannel.resize(length);
        break; 
    }
}

const unsigned int RGBImage::pos(unsigned int i, unsigned int j)
{
    return (i * width) + j;
}