#ifndef RGB_IMAGE_HPP
#define RGB_IMAGE_HPP

#include <vector>
#include <string>

class RGBImage
{
public:
    enum Format
    {
        PGM = 0,
        PPM = 1
    };

    enum Channel
    {
        RED     = 0,
        GREEN   = 1,
        BLUE    = 2,
        INVALID = 3
    };

private:
    const unsigned int pos(const unsigned int i, const unsigned int j);
    void clearChannels();
    void copyChannels(const RGBImage &image);
    void resizeChannels(const unsigned int length);

    std::vector<unsigned int> redChannel;
    std::vector<unsigned int> greenChannel;
    std::vector<unsigned int> blueChannel;
    std::vector<unsigned int> greyChannel;
    unsigned int height;
    unsigned int width;
    Format format;

public:
    RGBImage() : height(0), width(0), format(PGM) { }
    RGBImage(const std::string &fileName);
    RGBImage(const unsigned int height, const unsigned int width, const Format format);
    RGBImage(const RGBImage &image, Channel rgb);
    RGBImage(const RGBImage &image);

    void deleteImage();
    void copyImage(const RGBImage &image);
    void resize(const unsigned int height, const unsigned int width);

    bool validPixel(int i, int j);

    void setPixel(const unsigned int i, const unsigned int j, const unsigned int value);
    void setPixel(const unsigned int i, const unsigned int j, const unsigned int rgb, const unsigned int value);

    const unsigned int getHeight();
    const unsigned int getWidth();
    const Format getFormat();
    unsigned int getPixel(int i, int j);
    unsigned int getPixel(int i, int j, const int rgb);

    void outputImage(const std::string &fileName);
    void initImage(const std::string &fileName);
};

#endif // RGB_IMAGE_HPP