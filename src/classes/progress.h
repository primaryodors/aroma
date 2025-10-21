
#ifndef _PROGRESS
#define _PROGRESS

class Progressbar
{
    public:
    float minimum = 0;
    float maximum = 100;
    int width = 100;

    void set_color(int r, int g, int b);
    void update(float value);
    void erase();

    protected:
    int rbase = 98, gbase = 176, bbase = 224;       // PrimaryOdors blue.
    int spinchr = 0;
    float hueoffset = 0;
};

#endif