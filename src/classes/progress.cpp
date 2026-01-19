
#include <iostream>
#include "misc.h"
#include "progress.h"

using namespace std;

void Progressbar::set_color(int r, int g, int b)
{
    rbase = max(0, min(r, 255));
    gbase = max(0, min(g, 255));
    bbase = max(0, min(b, 255));
}

void Progressbar::update(float value)
{
    // value = value/(poses) + (float)(pose-1)*100.0/poses;
    float percentage = (value-minimum)/(maximum-minimum)*width;
    if (percentage > width) percentage = width;
    cout << "\033[A|";
    int i;
    bool grayyet = false;
    for (i=0; i<width; i++)
    {
        float cmpi = i;
        if (cmpi <= percentage)
        {
            float h = M_PI*2 * cmpi / 81 + hueoffset;
            int r, g, b;
            r = rbase +  24 * sin(h-0.333);
            g = gbase +  26 * sin(h+0.333);
            b = bbase +  31 * sin(h);
            colorrgb(r, g, b);
            cout << "\u2593";
            colorless();
        }
        else
        {
            if (!grayyet) colorrgb(64, 68, 70);
            cout << "\u2591";
            grayyet = true;
        }
    }
    colorless();
    cout << ("|/-\\")[spinchr/4] << " " << (int)percentage << "%.               " << endl;
    spinchr++;
    if (spinchr >= 16) spinchr = 0;
    hueoffset += 0.1;
}

void Progressbar::erase()
{
    cout << "\033[A\033[K";
}