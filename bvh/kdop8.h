#ifndef KDOP24_H
#define KDOP24_H

namespace icy { class kDOP8; }

class icy::kDOP8
{
public:
    float d[4], g[4]; // 0-3 lower boundaries; 4-7 higher boundaries
    float ctrX, ctrY;
    float dX, dY;

    void Reset();
    bool Overlaps(kDOP8 &b);
    void Expand(float x, float y);
    void Expand(kDOP8 &b);
    void ExpandBy(float radius);

private:
    inline void MinMax(float p, float &mi, float &ma);
};

#endif // KDOP24_H
