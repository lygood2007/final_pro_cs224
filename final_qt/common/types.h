/** types.h
 ** Brief: This is the header file of all the necessary types in the project
 ** Project: large-scale fluids
 ** Date: 04/10/2013
 ** Member: Scott, Hobarts, Yan Li
 **/
#ifndef TYPES_H
#define TYPES_H

struct Colorf
{
    Colorf() {r = g = b = 0.f; a= 1.f;}
    Colorf( float pr, float pg, float pb, float pa) {r = pr; g = pg; b = pb; a = pa;}
    float r;
    float g;
    float b;
    float a;
};

struct Index2D
{
    int indRow;
    int indCol;
};

struct Index1D
{
    unsigned int ind;
};

union Tri
{
    struct
    {
        // 2D index
        Index2D a2D;
        Index2D b2D;
        Index2D c2D;
    };

    struct
    {
        // 1D index
        Index1D a1D;
        Index1D b1D;
        Index1D c1D;
    };
};

#endif // TYPES_H
