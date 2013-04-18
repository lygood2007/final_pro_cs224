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

#endif // TYPES_H
