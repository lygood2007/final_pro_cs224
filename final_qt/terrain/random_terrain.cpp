/** random_terrain.cpp
 ** Brief: The source file of RandomTerrain Class
 ** Project: large-scale fluids
 ** Date: 04/15/2013
 ** Member: Scott, Hobarts, Yan Li
 **/

#include "random_terrain.h"

RandomTerrain::RandomTerrain(const int decay, const float roughness) : Terrain()
{
   m_decay = decay;
   m_roughness  = roughness;
}

RandomTerrain::~RandomTerrain()
{
    // No heap. empty
}

/**
 * Computes the amount to perturb the height of the vertex currently being processed.
 * Feel free to modify this.
 *
 * @param depth The current recursion depth
 */
double RandomTerrain::getPerturb(const int cur_depth) const
{
    /** roughness = (cur_depth/m_depth).^m_decay*k, where k is between -1 and 1*/
    return m_roughness
           * pow((double)cur_depth / m_depth, m_decay)
           * ((rand() % 200-100) / 100.0);
}

/**
* Subdivides a square by finding the vertices at its corners, the midpoints of each side, and
* the center (as the algorithm describes). Then recurs on each of the four sub-squares created.
*
* @param topLeft The grid coordinate of the top-left corner of the square to subdivide
* @param bottomRight The grid coordinate of the bottom-right corner of the square to subdivide
* @param depth The current recursion depth, decreasing as this function recurses deeper. The
*              function should stop recurring when this value reaches zero.
*/
void RandomTerrain::subdivideSquare(GridIndex topleft, GridIndex botright, GLint curDepth)
{

    // TL--TM--TR    +---> x
    // |   |   |     |
    // ML--MM--MR    V
    // |   |   |     y
    // BL--BM--BR

    // corner coordinates (in the grid space [x,y])
    GridIndex TL = GridIndex(topleft.x, topleft.y);
    GridIndex TR = GridIndex(botright.x, topleft.y);
    GridIndex BL = GridIndex(topleft.x, botright.y);
    GridIndex BR = GridIndex(botright.x, botright.y);

    // corner vertices on the terrain (in the grid space [x,y,z])
    Vector3 &vTL = m_vertices[getIndex(TL)];
    Vector3 &vTR = m_vertices[getIndex(TR)];
    Vector3 &vBL = m_vertices[getIndex(BL)];
    Vector3 &vBR = m_vertices[getIndex(BR)];

    if( curDepth  == 0 )
        return ;

    GridIndex ML = GridIndex( topleft.x, (topleft.y + botright.y)/2.0 );
    GridIndex TM = GridIndex( (topleft.x + botright.x)/2.0, topleft.y );
    GridIndex MR = GridIndex( botright.x, (topleft.y + botright.y )/2.0 );
    GridIndex BM = GridIndex( (topleft.x + botright.x )/2.0, botright.y );
    GridIndex MM = GridIndex( (topleft.x + botright.x )/2.0, (topleft.y + botright.y )/2.0 );


    Vector3 &vML = m_vertices[getIndex(ML)];
    Vector3 &vTM = m_vertices[getIndex(TM)];
    Vector3 &vMR = m_vertices[getIndex(MR)];
    Vector3 &vBM = m_vertices[getIndex(BM)];
    Vector3 &vMM = m_vertices[getIndex(MM)];

    vML.x = 0.5*(vTL.x + vBL.x);
    vML.z = 0.5*(vTL.z + vBL.z);
    vML.y = 0.5*(vTL.y + vBL.y );

    vTM.x = 0.5*(vTL.x + vTR.x);
    vTM.z = 0.5*(vTL.z + vTR.z);
    vTM.y = 0.5*( vTL.y + vTR.y );

    vMR.x = 0.5*(vTR.x + vBR.x);
    vMR.z = 0.5*(vTR.z + vBR.z);
    vMR.y = 0.5*( vTR.y + vBR.y );

    vBM.x = 0.5*(vBL.x + vBR.x);
    vBM.z = 0.5*(vBL.z + vBR.z);
    vBM.y = 0.5*( vBL.y + vBR.y );

    vMM.x = 0.25*( vTL.x + vTR.x + vBL.x + vBR.x);
    vMM.z = 0.25*( vTL.z + vTR.z + vBL.z + vBR.z);
    vMM.y = 0.25*( vTL.y + vTR.y + vBL.y + vBR.y ) + getPerturb( curDepth);


    subdivideSquare( TL, MM, curDepth-1 );
    subdivideSquare( TM, MR, curDepth-1 );
    subdivideSquare( MM, BR, curDepth-1 );
    subdivideSquare( ML, BM, curDepth-1 );
}

void RandomTerrain::populateTerrain()
{
    int avgHeight = (TERRAIN_MAX_HEIGHT + TERRAIN_MIN_HEIGHT)*0.5;
    Vector3 tl(-m_bound, avgHeight, -m_bound);
        Vector3 tr(m_bound, avgHeight, -m_bound);
        Vector3 bl(-m_bound, avgHeight, m_bound);
        Vector3 br(m_bound, avgHeight, m_bound);
        GridIndex tlg(0,0);
        GridIndex trg(0,m_gridLength-1);
        GridIndex blg(m_gridLength-1, 0);
        GridIndex brg(m_gridLength-1, m_gridLength-1);
        m_vertices[getIndex(tlg)] = tl;
        m_vertices[getIndex(trg)] = tr;
        m_vertices[getIndex(blg)] = bl;
        m_vertices[getIndex(brg)] = br;
        subdivideSquare(tlg, brg, m_depth);
}
