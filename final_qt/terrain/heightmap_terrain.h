/** heightmap_terrain.h
 ** Brief: The header file of HeightmapTerrain Class
 ** Project: large-scale fluids
 ** Date: 04/15/2013
 ** Member: Scott, Hobarts, Yan Li
 **/
#ifndef HEIGHTMAPTERRAIN_H
#define HEIGHTMAPTERRAIN_H

#include "terrain.h"

//added by hcreynol
#include <QFile>

#define HEIGHTMAP_FILENAME "./resource/h_map_1.jpg"
//added by hcreynol
#define DEFAULT_GRID_BOUNDS 10
#define DEFAULT_MIN_HEIGHT -3
#define DEFAULT_MAX_HEIGHT 5

class HeightmapTerrain : public Terrain
{
public:
    HeightmapTerrain();
    HeightmapTerrain(QString filename);
    virtual ~HeightmapTerrain();

private:

protected:
    virtual void populateTerrain();
    double determinePosition(double gridPosition); //added by hcreynol
    double interpolateHeight(QImage heightMap, double x, double y); //added by hcreynol
    // added by hcreynol

protected:
    QString m_filename; // filename of the heightmap file
};

#endif // HEIGHTMAPTERRAIN_H
