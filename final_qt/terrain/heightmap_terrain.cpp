/** heightmap_terrain.cpp
 ** Brief: The source file of HeightmapTerrain Class
 ** Project: large-scale fluids
 ** Date: 04/10/2013
 ** Member: Scott, Hobarts, Yan Li
 **/
#include "heightmap_terrain.h"

HeightmapTerrain::HeightmapTerrain()
    :Terrain(), m_filename(HEIGHTMAP_FILENAME)
{
}

HeightmapTerrain::HeightmapTerrain(QString filename)
    :Terrain(), m_filename(filename)
{
}

HeightmapTerrain::~HeightmapTerrain()
{
    // No heap. empty
}

void HeightmapTerrain::populateTerrain()
{
    //get the file
    QFile file(m_filename);
    assert(file.exists() && "The height does not exist");

    //load the file
    QImage image;
    image.load(file.fileName());

    //iterate through grid
    for(int i = 0; i < m_gridLength; i++){
        for(int j = 0; j < m_gridLength; j++){
            //interpolate the x, y, z positions
            double iProp = (double)i / (double)(m_gridLength - 1);
            double jProp = (double)j / (double)(m_gridLength - 1);

            double currX = determinePosition(jProp);
            double currZ = determinePosition(iProp);
            double currY = interpolateHeight(image, jProp, iProp);

            //set the values
            Vector3 &currVertex = m_vertices[getIndex(i, j)];
            currVertex.x = currX;
            currVertex.y = currY;
            currVertex.z = currZ;
        }
    }
}

double HeightmapTerrain::determinePosition(double gridPosition)
{
    return -(TERRAIN_BOUND - (2.0 * gridPosition * TERRAIN_BOUND));
}

double HeightmapTerrain::interpolateHeight(QImage heightMap, double x, double y)
{
    //get width and height of the map
    int width = heightMap.width();
    int height = heightMap.height();

    //get the current x and y positions
    double currX = x * (double)width;
    double currY = y * (double)height;

    //get the indices of the map
    double lowX = min(max(floor(currX), 0), width - 1);
    double highX = min(max(ceil(currX), 0), width - 1);
    double lowY = min(max(floor(currY), 0), height - 1);
    double highY = min(max(ceil(currY), 0), height - 1);

    //get the proportions
    double xDist = 1 - (currX - lowX);
    double yDist = 1 - (currY - lowY);

    int lowXlowYgray = qGray(heightMap.pixel((int)lowX, (int)lowY));
    int lowXhighYgray = qGray(heightMap.pixel((int)lowX, (int)highY));
    int highXlowYgray = qGray(heightMap.pixel((int)highX, (int)lowY));
    int highXhighYgray = qGray(heightMap.pixel((int)highX, (int)highY));

    double grayProp = (xDist * yDist * (double)lowXlowYgray) +
            (xDist * (1 - yDist) * (double)lowXhighYgray) +
            ((1 - xDist) * yDist * (double)highXlowYgray) +
            ((1 - xDist) * (1 - yDist) * (double)highXhighYgray);
    grayProp = grayProp / 255.0;

    return TERRAIN_MIN_HEIGHT + (grayProp * (TERRAIN_MAX_HEIGHT - TERRAIN_MIN_HEIGHT));
}
