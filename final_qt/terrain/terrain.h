/** terrain.h
 ** Brief: The header file of terrain Class
 ** Project: large-scale fluids
 ** Date: 04/15/2013
 ** Member: Scott, Hobarts, Yan Li
 **/

#ifndef TERRAIN_H
#define TERRAIN_H

#include <QList>
#include <qgl.h>
#include "vector.h"
//debugging
#include <iostream>
#include <assert.h>

#define DEFAULT_DEPTH 6
#define TERRAIN_BOUND 50
#define TERRAIN_MIN_HEIGHT -10
#define TERRAIN_MAX_HEIGHT 20

#define TEXTURE_DIR "./resource/terrain.jpg"

typedef Vector2 GridIndex;

extern GLuint loadTexture(const QString &filename);

class Terrain{
public:

    Terrain();
    //Terrain(const int decay, const float roughness, const int depth, const bool renderNormals);
    //Terrain(QString filename); //added by hcreynol
    virtual ~Terrain();

    /**
     * Draws a line at each vertex showing the direction of that vertex's normal. You may find
     * this to be a useful tool if you're having trouble getting the lighting to look right.
     * By default, this function is called in paintGL(), but only renders anything if
     * m_renderNormals is true. You do not need to modify this function.
     */
    void drawNormals() const;

    /**
     * Render the terrain to screen.
     */
    void draw() const;

    /**
     * Generate the terrain.
     */
    void generate();

    /**
     * Enable drawing normals
     */
    inline void enableNormal(){m_renderNormals = true;}
    /**
     * Disable drawing normals
     */
    inline void disableNormal(){m_renderNormals = false;}
    /**
     * Get m_renderNormal
     */
    inline bool isRenderingNormal() const {return m_renderNormals;}

private:

protected:

    /**
     * Converts a grid coordinate (row, column) to an index into a 1-dimensional array.
     * Can be used to index into m_terrain or m_normalMap.
     * Returns -1 if the grid coordinate entered is not valid.
     */
   inline int getIndex(const GridIndex &c) const;
   /**
     * Converts a grid coordinate (row, column) to an index into a 1-dimensional array.
     * Can be used to index into m_terrain or m_normalMap.
     * Returns -1 if the grid coordinate entered is not valid.
     */
   inline int getIndex(const int row, const int column) const;

   /**
    * Retrieves the position of each neighbor of the given grid coordinate (i.e. all grid
    * coordinates that can be found one unit horizontally, vertically, or diagonally from
    * the specified grid coordinate).
    *
    * @param coordinate The grid coordinate whose neighbors are to be retrieved
    */
   QList<Vector3*> getSurroundingVertices(const GridIndex &coordinate) const;

   /**
    * Populate the terrain, need to be implemented in children classes
    */
   virtual  void populateTerrain() = 0;
    /**
    * Computes the normal vector of each terrain vertex and stores it in the corresponding vertex.
    */
   void computeNormals();

   /**
    ** Load the texture into memory.
    **/
   void loadTextureToTerrain();

protected:
    Vector3* m_vertices; // The vertices
    Vector3* m_normals; // The normals

    int m_depth; // The number of recursion levels to use to generate terrain. Can be used as a level-of-detail parameter.
    bool m_renderNormals; // Flag for rendering normals
    int m_gridLength; // The grid length
    GLuint m_textureId; // The texture id
};

#endif // TERRAIN_H
