/** texture_loader.h
 ** Brief: The header file of terrain Class
 ** Project: large-scale fluids
 ** Date: 04/10/2013
 ** Member: Scott, Hobarts, Yan Li
 **/

#ifndef TERRAIN_H
#define TERRAIN_H

#include <QList>
#include <qgl.h>
#include "vector.h"

//debugging
#include <iostream>

//added by hcreynol
#include <QFile>
#include <assert.h>

#define DEFAULT_DEPTH 8
#define DEFAULT_DECAY 2
#define DEFAULT_ROUGHNESS 2
#define TEXTURE_DIR "./resource/terrain.jpg"

typedef Vector2 GridIndex;

extern GLuint loadTexture(const QString &filename);

class Terrain{
public:

    Terrain();
    Terrain(const int decay, const float roughness, const int depth, const bool renderNormals);
    Terrain(QString filename); //added by hcreynol
    ~Terrain();

    /**
     * Initialize the variables in the class
     */
    void init(const int decay = DEFAULT_DEPTH, const int depth = DEFAULT_DECAY,
              const float roughness = DEFAULT_ROUGHNESS, const bool renderNormals = false);

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
     * Computes the amount to perturb the height of the vertex currently being processed.
     * Feel free to modify this.
     *
     * @param depth The current recursion depth
     */
    double getPerturb(const int cur_depth) const;
    /**
     * Retrieves the position of each neighbor of the given grid coordinate (i.e. all grid
     * coordinates that can be found one unit horizontally, vertically, or diagonally from
     * the specified grid coordinate).
     *
     * @param coordinate The grid coordinate whose neighbors are to be retrieved
     */
    QList<Vector3*> getSurroundingVertices(const GridIndex &coordinate) const;
    /**
    * Subdivides a square by finding the vertices at its corners, the midpoints of each side, and
    * the center (as the algorithm describes). Then recurs on each of the four sub-squares created.
    *
    * @param topLeft The grid coordinate of the top-left corner of the square to subdivide
    * @param bottomRight The grid coordinate of the bottom-right corner of the square to subdivide
    * @param depth The current recursion depth, decreasing as this function recurses deeper. The
    *              function should stop recurring when this value reaches zero.
    */
   void subdivideSquare(const GridIndex tlg, const GridIndex brg, GLint curDepth);
   /**
   * Computes the normal vector of each terrain vertex and stores it in the corresponding vertex.
   */
  void computeNormals();
    /**
     * Sets default values for the four corners of the terrain grid and calls subdivideSquare()
     * to begin the terrain generation process. You do not need to modify this function.
     */
    void populateTerrain();

    void populateTerrainFromHeightmap(); //added by hcreynol
    double interpolateHeight(QImage heightMap, double x, double y); //added by hcreynol

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
     ** Load the texture into memory.
     **/
    void loadTextureToTerrain();

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

protected:
private:

    Vector3* m_vertices; // The vertices
    Vector3* m_normals; // The normals

    int m_decay; // Controls how much heights can vary per recursion depth level. Higher values generate smoother terrain.
    float m_roughness; // Controls the height of the terrain. Large values generate higher peaks and lower valleys.
    int m_depth; // The number of recursion levels to use to generate terrain. Can be used as a level-of-detail parameter.
    bool m_renderNormals; // Flag for rendering normals
    int m_gridLength; // The grid length
    GLuint m_textureId; // The texture id

    // added by hcreynol
    QString m_filename; // filename of the heightmap file
};

#endif // TERRAIN_H
