/** terrain.cpp
 ** Brief: The source file of terrain.h
 ** Project: large-scale fluids
 ** Date: 04/10/2013
 ** Member: Scott, Hobarts, Yan Li
 **/
#include "terrain.h"

Terrain::Terrain()
{
    init();
}

Terrain::Terrain( const int decay, const float roughness, const int depth, const bool renderNormals)
{
    init(decay, roughness, depth, renderNormals);
}

Terrain::~Terrain()
{
    if( m_normals )
        delete []m_normals;
    if( m_vertices )
        delete []m_vertices;

    if( m_textureId != 0 )
        glDeleteTextures(1,&m_textureId);
}

/**
 * Initialize the variables in the class
 */
void Terrain::init(const int decay, const int depth,
                   const float roughness, const bool renderNormals)
{
    m_decay = 2;
    m_depth = 6;
    m_roughness  = 5;
    m_gridLength = (1<<m_depth)+1;
    int terrain_array_size = m_gridLength*m_gridLength;
    m_vertices = new Vector3[terrain_array_size];
    m_normals = new Vector3[terrain_array_size];
    m_renderNormals = renderNormals;
    m_textureId = 0;
}

/**
  * Converts a grid coordinate (row, column) to an index into a 1-dimensional array.
  * Can be used to index into m_vertices or m_normals.
  * Returns -1 if the grid coordinate entered is not valid.
  */
inline int Terrain::getIndex(const GridIndex &c) const
{
    return getIndex(c.x, c.y);
}

/**
  * Converts a grid coordinate (row, column) to an index into a 1-dimensional array.
  * Can be used to index into m_vertices or m_normals.
  * Returns -1 if the grid coordinate entered is not valid.
  */
inline int Terrain::getIndex(const int row, const int col) const
{
    if (row < 0 || row >= m_gridLength || col < 0 || col >= m_gridLength)
        return -1;

    return row * m_gridLength + col;
}

/**
 * Computes the amount to perturb the height of the vertex currently being processed.
 * Feel free to modify this.
 *
 * @param depth The current recursion depth
 */
double Terrain::getPerturb(const int cur_depth) const
{
    /** roughness = (cur_depth/m_depth).^m_decay*k, where k is between -1 and 1*/
    return m_roughness
           * pow((double)cur_depth / m_depth, m_decay)
           * ((rand() % 200-100) / 100.0);
}

/**
 * Retrieves the position of each neighbor of the given grid coordinate (i.e. all grid
 * coordinates that can be found one unit horizontally, vertically, or diagonally from
 * the specified grid coordinate).
 *
 * @param coordinate The grid coordinate whose neighbors are to be retrieved
 */
QList<Vector3*> Terrain::getSurroundingVertices(const GridIndex &coordinate) const
{
    GridIndex coords[8];
    coords[0] = GridIndex(coordinate.x,     coordinate.y - 1);
    coords[1] = GridIndex(coordinate.x + 1, coordinate.y - 1);
    coords[2] = GridIndex(coordinate.x + 1, coordinate.y);
    coords[3] = GridIndex(coordinate.x + 1, coordinate.y + 1);
    coords[4] = GridIndex(coordinate.x,     coordinate.y + 1);
    coords[5] = GridIndex(coordinate.x - 1, coordinate.y + 1);
    coords[6] = GridIndex(coordinate.x - 1, coordinate.y);
    coords[7] = GridIndex(coordinate.x - 1, coordinate.y - 1);

    int index;
    QList<Vector3*> vecs;

    for (int i = 0; i < 8; i++)
    {
        index = getIndex(coords[i]);
        if (index != -1)
            vecs.push_back(& m_vertices[index]);
    }

    return vecs;
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
void Terrain::subdivideSquare(GridIndex topleft, GridIndex botright, GLint curDepth)
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

 void Terrain::computeNormals()
 {
     // For each vertex in the 2D grid...
     for (int row = 0; row < m_gridLength; row++)
     {
         for (int column = 0; column < m_gridLength; column++)
         {
             const GridIndex gridPosition(row, column);                // 2D coordinate of the vertex on the terrain grid
             const int terrainIndex = getIndex(gridPosition);          // Index into the 1D position and normal arrays
             const Vector3& vertexPosition  = m_vertices[terrainIndex]; // Position of the vertex

             // Get the neighbors of the vertex at (a,b)
             const QList<Vector3*>& neighbors = getSurroundingVertices(gridPosition);
             int numNeighbors = neighbors.size();

             // Compute a list of vectors from vertexPosition to each neighbor in neighbors
             Vector3 *offsets = new Vector3[numNeighbors];
             for (int i = 0; i < numNeighbors; ++i)
             {
                 offsets[i] = Vector3::zero(); // TODO
                 offsets[i].x = neighbors[i]->x - vertexPosition.x;
                 offsets[i].y = neighbors[i]->y - vertexPosition.y;
                 offsets[i].z = neighbors[i]->z - vertexPosition.z;
             }

             // Compute cross products for each neighbor
             Vector3 *normals = new Vector3[numNeighbors];
             for (int i = 0; i < numNeighbors; ++i)
             {
                 normals[i] = Vector3::zero(); // TODO
                 if( i+1 == numNeighbors )
                     normals[i] = offsets[i].cross(offsets[0]);
                 else
                     normals[i] = offsets[i].cross(offsets[i+1]);
             }

             // Average the normals and store the final value in the normal map
             Vector3 sum = Vector3::zero();
             for (int i = 0; i < numNeighbors; ++i)
                 sum += normals[i];
             m_normals[terrainIndex] = sum.getNormalized();

             delete[] offsets;
             delete[] normals;
         }
     }
 }

/**
 * Sets default values for the four corners of the terrain grid and calls subdivideSquare()
 * to begin the terrain generation process. You do not need to modify this function.
 */
void Terrain::populateTerrain()
{
    Vector3 tl(-10, 2, -10);
    Vector3 tr(10, 2, -10);
    Vector3 bl(-10, 2, 10);
    Vector3 br(10, 2, 10);
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

/**
 * Draws a line at each vertex showing the direction of that vertex's normal. You may find
 * this to be a useful tool if you're having trouble getting the lighting to look right.
 * By default, this function is called in paintGL(), but only renders anything if
 * m_renderNormals is true. You do not need to modify this function.
 */
void Terrain::drawNormals() const
{
    if (m_renderNormals)
    {
        glColor3f(1,1,1);

        for (int row = 0; row < m_gridLength; row++)
        {
            for (int column = 0; column < m_gridLength; column++)
            {
                glBegin(GL_LINES);

                Vector3 curVert = m_vertices[getIndex(row, column)];
                Vector3 curNorm = m_normals[getIndex(row, column)];

                glNormal3f(curNorm.x,curNorm.y,curNorm.z);
                glVertex3f(curVert.x, curVert.y, curVert.z);
                glVertex3f(curVert.x +curNorm.x,
                           curVert.y + curNorm.y,
                           curVert.z + curNorm.z);

                glEnd();
            }
        }
    }
}

/**
 * Render the terrain to screen
 */
void Terrain::draw() const
{
    // Clear the color buffer before you draw!

    glDisable( GL_TEXTURE_2D);


    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,m_textureId);
    glPushMatrix();

    int index = 0;
    for( int i = 0; i < m_gridLength-1; i++ )
    {

        glBegin( GL_TRIANGLE_STRIP );

        for( int j = 0; j < m_gridLength; j++ )
        {
            index = i*m_gridLength + j;
            glColor3f(1.0f,1.0f,1.0f);
            glNormal3f( m_normals[index].x, m_normals[index].y, m_normals[index].z );
            glTexCoord2f( ((float)(m_gridLength - i))/m_gridLength, ((float)j)/m_gridLength );
            glVertex3f(  m_vertices[index].x, m_vertices[index].y, m_vertices[index].z );

            index = (i+1)*m_gridLength+j;
            glColor3f(1.0f,1.0f,1.0f);
            glNormal3f( m_normals[index].x, m_normals[index].y, m_normals[index].z );
            glTexCoord2f( ((float)(m_gridLength - i-1))/m_gridLength, ((float)j)/m_gridLength );
            glVertex3f( m_vertices[index].x, m_vertices[index].y, m_vertices[index].z );
        }
        glEnd();
    }
    glDisable( GL_TEXTURE_2D);

    glEnable(GL_BLEND);
    glBlendFunc(GL_ONE,GL_ONE);
    glBegin(GL_QUADS);
    glNormal3f(0.f,1.f,0.f);
    glColor3f(0.3f,0.3f,1.f);
    glVertex3f(-10,4,-10);
    glVertex3f(-10,4,10);
    glVertex3f(10,4,10);
    glVertex3f(10,4,-10);
    glEnd();
    glDisable(GL_BLEND);
    drawNormals();

    glPopMatrix();
    // Force OpenGL to perform all pending operations -- usually a good idea to call this
    // We need that?
    glFlush();
    // We need that?
    // Swap the buffers to show what we have just drawn onto the screen
    //swapBuffers();
}

/**
 ** Load the texture into memory
 **/

void Terrain::loadTextureToTerrain()
{
    m_textureId = loadTexture(TEXTURE_DIR);
}

void Terrain::generate()
{
    loadTextureToTerrain();
    populateTerrain();
    computeNormals();
}
