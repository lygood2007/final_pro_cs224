/** terrain.cpp
 ** Brief: The source file of terrain.h
 ** Project: large-scale fluids
 ** Date: 04/10/2013
 ** Member: Scott, Hobarts, Yan Li
 **/
#include "terrain.h"

extern "C" {
    void glBindBuffer (GLenum target, GLuint buffer);
    void  glGenBuffers (GLsizei n, GLuint *buffers);
    void  glBufferData (GLenum target, GLsizeiptr size, const GLvoid *data, GLenum usage);
}

Terrain::Terrain()
{
    m_depth =  DEFAULT_TERRAIN_DEPTH;
    m_gridLength = (1<<m_depth)+1;
    int terrain_array_size = m_gridLength*m_gridLength;
    m_vertices = new Vector3[terrain_array_size];
    m_normals = new Vector3[terrain_array_size];
    m_uvs = new Vector2[terrain_array_size];
    m_renderNormals = false;
    m_textureId = 0;

    // The actual size of terrain is hard coded
    m_bound = TERRAIN_BOUND;
}

Terrain::~Terrain()
{
    if( m_normals )
        delete []m_normals;
    if( m_vertices )
        delete []m_vertices;
    if( m_uvs )
        delete []m_uvs;
    if( m_indices )
        delete []m_indices;

    if( m_textureId != 0 )
        glDeleteTextures(1,&m_textureId);
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

    /**
     * Old draw, without vbo
     */
/*
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
*/

    glBindBuffer( GL_ARRAY_BUFFER, m_vertexBuffer );
    glVertexPointer(3,GL_FLOAT,0,(char*)NULL);
    glEnableClientState( GL_VERTEX_ARRAY );

    glBindBuffer( GL_ARRAY_BUFFER, m_normalBuffer );
    glNormalPointer(GL_FLOAT,0,(char*)NULL);
    glEnableClientState( GL_NORMAL_ARRAY );

    glBindBuffer( GL_ARRAY_BUFFER, m_texBuffer );
    glTexCoordPointer(2,GL_FLOAT,0,(char*)NULL);
    glEnableClientState( GL_TEXTURE_COORD_ARRAY );

    glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, m_indexBuffer );

    int indexSize = (m_gridLength + 1)*(m_gridLength - 1)*2;
    glDrawElements( GL_TRIANGLE_STRIP, indexSize, GL_UNSIGNED_INT, 0 );

    glDisableClientState( GL_NORMAL_ARRAY );
     glDisableClientState( GL_TEXTURE_COORD_ARRAY );
    glDisableClientState( GL_VERTEX_ARRAY );
    glBindBuffer( GL_ARRAY_BUFFER, 0 );
     glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, 0 );
    glDisable( GL_TEXTURE_2D);
    drawNormals();
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
 * Draw a virtual transparent boundary for the terrain
 */
void Terrain::drawBoundary() const
{
    glDisable(GL_CULL_FACE);

//   glEnable(GL_BLEND);
    glColor3f(0.8f,0.8f,0.8f);
    glBegin( GL_QUADS );
    //glNormal3f( 1.f,0.f,0.f);
    glVertex3f( m_bound, TERRAIN_MAX_HEIGHT*0.6, m_bound );
    glVertex3f( m_bound, TERRAIN_MAX_HEIGHT*0.6, -m_bound );
    glVertex3f( m_bound, TERRAIN_MAX_HEIGHT*1.5, -m_bound );
    glVertex3f( m_bound, TERRAIN_MAX_HEIGHT*1.5, m_bound );

 //   glNormal3f(0.f, 0.f, 1.f );
    glVertex3f( -m_bound, TERRAIN_MAX_HEIGHT*0.6, m_bound );
    glVertex3f( m_bound, TERRAIN_MAX_HEIGHT*0.6, m_bound );
    glVertex3f( m_bound, TERRAIN_MAX_HEIGHT*1.5, m_bound );
    glVertex3f( -m_bound, TERRAIN_MAX_HEIGHT*1.5, m_bound );

   // glNormal3f(-1.f, 0.f, 0.f );
    glVertex3f( -m_bound, TERRAIN_MAX_HEIGHT*0.6, -m_bound );
    glVertex3f( -m_bound, TERRAIN_MAX_HEIGHT*0.6, m_bound );
    glVertex3f( -m_bound, TERRAIN_MAX_HEIGHT*1.5, m_bound );
    glVertex3f( -m_bound, TERRAIN_MAX_HEIGHT*1.5, -m_bound );

   // glNormal3f( 0.f, 0.f, -1.f );
    glVertex3f( m_bound, TERRAIN_MAX_HEIGHT*0.6, -m_bound );
    glVertex3f( -m_bound,TERRAIN_MAX_HEIGHT*0.6, -m_bound );
    glVertex3f( -m_bound, TERRAIN_MAX_HEIGHT*1.5, -m_bound );
    glVertex3f( m_bound, TERRAIN_MAX_HEIGHT*1.5, -m_bound );
    glEnd();
//    glDisable(GL_BLEND);
    glEnable(GL_CULL_FACE);
    //glDisable(GL_BLEND);
}

void Terrain::generate()
{
        loadTextureToTerrain();
        populateTerrain();
        computeNormals();
        generateUV();
        generateIndices();
        generateVBO();
}

/**
 * Generate the UV
 */
void Terrain::generateUV()
{
    int index = 0;
    for( int i = 0; i < m_gridLength; i++ )
    {
        for( int j = 0; j < m_gridLength; j++ )
        {
            index = i*m_gridLength + j;
            m_uvs[index].x = 1.f - i/(float)m_gridLength;
            m_uvs[index].y  = j/(float)m_gridLength;
        }
    }
}

/**
 * Generate indices for triangle strip
 */
void Terrain::generateIndices()
{
    int size = (m_gridLength+1)*(m_gridLength-1)*2;
    m_indices = new GLuint[size];
    int index = 0;
    for( int z = 0; z < m_gridLength-1; z++ )
    {
        int x;
        for( x =  0; x < m_gridLength; x++ )
        {

            m_indices[index] = x + z*m_gridLength;
            index++;
            m_indices[index] = x + (z+1)*m_gridLength;
            index++;
        }
        m_indices[index] = x-1 + (z+1)*m_gridLength;
        index++;
        m_indices[index] = 0 + (z+1)*m_gridLength;
        index++;
    }
}

/**
 * Generate the vertex buffer
 */
void Terrain::generateVBO()
{
    int terrainSize = m_gridLength*m_gridLength;
    glGenBuffers(1,&m_vertexBuffer );
    glBindBuffer( GL_ARRAY_BUFFER, m_vertexBuffer );
    glBufferData( GL_ARRAY_BUFFER, terrainSize*sizeof(Vector3), m_vertices, GL_STATIC_DRAW );

    glGenBuffers( 1, &m_normalBuffer );
    glBindBuffer( GL_ARRAY_BUFFER, m_normalBuffer );
    glBufferData( GL_ARRAY_BUFFER, terrainSize*sizeof(Vector3), m_normals, GL_STATIC_DRAW );

    glGenBuffers( 1, &m_texBuffer );
    glBindBuffer( GL_ARRAY_BUFFER, m_texBuffer );
    glBufferData( GL_ARRAY_BUFFER, terrainSize*sizeof(Vector2), m_uvs, GL_STATIC_DRAW );

    int indexSize = (m_gridLength+1)*(m_gridLength-1)*2;
    glGenBuffers(1,&m_indexBuffer );
    glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, m_indexBuffer );
    glBufferData( GL_ELEMENT_ARRAY_BUFFER, indexSize*sizeof(GLuint), m_indices, GL_STATIC_DRAW );
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
 ** Load the texture into memory
 **/

void Terrain::loadTextureToTerrain()
{
    m_textureId = loadTexture(TEXTURE_DIR);
}

