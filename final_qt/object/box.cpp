/** box.cpp
 ** Brief: This is the source file of the class Box
 ** Project: large-scale fluids
 ** Date: 04/10/2013
 ** Member: Scott, Hobarts, Yan Li
 **/

#include "box.h"

Box::Box():Object()
{
    m_length = 1;
    m_height = 1;
    m_width = 1;

}

Box::Box( FluidGPU* fluid,Vector3 position, float length, float height, float width, float dx, GLuint texID, const Colorf color )
    :Object( fluid, position, dx, texID, BOX_DENSITY, color )
{
    m_length = length;
    m_width = width;
    m_height = height;
}

Box::~Box()
{

}

/**
 * @brief compute the tessellation
 */
// This is different from the paper that the paper will use sub triangles.
// In Box, we can just make the tessellation larger without computing the area
void Box::computeTessel( )
{
    float dx = m_dx;
    for( int i = 1; i < MAX_TES_X; i++ )
    {
        if( m_length/(float)i < dx )
        {
            m_tessell[0] = i;
            break;
        }
    }
    if( m_tessell[0] == -1)
        m_tessell[0] = MAX_TES_X;

    for( int i = 1; i < MAX_TES_Y; i++ )
    {
        if( m_height/(float)i < dx )
        {
            m_tessell[1] = i;
            break;
        }
    }
    if( m_tessell[1] == -1)
        m_tessell[1] = MAX_TES_X;

    for( int i = 1; i < MAX_TES_Z; i++ )
    {
        if( m_width/float(i) < dx)
        {
            m_tessell[2] = i;
            break;
        }
    }
    if( m_tessell[2] == -1)
        m_tessell[2] = MAX_TES_X;
}

/**
 * @brief buildTriangleList build the triangle list
 */
void Box::buildTriangleList()
{
    /**
     * In this section we fix the rotation of vertices.
     * No more rotation will be taken into computation.
     * So this is an simplified version..
     **/
    // UV not supported now
    Vector3 origin;
    float unitL = m_length/m_tessell[0];
    float unitH = m_height/m_tessell[1];
    float unitW = m_width/m_tessell[2];
    int lsize = m_tessell[0];
    int hsize = m_tessell[1];
    int wsize = m_tessell[2];
    float lsizef = (float)lsize;
    float hsizef = (float)hsize;
    float wsizef = (float)wsize;
    // Push front and back faces
    origin = Vector3( -m_length/2, m_height/2,m_width/2 );
    for( int y = 0; y < hsize; y++ )
    {
        for( int x = 0; x < lsize; x++ )
        {
            Vector3 v01;
            Vector3 v02;
            // front face
            Tri t1;
            t1.verts[0] = origin + Vector3(x*unitL, -y*unitH,0);
            t1.verts[1] = origin + Vector3(x*unitL, -(y+1)*unitH,0);
            t1.verts[2] = origin + Vector3((x+1)*unitL, -(y+1)*unitH,0);
            t1.uvs[0] = Vector2(x/lsizef,1-y/hsizef);
            t1.uvs[1] = Vector2(x/lsizef,1-(y+1)/hsizef);
            t1.uvs[2] = Vector2((x+1)/lsizef,1-(y+1)/hsizef);

            v01 = t1.verts[1] - t1.verts[0];
            v02 = t1.verts[2] - t1.verts[1];
            Vector3 cro = v01.cross(v02);
            t1.area = cro.length()/2.f;
            t1.norms[0] = Vector3(0.f,0.f,1.f); t1.norms[1] = Vector3(0.f,0.f,1.f), t1.norms[2] = Vector3(0.f,0.f,1.f);
            t1.avgNorm = t1.norms[0];
            t1.avgPos = (t1.verts[0] + t1.verts[1] + t1.verts[2])/3.f+ m_position;
            m_tris.push_back( t1 );
            // back face
            Tri t2;
            t2.verts[0] = t1.verts[2]; t2.verts[0].z = -t2.verts[0].z;
            t2.verts[1] = t1.verts[1]; t2.verts[1].z = -t2.verts[1].z;
            t2.verts[2]= t1.verts[0];t2.verts[2].z = -t2.verts[2].z;
            t2.uvs[0] = t1.uvs[2];t2.uvs[0].x = 1-t2.uvs[0].x;
            t2.uvs[1] = t1.uvs[1];t2.uvs[1].x = 1-t2.uvs[1].x;
            t2.uvs[2] = t1.uvs[0];t2.uvs[2].x = 1-t2.uvs[2].x;
            t2.norms[0] = -t1.norms[2];
            t2.norms[1] = -t1.norms[1];
            t2.norms[2] = -t1.norms[0];
            t2.area = t1.area;
            t2.avgNorm = t1.norms[0];
            t2.avgPos = (t2.verts[0] + t2.verts[1] + t2.verts[2])*0.33f+ m_position;
            m_tris.push_back( t2 );
            //front face
            Tri t3;
            t3.verts[0] = origin + Vector3(x*unitL, -y*unitH,0);
            t3.verts[1] = origin + Vector3((x+1)*unitL, -(y+1)*unitH,0);
            t3.verts[2] = origin + Vector3((x+1)*unitL, -(y)*unitH,0);
            t3.uvs[0] = Vector2(x/lsizef,1-y/hsizef);
            t3.uvs[1] = Vector2((x+1)/lsizef,1-(y+1)/hsizef);
            t3.uvs[2] = Vector2((x+1)/lsizef,1-(y)/hsizef);
            v01 = t3.verts[1] - t3.verts[0];
            v02 = t3.verts[2] - t3.verts[0];
            t3.area = fabs(v01.cross(v02).length())/2.f;
            t3.norms[0] = Vector3(0.f,0.f,1.f); t3.norms[1] = Vector3(0.f,0.f,1.f); t3.norms[2] = Vector3(0.f,0.f,1.f);
            t3.avgNorm = t3.norms[0];
            t3.avgPos = (t3.verts[0] + t3.verts[1] + t3.verts[2])/3.f+ m_position;
            m_tris.push_back( t3 );

            // back face
            Tri t4;
            t4.verts[0] = t3.verts[2]; t4.verts[0].z = -t4.verts[0].z;
            t4.verts[1] = t3.verts[1]; t4.verts[1].z = -t4.verts[1].z;
            t4.verts[2] = t3.verts[0]; t4.verts[2].z = -t4.verts[2].z;
            t4.uvs[0] = t3.uvs[2],t4.uvs[0].x = 1.f-t4.uvs[0].x;
            t4.uvs[1] = t3.uvs[1],t4.uvs[1].x = 1.f-t4.uvs[1].x;
            t4.uvs[2] = t3.uvs[0]; t4.uvs[2].x = 1.f-t4.uvs[2].x;

            t4.norms[0] = -t3.norms[2];
            t4.norms[1] = -t3.norms[1];
            t4.norms[2] = -t3.norms[0];
            t4.area = t3.area;
            t4.avgNorm = t4.norms[0];
            t4.avgPos = (t4.verts[0] + t4.verts[1] + t4.verts[2])*0.33f+ m_position;
            m_tris.push_back( t4 );
        }
    }

    // Push left and right face
    origin = Vector3( -m_length/2, m_height/2, -m_width/2 );
    for( int y = 0; y < hsize; y++ )
    {
        for( int z = 0; z < wsize; z++ )
        {
            Vector3 v01;
            Vector3 v02;
            // left face
            Tri t1;
            t1.verts[0] = origin + Vector3(0, -y*unitH,z*unitW);
            t1.verts[1] = origin + Vector3(0, -(y+1)*unitH,z*unitW);
            t1.verts[2] = origin + Vector3(0, -(y+1)*unitH,(z+1)*unitW);
            t1.uvs[0] = Vector2(z/wsizef,1-y/hsizef);
            t1.uvs[1] = Vector2(z/wsizef,1-(y+1)/hsizef);
            t1.uvs[2] = Vector2((z+1)/wsizef,1-(y+1)/hsizef);
            t1.norms[0] = Vector3(-1.f,0.f,0.f); t1.norms[1] = Vector3(-1.f,0.f,0.f), t1.norms[2] = Vector3(-1.f,0.f,0.f);
            v01 = t1.verts[1] - t1.verts[0];
            v02 = t1.verts[2] - t1.verts[1];
            Vector3 cro = v01.cross(v02);
            t1.area = cro.length()/2.f;
            t1.avgNorm = t1.norms[0];
             t1.avgPos = (t1.verts[0] + t1.verts[1] + t1.verts[2])*0.33f+ m_position;
            m_tris.push_back( t1 );


            // right face
            Tri t2;
            t2.verts[0] = t1.verts[2]; t2.verts[0].x = -t2.verts[0].x;
            t2.verts[1] = t1.verts[1]; t2.verts[1].x = -t2.verts[1].x;
            t2.verts[2]= t1.verts[0];t2.verts[2].x = -t2.verts[2].x;
            t2.uvs[0] = t1.uvs[2];t2.uvs[0].x = 1-t2.uvs[0].x;
            t2.uvs[1] = t1.uvs[1];t2.uvs[1].x = 1-t2.uvs[1].x;
           t2.uvs[2] = t1.uvs[0];t2.uvs[2].x = 1-t2.uvs[2].x;

            t2.norms[0] = -t1.norms[2];
            t2.norms[1] = -t1.norms[1];
            t2.norms[2] = -t1.norms[0];
            t2.area = t1.area;
            t2.avgNorm = t2.norms[0];
            t2.avgPos = (t2.verts[0] + t2.verts[1] + t2.verts[2])*0.33f+ m_position;
            m_tris.push_back( t2 );


            Tri t3;
            t3.verts[0] = origin + Vector3(0, -y*unitH,z*unitW);
            t3.verts[1] = origin + Vector3(0, -(y+1)*unitH,(z+1)*unitW);
            t3.verts[2] = origin + Vector3(0, -(y)*unitH,(z+1)*unitW);
            t3.uvs[0] = Vector2(z/wsizef,1-y/hsizef);
            t3.uvs[1] = Vector2((z+1)/wsizef,1-(y+1)/hsizef);
            t3.uvs[2] = Vector2((z+1)/wsizef,1-(y)/hsizef);
            t3.norms[0] = Vector3(-1.f,0.f,0.f); t3.norms[1] = Vector3(-1.f,0.f,0.f); t3.norms[2] = Vector3(-1.f,0.f,0.f);
            v01 = t3.verts[1] - t3.verts[0];
            v02 = t3.verts[2] - t3.verts[0];
            t3.area = fabs(v01.cross(v02).length())/2.f;
            t3.avgNorm = t3.norms[0];
            t3.avgPos = (t3.verts[0] + t3.verts[1] + t3.verts[2])*0.33f+ m_position;
            m_tris.push_back( t3 );

            // left face
            Tri t4;
            t4.verts[0] = t3.verts[2]; t4.verts[0].x = -t4.verts[0].x;
            t4.verts[1] = t3.verts[1]; t4.verts[1].x = -t4.verts[1].x;
            t4.verts[2] = t3.verts[0]; t4.verts[2].x = -t4.verts[2].x;
            t4.uvs[0] = t3.uvs[2];t4.uvs[0].x = 1.f-t4.uvs[0].x;
            t4.uvs[1] = t3.uvs[1];t4.uvs[1].x = 1.f-t4.uvs[1].x;
            t4.uvs[2] = t3.uvs[0]; t4.uvs[2].x = 1.f-t4.uvs[2].x;
            t4.norms[0] = -t3.norms[2];
            t4.norms[1] = -t3.norms[1];
            t4.norms[2] = -t3.norms[0];
            t4.area = t3.area;
            t4.avgNorm = t4.norms[0];
            t4.avgPos = (t4.verts[0] + t4.verts[1] + t4.verts[2])*0.33f+ m_position;
            m_tris.push_back( t4 );
        }
    }

    // Push top and bottom faces
     origin = Vector3( -m_length/2, m_height/2, -m_width/2 );
     for( int z = 0; z < wsize; z++ )
     {
         for( int x = 0; x < lsize; x++ )
         {
             Vector3 v01;
             Vector3 v02;
             // top face
             Tri t1;
             t1.verts[0] = origin + Vector3(x*unitL, 0,z*unitW);
             t1.verts[1] = origin + Vector3((x+1)*unitL,0,(z+1)*unitW);
             t1.verts[2] = origin + Vector3((x+1)*unitL,0,z*unitW);
             t1.uvs[0] = Vector2(x/lsizef,z/wsizef);
             t1.uvs[1] = Vector2((x+1)/lsizef,(z+1)/wsizef);
             t1.uvs[2] = Vector2((x+1)/lsizef,(z)/wsizef);
             t1.norms[0] = Vector3(0.f,1.f,0.f); t1.norms[1] = Vector3(0.f,1.f,0.f), t1.norms[2] = Vector3(0.f,1.f,0.f);
             v01 = t1.verts[1] - t1.verts[0];
             v02 = t1.verts[2] - t1.verts[0];
             Vector3 cro = v01.cross(v02);
             t1.area = cro.length()/2.f;
             t1.avgNorm = t1.norms[0];
              t1.avgPos = (t1.verts[0] + t1.verts[1] + t1.verts[2])/3.f+ m_position;
             m_tris.push_back( t1 );
             //areaCheck1 += t1.area;
             // bottom face
             Tri t2;
             t2.verts[0] = t1.verts[2]; t2.verts[0].y = -t2.verts[0].y;
             t2.verts[1] = t1.verts[1]; t2.verts[1].y = -t2.verts[1].y;
             t2.verts[2]= t1.verts[0];t2.verts[2].y = -t2.verts[2].y;
             t2.uvs[0] = t1.uvs[2];t2.uvs[0].x = 1-t2.uvs[0].x;
             t2.uvs[1] = t1.uvs[1];t2.uvs[1].x = 1-t2.uvs[1].x;
             t2.uvs[2] = t1.uvs[0];t2.uvs[2].x = 1-t2.uvs[2].x;
             t2.norms[0] = -t1.norms[2];
             t2.norms[1] = -t1.norms[1];
             t2.norms[2] = -t1.norms[0];
              t2.area = t1.area;
              t2.avgNorm = t2.norms[0];
              t2.avgPos = (t2.verts[0] + t2.verts[1] + t2.verts[2])/3.f+ m_position;
             m_tris.push_back( t2 );
             // top face
             Tri t3;
             t3.verts[0] = origin + Vector3(x*unitL,0,z*unitW);
             t3.verts[1] = origin + Vector3(x*unitL,0,(z+1)*unitW);
             t3.verts[2] = origin + Vector3((x+1)*unitL,0,(z+1)*unitW);
             t3.uvs[0] = Vector2(x/lsizef,z/wsizef);
             t3.uvs[1] = Vector2((x)/lsizef,(z+1)/wsizef);
             t3.uvs[2] = Vector2((x+1)/lsizef,(z+1)/wsizef);
             t3.norms[0] = Vector3(0.f,1.f,0.f); t3.norms[1] = Vector3(0.f,1.f,0.f); t3.norms[2] = Vector3(0.f,1.f,0.f);
             v01 = t3.verts[1] - t3.verts[0];
             v02 = t3.verts[2] - t3.verts[0];
             t3.area = v01.cross(v02).length()/2.f;
              t3.avgNorm = t3.norms[0];
              t3.avgPos = (t3.verts[0] + t3.verts[1] + t3.verts[2])/3.f+ m_position;
             m_tris.push_back( t3 );

             // bottom face
             Tri t4;
             t4.verts[0] = t3.verts[2]; t4.verts[0].y = -t4.verts[0].y;
             t4.verts[1] = t3.verts[1]; t4.verts[1].y = -t4.verts[1].y;
             t4.verts[2] = t3.verts[0]; t4.verts[2].y = -t4.verts[2].y;
             t4.uvs[0] = t3.uvs[2];t4.uvs[0].x = 1.f-t4.uvs[0].x;
             t4.uvs[1] = t3.uvs[1];t4.uvs[1].x = 1.f-t4.uvs[1].x;
             t4.uvs[2] = t3.uvs[0]; t4.uvs[2].x = 1.f-t4.uvs[2].x;
             t4.norms[0] = -t3.norms[2];
             t4.norms[1] = -t3.norms[1];
             t4.norms[2] = -t3.norms[0];
             t4.area = t3.area;
              t4.avgNorm = t4.norms[0];
              t4.avgPos = (t4.verts[0] + t4.verts[1] + t4.verts[2])/3.f + m_position;
             m_tris.push_back( t4 );
         }
     }

     updatePosAndNorm();

     // Rotate
     /*for( int i = 0; i < m_tris.size(); i++ )
     {
         m_tris[i].verts[0] = getRotVec(m_tris[i].verts[0]);
         m_tris[i].verts[1] = getRotVec(m_tris[i].verts[1]);
         m_tris[i].verts[2] = getRotVec(m_tris[i].verts[2]);

         m_tris[i].norms[0] = getRotVec(m_tris[i].norms[0]);
         m_tris[i].norms[1] = getRotVec(m_tris[i].norms[1]);
         m_tris[i].norms[2] = getRotVec(m_tris[i].norms[2]);
         m_tris[i].avgNorm = 0.33*(m_tris[i].norms[0] + m_tris[i].norms[1] + m_tris[i].norms[2]);
         m_tris[i].avgPos = (m_tris[i].verts[0] + m_tris[i].verts[1] + m_tris[i].verts[2])*0.33f + m_position;
     }*/
    //Check the top and bottom

     /**
      * Test if the sum of the triangle areas are correct
      */

     /*float areaTest = m_length*m_height*2 + m_width*m_height*2 + m_width*m_length*2;
     float sumTriArea = 0;
     for( int i = 0; i < m_tris.size(); i++ )
     {
         sumTriArea += m_tris[i].area;
     }*/
     /*if( fabs(sumTriArea > areaTest) > EPSILON )
     {
         assert(0);
     }*/
}

/**
 * @brief computeMass Compute the mass of the object
 */
void Box::computeMass()
{
    m_mass = m_density*m_length*m_height*m_width;
    m_massInv = 1/m_mass;
}


/**
 * @brief computeBoundingRadius Compute the bounding radius
 */
void Box::computeBoundingRadius()
{
    m_boundingRadius = sqrt(m_length*m_length + m_width*m_width + m_height*m_height);
}
