/** object.cpp
 ** Brief: This is the source file of the class Object
 ** Project: large-scale fluids
 ** Date: 04/10/2013
 ** Member: Scott, Hobarts, Yan Li
 **/

#include "object.h"
#include <stdio.h>
#include <iostream>
#include "fluidGPU.h"

Object::Object()
{
    m_domainSize = TERRAIN_BOUND*2;
    m_position = Vector3(0.f,0.f,0.f);
    //m_fluiddx = 1.f;
    m_dx = 1.f;
    m_dxInv = 1.f/m_dx;
    m_color = Colorf(1.f,1.f,1.f,1.f);
    m_renderNormals = false;
    m_density = WATER_DENSITY;
    m_upwards = false;
    m_decay = false;
    m_dragCoeff = DRAG_COEFF;
    m_liftCoeff = LIFT_COEFF;
    m_w = DEFAULT_W;
    m_texID = 0;
    c1 = -WATER_DENSITY*m_dragCoeff*0.5;
    c2 = -WATER_DENSITY*m_liftCoeff*0.5;
    c3 = WATER_DENSITY*GRAVITY;
    memset( m_tessell, -1, 3*sizeof(int ) );
    computeOrigin();
}

Object::Object(FluidGPU *fluid, Vector3 position, float dx, GLuint texID, float density, const Colorf color)
{
    m_domainSize = fluid->getDomainSize();
    m_position = position;
    m_dx = dx;
    m_dxInv = 1.f/dx;
    m_texID = texID;
    m_color = color;
    m_renderNormals = false;
    m_dragCoeff = DRAG_COEFF;
    m_liftCoeff = LIFT_COEFF;
    m_w = DEFAULT_W;
    m_upwards = false;
    m_decay = false;
    m_density = density;
    c1 = -WATER_DENSITY*m_dragCoeff*0.5;
    c2 = -WATER_DENSITY*m_liftCoeff*0.5;
    c3 = WATER_DENSITY*GRAVITY;
    memset( m_tessell, -1, 3*sizeof(int ) );
    computeOrigin();
    initFluidInfo( fluid );
}

Object::~Object()
{

}


void Object::setFluid( FluidGPU* fluid )
{
    initFluidInfo( fluid );
}

void Object::draw()
{
    glEnable(GL_DEPTH_TEST);
    glMatrixMode( GL_MODELVIEW );
     glPushMatrix();
    glTranslatef( m_position.x, m_position.y, m_position.z );
    glRotatef( m_angle[4],0,0,1);
    glRotatef( m_angle[5],0,1,0);
    glRotatef( m_angle[6],1,0,0 );
    glEnable(GL_TEXTURE_2D );
    glBindTexture( GL_TEXTURE_2D, m_texID );
     glBegin(GL_TRIANGLES);
     glColor3f( 1.0f, 1.0f, 1.0f );
    for( int i = 0; i < m_tris.size(); i++ )
    {
        glNormal3fv( m_tris[i].norms[0].xyz );
        glTexCoord2f( m_tris[i].uvs[0].x,m_tris[i].uvs[0].y);
        glVertex3fv( m_tris[i].verts[0].xyz );
        glNormal3fv( m_tris[i].norms[1].xyz );
        glTexCoord2f( m_tris[i].uvs[1].x,m_tris[i].uvs[1].y);
        glVertex3fv( m_tris[i].verts[1].xyz );
        glNormal3fv( m_tris[i].norms[2].xyz );
        glTexCoord2f( m_tris[i].uvs[2].x,m_tris[i].uvs[2].y);
        glVertex3fv( m_tris[i].verts[2].xyz );
    }
    glEnd();
    glPopMatrix();
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_TEXTURE_2D );
    if( m_renderNormals )
        drawNormal();
}

void Object::update( float dt, FluidGPU* fluid )
{
    /**
     * Compute the buoyance
     */

  //  int buf;
    //printf("%f\n",fluid->getFieldArray(VEL_U,buf)[750]);
    //Vector3 decay = Vector3(1.f,1.f,1.f);
    Vector3 yUnit = Vector3(0.f,1.f,0.f );
    Vector3 buoyTotal = Vector3(0.f,0.f,0.f);
    Vector3 buoyAccTotal = Vector3(0.f,0.f,0.f);
    Vector3 liftTotal = Vector3(0.f,0.f,0.f);
    Vector3 liftAccTotal = Vector3(0.f,0.f,0.f);
    Vector3 dragTotal = Vector3(0.f,0.f,0.f);
    Vector3 dragAccTotal = Vector3(0.f,0.f,0.f);
    // Loop through all the triangles
    for( int i = 0; i < m_tris.size(); i++ )
    {
        Vector3 pos = m_tris[i].avgPos;
        float h;
        Vector3 fluidVelocity;
        getInterpVelocityAndHeight( pos, fluidVelocity,h);

        if( pos.y > h )
            continue;

        Vector3 relVelocity = m_velocity - fluidVelocity;
        float cos = m_tris[i].avgNorm.dot(relVelocity.getNormalized());
        float Aef;

        if( cos <= 0.001 )
            Aef = 0;
        else
        {
         // No w currently
            Aef =m_tris[i].area*cos;
        }
        float l = relVelocity.length();
        Vector3 drag = c1*l*relVelocity*Aef;
        Vector3 tmp1 = m_tris[i].avgNorm.cross( relVelocity);
        Vector3 lift;
        if( tmp1.length() <= 0.00001)
        {
            lift = Vector3::zero();
        }
        else
        {
        Vector3 tmp2 = tmp1.getNormalized();
         lift = c2*l
                *( relVelocity.cross( tmp2 ) )*Aef;
        }
        float buoy = c3*m_tris[i].area*(h - pos.y )*m_tris[i].avgNorm.dot(yUnit);
        buoyTotal += Vector3( 0.f, buoy, 0.f );
        dragTotal += drag;
        liftTotal += lift;
    }
    buoyAccTotal = buoyTotal*m_massInv;
    dragAccTotal = dragTotal*m_massInv;
    liftAccTotal = liftTotal*m_massInv;
   //printf("Drag a: %f, %f, %f\n",dragAccTotal.x, dragAccTotal.y, dragAccTotal.z );
  /* printf("Lift a: %f, %f, %f\n",liftAccTotal.x, liftAccTotal.y, liftAccTotal.z );
    if( m_upwards &&fabs(m_lastBuoAcc.y) > fabs(GRAVITY) && fabs(buoyAccTotal.y) <= fabs(GRAVITY) )
    {
        // In this case, it means the object pass through the point where Buoyancy
        // Equals GRAVITY
        //m_decay = false;
    }
    if( m_decay )
    {
        // Decay faster for y direction, slower for x,z direction
        //decay.x = 0.99f;
        decay.y = m_density/WATER_DENSITY;
        //decay.z = 0.99f;
    }
    m_velocity.x = m_velocity.x*decay.x;m_velocity.y = m_velocity.y*decay.y;m_velocity.z = m_velocity.z*decay.z;
*/
    m_velocity += 2*(buoyAccTotal*dt + Vector3(0.f,GRAVITY,0.f)*dt +dragAccTotal*dt +liftAccTotal*dt);
    //m_velocity.y += 2*liftAccTotal.y*dt;
    // Backup the buoyancy in this frame
    /*m_lastBuoAcc = buoyAccTotal;
    if( m_velocity.y < 0 )
    {
        m_decay = false;
        m_upwards = false;
    }
    else
        m_upwards = true;
*/
    // Check if the next position is below terrain, if it is, we set the velocity to zero
    // Actually if one of the box's triangle is below terrain
    bool hitBottom = false;
    int buff;
    float* terrainHeight = fluid->getFieldArray( TERRAINH, buff );
     float gridSize =fluid->getGridSize();
     Vector3 next = m_position + dt*(m_velocity);
    for( int i = 0; i < m_tris.size(); i++ )
    {
       // for( int j = 0; j < 3; j++ )
     //   {
            Vector3 tmpPos =  next + m_tris[i].verts[0];
            float th = bilinearInterpReal( terrainHeight, tmpPos,
                                           m_origX,m_origZ,m_dx,
                                           gridSize, gridSize
                                           );
            if( tmpPos.y < th )
            {
                hitBottom = true;
                break;
            }
     //   }
    }
    if( hitBottom )
    {
        // Hack
        const float hitDecay = 0.3;
        m_velocity.x = m_velocity.x*hitDecay;
        m_velocity.y = 0;
        m_velocity.z = m_velocity.z*hitDecay;

    }

    updatePosWithBoundaryCheck(dt);

    updatePosAndNorm();
}

/**
 * @brief Compute mass, tesselation,etc
 */
void Object::initPhysics()
{

#ifdef JITTER_ORIGIN
        m_velocity.x = -MAX_INIT_U + rand()%(2*MAX_INIT_U);
        m_velocity.y = 0;
        m_velocity.z = -MAX_INIT_W + rand()%(2*MAX_INIT_W);
        generateRotation();
#else
    m_velocity = Vector3(0.f,0.f,0.f);
    m_rotMat[0] = Matrix4x4::identity();
    m_rotMat[1] = Matrix4x4::identity();
    m_rotMat[2] = Matrix4x4::identity();
#endif
    computeTessel();
    computeMass();
    buildTriangleList();
}

/**
 * @brief initFluidInfo Initialize the data from fluid
 */
void Object::initFluidInfo( FluidGPU* fluid )
{
    int buff;
    m_curH  = fluid->getFieldArray( HEIGHT, buff );
    m_curU  = fluid->getFieldArray( VEL_U, buff );
    m_curW  = fluid->getFieldArray( VEL_W, buff );

    m_gridSize = fluid->getGridSize();
    m_uwidth = m_gridSize+1;
    m_wheight = m_gridSize+1;
}

/**
 * @brief drawNormal Draw the normals of this box
 */
void Object::drawNormal() const
{
    if( m_renderNormals )
    {
        const float magn = 0.4;
        glMatrixMode(GL_MODELVIEW);

        glPushMatrix();
         glTranslatef( m_position.x, m_position.y, m_position.z );
         glRotatef( m_angle[4],0,0,1);
         glRotatef( m_angle[5],0,1,0);
         glRotatef( m_angle[6],1,0,0 );
        for( int i = 0; i < m_tris.size(); i++ )
        {
            for( int j = 0; j < 3; j++ )
            {
                glBegin(GL_LINES);
                glNormal3fv( m_tris[i].norms[j].xyz );
                glVertex3fv(m_tris[i].verts[j].xyz );
                glVertex3f( m_tris[i].verts[j].x + magn*m_tris[i].norms[j].x,
                            m_tris[i].verts[j].y + magn*m_tris[i].norms[j].y,
                            m_tris[i].verts[j].z + magn*m_tris[i].norms[j].z
                            );
                glEnd();
            }
        }
        glPopMatrix();
    }
}

/**
     *  Update the triangles' postion and norms (average position and normal)
     *  The are used for buoyance computation
     */
void Object::updatePosAndNorm()
{
    // Currently, we don't update the normals because no rotation is included
    for( int i = 0; i < m_tris.size(); i++ )
    {
        m_tris[i].avgPos = (m_tris[i].verts[0] + m_tris[i].verts[1] + m_tris[i].verts[2])*0.33f + m_position;
    }
}

/**
 * @brief updatePosWithBoundaryCheck update the next position, avoiding moving out of boundary
 */
void Object::updatePosWithBoundaryCheck( float dt )
{
    Vector3 nextPos = m_position+m_velocity*dt;
    float halfDomain = m_domainSize/2.f;
    // To adjust the size of the boundary
    const float ratio = 0.9;
    const float decayRatio = 0.7;
    if( nextPos.z > ratio*halfDomain || nextPos.z < -ratio*halfDomain )
    {
        m_velocity.z = -m_velocity.z*decayRatio;
    }
    if( nextPos.x > ratio*halfDomain || nextPos.x < -ratio*halfDomain )
    {
        m_velocity.x = -m_velocity.x*decayRatio;
    }
    m_position += m_velocity*dt;
}

/**
 * @brief bilinearInterpReal Get
 * @brief This function will do bilinear intepolation based on x,z value of the vector
 * @param fluid The fluid
 * @param pos The position
 * @param fv the fluid velocity
 * @param the interpolated height
 */
 void Object::getInterpVelocityAndHeight(/* const FluidGPU* fluid,*/ Vector3 pos,
                                    Vector3& fv, float&h )
 {
  /*   int buffLength;
     // Get the field
     float* hf = fluid->getFieldArray( HEIGHT, buffLength );
     float* uf  =fluid->getFieldArray( VEL_U, buffLength );
     float* wf  =fluid->getFieldArray( VEL_W, buffLength );
*/
  //   float a = wf[2];

     /**
      * Get the size. width and height are the size of the height field
      */
   /*  int width = fluid->getGridSize();
     int height = width;
     // width of velocity u
     int uwidth = fluid->getGridSize()+1;
     // height of velocity v
     int wheight = fluid->getGridSize()+1;
*/
     float x = pos.x - m_origX;
     float z = pos.z - m_origZ;
//     float dx = m_fluiddx

     if( x < 0 )
         x = 0.f;
     if( z < 0 )
         z = 0.f;
     float xx = (m_gridSize - 1)*m_dx;
     if( x > xx )
         x = xx;
     float zz = (m_gridSize - 1)*m_dx;
     if( z > zz )
         z = zz;

     // Get the index
     const int X = (int)(x*m_dxInv);
     const int Y = (int)(z*m_dxInv);
     const float s1 = x - X*m_dx;
     const float s0 = 1.f - s1;
     const float t1 = z - Y*m_dx;
     const float t0 = 1.f-t1;
     float e1, e2, e3,e4;
     e1 = e2 = e3 = e4 = 0;
     e1 = m_curH[Y*m_gridSize+X];
     if( Y+1 <= m_gridSize- 1 )
     {
         e2 = m_curH[(Y+1)*m_gridSize+ X];
     }
     if( X +1 <= m_gridSize-1 )
     {
         e3 = m_curH[Y*m_gridSize + X+1];
     }
     if( Y+1 <= m_gridSize - 1 && X + 1 <= m_gridSize - 1)
     {
         e4 = m_curH[(Y+1)*m_gridSize + X+1];
     }

     h= s0*(t0*e1 + t1*e2 )+
             s1*(t0*e3  + t1*e4 );


     // Deal with the fluid's velocity
     fv = Vector3(0.f,0.f,0.f);
     if( X+1 <= m_uwidth -1 )
     {
         // The u component,average the two
         fv.x = 0.5*(m_curU[Y*m_uwidth + X+1 ] + m_curU[Y*m_uwidth +X]);
     }
     else
     {
         fv.x = m_curU[Y*m_uwidth +X];
     }

     if( Y+1 <= m_wheight - 1)
     {
         fv.z = 0.5*(m_curW[Y*m_gridSize + X] + m_curW[(Y+1)*m_gridSize+ X]);
     }
     else
         fv.z = m_curW[Y*m_gridSize + X ];

     float c1 = 0.f;
     float c2 = 0.f;
     if( X-1 >= 0 && X+1 <= m_gridSize - 1)
     {
         c1 = (m_curH[Y*m_gridSize+ X+1] - m_curH[Y*m_gridSize+ X - 1])*0.5*m_dxInv;
     }
     if( Y-1 >= 0 && Y + 1 <= m_gridSize -1 )
     {
         c2 = (m_curH[(Y+1)*m_gridSize + X ] - m_curH[(Y-1)*m_gridSize+X])*0.5*m_dxInv;
     }
    fv.y  = c1*fv.x + c2*fv.z;
     fv.x *= MAG_U;
    fv.z *= MAG_W;
 }

/**
 * @brief computeOrigin compute the origin position of x and z
 */
void Object::computeOrigin()
{
    m_origX = -m_domainSize/2.f;
    m_origZ = -m_domainSize/2.f;
}

void Object::generateRotation()
{
    float angle = 2*M_PI*(rand()/(float)RAND_MAX);
    m_angle[0] = angle;

    angle = 2*M_PI*(rand()/(float)RAND_MAX);
    m_angle[1] = angle;
    angle = 2*M_PI*(rand()/(float)RAND_MAX);
    m_angle[2] = angle;
    m_angle[3] = m_angle[2]/(2*M_PI)*360;
    m_angle[4] = m_angle[1]/(2*M_PI)*360;
    m_angle[5] = m_angle[0]/(2*M_PI)*360;


    m_rotMat[0] = Matrix4x4::identity();
    m_rotMat[1] = Matrix4x4::identity();
    m_rotMat[2] = Matrix4x4::identity();

    // Rotation along x-axis
    m_rotMat[0].data[5] = cos(m_angle[0]);
    m_rotMat[0].data[6] = -sin(m_angle[0]);
    m_rotMat[0].data[9] = sin(m_angle[0]);
    m_rotMat[0].data[10] = cos(m_angle[0]);

    // Rotation along y-axis
    m_rotMat[1].data[0] = cos(m_angle[1]);
    m_rotMat[1].data[2] = sin(m_angle[1]);
    m_rotMat[1].data[8] = -sin(m_angle[1]);
    m_rotMat[1].data[10] = cos(m_angle[1]);

    // Rotation along z-axis
    m_rotMat[2].data[0] = cos(m_angle[2]);
    m_rotMat[2].data[1] = -sin(m_angle[2]);
    m_rotMat[2].data[4] = sin(m_angle[2]);
    m_rotMat[2].data[5] = cos(m_angle[2]);
}

/**
 * Get the rotated position using m_angle which defines the angle of ration along x,y,z     *
 */
Vector3 Object::getRotVec( Vector3 vec )
{
    Vector4 tmp;
    tmp.x = vec.x;
    tmp.y = vec.y;
    tmp.z = vec.z;
    tmp.w  =1;
    Vector4 out;
    m_rotMat[0].mulVec4( tmp, out );
    m_rotMat[1].mulVec4( out, tmp );
    m_rotMat[2].mulVec4( tmp, out );
    Vector3 out3;
    out3.x  = out.x;
    out3.y = out.y;
    out3.z = out.z;
    return out3;
}
