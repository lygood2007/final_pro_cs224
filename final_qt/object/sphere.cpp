/** box.cpp
 ** Brief: This is the source file of the class Sphere
 ** Project: large-scale fluids
 ** Date: 04/10/2013
 ** Member: Scott, Hobarts, Yan Li
 **/

#include "sphere.h"
#include "CS123Common.h"

Sphere::Sphere()
{
    m_radius = 1.f;
}

Sphere::Sphere(FluidGPU *fluid, Vector3 position, float radius, float dx, GLuint texID,
               const Colorf color)
    :Object(fluid,position,dx,texID,SPHERE_DENSITY,color)
{
    m_radius = radius;
}

Sphere::~Sphere()
{

}

/**
 * @brief compute the tessellation
 */
void Sphere::computeTessel()
{
    float dx = m_dx;
    float circum = 2*M_PI*m_radius;

    /*****************************
     * m_tessel[0] controls the latitude
     * m_tessel[1] controls the longitude
     *****************************/
    // HACK
    for( int i = 1; i < MAX_TES_X; i++ )
    {
        if( (0.5*circum)/(float)i < dx )
        {
            m_tessell[0] = i;
            break;
        }
    }

    for( int i = 1; i < MAX_TES_Y; i++ )
    {
        if(circum/(float)i < dx )
        {
            m_tessell[1] = i;
            break;
        }
    }
}

/**
 * @brief Sphere::buildTriangleList Build the triangle list
 */
void Sphere::buildTriangleList()
{
 /*   float theta = 0;
    float phi = M_PI/m_tessell[0];
    float unitTheta = 2*M_PI/m_tessell[1];

    Tri tri;
    // push top points to the triangle list
    for( int i = 0; i < m_parameter2; i++ )
    {
        theta = 2*M_PI-unitTheta*i;
        tri.p1.position = vec3<float>( 0, 0.5, 0 );
        tri.p1.normal = vec3<float>( 0, 1, 0 );
        tri.p1.u = (2*M_PI - theta)/(2*M_PI); tri.p1.v = 1;
        tri.p2.position = vec3<float>( m_radius*sin(phi)*cos(theta) ,
                                       m_radius*cos(phi), m_radius*sin(phi)*sin(theta));
        tri.p2.normal = vec3<float>( sin(phi)*cos(theta) ,
                                     cos(phi), sin(phi)*sin(theta));
        tri.p2.u = (2*M_PI - theta)/(2*M_PI); tri.p2.v = 1 - phi/M_PI;

        theta = 2*M_PI-unitTheta*(i+1);
        tri.p3.position = vec3<float>( m_radius*sin(phi)*cos(theta) ,
                                       m_radius*cos(phi), m_radius*sin(phi)*sin(theta));
        tri.p3.normal = vec3<float>( sin(phi)*cos(theta) ,
                                     cos(phi), sin(phi)*sin(theta));
        tri.p3.u = (2*M_PI - theta)/(2*M_PI); tri.p3.v = 1 - phi/M_PI;
        m_triangles.push_back(tri);
    }

    // push bottom points to the triangle list
    phi = M_PI - phi;
    for( int i = 0; i < m_parameter2; i++ )
    {
        theta = unitTheta*i;
        tri.p1.position = vec3<float>( 0,-0.5, 0 );
        tri.p1.normal = vec3<float>( 0, -1, 0 );
        tri.p1.u = (2*M_PI - theta)/(2*M_PI); tri.p1.v = 0;
        tri.p2.position = vec3<float>( m_radius*sin(phi)*cos(theta) ,
                                       m_radius*cos(phi), m_radius*sin(phi)*sin(theta));
        tri.p2.normal = vec3<float>( sin(phi)*cos(theta) ,
                                     cos(phi),sin(phi)*sin(theta));
        tri.p2.u = (2*M_PI - theta)/(2*M_PI); tri.p2.v = 1 - phi/M_PI;
        theta = unitTheta*(i+1);

        tri.p3.position = vec3<float>( m_radius*sin(phi)*cos(theta) ,
                                       m_radius*cos(phi), m_radius*sin(phi)*sin(theta));
        tri.p3.normal = vec3<float>( sin(phi)*cos(theta) ,
                                     cos(phi), sin(phi)*sin(theta));
        tri.p3.u = (2*M_PI - theta)/(2*M_PI); tri.p3.v = 1 - phi/M_PI;
        m_triangles.push_back(tri);
    }

    // push side points to the triangle list
    for( int i = 1; i < m_parameter1-1; i++ )
    {
        for( int j = 0; j < m_parameter2; j++ )
        {
            for( int m = 0; m < 2; m++ )
            {
                if( m == 0)
                {
                    phi = (M_PI/m_parameter1)*i;
                    theta = 2*M_PI-unitTheta*j;
                    tri.p1.position = vec3<float>( m_radius*sin(phi)*cos(theta) ,
                                                   m_radius*cos(phi), m_radius*sin(phi)*sin(theta) );
                    tri.p1.normal = vec3<float>( sin(phi)*cos( theta ),
                                                 cos(phi), sin( phi )*sin( theta ) );
                    tri.p1.u = (2*M_PI - theta)/(2*M_PI); tri.p1.v = 1-phi/M_PI;

                    phi = (M_PI/m_parameter1)*(i+1);
                    theta = 2*M_PI-unitTheta*(j);
                    tri.p2.position = vec3<float>( m_radius*sin(phi)*cos(theta) ,
                                                   m_radius*cos(phi), m_radius*sin(phi)*sin(theta) );
                    tri.p2.normal = vec3<float>( sin(phi)*cos( theta ),
                                                 cos(phi), sin( phi )*sin( theta ) );
                    tri.p2.u = (2*M_PI - theta)/(2*M_PI); tri.p2.v = 1-phi/M_PI;
                    phi = (M_PI/m_parameter1)*(i+1);
                    theta = 2*M_PI-unitTheta*(j+1);

                    tri.p3.position = vec3<float>( m_radius*sin(phi)*cos(theta) ,
                                                   m_radius*cos(phi), m_radius*sin(phi)*sin(theta) );
                    tri.p3.normal = vec3<float>( sin(phi)*cos( theta ),
                                                 cos(phi), sin( phi )*sin( theta ) );
                    tri.p3.u = (2*M_PI - theta)/(2*M_PI); tri.p3.v = 1-phi/M_PI;
                }
                else
                {
                    phi = (M_PI/m_parameter1)*i;
                    theta = 2*M_PI-unitTheta*j;
                    tri.p1.position = vec3<float>( m_radius*sin(phi)*cos(theta) ,
                                                   m_radius*cos(phi), m_radius*sin(phi)*sin(theta) );
                    tri.p1.normal = vec3<float>( sin(phi)*cos( theta ),
                                                 cos(phi), sin( phi )*sin( theta ) );
                    tri.p1.u = (2*M_PI - theta)/(2*M_PI); tri.p1.v = 1-phi/M_PI;
                    phi = (M_PI/m_parameter1)*(i+1);
                    theta = 2*M_PI-unitTheta*(j+1);
                    tri.p2.position = vec3<float>( m_radius*sin(phi)*cos(theta) ,
                                                   m_radius*cos(phi), m_radius*sin(phi)*sin(theta) );
                    tri.p2.normal = vec3<float>( sin(phi)*cos( theta ),
                                                 cos(phi), sin( phi )*sin( theta ) );
                    tri.p2.u = (2*M_PI - theta)/(2*M_PI); tri.p2.v = 1-phi/M_PI;
                    phi = (M_PI/m_parameter1)*i;
                    theta = 2*M_PI-unitTheta*(j+1);

                    tri.p3.position = vec3<float>( m_radius*sin(phi)*cos(theta) ,
                                                   m_radius*cos(phi), m_radius*sin(phi)*sin(theta) );
                    tri.p3.normal = vec3<float>( sin(phi)*cos( theta ),
                                                 cos(phi), sin( phi )*sin( theta ) );
                    tri.p3.u = (2*M_PI - theta)/(2*M_PI); tri.p3.v = 1-phi/M_PI;
                }
                m_triangles.push_back(tri);
            }
        }
    }*/
}

/**
 * @brief computeMass Compute the mass of the object
 */
void Sphere::computeMass()
{
    m_mass = 4.f/3.f*M_PI*m_radius*m_radius*m_radius;
}
