/*
 *@NOTE the following assumptions:
 *  all distances in meters
 *  all time steps in seconds
 *  gravity is in the -y plane
 *  2D height field plane is the xz plane
 *
 *Comment history:
 *  04/14 - added comments for the algo structure from the paper - SH
 *
 *
**/



#include "fluid.h"

fluid::fluid()
{
}

//where the magic happens
void fluid::draw()
{

//Main algo loop - once per time step

    //Height field fluid simulation - Sec 2.1
        //discretize the simulation domain where the heights are stored at the cell centers
        //and the velocity components on faces

        //employing time splitting by first solving for self advection of the velocity field
            //advection is the unconditionally stable modified MacCormack method
            //but we fall back onto the semi-Lagrangian method if the first equation result is out of bounds

        //then integrating the height field and velocity field forward in time
            //explicity intergrate the height of the field using equation (7) to guarentee mass preservation
            //@TODO - modify this step to take waterfall discontinuties into account for Sec 2.4.3

            //update face velocities taking the gradient of the water height into account
            //@TODO - this also needs to be modified to account for waterfalls, Sec 2.4.3

        //boundary conditions
            //dealing with reflective and static surfaces by setting face velocity to zero
            //the entire "wet" cell needs to be higher than the terrain level in order to flow, if not treat as stopped

            //for open water scenes implement Perfectly Matched Layers to dampen out the waves when the get near the edge
            //the width of the dampening region is 10 cells

            //@NOTE: for stability, be sure to clamp h_i,j to always be >= zero

            //@NOTE:  for violent wave stabilty clamp u_i +1/2, j and v_i, j+1/2 to alpha dx/dt = 0.5

            //@NOTE: limit the water depth used for the height intergration, beta = 2

        //Overshooting reduction - Sec 2.2
        //to prevent triangle overlap when moving to shallow regions from deeper water
        //to detect edges of the waves and reduce the magnitude - need to fix in x and z planes

    //Solids simulation
    //we currently have the terrain right in glwidget - will either need to move that here to pass a pointer
    //for the next step of the algo
    //Two-way coupling of height field and solids - Sec 2.3

        //recursively divide each triangle in fluid and solids into sub triangles smaller than kdx^2 where k = 1
        //let p = the position at this time step and v = the velocity at this time step of the centroid of the sub triangle
        //A is the area of the subtriangle with  p and v obtained by baycentric interpolation from the original triangle
        //with normal of n from the vector triangle

        //calculate the buoyancy of the triangle
        //calculate the lift of the triangle
        //calculate the drag of the triangle

        //add these forces to the solid the subtriangle belongs to, for a fluid weight the forces to the three vertices of the original tri


        //modify the height and velocity of the fluid due to solids using algo 2


    //Particles generation and simulation - Sec 2.4

    //spray/splash and rform





}
