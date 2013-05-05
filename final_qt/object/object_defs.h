#ifndef OBJECT_DEFS_H
#define OBJECT_DEFS_H


#define DRAG_COEFF 1
#define LIFT_COEFF 0.1
#define DEFAULT_W 1.0

#define MAX_TES_X 50
#define MAX_TES_Y 50
#define MAX_TES_Z 50

#define JITTER_ORIGIN

#define MAX_INIT_U 20
#define MAX_INIT_W 20
#define MAG_U 2.2
#define MAG_W 2.2
#define MIN_DENSITY 300

enum ObjectType
{
    BOX,
    SPHERE
};

#endif // OBJECT_DEFS_H
