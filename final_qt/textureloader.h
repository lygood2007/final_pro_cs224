/** texture_loader.h
 ** Brief: The decalaration of the texture loader function for loading image into OpenGL index
 ** Project: large-scale fluids
 ** Date: 04/10/2013
 ** Member: Scott, Hobarts, Yan Li
 **/
#ifndef TEXTURE_LOADER_H
#define TEXTURE_LOADER_H

#include "qgl.h"
#include "GL/glut.h"

/**
 * @brief loadTexture: Load the texture into memory
 * @return: The index for new texture
 */
GLuint loadTexture(const QString &filename);

#endif // TEXTURE_LOADER_H
