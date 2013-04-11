/** texture_loader.cpp
 ** Brief: The source file of textureloader.h
 ** Project: large-scale fluids
 ** Date: 04/10/2013
 ** Member: Scott, Hobarts, Yan Li
 **/

#include "textureloader.h"
#include <QFile>
#include <assert.h>
/**
 * @brief loadTexture: Load the texture into memory
 * @return: The index for new texture
 */
GLuint loadTexture(const QString &filename)
{
    // Make sure the image file exists
    QFile file(filename);
    assert(file.exists() && "The texture does not exist");

    // Load the file into memory
    QImage image;
    image.load(file.fileName());
  //  image.save("aaa.jpg");
    image = image.mirrored(false, true);
    QImage texture = QGLWidget::convertToGLFormat(image);
    //glEnable(GL_TEXTURE_2D);
    // Generate a new OpenGL texture ID to put our image into
    GLuint id =0;
    glGenTextures(1, &id);

    // Make the texture we just created the new active texture
    glBindTexture(GL_TEXTURE_2D, id);

    // Copy the image data into the OpenGL texture
 //   gluBuild2DMipmaps(GL_TEXTURE_2D, 3, texture.width(), texture.height(), GL_RGBA, GL_UNSIGNED_BYTE, texture.bits());
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA,texture.width(), texture.height(), 0, GL_RGBA, GL_UNSIGNED_BYTE, texture.bits());
    // Set filtering options
    //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

    // Set coordinate wrapping options
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);

    glBindTexture(GL_TEXTURE_2D,0);

    return id;
}
