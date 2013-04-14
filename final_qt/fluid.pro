#-------------------------------------------------
#
# Project created by QtCreator 2012-03-29T17:53:03
#
#-------------------------------------------------

QT       += core gui opengl

TARGET = fluid
TEMPLATE = app

INCLUDEPATH += common \

DEPENDPATH += common \

SOURCES += main.cpp\
        mainwindow.cpp \
    glwidget.cpp \
    textureloader.cpp \
    terrain.cpp \
    camera.cpp \
    common/CS123Matrix.inl \
    common/CS123Vector.inl \
    fluid.cpp \
    resourceloader.cpp \
    glm.cpp \
    targa.cpp

HEADERS  += mainwindow.h \
    glwidget.h \
    textureloader.h \
    terrain.h \
    camera.h \
    common/vector.h \
    common/CS123Algebra.h \
    common/vector.h \
    common/CS123Common.h \
    fluid.h \
    resourceloader.h \
    glm.h \
    targa.h

OTHER_FILES += cuda/test.cu \
    shaders/refract.vert \
    shaders/refract.frag \
    shaders/reflect.vert \
    shaders/reflect.frag \
    shaders/brightpass.frag \
    shaders/blur.frag
CUDA_SOURCES += cuda/test.cu

# you shouldn't have to change anything under this line

BUILD_DIR = build
OBJECTS_DIR = $$BUILD_DIR/obj
DESTDIR = .

# paths to cuda sdk on filesystem
CUDA_SDK = /contrib/projects/cuda-sdk/C
CUDA_DIR = /contrib/projects/cuda-toolkit/cuda
# cuda architecture, we are 2.1 which is sm_21
CUDA_ARCH = sm_21
# flags for the cuda compiler, in particular, verbosity about what ptx assembly is doing
NVCC_FLAGS = --compiler-options -fno-strict-aliasing -use_fast_math --ptxas-options=-v

#include paths for cuda
INCLUDEPATH += $$CUDA_DIR/include \
                $$CUDA_SDK/common/inc \
                $$CUDA_SDK/../shared/inc

#libs
LIBS += -L$$CUDA_DIR/lib \
        -L$$CUDA_SDK/lib \
        -L$$CUDA_SDK/common/lib/linux \
        -L$$CUDA_SDK/../shared/lib

LIBS += -lcudart -lcutil_i386
CUDA_INC = $$join(INCLUDEPATH, ' -I', '-I', ' ')

cuda.input = CUDA_SOURCES
cuda.output = ${OBJECTS_DIR}${QMAKE_FILE_BASE}_cuda.o

# tell it to use cuda compilers
cuda.commands = $$CUDA_DIR/bin/nvcc -g -G -arch=$$CUDA_ARCH -c $$NVCC_FLAGS $$CUDA_INC $$LIBS ${QMAKE_FILE_NAME} -o ${QMAKE_FILE_OUT}

cuda.dependency_type = TYPE_C
cuda.depend_command = $$CUDA_DIR/bin/nvcc -g -G -M $$CUDA_INC $$NVCCFLAGS ${QMAKE_FILE_NAME}

QMAKE_EXTRA_COMPILERS += cuda

FORMS    += mainwindow.ui

#unix|win32: LIBS += -lGLU -lglut

unix:!macx:!symbian: LIBS += -lGLU -lglut
