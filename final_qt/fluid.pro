#-------------------------------------------------
#
# Project created by QtCreator 2012-03-29T17:53:03
#
#-------------------------------------------------

QT       += core gui opengl

TARGET = fluid
TEMPLATE = app

INCLUDEPATH += common \
    support \
    terrain \
    water \
    camera \
    loader \
    object

DEPENDPATH += common \
    support \
    terrain \
    water \
    camera  \
    loader \
    object

SOURCES += \
    camera/camera.cpp \
    terrain/terrain.cpp \
    support/glwidget.cpp \
    support/main.cpp \
    support/mainwindow.cpp \
    loader/glm.cpp \
    loader/resourceloader.cpp \
    loader/targa.cpp \
    loader/textureloader.cpp \
    common/CS123Matrix.cpp \
    common/CS123Matrix.inl \
    common/CS123Vector.inl \
    terrain/random_terrain.cpp \
    terrain/heightmap_terrain.cpp \
    common/utils.cpp \
    water/fluidCPU.cpp \
    water/fluidGPU.cpp \
    water/particle.cpp \
    water/particlesource.cpp \
    object/box.cpp \
    object/object.cpp

HEADERS  += \
    camera/camera.h \
    terrain/terrain.h \
    support/glwidget.h \
    support/mainwindow.h \
    support/ui_mainwindow.h \
    loader/glm.h \
    loader/resourceloader.h \
    loader/targa.h \
    loader/textureloader.h \
    common/CS123Algebra.h \
    common/CS123Common.h \
    common/types.h \
    common/vector.h \
    terrain/random_terrain.h \
    terrain/heightmap_terrain.h \
    common/debug_marco.h \
    common/utils.h \
    common/grid.h \
    water/fluidCPU.h \
    water/fluid_global.h \
    water/fluidGPU.h \
    water/particle.h \
    water/particlesource.h \
    object/box.h \
    object/object.h \
    object/object_defs.h

OTHER_FILES += cuda/fluid_compute.cu \
    cuda/test.cu \
    shaders/refract.vert \
    shaders/refract.frag \
    shaders/reflect.vert \
    shaders/reflect.frag \
    shaders/brightpass.frag \
    shaders/blur.frag \
    shaders/fresnel.vert \
    shaders/fresnel.vars \
    shaders/fresnel.frag \
    shaders/f2.frag \
    shaders/f2.vert \
    shaders/point.frag \
    shaders/point.vert \
    shaders/splash.vert \
    shaders/splash.frag \
    shaders/foam.vert \
    shaders/foam.frag

CUDA_SOURCES += cuda/test.cu cuda/fluid_compute.cu

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

FORMS    += \
    support/mainwindow.ui

#unix|win32: LIBS += -lGLU -lglut

unix:!macx:!symbian: LIBS += -lGLU -lglut
