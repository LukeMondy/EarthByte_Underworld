####################################################################################################################### 
# Omega Lib build for gLucifer
#---------------------------------------------------------------------------------------------------------------------
# Copyright 2013							Monash University
# Authors:										
#  Owen Kaluza							owen.kaluza@monash.edu
#######################################################################################################################
# Source files
SET( srcs 
        OpenGLViewer.cpp
          base64.cpp
          ColourMap.cpp
          DrawingObject.cpp
          Extensions.cpp
          FontSans.cpp
          Geometry.cpp
          GraphicsUtil.cpp
          Lines.cpp
          Model.cpp
          Points.cpp
          QuadSurfaces.cpp
          Shaders.cpp
          Shapes.cpp
          Tracers.cpp
          TriSurfaces.cpp
          Vectors.cpp
          VideoEncoder.cpp
          View.cpp
          Win.cpp
          jpeg/jpge.cpp
          sqlite3/src/sqlite3.c
        )
    
#######################################################################################################################
# Headers
SET( headers 
        ViewerApp.h
        OpenGLViewer.h
        ) 
        
#######################################################################################################################
# Setup compile info

add_definitions(-DUSE_OMEGALIB)
add_library(common OBJECT ${srcs})

add_executable(gLucifer GLucifer.cpp GLuciferViewer.cpp InteractiveViewer.cpp ${headers} $<TARGET_OBJECTS:common>)

add_executable(FractalView FractalView.cpp FractalViewer.cpp ${headers} $<TARGET_OBJECTS:common>)

set_target_properties(gLucifer PROPERTIES FOLDER examples)
target_link_libraries(gLucifer dl	omega)

set_target_properties(FractalView PROPERTIES FOLDER examples)
target_link_libraries(FractalView dl	omega)

