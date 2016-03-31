#-------------------------------------------------
#
# Project created by QtCreator 2015-12-14T18:02:44
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = ErrorPropagation
TEMPLATE = app
DESTDIR = ./bin
OBJECTS_DIR = ./obj
MOC_DIR = ./moc
RCC_DIR = ./rcc
UI_DIR = ./ui

SOURCES += main.cpp\
        errorprop.cpp \
    propagationform.cpp \
    firststep.cpp \
    thirdstep.cpp \
    geomtransform.cpp \
    allsteps.cpp \
    enhancedallsteps.cpp

HEADERS  += errorprop.h \
    propagationform.h \
    firststep.h \
    thirdstep.h \
    geomtransform.h \
    allsteps.h \
    enhancedallsteps.h

FORMS    += errorprop.ui \
    propagationform.ui
