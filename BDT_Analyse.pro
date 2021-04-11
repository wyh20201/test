TEMPLATE = app
CONFIG += console
CONFIG += c++11
QT += core

SOURCES += main.cpp \
    source.cpp \
    matrix.cpp

include(deployment.pri)
qtcAddDeployment()

HEADERS += \
    source.h \
    matrix.h

