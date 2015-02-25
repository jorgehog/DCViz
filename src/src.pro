CONFIG -= app_bundle
CONFIG -= qt

COMMON_CXXFLAGS = -std=c++11

QMAKE_CXXFLAGS += \
    $$COMMON_CXXFLAGS

QMAKE_CXXFLAGS_DEBUG += \
    $$COMMON_CXXFLAGS

QMAKE_CXXFLAGS_RELEASE += \
    $$COMMON_CXXFLAGS \
    -O3 \
    -DNDEBUG

QMAKE_LFLAGS_RELEASE -= -O1
QMAKE_LFLAGS_RELEASE += -O3

INCLUDEPATH += /usr/include/python2.7
LIBS += -lpython2.7



TOP_PWD = $$PWD/../

TEMPLATE = lib

TARGET = ../lib/DCViz

QMAKE_LFLAGS_DEBUG += -g

HEADERS = DCViz.h

SOURCES = DCViz.cpp

QMAKE_PRE_LINK += $(MKDIR) $$PWD/../lib $$shadowed($$PWD)/../lib

!equals(PWD, $${OUT_PWD}) {
    QMAKE_POST_LINK += $(COPY_DIR) $$OUT_PWD/../lib $$TOP_PWD
}



