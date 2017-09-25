
QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets
#QT += widgets
#qtHaveModule(printsupport): QT += printsupport

TARGET = apTool
TEMPLATE = apTool
#LIBS += -L/usr/lib/x86_64-linux-gnu -L/usr/local/lib -lopencv_highgui -lopencv_imgproc -lopencv_core -lQtGui -lQtCore -lpthread


#CONFIG += debug
LIBS += `pkg-config opencv --libs`

LIBS += -lalglib

#LIBS += -L/usr/local/lib -lopencv_highgui -lopencv_core -lopencv_imgcodecs -lopencv_video
INCLUDEPATH += /usr/local/include/opencv2/
INCLUDEPATH += /usr/include/libalglib/
#INCLUDEPATH += /usr/include/opencv2/

TARGET = apTool
TEMPLATE = app


SOURCES += main.cpp\
        aptool.cpp \
    s_hull_pro.cpp \
    imageview.cpp \
    lwidget.cpp

HEADERS  += aptool.h \
    s_hull_pro.h \
    imageview.h \
    lwidget.h

FORMS    += aptool.ui \
  imageview.ui
