##################################################
## test project file for unit tests
##################################################

TEMPLATE = app

QT += testlib widgets printsupport

TARGET = MieSimulatorGUI_test
CONFIG += warn_on console
DEFINES += QT_DEPRECATED_WARNINGS

INCLUDEPATH += $$PWD/../app \
               $$PWD/../app/calc \
               $$PWD/../app/dialog \
               $$PWD/../app/lib

FORMS += \
    ../app/dialog/mainwindow.ui

SOURCES += \
    test_main.cpp \
    test_parameters.cpp \
    test_calculate.cpp \
    test_miesimulation.cpp \
    test_utilities.cpp \
    test_plotdata.cpp \
    test_mainwindowsupport.cpp \
    ../app/parameters.cpp \
    ../app/calc/calculate.cpp \
    ../app/calc/miesimulation.cpp \
    ../app/calc/utilities.cpp \
    ../app/dialog/plotdata.cpp \
    ../app/dialog/mainwindowsupport.cpp \
    ../app/lib/qcustomplot.cpp

HEADERS += \
    test_parameters.h \
    test_calculate.h \
    test_miesimulation.h \
    test_utilities.h \
    test_plotdata.h \
    test_mainwindowsupport.h \
    ../app/parameters.h \
    ../app/calc/calculate.h \
    ../app/calc/miesimulation.h \
    ../app/calc/utilities.h \
    ../app/dialog/plotdata.h \
    ../app/dialog/mainwindowsupport.h \
    ../app/lib/qcustomplot.h