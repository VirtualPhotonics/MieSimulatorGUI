##################################################
## test project file for unit tests
##################################################

TEMPLATE = app

QT += testlib widgets

TARGET = MieSimulatorGUI_test
CONFIG += warn_on console
DEFINES += QT_DEPRECATED_WARNINGS

INCLUDEPATH += ../app \
               ../app/calc

SOURCES += \
    test_main.cpp \
    test_parameters.cpp \
    test_calculate.cpp \
    test_miesimulation.cpp \
    test_utilities.cpp \
    ../app/parameters.cpp \
    ../app/calc/calculate.cpp \
    ../app/calc/miesimulation.cpp \
    ../app/calc/utilities.cpp

HEADERS += \    
    test_parameters.h \
    test_calculate.h \
    test_miesimulation.h \
    test_utilities.h \
    ../app/parameters.h \
    ../app/calc/calculate.h \
    ../app/calc/miesimulation.h \
    ../app/calc/utilities.h

