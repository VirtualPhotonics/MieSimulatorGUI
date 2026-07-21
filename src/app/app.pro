##################################################
## app project file
##################################################

QT += core gui widgets printsupport

greaterThan(QT_MAJOR_VERSION, 6): QT += widgets

TARGET = MieSimulatorGUI_v2_1
CONFIG -= -qt-freetype
TEMPLATE = app
DEFINES += QT_DEPRECATED_WARNINGS
qcustomplot_compiler_flags.input = SOURCES
QMAKE_CXXFLAGS += -Wno-deprecated-declarations

win32 {
    RC_ICONS = misc/MieSimulatorGUI.ico
}

macos {
    ICON = misc/MieSimulatorGUI.icns
}

unix {
    isEmpty(PREFIX) {
        PREFIX = /usr
    }

    target.path = $$PREFIX/bin

    shortcutfiles.files = misc/MieSimulatorGUI.desktop
    shortcutfiles.path = $$PREFIX/share/applications/
    data.files += misc/MieSimulatorGUI.png
    data.path = $$PREFIX/share/icons/hicolor/256x256/

    INSTALLS += shortcutfiles
    INSTALLS += data
}

INSTALLS += target

RESOURCES += MieResource.qrc

FORMS += \
    dialog/displaydialog.ui \
    dialog/mainwindow.ui \
    dialog/optionsdialog.ui

DISTFILES += \
    misc/MieSimulatorGUI.desktop \
    misc/MieSimulatorGUI.ico \
    misc/MieSimulatorGUI.png

HEADERS += \
    lib/qcustomplot.h \
    calc/calculate.h \
    calc/miesimulation.h \
    calc/utilities.h \
    dialog/displaydialog.h \
    dialog/mainwindow.h \
    dialog/mainwindowsupport.h \
    dialog/optionsdialog.h \
    dialog/plotdata.h \
    parameters.h

SOURCES += \    
    lib/qcustomplot.cpp \
    calc/calculate.cpp \
    calc/miesimulation.cpp \
    calc/utilities.cpp \
    dialog/displaydialog.cpp \
    dialog/mainwindow.cpp \
    dialog/mainwindowsupport.cpp \
    dialog/optionsdialog.cpp \
    dialog/plotdata.cpp \
    parameters.cpp \
    main.cpp
