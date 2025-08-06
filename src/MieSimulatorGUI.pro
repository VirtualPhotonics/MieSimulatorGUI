#-------------------------------------------------
#
# Project created by QtCreator 2014-06-06T12:06:46
#
#-------------------------------------------------
QT       += core gui \
        core gui widgets \
        printsupport


greaterThan(QT_MAJOR_VERSION, 5): QT += widgets

TARGET = MieSimulatorGUI_continuous
CONFIG -= -qt-freetype
TEMPLATE = app
DEFINES += QT_DEPRECATED_WARNINGS

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

RESOURCES += \
    MieResource.qrc

FORMS += \
    displaydialog.ui \
    mainwindow.ui \
    optionsdialog.ui

DISTFILES += \
    misc/MieSimulatorGUI.desktop \
    misc/MieSimulatorGUI.ico \
    misc/MieSimulatorGUI.png

HEADERS += \
    lib/qwt/qtconcurrentrun.h \
    lib/qwt/qwt.h \
    lib/qwt/qwt_abstract_legend.h \
    lib/qwt/qwt_abstract_scale.h \
    lib/qwt/qwt_abstract_scale_draw.h \
    lib/qwt/qwt_abstract_slider.h \
    lib/qwt/qwt_clipper.h \
    lib/qwt/qwt_color_map.h \
    lib/qwt/qwt_curve_fitter.h \
    lib/qwt/qwt_dyngrid_layout.h \
    lib/qwt/qwt_global.h \
    lib/qwt/qwt_graphic.h \
    lib/qwt/qwt_interval.h \
    lib/qwt/qwt_legend.h \
    lib/qwt/qwt_legend_data.h \
    lib/qwt/qwt_legend_label.h \
    lib/qwt/qwt_math.h \
    lib/qwt/qwt_null_paintdevice.h \
    lib/qwt/qwt_painter.h \
    lib/qwt/qwt_painter_command.h \
    lib/qwt/qwt_pixel_matrix.h \
    lib/qwt/qwt_plot.h \
    lib/qwt/qwt_plot_canvas.h \
    lib/qwt/qwt_plot_curve.h \
    lib/qwt/qwt_plot_dict.h \
    lib/qwt/qwt_plot_item.h \
    lib/qwt/qwt_plot_layout.h \
    lib/qwt/qwt_plot_seriesitem.h \
    lib/qwt/qwt_point_3d.h \
    lib/qwt/qwt_point_data.h \
    lib/qwt/qwt_point_mapper.h \
    lib/qwt/qwt_point_polar.h \
    lib/qwt/qwt_polar.h \
    lib/qwt/qwt_polar_canvas.h \
    lib/qwt/qwt_polar_curve.h \
    lib/qwt/qwt_polar_global.h \
    lib/qwt/qwt_polar_grid.h \
    lib/qwt/qwt_polar_item.h \
    lib/qwt/qwt_polar_itemdict.h \
    lib/qwt/qwt_polar_layout.h \
    lib/qwt/qwt_polar_plot.h \
    lib/qwt/qwt_round_scale_draw.h \
    lib/qwt/qwt_samples.h \
    lib/qwt/qwt_scale_div.h \
    lib/qwt/qwt_scale_draw.h \
    lib/qwt/qwt_scale_engine.h \
    lib/qwt/qwt_scale_map.h \
    lib/qwt/qwt_scale_widget.h \
    lib/qwt/qwt_series_data.h \
    lib/qwt/qwt_series_store.h \
    lib/qwt/qwt_slider.h \
    lib/qwt/qwt_spline.h \
    lib/qwt/qwt_symbol.h \
    lib/qwt/qwt_text.h \
    lib/qwt/qwt_text_engine.h \
    lib/qwt/qwt_text_label.h \
    lib/qwt/qwt_transform.h \
    lib/qwt/qwt_widget_overlay.h \
    lib/qcustomplot.h \
    calculate.h \
    displaydialog.h \
    mainwindow.h \
    mainwindowsupport.h \
    miesimulation.h \
    optionsdialog.h \
    parameters.h \
    plotdata.h \

    utilities.h

SOURCES += \
    lib/qwt/qwt_abstract_legend.cpp \
    lib/qwt/qwt_abstract_scale.cpp \
    lib/qwt/qwt_abstract_scale_draw.cpp \
    lib/qwt/qwt_abstract_slider.cpp \
    lib/qwt/qwt_clipper.cpp \
    lib/qwt/qwt_color_map.cpp \
    lib/qwt/qwt_curve_fitter.cpp \
    lib/qwt/qwt_dyngrid_layout.cpp \
    lib/qwt/qwt_graphic.cpp \
    lib/qwt/qwt_interval.cpp \
    lib/qwt/qwt_legend.cpp \
    lib/qwt/qwt_legend_data.cpp \
    lib/qwt/qwt_legend_label.cpp \
    lib/qwt/qwt_math.cpp \
    lib/qwt/qwt_null_paintdevice.cpp \
    lib/qwt/qwt_painter.cpp \
    lib/qwt/qwt_painter_command.cpp \
    lib/qwt/qwt_pixel_matrix.cpp \
    lib/qwt/qwt_plot.cpp \
    lib/qwt/qwt_plot_axis.cpp \
    lib/qwt/qwt_plot_canvas.cpp \
    lib/qwt/qwt_plot_curve.cpp \
    lib/qwt/qwt_plot_dict.cpp \
    lib/qwt/qwt_plot_item.cpp \
    lib/qwt/qwt_plot_layout.cpp \
    lib/qwt/qwt_plot_seriesitem.cpp \
    lib/qwt/qwt_point_3d.cpp \
    lib/qwt/qwt_point_data.cpp \
    lib/qwt/qwt_point_mapper.cpp \
    lib/qwt/qwt_point_polar.cpp \
    lib/qwt/qwt_polar_canvas.cpp \
    lib/qwt/qwt_polar_curve.cpp \
    lib/qwt/qwt_polar_grid.cpp \
    lib/qwt/qwt_polar_item.cpp \
    lib/qwt/qwt_polar_itemdict.cpp \
    lib/qwt/qwt_polar_layout.cpp \
    lib/qwt/qwt_polar_plot.cpp \
    lib/qwt/qwt_round_scale_draw.cpp \
    lib/qwt/qwt_scale_div.cpp \
    lib/qwt/qwt_scale_draw.cpp \
    lib/qwt/qwt_scale_engine.cpp \
    lib/qwt/qwt_scale_map.cpp \
    lib/qwt/qwt_scale_widget.cpp \
    lib/qwt/qwt_series_data.cpp \
    lib/qwt/qwt_slider.cpp \
    lib/qwt/qwt_spline.cpp \
    lib/qwt/qwt_symbol.cpp \
    lib/qwt/qwt_text.cpp \
    lib/qwt/qwt_text_engine.cpp \
    lib/qwt/qwt_text_label.cpp \
    lib/qwt/qwt_transform.cpp \
    lib/qwt/qwt_widget_overlay.cpp \
    lib/qcustomplot.cpp \
    calculate.cpp \
    displaydialog.cpp \
    main.cpp \
    mainwindow.cpp \
    mainwindowsupport.cpp \
    miesimulation.cpp \
    optionsdialog.cpp \
    parameters.cpp \
    plotdata.cpp \    
    utilities.cpp

