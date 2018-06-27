#-------------------------------------------------
#
# Project created by QtCreator 2014-06-06T12:06:46
#
#-------------------------------------------------
QT       += core gui \
        core gui widgets \
        printsupport


greaterThan(QT_MAJOR_VERSION, 5): QT += widgets

TARGET = MieSimulatorGUI_v1_1
CONFIG -= -qt-freetype
TEMPLATE = app
DEFINES += QT_DEPRECATED_WARNINGS

RESOURCES += \
    MieRes.qrc

FORMS += \
    displaydialog.ui \
    mainwindow.ui \
    optionsdialog.ui

DISTFILES += \
    icon/Banner.png \
    icon/vpIcon.png \
    vpIcon.rc

HEADERS += \
    qwt/qtconcurrentrun.h \
    qwt/qwt.h \
    qwt/qwt_abstract_legend.h \
    qwt/qwt_abstract_scale.h \
    qwt/qwt_abstract_scale_draw.h \
    qwt/qwt_abstract_slider.h \
    qwt/qwt_clipper.h \
    qwt/qwt_color_map.h \
    qwt/qwt_curve_fitter.h \
    qwt/qwt_dyngrid_layout.h \
    qwt/qwt_global.h \
    qwt/qwt_graphic.h \
    qwt/qwt_interval.h \
    qwt/qwt_legend.h \
    qwt/qwt_legend_data.h \
    qwt/qwt_legend_label.h \
    qwt/qwt_math.h \
    qwt/qwt_null_paintdevice.h \
    qwt/qwt_painter.h \
    qwt/qwt_painter_command.h \
    qwt/qwt_pixel_matrix.h \
    qwt/qwt_plot.h \
    qwt/qwt_plot_canvas.h \
    qwt/qwt_plot_curve.h \
    qwt/qwt_plot_dict.h \
    qwt/qwt_plot_item.h \
    qwt/qwt_plot_layout.h \
    qwt/qwt_plot_seriesitem.h \
    qwt/qwt_point_3d.h \
    qwt/qwt_point_data.h \
    qwt/qwt_point_mapper.h \
    qwt/qwt_point_polar.h \
    qwt/qwt_polar.h \
    qwt/qwt_polar_canvas.h \
    qwt/qwt_polar_curve.h \
    qwt/qwt_polar_global.h \
    qwt/qwt_polar_grid.h \
    qwt/qwt_polar_item.h \
    qwt/qwt_polar_itemdict.h \
    qwt/qwt_polar_layout.h \
    qwt/qwt_polar_plot.h \
    qwt/qwt_round_scale_draw.h \
    qwt/qwt_samples.h \
    qwt/qwt_scale_div.h \
    qwt/qwt_scale_draw.h \
    qwt/qwt_scale_engine.h \
    qwt/qwt_scale_map.h \
    qwt/qwt_scale_widget.h \
    qwt/qwt_series_data.h \
    qwt/qwt_series_store.h \
    qwt/qwt_slider.h \
    qwt/qwt_spline.h \
    qwt/qwt_symbol.h \
    qwt/qwt_text.h \
    qwt/qwt_text_engine.h \
    qwt/qwt_text_label.h \
    qwt/qwt_transform.h \
    qwt/qwt_widget_overlay.h \
    calculate.h \
    displaydialog.h \
    mainwindow.h \
    mainwindowsupport.h \
    miesimulation.h \
    optionsdialog.h \
    parameters.h \
    plotdata.h \
    qcustomplot.h \
    utilities.h

SOURCES += \
    qwt/qwt_abstract_legend.cpp \
    qwt/qwt_abstract_scale.cpp \
    qwt/qwt_abstract_scale_draw.cpp \
    qwt/qwt_abstract_slider.cpp \
    qwt/qwt_clipper.cpp \
    qwt/qwt_color_map.cpp \
    qwt/qwt_curve_fitter.cpp \
    qwt/qwt_dyngrid_layout.cpp \
    qwt/qwt_graphic.cpp \
    qwt/qwt_interval.cpp \
    qwt/qwt_legend.cpp \
    qwt/qwt_legend_data.cpp \
    qwt/qwt_legend_label.cpp \
    qwt/qwt_math.cpp \
    qwt/qwt_null_paintdevice.cpp \
    qwt/qwt_painter.cpp \
    qwt/qwt_painter_command.cpp \
    qwt/qwt_pixel_matrix.cpp \
    qwt/qwt_plot.cpp \
    qwt/qwt_plot_axis.cpp \
    qwt/qwt_plot_canvas.cpp \
    qwt/qwt_plot_curve.cpp \
    qwt/qwt_plot_dict.cpp \
    qwt/qwt_plot_item.cpp \
    qwt/qwt_plot_layout.cpp \
    qwt/qwt_plot_seriesitem.cpp \
    qwt/qwt_point_3d.cpp \
    qwt/qwt_point_data.cpp \
    qwt/qwt_point_mapper.cpp \
    qwt/qwt_point_polar.cpp \
    qwt/qwt_polar_canvas.cpp \
    qwt/qwt_polar_curve.cpp \
    qwt/qwt_polar_grid.cpp \
    qwt/qwt_polar_item.cpp \
    qwt/qwt_polar_itemdict.cpp \
    qwt/qwt_polar_layout.cpp \
    qwt/qwt_polar_plot.cpp \
    qwt/qwt_round_scale_draw.cpp \
    qwt/qwt_scale_div.cpp \
    qwt/qwt_scale_draw.cpp \
    qwt/qwt_scale_engine.cpp \
    qwt/qwt_scale_map.cpp \
    qwt/qwt_scale_widget.cpp \
    qwt/qwt_series_data.cpp \
    qwt/qwt_slider.cpp \
    qwt/qwt_spline.cpp \
    qwt/qwt_symbol.cpp \
    qwt/qwt_text.cpp \
    qwt/qwt_text_engine.cpp \
    qwt/qwt_text_label.cpp \
    qwt/qwt_transform.cpp \
    qwt/qwt_widget_overlay.cpp \
    calculate.cpp \
    displaydialog.cpp \
    main.cpp \
    mainwindow.cpp \
    mainwindowsupport.cpp \
    miesimulation.cpp \
    optionsdialog.cpp \
    parameters.cpp \
    plotdata.cpp \
    qcustomplot.cpp \
    utilities.cpp

