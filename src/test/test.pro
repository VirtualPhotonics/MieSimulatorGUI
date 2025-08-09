# This is the project file for the unit test.

# Set the template to create an application executable.
TEMPLATE = app

# Add the Qt Test library, which is required for unit testing.
# We also include 'widgets' as the app under test is a GUI app,
# and some functions might implicitly require it (though not directly in 'calculate').
QT += testlib widgets

# The target name for our test executable.
TARGET = test_calculate

# Configure the build to be in debug and release modes.
CONFIG += warn_on console

# Include the headers from the main application to access the classes we want to test.
INCLUDEPATH += ../app \
               ../app/calc

# Add the sources for our test suite and the classes we are testing.
# We include the calculate.cpp directly, as well as mock implementations for its dependencies.
# The `mock_` files are simple stubs that allow us to test the `calculate` class in isolation.
SOURCES += \
    test_calculate.cpp \
    ../app/parameters.cpp \
    ../app/calc/calculate.cpp \
    ../app/calc/miesimulation.cpp \
    ../app/calc/utilities.cpp

# Add the headers for the test suite and the classes we are testing.
HEADERS += \
    test_calculate.h \
    ../app/parameters.h \
    ../app/calc/calculate.h \
    ../app/calc/miesimulation.h \
    ../app/calc/utilities.h

# Define some preprocessor macros if needed.
DEFINES += QT_TESTLIB_WIDGETS_LIB
