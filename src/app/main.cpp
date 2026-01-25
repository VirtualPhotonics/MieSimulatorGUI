/***********************************************************************
** This program provides a GUI tool for calculating the characteristics
** of Mie scatterers. The tool calculates the spectral dependence of the
** scattering coefficient, scattering cross section, reduced scattering
** coefficient, phase function, average cosine of the phase function,
** and scattering matrix entries for a single or series of wavelengths.
** These calculations can be run for mono-disperse and poly-disperse
** distributions based on Lorentz-Mie solution. The tool is also capable
** of estimating fitting parameters from the reduced scattering
** coefficient curve and is limited to the characteristics of the
** far-field scattering pattern. The Mie Simulator algorithm is based on
** the BHMIE code available from Bohren and Huffman. It was originally
** developed for National Short Course in Computational Biophotonics
** at Beckman Laser Institute (BLI) by Janaka Ranasinghesagara.
**
** This library is distributed in the hope that it will be useful, but
** WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
** Lesser General Public License and QT License for more details
**
** $QT_BEGIN_LICENSE:GPL$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing
** terms and conditions see https://www.qt.io/terms-conditions. For
** further information use the contact form at
** https://www.qt.io/contact-us.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 3 or (at your option) any later version
** approved by the KDE Free Qt Foundation. The licenses are as published
** by the Free Software Foundation and appearing in the file LICENSE.GPL3
** included in the packaging of this file. Please review the following
** information to ensure the GNU General Public License requirements will
** be met: https://www.gnu.org/licenses/gpl-3.0.html.**
** $QT_END_LICENSE$//
**
************************************************************************/

#include <QApplication>
#include <QtPlugin>
#include <QScreen>
#include <QFont>
#include <QStyleFactory>
#include <QGuiApplication>
#include "dialog/mainwindow.h"

int main(int argc, char *argv[])
{
    QApplication::setHighDpiScaleFactorRoundingPolicy(
        Qt::HighDpiScaleFactorRoundingPolicy::PassThrough);

    QApplication app(argc, argv);

    // Use this for cross-platform font consistency
    #ifdef Q_OS_WIN
        QString family = "Arial";
        QFont appFont(family, 9);
    #elif defined(Q_OS_MAC)
        QString family = "Helvetica";
        QFont appFont(family, 11);
    #else
        QString family = "Liberation Sans";
        QFont appFont(family, 10);
    #endif
    QApplication::setFont(appFont);
    app.setStyle(QStyleFactory::create("Fusion"));

    MainWindow mainWindow;

    // Use the primary window setting
    QScreen *primaryScreen = QGuiApplication::primaryScreen();
    if (primaryScreen)
    {
        // Get availableSize()
        QSize availableSize = primaryScreen->availableSize();

        // Use a reference resolution
        QSize referenceSize(1920, 1080);

        // Calculate scale factor
        double scaleWidth = (double)availableSize.width() / referenceSize.width();
        double scaleHeight = (double)availableSize.height() / referenceSize.height();

        // Use the smaller ratio to ensure the window always fits the screen
        double scaleFactor = qMin(scaleWidth, scaleHeight);

        // Define your base design size
        QSize baseSize(1120, 760);

        // Apply scale factor
        int newWidth = static_cast<int>(baseSize.width() * scaleFactor);
        int newHeight = static_cast<int>(baseSize.height() * scaleFactor);

        // Subtract a small margin for window decorations/taskbar padding
        newWidth = qMin(newWidth, availableSize.width() - 50);
        newHeight = qMin(newHeight, availableSize.height() - 50);

        mainWindow.resize(newWidth, newHeight);
    }

    mainWindow.show();
    return app.exec();
}
