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
    QFont defaultFont;
    defaultFont.setPointSize(9);
    defaultFont.setFamily("Arial");
    app.setFont(defaultFont);
    app.setStyle(QStyleFactory::create("Fusion"));

    MainWindow mainWindow;

    // Use the primary window setting
    QScreen *primaryScreen = QGuiApplication::primaryScreen();
    if (primaryScreen)
    {
        QSize currentScreenSize = primaryScreen->availableSize();

        QSize referenceScreenSize(1920, 1080);

        double scaleFactor = qMin(
            (double)currentScreenSize.width() / referenceScreenSize.width(),
            (double)currentScreenSize.height() / referenceScreenSize.height()
            );

        QSize baseSize(1120, 768);

        int newWidth = baseSize.width() * scaleFactor;
        int newHeight = baseSize.height() * scaleFactor;

        // Resize the window
        mainWindow.resize(newWidth, newHeight);
    }

    mainWindow.show();
    return app.exec();
}
