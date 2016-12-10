//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2014 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: main.cc
//    - Description: Main file for running FCST GUI
//    - Developers: Phil Wardlaw, 2014
//
//---------------------------------------------------------------------------


#include <QApplication>
#include <QSplashScreen>
#include <QTimer>

#include <main_window.h>


/**
 * The main file forFCST parameterGUI.
 * Modified from the deal.II parameterGUI, which was developed by Martin Steigemann, Wolfgang Bangerth, 2010.
 *
 * The majority of the implementation lies in main_window.h, main_window.cc
 *
 * @ingroup ParameterGui
 * 
 * @author Philip Wardlaw 2014
 */
int main(int argc, char *argv[])
{
  // init resources such as icons or graphics
  Q_INIT_RESOURCE(application);
  QApplication app(argc, argv);

  // setup a splash screen
  QSplashScreen* splash = new QSplashScreen;
  splash->setPixmap(QPixmap(":/images/logo_fcst_64.png"));
  splash->show();

  // and close it after 3000 ms
  QTimer::singleShot(3000, splash, SLOT(close()));

  // setup the application name
  app.setApplicationName("OpenFCST: Fuel Cell Simulation Toolbox");

  // give command line arguments to main_win, which it actually just ignores
  FCSTGUI::MainWindow* main_win =    new FCSTGUI::MainWindow (argv[1]);

  // show the main window with a short delay
  QTimer::singleShot(1500, main_win, SLOT(showMaximized()));
  
  // so we can see the splash screen
  return app.exec();
}

