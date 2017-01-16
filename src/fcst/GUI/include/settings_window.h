//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2014 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: settings_window.h
//    - Description: GUI for showing FCST running settings
//    - Developers: Simon Mattern, 2016
//
//---------------------------------------------------------------------------
#ifndef SETTINGS_WINDOW_H
#define SETTINGS_WINDOW_H

//QT
#include <QtGui>

//STD
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <string>
#include <ctime>

namespace FCSTGUI
{
    /**
     * 
     * \brief Class used to create simulation settings window.
     * 
     * This class is used to show the user a window with simulation settings. The window is opened when the user selects either
     * Options > Advanced Settings or the tools button on the toolbar on the main window.
     * 
     * @author Simon Matterns, 2016
     * 
     */
    class settingsWindow : public QMainWindow
    {
        Q_OBJECT // This is a Qt object necessary to give Qt signal and slots information. Otherwise button in the toolbar would not work.
    public:
        /**
         * Destructor
         */
        ~settingsWindow();
        
        /**
         * second Constructor.
         */
        settingsWindow(QSettings *guiSettings =0);
        
    private slots:
        /**
         * Close the settings window
         */       
        void closeSettings();
        /**
         * Save the changes in settings.ini file
         */
        void saveSettings();
        /**
         * Slot for browsing a new bin file
         */
        void browse();
        
    private:
        /**
         * Creates all Widgets and Widgets settings
         */
        void createWidgets();
        /**
         * Ceates Layout of the window
         */
        void createLayout();
        /**
         * Loads the values from settings.ini into the widgets
         */
        void LoadSettings();
        /**
         * Function to check if a file exists
         */
        bool fileExists(QString path);
        /**
         * Create QLineEdit
         */    
        QLineEdit* createQLineEdit();
        
        /**
         * Hold the reviewed data from settings.ini file
         */
        QPointer<QSettings> gui_setting;
        /**
         * Shown Window
         */
        QPointer<QDialog> settings_dialog;
        /**
         * Spacer to hold QPushButton on the right position
         */   
        QSpacerItem * spacer;
        /**
         * ComboBox shows selected binary
         */       
        QPointer<QComboBox> used_binary_box;
        /**
         * QLineEdit for showing argument
         */
        QPointer<QLineEdit> sim_argument_edit;
        /**
         * QLineEdit for showing main file name
         */
        QPointer<QLineEdit> main_file_edit;
        /**
         * QLineEdit for showing data file name
         */
        QPointer<QLineEdit> data_file_edit;
        /**
         * QLineEdit for showing opt. file name
         */
        QPointer<QLineEdit> opt_file_edit;
        /**
         * QSpinBox for setting number of CPU cores
         */
        QPointer<QSpinBox> core_number_spinbox;
        /**
         * Showing 2D bin path
         */
        QPointer<QComboBox> bin_2D_path_box;
        /**
         * Showing 3D bin path
         */
        QPointer<QComboBox> bin_3D_path_box;
        /**
         * Button for closing the window
         */
        QPointer<QPushButton> cancel_button;
        /**
         * Button for saving the changes
         */
        QPointer<QPushButton> save_button;
        /**
         * Button for browsing 
         */
        QPointer<QPushButton> browse_2Dpath_button;
        /**
         * Button for browsing 
         */
        QPointer<QPushButton> browse_3Dpath_button; 
        
        
    };
}

#endif // SETTINGS_WINDOW_H
