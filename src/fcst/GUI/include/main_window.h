//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2014 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: main_window.h
//    - Description: GUI for running FCST
//    - Developers: Phil Wardlaw, Simon Mattern, Marc Secanell, 2014-16
//
//---------------------------------------------------------------------------


#ifndef MAINWINDOW_H
#define MAINWINDOW_H


//QT
#include <QApplication>
#include <QMainWindow>
#include <QTreeWidget>
#include <QDialog>
#include <QSettings>
#include <QList>
#include <QtGui>
#include <QTextEdit>
#include <QPushButton>
#include <QDir>
#include <QStringList>
#include <QLabel>
#include <QFile>
#include <QFileInfo>
#include <QListWidget>
#include <QListWidgetItem>
#include <QFileSystemWatcher>
#include <QFont>
//STD
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <string>
#include <ctime>
//Linux Sleep mechanism
#include <unistd.h>
//Deal.ii contribs developed by Martin Steigemann, Wolfgang Bangerth, 2010.
#include "parameter_delegate.h"
#include "xml_parameter_reader.h"
#include "xml_parameter_writer.h"
#include "settings_window.h"


namespace FCSTGUI
{
/**
 * \brief The MainWindow class of the the FCST parameterGUI.
 *
 * Modified from the deal.II parameterGUI, which was developed by Martin Steigemann, Wolfgang Bangerth, 2010.
 *
 * The FCST parameter GUI allows users to start a new "FCST Project" whereby fcst is called to generate
 * default parameter files in a logical sequence (main >> data & opt), allowing the user to manually edit
 * these files, and then eventually run a simulation using FCST and the edited parameter files.
 * Output streams and files are displayed within the user interface, providing the user
 * with a integrated environment in which they can setup and run a simulation from start to finish.
 *
 * <h2>Developer Guide </h2>
 * This application inherits the Qt QMainWindow to provide the functionality of a "main" application,
 * using many of its interfaces.
 * (For developers who working on the code: Use "textOut->insertPlanText(VariableYouWantToSee)"
 * too see the content of any variable in right top corner of the interface (QTextEdit ouput).
 * Note to delete the output code before you create a request or publishing the GUI code).
 *
 * <h3>Visual Elements </h3>
 * The body of the application is composed of several widgets. Each widget is associated with a
 * Box layout element so that it can dynamically expand as the QMainWindow size changes.
 *
 *
 * Top left: A QTabWidget is used to house several QTreeWidgets which display the loaded parameter files.
 * Bottom Left: A QWidget supports a QLabel (used to display the logo image) and QPushButton
 * Top Right: A read only QTextEdit is used to display fcst output
 * Bottom Right: A QListWidget widget is used to represent a list of files created by fcst
 *
 * Menus: The menuBar, a member of QMainWindow is populated by two other QMenu's
 * Status Bar: The statusBar is a member of QMainWindow
 *
 * <h3>Maintaining and modifying</h3>
 * Private functions are used to load and set up aspects of the GUI, whilst private slots are
 * used to handle events invoked by the GUI (i.e. a button click). To modify the behaviour of
 * the how the GUI reacts to user interactions please review the slot functions.
 *
 * The main driving function of the program is the next_click() slot. When the user presses the
 * next button/starts a new project/opens an existing project next_click() will try progress
 * the program state "currentState" along the GUI's enumeration of states "States". The operations
 * associated with each state are as follows:
 *
 *
 * Start:
 *  - The program has just started
 *  - Working file directory must be set
 *  - A main file must be prepared by calling FCST
 *  - Load the main file to the GUI
 *  - Progress to next state
 *
 * Main:
 *  - Assumption: FCST has been used to generate a main file and the user has had time to edit it
 *  - Save the main file and use it to generate data and optimization files
 *  - Load these files to GUI
 *  - Progress to next state
 *
 * Ready:
 *  - Assumption: main and data files have been loaded and the user has had time to edit them
 *  - The GUI can get to this state if the user loads an existing project, by progressing from the Main state, or by retreating from the Running state
 *  - Save these files and start a simulation using FCST
 *  - Begin watching for new files in working directory (detect as FCST output)
 *  - Progress to next state
 *
 * Running:
 *  - FCST has been called, output is being piped to text element
 *  - The user has the option to end current simulation
 *  - Retreat to previous state
 *
 * <h3>Assumptions of FCST</h3>
 * This program is written with the following assumption of FCST binary's (BIN) interfaces
 *
 * BIN -p:
 *  - Result: FCST will output a default main file called "main.xml"
 *
 * BIN -p file_name:
 *  - FCST will assume that file_name is the name of a main file
 *  - FCST will read the main file and generate data parameter file "data.xml"
 *  - FCST may create and  optimization parameter files "opt.xml", but not required to
 *
 * <h3>Notes for Qt beginners</h3>
 * Qt online class documentation is very good! If you are wondering about an API just Google any class name.
 *
 * Most objects in Qt are handled as naked pointers. Object constructors take a pointer to a parent
 * object. The parent will manage their destruction. When the MainWindow is destroyed by the
 *
 * \note It is important that all developers learn and understand the Signal/Slots mechanism.
 *
 *
 * @ingroup ParameterGui
 * @author Phil Wardlaw, Simon Mattern, Marc Secanell, 2014-16
 */
class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    /**
     * Constructor.
     * If a @p filename is given,
     * the MainWindow tries to open
     * and parse the file.
     */
    MainWindow(const QString  &filename = "");

protected:
    /**
     * Reimplemented from QMainWindow.
     * We ask, if changes should be saved.
     */
    void closeEvent(QCloseEvent *event);

private slots:
    /**
     * Exports the current opened project into prm
     */
    void export_project_to_prm();
    /**
     * Opens a new window with all the settings
     */
    void openSettingsWindow();
    /**
     * Create a new window with spinbox
     */
    void createSpinBoxWindow();
    /**
     * Sets the amount of CPU cores for running the simulation
     */
    void setCoreNumber();
    /**
     * Calls function to open a project.
     */
    void openProjectSlot();
    /**
     * Save the parameter file.
     */
    bool saveAll();
    /**
     * Open a file dialog to save the parameter file.
     */
    bool save_as();
    /**
     * Show some information on the parameterGUI
     */
    void about();
    /**
     * Show some information on the OpenFCST
     */
    void aboutOpenFCST();
    /**
     * A <tt>slot</tt> that should be always called,
     * if parameter values were changed.
     */
    void tree_was_modified();
    /**
     * A <tt>slot</tt>  for the next button, progresses the GUI's state.
     */
    void next_click();
    /**
     * A <tt>slot</tt>  that should be called to take output from the current process's
     * std and err streams.
     */
    void printOutput();
    /**
     * A <tt>slot</tt>  that is called when FCST ends, updates GUI state.
     */
    void finishedFCST(int status ,QProcess::ExitStatus qStatus);
    /**
     * A <tt>slot</tt>  attached to a menu action, to start a new project.
     *
     * If no files are currently loaded then next_click() is called, otherwise
     * a new instance of FCSTGUI::MainWindow is called by starting a separate process.
     */
    void newProject();
    /**
     * A <tt>slot</tt>  that is called when fileWatcher reports new files in the  workingDir.
     */
    void updateFolderView(const QString & path);
    /**
     * A <tt>slot</tt>  to manage user interaction when a file name is double clicked in the file
     * browser.
     */
    void fileBrowserClicked( QListWidgetItem * item);
    /**
     * A <tt>slot</tt>  to manage user interaction when the menu "Save Log File" action is triggered.
     * Saves the OpenFCST simulation log stored in textOut class member.
     */
    void save_log();
    /**
     * A <tt>slot</tt>  that is called when the fileWatcher wants to close the open Projectfiles.
     * (tabIndex to close just a single tab not used yet, but already included for further use)
     */
    void closeProject();
    /**
     * Sets the dimension: 2D or 3D
     */
    void setDimension(bool isTriggered);
    /**
     * Starts a new instance of FCST
     */
    void newInstance();
    /**
     * This slot calls the conversion function from prm to xml 
     * by pressing the button
     */
    void convert_prm_xml_slot();
    /**
     * This slot calls the conversion function from xml to prm
     * by pressing the button
     */
    void convert_xml_prm_slot();

private:    
    /**
     * Starts converting user-selected file from .xml to .prm
     */
    void convert_xml_prm(QString FilePath);
    /**
     * Starts converting user-selected file from .prm to .xml
     */
    void convert_prm_xml(QString FilePath);
    /**
     * Opens a project
     */
    void openProject();
    /**
     * This funtions creates the toolbar
     */
    void create_toolbar();
    /**
     * This function creates settings.ini file an its contant.
     */
    void create_settings();  
    /**
     * This function creates all actions.
     */
    void create_actions();
    /**
     * This function creates all menus.
     */
    void create_menus();
    /**
     * This function checks, if parameters were changed
     * and show a dialog, if changes should be saved.
     * This function should be always called,
     * before open a new parameter file or before closing the GUI
     */
    bool maybe_save ();
    /**
     * Save parameters to @p filename in XML format.
     */
    bool save_file (const QString &filename, const QString &newfilename = "");
    /**
     * Load parameters from @p filename in XML format.
     */
    void load_file (const QString &filename, bool showMsg = true);
    /**
     * Function for setting up GUI, used by constructor
     */
    void create_widgets();
    /**
     * Function for creating tab given a file in the local working directory
     * Throws an std::runtime_error if file is not found.
     */
    void create_tab(QString name_);

    /**
     * Functions for calling and killing FCST, if pipeMsg is true then FCST output
     * will be displayed on screen. WorkDirection set the direction. 
     */
    bool callFCST(QStringList arguments, bool pipeMsg, QString WorkDirection, QString binaryAdd);
    bool callQuit; //Bool associated with the killFCST command
    void killFCST();

    /**
     * Function for getting working directory from user
     */
    bool setWorkingDir();
    /**
     * Function to check if a file exists
     */
    bool fileExists(QString path);
    /**
     * Logging function, appends to gui_log.txt
     */
    void log(const std::string &input);

    /**
     * ofstream object for logging
     */
    std::ofstream flog;
    /**
     * This is the trees  which we display all parameters.
     */
    std::map<QString ,QTreeWidget *> tree_widget;
    /**
     * The status of the individual trees, true if they do not need to be saved.
     */
    std::map<QString , bool> treeStatus;
    /**
     * This menu provides all file actions as <tt>open</tt>, <tt>save</tt>, <tt>save as</tt>
     * and <tt>exit</tt>
     */
    QPointer<QMenu> file_menu;
    /**
     * This menu provides some informations <tt>about</tt> the parameterGUI
     * and <tt>about Qt</tt>
     */
    QPointer<QMenu> help_menu;
    /**
     * Converting sub menu
     */
    QPointer<QMenu> convert_submenu;

    /**
     * QAction <tt>Start</tt> a new project.
     */
    QPointer<QAction> save_log_act;
    /**
     * QAction <tt>Start</tt> a new project.
     */
    QPointer<QAction> new_act;
    /**
     * QAction <tt>open</tt> a file.
     */
    QPointer<QAction> open_act;
    /**
     * QAction <tt>save</tt> a file.
     */
    QPointer<QAction> save_act;
    /**
     * QAction <tt>save as</tt> a file.
     */
    QPointer<QAction> save_as_act;
    /**
     * QAction <tt>save as prm</tt> a file.
     */
    QPointer<QAction> export_project_as_prm_act;
    /**
     * QAction <tt>exit</tt> the GUI.
     */
    QPointer<QAction> exit_act;
    /**
     * QAction <tt>about</tt> the parameterGUI.
     */
    QPointer<QAction> about_act;
    /**
     * QAction <tt>about</tt> the parameterGUI.
     */
    QPointer<QAction> aboutFCST_act;
    /**
     * QAction <tt>about</tt> Qt.
     */
    QPointer<QAction> about_qt_act;
    /**
     * This value stores the current <tt>filename</tt> we work on.
     */
    QStringList   current_files;
    /**
     * An object for storing user settings.
     */
    QPointer<QSettings> gui_settings;
    /**
     * Button for driving user events forward.
     */
    QPointer<QPushButton> nextButton;
    /**
     * QWidget elements, essentially panels, to house other elements of interest.
     */
    QPointer<QWidget> body , buttonPanel;
    /**
     * QTextEdit for displaying FCST output to user.
     */
    QPointer<QTextEdit> textOut;
    /**
     * Elements for displaying FCST logo in bottom left corner.
     */
    QPointer<QLabel> logoLabel;
    /**
     * Dynamic sizers for managing the size and position of visual elements.
     */
    QPointer<QHBoxLayout> hBox ;
    QPointer<QHBoxLayout> buttonBox;
    QPointer<QVBoxLayout> vBoxLeft;
    QPointer<QVBoxLayout> vBoxRight;
    QPointer<QTabWidget> tabWidget;
    /**
     * Enumeration describing states of project creation/simulation execution.
     */
    enum States {
        Start = 0,
        Main,   //preparing main
        Ready,  //ready to run
        Running //running
    };
    /**
     * The GUI's current state
     */
    States currentState;
    /**
     * Elements for displaying files output by FCST.
     */
    QPointer<QListWidget> outputList;
    QPointer<QFileSystemWatcher> fileWatcher;
    QFont  font;
    QPointer<QLabel>  fileLabel;
    QPointer<QLabel>FCSToutputLabel;
    /**
     * Object for calling FCST.
     */
    QProcess *process;
    /**
     * Strings for important locations.
     */
    QString workingDir;
    QString FCSTBinAddr;

    /**
     * OpenFCST binary arguments and file names produced by OpenFCST
     */
    QString paramArg;
    QString mainFileName, dataFileName, optFileName, currentProject;
    /**
     * Sets the numbers of CPU's used for the simulation
     */
    QString NumberOfCores;
    /**
     * Action to set fuel_chell-2d.bin
     * and run OpenFCST in 2D simulation
     */
    QPointer<QAction> set2Dsimulation_act;
    /**
     * Action to set fuel_chell-3d.bin
     * and run OpenFCST in 3D simulation
     */
    QPointer<QAction> set3Dsimulation_act;
    /**
     * This menu provides checkboxes about 2D, 3D or parallel running settings
     */
    QPointer<QMenu> settings_menu;
    /**
     * Stores the 2D binary path 
     */
    QString binary_path_2D;
    /**
     * Stores the 3D binary path
     */
    QString binary_path_3D;
    /**
     * Action for converting .prm to .xml file
     */
    QPointer<QAction> convert_prm_xml_act;
    /**
     * Action for CPU amount setting
     */
    QPointer<QAction> parallel_sim_act;
    /**
     * Action to close open project
     */
    QPointer<QAction> close_project_act;
    /**
     * Action to call settings window
     */
    QPointer<QAction> settings_act;
    /**
     * Action to call settings window
     */
    QPointer<QAction> convert_xml_prm_act;
    /**
     * Diablog window for setting the CPU's amount 
     */
    QPointer<QDialog> parallel_running_dialog;
    /**
     * SpinBox to set numbers of CPU's 
     */
    QPointer<QSpinBox> CoreNumber_SpinBox;
    /**
     * Function for creating QActions, returns QAction
     */
    QAction * createQAction(const QString &text,const char* member);
    
    
};
}

#endif
