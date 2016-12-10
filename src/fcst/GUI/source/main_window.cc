//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2006-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: main_window.cc
//    - Description: GUI for running FCST
//    - Developers: Phil Wardlaw, Simon Mattern, Marc Secanell, 2014-16
//
//---------------------------------------------------------------------------


#include <QtGui>

#include "main_window.h"

using namespace dealii::ParameterGui;

namespace FCSTGUI
{

//------------------------------------------------------------------------------------------------//
//----------------------------------Constructors and GUI Setup------------------------------------//
//------------------------------------------------------------------------------------------------//

MainWindow::MainWindow(const QString  &filename)
{
    log("Starting OpenFCST GUI \n");
    // Initialize a state
    currentState = Start;
    
    create_settings();
    create_widgets();
    create_actions();
    create_toolbar();
    create_menus();
    setWindowTitle(tr("[*] OpenFCST: Fuel Cell Simulation Toolbox "));
    this->statusBar()->show();
}

//------------------------------------------------------------------------------------------------//
void MainWindow::create_settings()
{
    log("Creating settings\n");
    
    // Load settings
    QString  settings_file = QDir::currentPath() + "/settings.ini";
    gui_settings = new QSettings (settings_file, QSettings::IniFormat);
    
    
    FCSTBinAddr = gui_settings->value("OpenFCSTbin").toString();
    if(FCSTBinAddr == "")
    {
        //Set default
        FCSTBinAddr = QDir::currentPath() + "/fuel_cell-2d.bin";
    }
    
    binary_path_2D = gui_settings->value("bin2DPath").toString();
    if(binary_path_2D == "")
    {
        //Set default
        binary_path_2D = QDir::currentPath() + "/fuel_cell-2d.bin";
    }
    
    binary_path_3D = gui_settings->value("bin3DPath").toString();
    if(binary_path_3D == "")
    {
        //Set default
        binary_path_3D = QDir::currentPath() + "/fuel_cell-3d.bin";
    }
    //Check if default or selected paths existed 
    if (fileExists(binary_path_2D) == false || fileExists(binary_path_3D) == false)
    {
        log("Getting OpenFCST binary from user... \n");
        QMessageBox::question(this, "OpenFCST Binary File",
                "OpenFCST binary not found, please select existing binary.",
                QMessageBox::Ok);
        //select 2D binary
        binary_path_2D = QFileDialog::getOpenFileName(this, tr("Select OpenFCST 2D simulation Binary File"),
                QDir::currentPath(),
                tr("2D Bin File(*2d.bin) (*2d.bin)"));
        
        //select 3D binary
        binary_path_3D = QFileDialog::getOpenFileName(this, tr("Select OpenFCST 3D simulation Binary File"),
                QDir::currentPath(),
                tr("3D Bin File(*3d.bin) (*3d.bin)"));
    }
    //Set working bin path
    gui_settings->setValue("OpenFCSTbin", binary_path_2D);
    
    //Load arguments for calling fcst
    paramArg = gui_settings->value("OpenFCSTparamArg").toString();
    if(paramArg == "")
    {
        //Set default
        paramArg = "-p";
        gui_settings->setValue("OpenFCSTparamArg", paramArg);
    }

    //Loaded expected file names to be produced by OpenFCST
    mainFileName = gui_settings->value("mainFileName").toString();
    if(mainFileName == "")
    {
        //Set default
        mainFileName = "main.xml";
        gui_settings->setValue("mainFileName", mainFileName);
    }

    dataFileName = gui_settings->value("dataFileName").toString();
    if(dataFileName == ""){
        //Set default
        dataFileName = "data.xml";
        gui_settings->setValue("dataFileName", dataFileName);
    }

    optFileName = gui_settings->value("optFileName").toString();
    if(optFileName == ""){
        //Set default
        optFileName = "opt.xml";
        gui_settings->setValue("optFileName", optFileName);
    }
    //sets the amount of the CPU's 
    CpuNumbers = gui_settings->value("CpuAmount").toString();
    if(CpuNumbers == ""){
        //Set default
        CpuNumbers = "1";
        gui_settings->setValue("CpuAmount", CpuNumbers);
    }
    
    workingDir = gui_settings->value("WorkDir").toString();
    if(workingDir == ""){
        //Set default
        workingDir = QDir::currentPath();
        gui_settings->setValue("WorkDir", workingDir);
    }
    //set and reset always the current project
    currentProject = gui_settings->value("CurrentProject").toString();
    currentProject = "";
    gui_settings->setValue("CurrentProject", currentProject);
    
    //set paths
    gui_settings->setValue("bin2DPath", binary_path_2D);
    gui_settings->setValue("bin3DPath", binary_path_3D);

}
//------------------------------------------------------------------------------------------------//

void MainWindow::create_widgets()
{
    log("Creating widgets\n");
    //Create body widget
    body = new QWidget(this);

    
    //Tab creation
    tabWidget = new QTabWidget(this);
    // set the tabs closeable
    tabWidget->setTabsClosable(false);

    //Button creation
    buttonPanel = new QWidget(this);
    nextButton = new QPushButton(tr("Start"));
    nextButton->setStyleSheet(QString::fromUtf8("min-height: 30px"));
    nextButton->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    //Connect nextButton click event to function
    connect(nextButton, SIGNAL(clicked(bool)), this, SLOT(next_click()));

    //Logo image creation
    logoLabel = new QLabel(this);
    logoLabel->setPixmap( QPixmap(":/images/logo_fcst_64.png"));
    //Set chamfered borders to logoLabel
    logoLabel->setStyleSheet(QString::fromUtf8("border-top: 2px solid qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1,"
            "stop:0 rgba(192, 192, 192, 255), stop:1 rgba(64, 64, 64, 255));"
            "border-left: 2px solid qlineargradient(spread:pad, x1:0, y1:0, x2:1, y2:0,"
            "stop:0 rgba(192, 192, 192, 255), stop:1 rgba(64, 64, 64, 255));"
            "border-right: 2px solid qlineargradient(spread:pad, x1:0, y1:0, x2:1, y2:0,"
            "stop:0 rgba(192, 192, 192, 255), stop:1 rgba(255, 255, 255, 255));"
            "border-bottom: 2px solid qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1,"
            "stop:0 rgba(192, 192, 192, 255), stop:1 rgba(255, 255, 255, 255));"));

    //Add items to horizontal sizer
    buttonBox = new QHBoxLayout(this);
    buttonBox->addWidget(logoLabel,2);
    buttonBox->addWidget(nextButton,1);
    buttonPanel->setLayout(buttonBox);

    //Text output creation
    FCSToutputLabel = new QLabel("Application Output:",this);
    FCSToutputLabel->setAlignment(Qt::AlignBottom);

    //+Font
    font = FCSToutputLabel->font();
    font.setPointSize(12);
    font.setBold(true);
    FCSToutputLabel->setFont(font);
    textOut = new QTextEdit(this);
    //Set textOut black background & grey text
    textOut->setStyleSheet(QString::fromUtf8("background-color: rgb(0, 0, 0);  color:  rgb(200, 200, 200);"));
    textOut->setReadOnly(true);

    //Button also uses same font
    nextButton->setFont(font);

    //File list creation
    fileLabel = new QLabel("Output Files:",this);
    fileLabel->setAlignment(Qt::AlignBottom);
    fileLabel->setFont(font);
    outputList = new QListWidget(this);
    outputList->setToolTip(tr("Double click to open file with operating system default. "
            "To change default program for a given file type use 'Open with' dialogue in standard file browser."));
    //Connect outputList double click event to function
    connect(outputList, SIGNAL(itemDoubleClicked (QListWidgetItem *)), this, SLOT(fileBrowserClicked(QListWidgetItem *)));


    //Set up dynamic sizers to manage alignment and proportions of all visual elements
    vBoxLeft  = new QVBoxLayout(this);
    vBoxRight = new QVBoxLayout(this);
    hBox      = new QHBoxLayout(this);

    vBoxLeft->addWidget(tabWidget, 8, 0);
    vBoxLeft->addWidget(buttonPanel, 1, 0);

    vBoxRight->addWidget(FCSToutputLabel, 1, 0);
    vBoxRight->addWidget(textOut, 15, 0);
    vBoxRight->addWidget(fileLabel, 1, 0);
    vBoxRight->addWidget(outputList, 7, 0);

    hBox->addLayout(vBoxLeft, 4);
    hBox->addLayout(vBoxRight, 3);

    body->setLayout(hBox);
    setCentralWidget(body);

}

//------------------------------------------------------------------------------------------------//
void MainWindow::create_actions()
{   
    log("Creating actions \n");
    new_act = createQAction(tr("&New Project..."), SLOT(newProject()));  // create actions
    new_act->setShortcut(Qt::CTRL + Qt::Key_N);                         // set a short cut
    new_act->setIcon(QIcon(":/images/new_icon.png"));
    new_act->setStatusTip(tr("Start a new project"));                   // set a status tip

    save_log_act = createQAction(tr("&Save Log..."), SLOT(save_log()));    // create actions
    save_log_act->setShortcut(Qt::CTRL + Qt::Key_L);                    // set a short cut
    save_log_act->setStatusTip(tr("Save the text log."));               // set a status tip

    open_act = createQAction(tr("&Open Project..."), SLOT(openProjectSlot()));   // create actions
    open_act->setShortcut(Qt::CTRL + Qt::Key_O);                        // set a short cut
    open_act->setIcon(QIcon(":/images/open_icon.png"));
    open_act->setStatusTip(tr("Open a XML file"));                      // set a status tip
    
    save_act = createQAction(tr("&Save All..."), SLOT(saveAll()));
    save_act->setShortcut(Qt::CTRL + Qt::SHIFT + Qt::Key_S );
    save_act->setIcon(QIcon(":/images/save_icon.png"));
    save_act->setStatusTip(tr("Save the current XML files."));
        
    save_as_act = createQAction(tr("&Save Current Tab As..."), SLOT(save_as()));
    save_as_act->setShortcut(Qt::CTRL + Qt::Key_S );
    save_as_act->setIcon(QIcon(":/images/save_as_icon.png"));
    save_as_act->setStatusTip(tr("Save the current XML file as"));

    exit_act = createQAction(tr("E&xit"),SLOT(close()));
    exit_act->setShortcut(Qt::CTRL + Qt::Key_Q);
    exit_act->setStatusTip(tr("Exit the parameterGUI application"));

    about_act = createQAction(tr("&About"), SLOT(about()));
    about_act->setStatusTip(tr("Show the parameterGUI About box"));

    aboutFCST_act = createQAction(tr("&About OpenFCST"), SLOT(aboutOpenFCST()));
    aboutFCST_act->setStatusTip(tr("Show the OpenFCST About box"));

    about_qt_act = new QAction(tr("About &Qt"), this);
    about_qt_act->setStatusTip(tr("Show the Qt library's About box"));
    connect(about_qt_act, SIGNAL(triggered()), qApp, SLOT(aboutQt()));
    
    set3Dsimulation_act = new QAction(tr("&3D Simulation"),this);
    set3Dsimulation_act->setCheckable(true); //set the Action checkable
    set3Dsimulation_act->setChecked(false);
    set3Dsimulation_act->setStatusTip(tr("Set the simulation to run in 3D"));    
    set3Dsimulation_act->setIcon(QIcon(":/images/3D_icon.png"));
    connect(set3Dsimulation_act, SIGNAL(triggered(bool)), this, SLOT(setDemensionenType(bool)));

    set2Dsimulation_act = new QAction(tr("&2D Simulation"),this);
    set2Dsimulation_act->setCheckable(true); //set the Action checkable
    set2Dsimulation_act->setChecked(true);
    set2Dsimulation_act->setStatusTip(tr("Set the simulation to run in 2D"));
    set2Dsimulation_act->setIcon(QIcon(":/images/2D_icon.png"));
    connect(set2Dsimulation_act, SIGNAL(triggered(bool)), this, SLOT(setDemensionenType(bool)));
     
    convert_prm_xml_act = createQAction(tr("Prm to Xml"),SLOT(convert_prm_xml_slot()));
    convert_prm_xml_act->setStatusTip(tr("Converts selected prm file into xml in same path"));
	
	convert_xml_prm_act = createQAction(tr("Xml to Prm"),SLOT(convert_xml_prm_slot()));
	convert_xml_prm_act->setStatusTip(tr("Converts selected xml file into prm in same path"));
    
    parallel_sim_act = createQAction(tr("Parallel Running"),SLOT(createSpinBoxWindow()));
    parallel_sim_act->setStatusTip(tr("Sets the amount of CPU's which are used for the simulation"));

    close_project_act= createQAction(tr("&Close Project"),SLOT(closeProject()));
    close_project_act->setIcon(QIcon(":/images/close_icon.png"));
    close_project_act->setStatusTip(tr("Close the open project"));
    
    settings_act = createQAction(tr("&Advanced Settings"), SLOT(openSettingsWindow()));
    settings_act->setIcon(QIcon(":/images/settings_icon.png"));
    settings_act->setStatusTip(tr("Open the advance settings window"));
    
    export_project_as_prm_act = createQAction(tr("&Export to prm"),SLOT(export_project_to_prm()));
    export_project_as_prm_act->setStatusTip(tr("Export the current open project in PRM files"));
}
//------------------------------------------------------------------------------------------------//
void MainWindow::create_toolbar()
{
    //create maintoolbar for icons
    QToolBar * mainToolBar = new QToolBar("Main Tool Bar");
    //add toolbar to the mainwindow
    this->addToolBar(mainToolBar);
    mainToolBar->setMovable(false);
    mainToolBar->addAction(new_act);
    mainToolBar->addAction(open_act);
    mainToolBar->addAction(save_act);
    mainToolBar->addAction(close_project_act);
    mainToolBar->addAction(settings_act);
    mainToolBar->addSeparator();
    mainToolBar->addAction(set2Dsimulation_act);
    mainToolBar->addAction(set3Dsimulation_act);
}
//------------------------------------------------------------------------------------------------//
void MainWindow::create_menus()
{
    log("Creating menus\n");
    file_menu = menuBar()->addMenu(tr("&File"));    // create a file menu
    file_menu->addAction(new_act);                  // and add actions
    file_menu->addAction(open_act);
    file_menu->addAction(save_act);
    file_menu->addAction(save_as_act);
    file_menu->addAction(save_log_act);
    file_menu->addAction(export_project_as_prm_act);
    file_menu->addAction(close_project_act);
    file_menu->addAction(exit_act);

    menuBar()->addSeparator();
    
    settings_menu = menuBar()->addMenu(tr("&Options")); //create a settingsmenu
    settings_menu->addAction(set2Dsimulation_act);
    settings_menu->addAction(set3Dsimulation_act);
    settings_menu->addAction(parallel_sim_act);
    convert_submenu = settings_menu->addMenu(tr("&Convert project..."));
    convert_submenu->addAction(convert_prm_xml_act);
    convert_submenu->addAction(convert_xml_prm_act);
    settings_menu->addAction(settings_act);
    
    menuBar()->addSeparator();

    help_menu = menuBar()->addMenu(tr("&Help"));    // create a help menu
    help_menu->addAction(about_act);
    help_menu->addAction(aboutFCST_act);
    help_menu->addAction(about_qt_act);
}


//------------------------------------------------------------------------------------------------//
void  MainWindow::create_tab(QString name_)
{
    QDir qdir;

    //name_ is just the file name, get full path
    QString absName = workingDir + "/" +  name_;
    log("Loading file to tab '" + absName.toStdString() + "'\n");

    //Test if the file exist, throw an exception
    if (qdir.exists(absName))
    {
        //Tree for showing XML tags, map stores tree widget for each tab
        // Setup a new tree
        tree_widget[absName] = new QTreeWidget;
        tree_widget[absName]->header()->setResizeMode(QHeaderView::ResizeToContents);

        //Setup titles
        tree_widget[absName]->setHeaderLabels(QStringList() << tr("(Sub)Sections/Parameters")
                << tr("Value"));

        // enables mouse events e.g. showing ToolTips
        tree_widget[absName]->setMouseTracking(true);

        // Set which actions will initiate item editing: Editing starts when:
        // DoubleClicked: an item is double clicked
        // EditKeyPressed: the platform edit key has been pressed over an item
        tree_widget[absName]->setEditTriggers(QAbstractItemView::DoubleClicked|
                QAbstractItemView::EditKeyPressed);

        // set the delegate for editing items
        tree_widget[absName]->setItemDelegate(new ParameterDelegate(1));

        //Load the file at this absolute path
        //File will be loaded to tree_widget[absName]
        load_file(absName, false);

        //Go to our new loaded file, name it after name_ (shorthand)
        tabWidget->setCurrentIndex(tabWidget->addTab(tree_widget[absName], name_));

        // This tree element hasn't changed and doesn't be saved.
        treeStatus[name_] = true;

        //connect slot: if the tree changes, the window will know,
        //it is important that this is done AFTER the file has been loaded
        connect(tree_widget[absName], SIGNAL(itemChanged(QTreeWidgetItem*, int)), this, SLOT(tree_was_modified()));

        //Add to list the current files we have loaded
        current_files << absName;
    }
    else
        throw std::runtime_error("Error loading file " +name_.toStdString() + ". File not found!");
}

//------------------------------------------------------------------------------------------------//
//--------------------------------Protected and Private Functions---------------------------------//
//------------------------------------------------------------------------------------------------//
void MainWindow::convert_xml_prm(QString FilePath)
{
    if(fileExists(FilePath))
    {
        //get the absolutePath of the FilePath to opperate in the same folder as the prm 
        QFileInfo fi;
        fi.setFile(FilePath);
        QString workDir =fi.absolutePath();

        QString binPath;
        binPath=gui_settings->value("OpenFCSTbin").toString();
        callFCST(QStringList() << "-r" << FilePath, true, workDir,binPath); // ./fuel-cell.bin -c

        //sleeps while fcst converts file
        usleep(5000);
    }
    else return;
}
//------------------------------------------------------------------------------------------------//
void MainWindow::convert_prm_xml(QString FilePath)
{
    //get the absolutePath of the FilePath to opperate in the same folder as the prm 
    if(fileExists(FilePath))
    {
        /// @todo check if user selected the correct file (main.prm) and throw a warning if data.prm was selected
        QFileInfo fi;
        fi.setFile(FilePath);
        QString workDir =fi.absolutePath();
        QString binPath;
        //Ask user if converting to 2D or 3D
        QMessageBox msgBox;
        msgBox.setText(tr("Choose convertion dimension"));
        QPushButton *convert2D = msgBox.addButton(tr("2D"), QMessageBox::ActionRole);
        QPushButton *convert3D = msgBox.addButton(tr("3D"), QMessageBox::ActionRole);
        msgBox.exec();
        
        if(msgBox.clickedButton()==convert2D)
        {
            binPath=gui_settings->value("bin2DPath").toString();
            callFCST(QStringList() << "-c" << FilePath, true, workDir,binPath); // ./fuel-cell.bin -c
        }
        else if (msgBox.clickedButton()==convert3D)
        {
            binPath=gui_settings->value("bin3DPath").toString();
            callFCST(QStringList() << "-c" << FilePath, true, workDir,binPath); // ./fuel-cell.bin -c
        }
        //sleeps while fcst converts file
        usleep(5000);
        //Ask user for open a project
        QMessageBox::StandardButton reply;
        reply= QMessageBox::question(this, "Open Project","Would you like to open a project?",
                    QMessageBox::Yes | QMessageBox::No);
        if (reply == QMessageBox::Yes)
        {
            workingDir = workDir;
            gui_settings->setValue("WorkDir", binary_path_2D);
            openProject();
        }
        else
        {
            return;
        }
    }
    else return;
    
}

//------------------------------------------------------------------------------------------------//
bool MainWindow::fileExists(QString path)
{
    QFileInfo checkFile(path);
    if(checkFile.exists() && checkFile.isFile())
    {
        return true;
    }
    else 
    {
        return false;
    }
}
//------------------------------------------------------------------------------------------------//
bool MainWindow::callFCST(QStringList arguments, bool pipeMsg, QString WorkDirection,QString binaryAdd)
{
    callQuit = false;
    log("Starting OpenFCST\n");
    process = new QProcess(this);
    //If we desire connect the output of this process to the GUI output screen
    if(pipeMsg){
        connect(process, SIGNAL(readyReadStandardError()), this, SLOT(printOutput()));
        connect(process, SIGNAL(readyReadStandardOutput()), this, SLOT(printOutput()));
    }
    //Connect the finished single to slot function so we can note when OpenFCST finishes
    connect(process, SIGNAL(finished(int,QProcess::ExitStatus)), this, SLOT(finishedFCST(int,QProcess::ExitStatus)));
    //Run fcst from WorkDirection
    process->setWorkingDirectory(WorkDirection);
    //get the amount of CPU's
    QString CpuNum = gui_settings->value("CpuAmount").toString();
    
    //     add mpirun argument to run
    //     @caution "-quiet" is called because of an error during the conversion from PRM into XML by calling mpirun
    //     source of the error not found yet: error output: 
    //      ------------------------------------------------------------------------
    //      mpirun has exited due to process rank 0 with PID 4694 on
    //      n ode seca*nell-grad15 exiting improperly. There are three reasons this could occur:
    //      
    //      1. this process did not call "init" before exiting, but others in
    //      the job did. This can cause a job to hang indefinitely while it waits
    //      for all processes to call "init". By rule, if one process calls "init",
    //          then ALL processes must call "init" prior to termination.
    //          
    //      2. this process called "init", but exited without calling "finalize".
    //      By rule, all processes that call "init" MUST call "finalize" prior to
    //      exiting or it will be considered an "abnormal termination"
    //          
    //      3. this process called "MPI_Abort" or "orte_abort" and the mca parameter
    //      orte_create_session_dirs is set to false. In this case, the run-time cannot
    //      detect that the abort call was an abnormal termination. Hence, the only
    //      error message you will receive is this one.
    //          
    //      This may have caused other processes in the application to be
    //      terminated by signals sent by mpirun (as reported here).
    //          
    //      You can avoid this message by specifying -quiet on the mpirun command line.
    QString runAdd = "mpirun -quiet -np "+CpuNum+" "+ binaryAdd +" ";
    for (QStringList::iterator it = arguments.begin(); it != arguments.end(); ++it)
    {
        QString current = *it;
        runAdd.push_back(current);
    }
    process->start(runAdd,QIODevice::ReadWrite | QIODevice::Text);
    return true;
}
//------------------------------------------------------------------------------------------------//
void MainWindow::killFCST()
{
    log("Call to terminate OpenFCST \n");
    callQuit = true;
    if(process != NULL)
        process->terminate();
}

//------------------------------------------------------------------------------------------------//
bool MainWindow::setWorkingDir()
{
    bool answer = false;
    log("Setting working directory \n");
    //Ask the user for a workign directory
    workingDir =  QFileDialog::getExistingDirectory(this, tr("Select directory to launch OpenFCST from."),
            "~", QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);

    if(not workingDir.isEmpty()){
        log("Working directory set to:" + workingDir.toStdString() +"\n");

        QDir qdir(workingDir);
        QStringList files = qdir.entryList(QDir::Files);


        //Check that no project files already exist in the folder

        if((files.indexOf(mainFileName) == -1) and
           (files.indexOf(dataFileName)  == -1)  and
           (files.indexOf(optFileName) == -1)){
            answer = true; //no files in folder
        }
        else{
            log("XML files found in working directory \n");
            //Give the user the option to overwrite or not
            answer = (QMessageBox::warning(this, "Existing project", "Some project files already exist in " +
                    workingDir +". Overwrite?", QMessageBox::Yes|QMessageBox::No) == QMessageBox::Yes);
        }

    }
    //Set new work dir
    gui_settings->setValue("WorkDir", workingDir);
    return answer;
}

//------------------------------------------------------------------------------------------------//
void MainWindow::closeEvent(QCloseEvent *event)
{
    if (maybe_save())								// reimplement the closeEvent from the QMainWindow class
        event->accept();							// check, if we have to save modified content,
    else									        // if content was saved, accept the event,
        event->ignore();							// otherwise ignore it
}
//------------------------------------------------------------------------------------------------//
bool MainWindow::maybe_save()
{
    log("Asking user if they want to save project.\n");
    if (isWindowModified())							// if content was modified
    {
        //Get unsaved file names to display to user
        QString unsavedFiles;

        for (std::map<QString, bool>::iterator x = treeStatus.begin(); x != treeStatus.end(); x++)
            if(x->second == false)
                unsavedFiles += x->first + " ";


        QString msg =  "The contents of files (" + unsavedFiles + ") has been modified.\n" +
                "Do you wish to close your files without saving your changes?";

        //Ask if content should be saved
        QMessageBox::StandardButton ret;
        ret = QMessageBox::warning(this, tr("OpenFCST parameterGUI"), msg,
                QMessageBox::Yes |QMessageBox::No);

        if (ret == QMessageBox::Yes)
            return true;
        else if (ret == QMessageBox::No)
            return false;
    };
    return true;
}
//------------------------------------------------------------------------------------------------//
void MainWindow::load_file(const QString &filename, bool showMsg)
{
    log("Loading file \n");
    QFile  file(filename);

    if (!file.open(QFile::ReadOnly | QFile::Text))              // open the file
    {
        log("Cannot read file '"+ filename.toStdString() + "'\n");
        QMessageBox::warning(this, tr("OpenFCST parameterGUI"),
                tr("Cannot read file %1:\n%2.")
                .arg(filename)
                .arg(file.errorString()));
        log(file.errorString().toStdString() + "\n");
        return;
    };

    tree_widget[filename]->clear();                             // clear the tree

    XMLParameterReader  xml_reader(tree_widget[filename]);      // and read the xml file

    if (!xml_reader.read_xml_file(&file))
    {
        log("Parse error in file '"+ filename.toStdString() + "'\n");
        QMessageBox::warning(this, tr("OpenFCST parameterGUI"),
                tr("Parse error in file %1:\n\n%2")
                .arg(filename)
                .arg(xml_reader.error_string()));
        log(xml_reader.error_string().toStdString() + "\n");
    }

}

//------------------------------------------------------------------------------------------------//
bool MainWindow::save_file(const QString &filename, const QString &newfilename)
{
    QString savePath;

    if(newfilename == "")
        savePath = filename;
    else
        savePath = newfilename;


    QFile  file(savePath);

    if (!file.open(QFile::WriteOnly | QFile::Text))               // open a file dialog
    {
        log("Cannot write file '"+ filename.toStdString() + "'\n");
        QMessageBox::warning(this, tr("OpenFCST parameterGUI"),
                tr("Cannot write file %1:\n%2.")
                .arg(savePath)
                .arg(file.errorString()));
        log(file.errorString().toStdString() + "\n");
        return false;
    };

    XMLParameterWriter  xml_writer(tree_widget[filename]);        // create a xml_writer

    if (!xml_writer.write_xml_file(&file))                        // and read the xml file
        return false;

    statusBar()->showMessage(tr("File saved"), 2000);             // if we succeed, show a message

    treeStatus[QFileInfo(filename).fileName()] = true;
    return true;
}

//------------------------------------------------------------------------------------------------//
void MainWindow::openProject()
{
    
    log("Opening project \n");
    // check, if the content was modified, ask the user if they want to save changes
    if (maybe_save())
    {
        
        //Content was not modified, but files are open, ask user if they wish to close these files
        if(not isWindowModified() and current_files.size() != 0){
            
            QString msg = "Close existing project (no changes to save)?";
            
            if(QMessageBox::question(this, "Close Project?",
                msg, QMessageBox::Yes|QMessageBox::No)== QMessageBox::No)
            {
                return;
            }
            
        }
        //Clear and close current files
        current_files.clear();
        tree_widget.clear();
        tabWidget->clear();
        
        QStringList tempList;
        QStringList additionalFiles;
        //Get main file
        QString answer = QFileDialog::getOpenFileName(this, tr("Open Main Parameter File"),
                                                        workingDir,
                                                        tr("Main XML File (*main*.xml) / Main PRM File (*main*.prm) (*main*.xml *main*.prm);;"
                                                        "Main XML File (*.xml) (*.xml);;"
                                                        "(Main PRM File (*main*.prm) (*main*.prm)"));
        
        
        //Check if main file is .prm
        if(answer.endsWith(".prm", Qt::CaseSensitive))
        {
            //ask user for conversion
            QString msg = "You select a .prm file.\n\n"
                          "Would you like to convert the .prm files into .xml in the same directory?";
            QMessageBox::StandardButton reply;
            reply= QMessageBox::question(this, "Convert Project",msg,QMessageBox::Yes | QMessageBox::No);
            if (reply == QMessageBox::Yes)
            {
                convert_prm_xml(answer);
            }
            else
            {
                return;
            }
        }
        else
        {      
            if(answer != "")
            {
                tempList << answer;
                //set current opened main file
                gui_settings->setValue("CurrentProject", answer);
            }
            else
                return; //They did not select a necessary file, return
                
            //Get directory
            workingDir = QFileInfo(tempList.at(0)).path();
            
            //Get Additional files
            answer = QFileDialog::getOpenFileName(this, tr("Select Data Parameter File"),
                                                    workingDir, tr("XML Files (*.xml)"));
            
            if(answer != "")
                additionalFiles << answer;
            else
            {
                return; //They did not select necessary a file, return
                //Reset current project if user didn't select data file
                QString cp = "";
                gui_settings->setValue("CurrentProject", cp);
            }

            //Main and Data now loaded, extra files are optional
            //Add as many more files as they desire
            QString msg = "Would you like to load additional parameter files?";
            while(QMessageBox::question(this, "Load more files...", msg, QMessageBox::Yes|QMessageBox::No) == QMessageBox::Yes){
                QString answer = QFileDialog::getOpenFileName(this, tr("Select Additional Parameter File"),
                                                            workingDir, tr("XML Files (*.xml)"));
                
                if(answer != "")
                    additionalFiles << answer;
            }
            //Check that the additional files share the same path, else reject
            for(unsigned int i = 0; i < additionalFiles.size(); i++){
                if (QFileInfo(additionalFiles.at(i)).path() != workingDir){
                    //File will be skipped, ask the user if they would like to select an alternative
                    //The alternative will be checked at a later iteration of this for loop
                    
                    QString msg = "Parameter file " + QFileInfo(additionalFiles.at(i)).fileName() + " must reside in the same folder as the main parameter file (" +
                    workingDir +  "). Would you like to select an alternative file?";
                    
                    
                    if (QMessageBox::warning(this, tr("Error"), msg, QMessageBox::Yes |QMessageBox::No) == QMessageBox::Yes){
                        QString answer = QFileDialog::getOpenFileName(this, tr("Select additional Parameter File"),
                                                                    workingDir,
                                                                    tr("XML Files (*.xml)"));
                        if(answer != "")
                            additionalFiles << answer;
                    }
                    
                }
                else {
                    tempList << additionalFiles.at(i);
                }
            }
            //Potential bug: the user may have supplied a data file from an invalid location
            //which was not replaced correctly.
            
            //Files loaded and checked, display files in tabs
            for(unsigned int i = 0; i < tempList.size(); i++)
            {
                if(tempList.at(i) !=""){
                    try{
                        create_tab(QFileInfo(tempList.at(i)).fileName());
                    }
                    catch(std::runtime_error &e){
                        log("Error opening project files. "+ std::string(e.what()) + "\n");
                        QMessageBox::warning(this, tr("Error"), tr(e.what()), QMessageBox::Ok );
                    }
                }
            }
            setWindowModified(false);
            //Necessary files loaded, allow user to run files.
            nextButton->setText("Run");
            currentState = Ready;
            //Set new work dir
            gui_settings->setValue("WorkDir", workingDir);
        }
    };
}
//-----------------------------------------------------------------------------
QAction* MainWindow::createQAction(const QString& text, const char* member)
{
    QAction * Act = new QAction(text,this);
    connect(Act,SIGNAL(triggered()),this,member);
    return Act;
}
//-----------------------------------------------------------------------------
void  MainWindow::log(const std::string &input)
{
    flog.open("gui_log.txt", std::ios_base::app);
    std::time_t result = std::time(nullptr);
    flog <<  std::asctime(std::localtime(&result)) << "   " << input;
    flog.close();
}
//------------------------------------------------------------------------------------------------//
//-------------------------------------------SLOTS------------------------------------------------//
//------------------------------------------------------------------------------------------------//
void MainWindow::export_project_to_prm()
{
    //get current path of the project
    QString currentProjectPath = gui_settings->value("CurrentProject").toString();
    if(currentProjectPath != "" )
    {
        QMessageBox::StandardButton reply;
        reply= QMessageBox::question(this, "Export project",
                                     tr("<p align='center'> Would you like to export "
                                            "the current project to prm?<br><br>"
                                            "The original .prm project files will be overwritten!"),
                                        QMessageBox::Yes | QMessageBox::No);
        if (reply == QMessageBox::Yes)
        {
            //save actuall project 
            saveAll();
            //call function to save file 
            convert_xml_prm(currentProjectPath);
        }
        else
        {
            return;
        }
    }
    else
    {
        return ;
    }
    
}

//------------------------------------------------------------------------------------------------//
void MainWindow::openSettingsWindow()
{
  settingsWindow * sw =  new settingsWindow(gui_settings);
}
//------------------------------------------------------------------------------------------------//
void MainWindow::createSpinBoxWindow()
{
    QVBoxLayout * layo = new QVBoxLayout();
    QPushButton * setButton = new QPushButton(tr("Ok"));
    QLabel * infoLabel = new QLabel(tr("Set amount usage of CPU's \n for running the simulation"));
    //connectiopn to set user value into settings.ini file
    connect(setButton, SIGNAL(pressed()),this, SLOT(setCpuAmount()));
    CpuAmound_SpinBox = new QSpinBox();
    //range only between 1-4 cores
    CpuAmound_SpinBox->setRange(1,4);
    int value = gui_settings->value("CpuAmount").toInt();
    CpuAmound_SpinBox->setValue(value);
    layo->addWidget(infoLabel);
    layo->addWidget(CpuAmound_SpinBox);
    layo->addWidget(setButton);
    
    parallel_running_dialog = new QDialog();
    parallel_running_dialog->setLayout(layo);
    parallel_running_dialog->setAttribute(Qt::WA_DeleteOnClose);
    parallel_running_dialog->exec();
}
//------------------------------------------------------------------------------------------------//

void MainWindow::setCpuAmount()
{
    int i_value = CpuAmound_SpinBox->value();
    QString s_val = QString::number(i_value);
    gui_settings->setValue("CpuAmount", s_val);
    parallel_running_dialog->close();
}
//------------------------------------------------------------------------------------------------//
void MainWindow::convert_xml_prm_slot()
{
    //Warning about the overwrite 
    QMessageBox::warning(this, 
                         tr("Warning"),
                         tr("<p align='center'> If you choose to convert a "
                                               "project from .xml into .prm, "
                                               "the original project files (*.prm) "
                                               "in the selected working folder will be overwritten!"
                                               "<br><br>Using current selected binary for conversion"), 
                         QMessageBox::Ok );

    //get converting file path from user
    QString FilePath = QFileDialog::getOpenFileName(this, tr("Select file to convert"),
                QDir::currentPath(),
                tr("main xml file(*main*.xml) (*main*.xml)"));
    convert_xml_prm(FilePath);
}
//------------------------------------------------------------------------------------------------//

void MainWindow::convert_prm_xml_slot()
{
    //get converting file path from user
    QString FilePath = QFileDialog::getOpenFileName(this, tr("Select main file to convert"),
                workingDir,
                tr("main prm file(*main*.prm) (*main*.prm)"));
    convert_prm_xml(FilePath);
}
//------------------------------------------------------------------------------------------------//

void MainWindow::newInstance()
{
    //Statring a new Instance
    log("Starting a new process: " + QCoreApplication::applicationFilePath().toStdString() +"\n");
    QProcess *newProcess = new QProcess(this);
    newProcess->start(QCoreApplication::applicationFilePath());
}
//------------------------------------------------------------------------------------------------//
void MainWindow::openProjectSlot()
{
    openProject();
}
//------------------------------------------------------------------------------------------------//
bool MainWindow::saveAll()
{
    log("Saving all xml files \n");
    for (int i = 0; i < current_files.size(); i++){
        save_file(current_files.at(i));
    }

    setWindowModified(false);
    return true;

}

//------------------------------------------------------------------------------------------------//
bool MainWindow::save_as()
{
    log("Saving individual file \n");
    if(current_files.size() == 0){
        statusBar()->showMessage(tr("Save as: Nothing to save"), 5000);
        return false; //Nothing to save
    }
    
    QString currentFilePath = workingDir + "/" + tabWidget->tabText(tabWidget->currentIndex());
    // open a file dialog
    QString  file_name = QFileDialog::getSaveFileName(this, tr("Save XML Parameter File"),
                    currentFilePath,
                    tr("XML Files (*.xml)"));

    if (file_name.isEmpty())                               // if no file was selected
        return false;                                      // return false
    else
        return save_file(currentFilePath, file_name);      // otherwise save content to file
}

//------------------------------------------------------------------------------------------------//
void MainWindow::tree_was_modified()
{
    log("Modified: " + tabWidget->tabText(tabWidget->currentIndex()).toStdString() + "\n");
    treeStatus[tabWidget->tabText(tabWidget->currentIndex())] = false; //Note which tab needs to be saved.
    setWindowModified(true);                                           // store, that the window was modified
    // this is a function from the QMainWindow class
    // and we use the windowModified mechanism to show a "*"
    // in the window title, if content was modified
}

//------------------------------------------------------------------------------------------------//
void MainWindow::about()
{
#ifdef Q_WS_MAC
    static QPointer<QMessageBox> old_msg_box;

    if (old_msg_box)
    {
        old_msg_box->show();
        old_msg_box->raise();
        old_msg_box->activateWindow();
        return;
    };
#endif

    QString title = "About OpenFCST parameterGUI";

    QString trAboutparameterGUIcaption;
    trAboutparameterGUIcaption = QMessageBox::tr(
            "<h3>FCST parameterGUI version 1.0</h3>"
            "<p>This program uses Qt version %1.</p>"
    ).arg(QLatin1String(QT_VERSION_STR));

    QString trAboutparameterGUItext;
    trAboutparameterGUItext = QMessageBox::tr(
            "<p>The <a href=\"http://openfcst.mece.ualberta.ca\">OpenFCST</a> parameterGUI is a graphical user interface for editing XML parameter files "
            "adapted from the deal.II parameterGUI(<a href=\"http://www.dealii.org/7.0.0/doxygen/deal.II/classParameterHandler.html\">dealii.org/doc</a>). </p>"
            "<p> Developed by Philip Wardlaw (<a href=\"https://www.linkedin.com/pub/phil-wardlaw/94/133/b09\">View  profile</a>) </p>"
    );

    QPointer<QMessageBox> msg_box = new QMessageBox;
    msg_box->setAttribute(Qt::WA_DeleteOnClose);
    msg_box->setWindowTitle(title);
    msg_box->setText(trAboutparameterGUIcaption);
    msg_box->setInformativeText(trAboutparameterGUItext);

#ifdef Q_WS_MAC
    old_msg_box = msg_box;
    msg_box->show();
#else
    msg_box->exec();
#endif
}
//------------------------------------------------------------------------------------------------//
void MainWindow::aboutOpenFCST()
{
#ifdef Q_WS_MAC
    static QPointer<QMessageBox> old_msg_box;

    if (old_msg_box)
    {
        old_msg_box->show();
        old_msg_box->raise();
        old_msg_box->activateWindow();
        return;
    };
#endif

    QString title = "About OpenFCST";

    QString trAboutparameterGUIcaption;
    trAboutparameterGUIcaption = QMessageBox::tr(
            "<h3>Authors:</h3>"
           "<p><a href=\"mailto:secanell@ualberta.ca\">Marc Secanell</a>, Valentin N. Zingan, Andreas Putz,"
           " Madhur Bhaiya, Michael Moore, Peter Dobson,"
           " Philip Wardlaw, Chad Balen,"
           " Jie Zhou, Kailyn Domican, Aslan Kosakian and Mayank Sabharwal.</p>"
    ).arg(QLatin1String(QT_VERSION_STR));


    QString trAboutparameterGUItext;
    trAboutparameterGUItext = QMessageBox::tr(
            "<h3>About OpenFCST:</h3>" 
            " <p><a href=\"http://openfcst.mece.ualberta.ca\">OpenFCST</a> (open-source fuel cell simulation toolbox) is an open-"
            " source, finite element method based, multi-dimensional mathematical modeling software for polymer electrolyte fuel cells.<\p>"
            " <p>The aim of the software is to develop a platform for collaborative development of fuel cell mathematical models. Detailed information"
            " regarding the latest updates on the software can be found at the <a href=\"http://openfcst.mece.ualberta.ca\">OpenFCST</a> website "
            " including the User Guide for the software. </p>"
            " <p>OpenFCST currently includes physical models for gas, electron, ion, ionomer-"
            " bound water and heat transport. It also contains effective transport media relations to estimate transport properties for gas diffusion"
            " layers, micro-porous layers and catalyst layers as well as several kinetic"
            " models for the fuel cell electrochemical reactions. OpenFCST has been structured as a toolbox such that it is easier for new users to"
            " integrate new physical models with existing framework.</p>"
            ""
            " <h3>Copyright:</h3>"
            " <p> FCST: Fuel Cell Simulation Toolbox is distributed under the MIT License.</p>"
            " <p> Copyright (C) 2013 Energy Systems Design Laboratory, University of Alberta. </p>"
            " <p> The MIT License (MIT)</p>"
            " <p>Permission is hereby granted, free of charge, to any person obtaining a copy of this software "
            "and associated documentation files (the \"Software\"), to deal in the Software without restriction, "
            "including without limitation the rights to use, copy, modify, merge, publish, distribute,"
            "sublicense, and/or sell copies of the Software, and to permit persons to whom the Software "
            "is furnished to do so, subject to the following conditions:"
            "The above copyright notice and this permission notice shall be included in all "
            "copies or substantial portions of the Software.</p>"
            " <p>THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, "
            "INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR "
            "PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE "
            "FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, "
            "ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.");

    QPointer<QMessageBox> msg_box = new QMessageBox;
    msg_box->setAttribute(Qt::WA_DeleteOnClose);
    msg_box->setWindowTitle(title);
    msg_box->setText(trAboutparameterGUIcaption);
    msg_box->setInformativeText(trAboutparameterGUItext);

#ifdef Q_WS_MAC
    old_msg_box = msg_box;
    msg_box->show();
#else
    msg_box->exec();
#endif
}
//------------------------------------------------------------------------------------------------//
void MainWindow::printOutput()
{
    //Output error and normal streams, as OpenFCST may pipe out on both
    QByteArray byteArrayE = process->readAllStandardError();
    QByteArray byteArrayS = process->readAllStandardOutput();
    textOut->moveCursor(QTextCursor::End);
    textOut->insertPlainText(byteArrayE);
    textOut->moveCursor(QTextCursor::End);
    textOut->insertPlainText(byteArrayS);
    textOut->moveCursor(QTextCursor::End);
}

//------------------------------------------------------------------------------------------------//
void MainWindow::finishedFCST(int status ,QProcess::ExitStatus qStatus)
{
    log("OpenFCST exited with code " + std::to_string(status) + "\n");
    if(currentState == Running){
        //FCST was running a simulation, default back to ready state
        currentState = Ready;
        nextButton->setText("Run");
    }

    if (qStatus == QProcess::CrashExit and not callQuit){
        log("OpenFCST Binary crashed \n");
        QMessageBox::information(this, "OpenFCST Crashed",
                       "OpenFCST appears to have crashed. For more information regarding the crash scenario you may wish to run "
                       "OpenFCST in debug mode - see the installation script for options.",
                       QMessageBox::Ok);
    }
}

//------------------------------------------------------------------------------------------------//
void MainWindow::updateFolderView(const QString & path)
{
    //Gets the list of files in path
    QDir qdir(path);
    QStringList files = qdir.entryList(QDir::Files);

    for(unsigned int i = 0; i < files.size(); i++){
        //Adds the files to the file list if they are not xml and aren't duplicates
        if(not files.at(i).endsWith(".xml")){
            if(outputList->findItems(files.at(i), Qt::MatchFixedString).size() == 0)
                outputList->addItem(new QListWidgetItem(files.at(i)));
        }
    }

}

//------------------------------------------------------------------------------------------------//
void MainWindow::save_log()
{
    //Saves the log produced by running OpenFCST simulation
    //Open a file dialog
    QString  file_name = QFileDialog::getSaveFileName(this, tr("Save Log File"),
                    workingDir,
                    tr("Log Files (*.log)"));

    //Write file
    QFile outfile;
    outfile.setFileName(file_name);
    outfile.open(QIODevice::Append);
    QTextStream out(&outfile);
    out << textOut->toPlainText() << endl;
}

//------------------------------------------------------------------------------------------------//
void MainWindow::fileBrowserClicked(QListWidgetItem * item)
{
    //Try open the file using operating system
    QString name = item->text();
    log("Call OS to open " + name.toStdString() + "\n");
    QDesktopServices::openUrl(QUrl(workingDir + "/" + name));
}
//------------------------------------------------------------------------------------------------//
void MainWindow::closeProject()
{
    // check, if the content was modified, ask the user if they want to save changes
    if(maybe_save())
    {
        //Content was not modified, but files are open, ask user if they wish to close these files
        if(not isWindowModified() and current_files.size() != 0)
        {
            QString msg = "Close existing project (no changes to save)?";
           
            if(QMessageBox::question(this, "Close Project?",
                    msg, QMessageBox::Yes|QMessageBox::No)== QMessageBox::No)
            {
                return;
            }
        }
        //Clear everything and close current files
        currentState = Start;
        current_files.clear();
        tree_widget.clear();
        tabWidget->clear();
        setWindowModified(false);
        nextButton->setText("Start");
        textOut->clear();
        outputList->clear();
        gui_settings->setValue("CurrentProject", "");
    }    
}
//------------------------------------------------------------------------------------------------//
void MainWindow::setDemensionenType(bool isTriggered)
{
    QObject * obj = sender();
    if (obj == set2Dsimulation_act)
    {
        //if user clicks already selected button again
        if(set2Dsimulation_act->isChecked()==false && isTriggered==false)
        {
            //let the button selected
            set2Dsimulation_act->setChecked(true);
        }
        else
        {
            set2Dsimulation_act->setChecked(true);
            set3Dsimulation_act->setChecked(false);
            
            //Load path into settings.ini
            gui_settings->setValue("OpenFCSTbin", binary_path_2D);
            
            QMessageBox::information(this, "Binary selected", "Simulation set to 2D", QMessageBox::Ok );
        }
    }
    else if (obj == set3Dsimulation_act)
    {
        if(set3Dsimulation_act->isChecked()==false && isTriggered==false)
        {
            set3Dsimulation_act->setChecked(true);
        }
        else
        {
            set3Dsimulation_act->setChecked(true);
            set2Dsimulation_act->setChecked(false);
            //Load bin path into settings.ini
            gui_settings->setValue("OpenFCSTbin", binary_path_3D);
          
            QMessageBox::information(this, "Binary selected", "Simulation set to 3D", QMessageBox::Ok );
        }   
    }
}
//------------------------------------------------------------------------------------------------//
void MainWindow::newProject()
{
    if(current_files.size() >0)
    {
        newInstance();
    }
    else
        next_click();
}

//------------------------------------------------------------------------------------------------//
void MainWindow::next_click()
{

    switch (currentState)
    {
    case Start:
    {
        //get a file list of current files
        if(setWorkingDir())
        {
            FCSTBinAddr = gui_settings->value("OpenFCSTbin").toString();
            log("Creating main param file.\n");
            callFCST(QStringList() << paramArg, false, workingDir,FCSTBinAddr); // ./fuel-cell.bin -p

            //if all went ok open a tab
            usleep(1000000); //Sleep while fcst produces the files and OS file system updates

            //mainFileName = gui_settings->setValue("WorkDir", workingDir);

            try
            {
                create_tab(mainFileName);
            }
            catch(std::runtime_error &e)
            {
                log( std::string(e.what())+"\n");
                QMessageBox::warning(this, tr("Error"), tr(e.what()), QMessageBox::Ok );
                return; //Do not proceed in state
            }
            //Proceed to next state
            currentState = Main;
            nextButton->setText("Next");

        }
    }
    break;
    case Main:
    {
        //Process the main
        saveAll();
        log("Creating other param files.\n");
        FCSTBinAddr = gui_settings->value("OpenFCSTbin").toString();
        callFCST(QStringList() << paramArg << mainFileName, false, workingDir, FCSTBinAddr); // ./fuel-cell.bin -p main.xml

        //Get default data/opt files, open tabs
        usleep(1000000); //Sleep while fcst produces the files and OS file system updates

        try{
            create_tab(dataFileName);
        }
        catch(std::runtime_error &e){
            log( std::string(e.what())+"\n");
            QMessageBox::warning(this, tr("Error"), tr(e.what()), QMessageBox::Ok );
            return; //Do not proceed in state
        }


        try{
            create_tab(optFileName);
        }
        catch(std::runtime_error &e){
            //We can live with this error. FCST may not be outputting a opt file
            log( std::string(e.what())+"\n");
        }

        setWindowModified(false);

        //Proceed to next state
        currentState = Ready;
        nextButton->setText("Run");

    }
    break;
    case Ready:
    {
        log("Getting ready to start simulation.\n");
        //We have main + other files, presuming user has edited them we may proceed to Running
        //First save files
        saveAll();

        //Set up the file watched tol show the user new OpenFCST outputfiles
        fileWatcher = new QFileSystemWatcher(this);
        fileWatcher->addPath(workingDir);
        connect(fileWatcher, SIGNAL(directoryChanged (QString)), this, SLOT(updateFolderView(QString)));
        updateFolderView(workingDir);
        FCSTBinAddr = gui_settings->value("OpenFCSTbin").toString();
        //Run OpenFCST
        callFCST(QStringList() << mainFileName, true, workingDir,FCSTBinAddr); // ./fuel-cell.bin  main.xml

        //Proceed to next state
        currentState = Running;

        //Button changes to "End Simulation" in case to user wishes to kill it
        nextButton->setText("End Simulation");
    }
    break;
    case Running:
    {
        //Give user the option to kill the running simulation
        if(QMessageBox::question(this, "Terminate simulation", ""
                "Do you really want to end the running simulation?",
                QMessageBox::Yes|QMessageBox::No) == QMessageBox::Yes)
        {
            //User chose to terminate
            killFCST();

            //Go back to Ready state
            currentState = Ready;
            nextButton->setText("Run");
        }

    }
    break;
    default:
    {
        //State gone out of scope
        currentState = Start;
    }
    }

}

//----------------------------------------------EOF-----------------------------------------------//

}
