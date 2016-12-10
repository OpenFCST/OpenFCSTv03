//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2014 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: settings_window.cc
//    - Description: GUI for showing FCST running settings
//    - Developers: Simon Mattern, 2016
//
//---------------------------------------------------------------------------

#include "settings_window.h"

using namespace FCSTGUI;

//------------------------------------------------------------------------------------------------//
settingsWindow::settingsWindow(QSettings *guiSettings)
{
    gui_setting = guiSettings;
    settings_dialog = new QDialog();
    settings_dialog->setWindowTitle(tr("Settings"));
    
    createWidgets();
    LoadSettings();
    createLayout();
    
}
//------------------------------------------------------------------------------------------------//
void settingsWindow::createWidgets()
{
    save_button = new QPushButton(tr("Save && Quit"));
    connect(save_button,SIGNAL(pressed()),this,SLOT(saveSettings()));
    
    cancel_button = new QPushButton(tr("Cancel"));
    connect(cancel_button,SIGNAL(pressed()),this,SLOT(closeSettings()));
    
    browse_2Dpath_button = new QPushButton("...");
    connect(browse_2Dpath_button,SIGNAL(pressed()),this, SLOT(browse()));
    
    browse_3Dpath_button =new QPushButton("...");
    connect(browse_3Dpath_button,SIGNAL(pressed()),this, SLOT(browse()));
    
    sim_argument_edit = createQLineEdit();
    main_file_edit  = createQLineEdit();
    data_file_edit  = createQLineEdit();
    opt_file_edit  = createQLineEdit();
    
    cpu_amount_spinbox = new QSpinBox();
    cpu_amount_spinbox->setRange(1,4);
    cpu_amount_spinbox->setMaximumSize(200,150);
    cpu_amount_spinbox->setEnabled(true);
    
    //path settings
    used_binary_box = new QComboBox();
    bin_2D_path_box = new QComboBox();
    bin_2D_path_box->setDuplicatesEnabled(false);
    bin_3D_path_box = new QComboBox();
    bin_3D_path_box->setDuplicatesEnabled(false);
    
    spacer = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);
}
//------------------------------------------------------------------------------------------------//
void settingsWindow::LoadSettings()
{
    cpu_amount_spinbox->setValue(gui_setting->value("CpuAmount").toInt());
    sim_argument_edit->setText(gui_setting->value("OpenFCSTparamArg").toString());
    main_file_edit->setText(gui_setting->value("mainFileName").toString());
    data_file_edit->setText(gui_setting->value("dataFileName").toString());
    opt_file_edit->setText(gui_setting->value("optFileName").toString());
    
    used_binary_box->addItem(gui_setting->value("OpenFCSTbin").toString());
    used_binary_box->setEditable(false);
    bin_2D_path_box->addItem(gui_setting->value("bin2DPath").toString());
    bin_3D_path_box->addItem(gui_setting->value("bin3DPath").toString());
}
//------------------------------------------------------------------------------------------------//
void settingsWindow::createLayout()
{
    QPointer<QLabel> bin2DLabel = new QLabel("2D binary path");
    QPointer<QLabel> bin3DLabel = new QLabel("3D binary path");
    QPointer<QLabel> binLabel =new QLabel("Used binary \nduring the simulation");
    QPointer<QLabel> parallelLabel =new QLabel("Parralel Running \n(Numbers of CPU used during the simulation)");
    QPointer<QLabel> mainfileLabel =new QLabel("Main File Name");
    QPointer<QLabel> datafileLabel =new QLabel("Data File Name");
    QPointer<QLabel> optfileLabel =new QLabel("Opt File Name");
    QPointer<QLabel> usedparaLabel =new QLabel("Used parameter during the simulation");
    
    //Gridlayout
    QPointer<QGridLayout> mainlayout = new QGridLayout();
    mainlayout->addWidget(bin3DLabel,0,0);
    mainlayout->addWidget(bin_3D_path_box,0,1);
    mainlayout->addWidget(browse_3Dpath_button,0,2);
    mainlayout->addWidget(bin2DLabel,1,0);
    mainlayout->addWidget(bin_2D_path_box,1,1);
    mainlayout->addWidget(browse_2Dpath_button,1,2);
    mainlayout->addWidget(binLabel,2,0 );
    mainlayout->addWidget(used_binary_box,2,1 );
    mainlayout->addWidget(parallelLabel,3,0);
    mainlayout->addWidget(cpu_amount_spinbox,3,1);
    mainlayout->addWidget(mainfileLabel,4,0);
    mainlayout->addWidget(main_file_edit,4,1);
    mainlayout->addWidget(datafileLabel,5,0);
    mainlayout->addWidget(data_file_edit,5,1);
    mainlayout->addWidget(optfileLabel,7,0);
    mainlayout->addWidget(opt_file_edit,7,1);
    mainlayout->addWidget(usedparaLabel,8,0);
    mainlayout->addWidget(sim_argument_edit,8,1);
    mainlayout->addItem(spacer,9,0);
    mainlayout->addWidget(save_button,9,1,Qt::AlignRight);
    mainlayout->addWidget(cancel_button,9,2);
    
    //mainwindow with mainlayout (QgridLayout) in it
    settings_dialog->setLayout(mainlayout);
    settings_dialog->setAttribute(Qt::WA_DeleteOnClose);
    settings_dialog->exec();
}
//--------------------------------------------------------------------------//
//--------------------------------------------------------------------------//
//---------------------------Protected and Private Functions----------------//
//--------------------------------------------------------------------------//
//--------------------------------------------------------------------------//
QLineEdit* settingsWindow::createQLineEdit()
{
    QPointer<QLineEdit> LineEd = new QLineEdit();
    LineEd->setReadOnly(true);
    LineEd->setMaximumSize(200,150);
    return LineEd;
}
//-----------------------------------------------------------------------------
bool settingsWindow::fileExists(QString path)
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
//--------------------------------------------------------------------------//
//--------------------------------------------------------------------------//
//---------------------------SLOTS------------------------------------------//
//--------------------------------------------------------------------------//
//--------------------------------------------------------------------------//
void settingsWindow::closeSettings()
{
    settings_dialog->close();
}

//--------------------------------------------------------------------------//
void settingsWindow::saveSettings()
{
    int iValue = cpu_amount_spinbox->value();
    QString sVal = QString::number(iValue);
    gui_setting->setValue("CpuAmount", sVal);
    gui_setting->setValue("bin2DPath", bin_2D_path_box->currentText());
    gui_setting->setValue("bin3DPath", bin_3D_path_box->currentText());
    settings_dialog->close();
}

//--------------------------------------------------------------------------//
void settingsWindow::browse()
{
    QPointer<QObject> obj = sender();
    QString binary_path_2D = bin_2D_path_box->currentText();
    QString binary_path_3D = bin_3D_path_box->currentText();
    do
    {
        if(obj==browse_2Dpath_button)
        {
            //select 2D binary
            binary_path_2D = QFileDialog::getOpenFileName(this, 
                                                          tr("Select OpenFCST 2D simulation Binary File"),
                                                          QDir::currentPath(),
                                                          tr("2D Bin File(*2d.bin) (*2d.bin)"));
        }
        else if (obj==browse_3Dpath_button)
        {
            //select 3D binary
            binary_path_3D = QFileDialog::getOpenFileName(this,
                                                          tr("Select OpenFCST 3D simulation Binary File"),
                                                          QDir::currentPath(),
                                                          tr("3D Bin File(*3d.bin) (*3d.bin)"));
        }
    }  while ((fileExists(binary_path_2D) == false && binary_path_2D != "") ||
    (fileExists(binary_path_3D) == false && binary_path_3D != ""));
    if(obj==browse_2Dpath_button && binary_path_2D != "") 
    { 
        bin_2D_path_box->addItem(binary_path_2D); 
    }
    else if (obj==browse_3Dpath_button && binary_path_3D != "") 
    {
        bin_3D_path_box->addItem(binary_path_3D);
    }
}
//--------------------------------------------------------------------------//
settingsWindow::~settingsWindow()
{
    
}

