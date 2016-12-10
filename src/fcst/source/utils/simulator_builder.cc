//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2009-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: simulation_builder.cc
//    - Description: Class that generates the simulations (Called by main.cc)
//    - Developers: M. Secanell, P Dobson, Valentin N. Zingan, M. Sabharwal
//
//---------------------------------------------------------------------------

#include "utils/simulator_builder.h"

using namespace boost;
namespace po = boost::program_options;

//---------------------------------------------------------------------------
template<int dim>
SimulatorBuilder<dim>::SimulatorBuilder()
:
dakota_use(false),
dakota_direct(false),
program_name("openFCST"),
program_version("Version 0.3 (beta), April 2015"),
data(shared_ptr<FuelCell::ApplicationCore::ApplicationData>(new FuelCell::ApplicationCore::ApplicationData()))
{
}

//---------------------------------------------------------------------------
template<int dim>
SimulatorBuilder<dim>::~SimulatorBuilder()
{
    data.reset();
}

//---------------------------------------------------------------------------
template<int dim>
void SimulatorBuilder<dim>::parse_inputs(int argc, char *argv[])
{
    // Create options description
    po::options_description desc("Allowed options");

    //-- Declare parameters:
    //
    //The desc.add_options() returns an object with operator().
    //This allows us to write the option description in the following way. Pretty neat!
    desc.add_options()
        ("help,h", "Output help information")
        ("about,a", "Output information about OpenFCST")
        ("version,v", "Output version information")

        ("default-param,p", po::value<std::string>()->implicit_value(""), "Output default main parameter file. \n"
                                                                          "If file name is supplied it will be interpreted as a main file"
                                                                          "and corresponding files, i.e. data files and optimizatin files"
                                                                          "(if openFCST compiled using Dakota) will be generated.")

        ("prm2xml,c", po::value<std::string>(),                "Convert .prm parameter file to .xml file. Expected argument is main file."
                                                               "Multiple .xml parameter files corresponding to main, data and optimization files"
                                                               "may be produced based on of the input file.")

        ("dakota-files,D", po::value< std::vector<std::string> >(),  "Specify Dakota parameter files, "
                                                                     "'parameter_file' to set optimization parameters,"
                                                                     "'results_file' to set responses, supply file "
                                                                     "names in that order. ")

        ("input-file", po::value<std::string>(),       "Main input file. Filename can be supplied without flag.")

        ("xml2prm,r", po::value<std::string>(),        "Re-Converts .xml paramter file to .prm. Expected argument is main file."
                                                       "Multiple .prm parameter files corresponding to main, data and optimization files"
                                                       "may be produced based on of the input file.");


    //Set input-file to be "positional" option
    //i.e. we do not have to supply the argument name
    //we can just call the binary with and supply the file name
    po::positional_options_description positional_option;
    positional_option.add("input-file", -1);


    //Variable map where the user inputs will be stored
    //variable_map inherits std::map<std::string, variable_value>
    po::variables_map vm;

    //Try parse the users input
    //User mistakes will cause exceptions to be thrown
    try{
        po::store(po::command_line_parser(argc, argv).options(desc).positional(positional_option).run(), vm);
        po::notify(vm);
    }
    catch(po::error& e)
    {
        //Report the user's mistake
        FcstUtilities::log << "ERROR: " << e.what() << std::endl << std::endl;
        FcstUtilities::log << desc << std::endl;
        exit(1);
    }

    //Begin evaluating user inputs
    if (vm.count("help")) {
        FcstUtilities::log << desc << "\n";
        exit(0);
    }

    if (vm.count("about")) {
        print_logo();
        exit(0);
    }

    if (vm.count("version")) {
        FcstUtilities::log << program_name << " :: "
                << program_version << std::endl;
        exit(0);
    }

    if (vm.count("default-param")) {

        std::string fileName = vm["default-param"].as<std::string>();

        if(fileName == ""){
            output_default_main();
        }
        else{
            output_default_other(fileName);
        }

        exit(0);
    }

    if (vm.count("prm2xml")) {
        
        std::string fileName = vm["prm2xml"].as<std::string>();
        convert_file(fileName);
        exit(0);
    }    
    
    if (vm.count("xml2prm")) {        
        
        std::string fileName = vm["xml2prm"].as<std::string>();
        convert_file(fileName);
        exit(0);
    }
    
    if (vm.count("dakota-files")) {

        dakota_use = true;
        std::vector<std::string> fileNames = vm["dakota-files"].as<std::vector<std::string>>();

        if(fileNames.size() != 2){
            FcstUtilities::log << "ERROR: 2 file names must be supplied for Dakota \n";
            exit(1);
        }

        dakota_parameters = fileNames[0];
        dakota_results = fileNames[1];
    }


    //Finally we check for an input parameter file
    if (vm.count("input-file")) {
        input_file =  vm["input-file"].as<std::string>();
    }
    else{
        FcstUtilities::log << "ERROR: No input file/arguments specified \n";
        FcstUtilities::log << desc << "\n";
        exit(1);
    }

}

//---------------------------------------------------------------------------
template<int dim>
void SimulatorBuilder<dim>::scan()
{
    sim_selector = shared_ptr<SimulationSelector<dim> > (new SimulationSelector<dim>(data));
    this->declare_parameters(param);

    FcstUtilities::read_parameter_files(param,
                                        input_file);

    initialize(param);
}

//---------------------------------------------------------------------------
template<int dim>
void SimulatorBuilder<dim>::run()
{

    FcstUtilities::log.push("MAIN");
    print_logo();

    FcstUtilities::log<<"Running "<<analysis_type<<" simulation"<<std::endl;
    timer.restart();
         
    if (analysis_type.compare("Analysis") == 0) {
        // Select linear, nonlinear and adaptive refinement objects:
        app_lin = sim_selector->select_application();
        newton = sim_selector->select_solver(app_lin.get());        
        solver = sim_selector->select_solver_method(app_lin.get(), newton.get());
        // Solve problem:
        solver->solve(simulator_parameter_file_name, param);
    }
    else if (analysis_type.compare("UnitTests") == 0) 
        run_test();
    
    else if (analysis_type.compare("ParametricStudy") == 0)
        param_study.run(param, simulator_parameter_file_name, sim_selector);
        
    else if (analysis_type.compare("PolarizationCurve") == 0)
        curve.run(param, simulator_parameter_file_name, sim_selector);    
    
    else if (analysis_type.compare("Optimization") == 0){
        run_optimization();
    }
    
    else {
        FcstUtilities::log << " == ERROR =="<<std::endl;
        FcstUtilities::log << " set Analysis type in the main file does not match any of the following: UnitTests | Analysis | ParametricStudy | PolarizationCurve | Optimization"<<std::endl;
        FcstUtilities::log << " Please modify either main.prm or main.xml file as appropriate."<<std::endl;
        FcstUtilities::log << " "<<std::endl;
        FcstUtilities::log << "Exiting program... "<<std::endl;
        FcstUtilities::log << " "<<std::endl;
        FcstUtilities::log << " ===="<<std::endl;
    }        
    
    timer.stop();
    print_timer_info();
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// PRIVATE
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
template<int dim>
void SimulatorBuilder<dim>::declare_parameters(ParameterHandler& param) const
{
    
    sim_selector->declare_parameters(param);
    
    param.enter_subsection("Simulator");
    {
        param.declare_entry("simulator parameter file name",
                            "",
                            Patterns::Anything(),
                            "Enter the file you would like to use to read data"
                            "for the application. This is the data.xml file that is generated"
                            "using the GUI");
                
        param.declare_entry("Analysis type",
                            "None",
                            Patterns::Selection("None | UnitTests | Analysis | ParametricStudy | PolarizationCurve | Optimization"),
                            "Specify the type of study you would like to perform. The options are: \n"
                            "Unit Tests: Performs internal tests to make sure openFCST is running correctly. \n"
                            "Analysis: Performs an analysis run, i.e. uses the model in simulator name and the provided here"
                            " to perform an analysis run. \n"
                            "ParametricStudy: Performs a parameteric study using same values as analysis run and "
                            " modifying parameter specified in Parametric study \n"
                            "PolarizationCurve: Computes a polarization curve using same values as analysis run. \n"
                            "Optimization: Solves an optimization problem using the values in the analysis files"
                            " and the data in \"Optimization study\" and within.");
        
        param.enter_subsection("Optimization");
        {
            param.declare_entry("optimization parameter file name",
                                "",
                                Patterns::Anything(),
                                "Enter the file you would like to use to read optimization"
                                "options from. This is the opt.xml file that is generated"
                                "using the GUI");
            param.declare_entry("Dakota direct",
                                "true",
                                Patterns::Bool(),
                                "Set to true if you would like OpenFCST to control Dakota. If set to false, "
                                "openFCST will simply run once and generate a file that Dakota can read.");
        }
        param.leave_subsection();
    }
    param.leave_subsection();
    
    param_study.declare_parameters(param);
    
    curve.declare_parameters(param);
    
    param.enter_subsection("Logfile");
    {
        param.declare_entry("Logfile name", 
                            "logfile.log", 
                            Patterns::Anything(),
                            "File name where all the output to screen is recoreded");
        param.declare_entry("FileDepth", "10000", Patterns::Integer());
        param.declare_entry("ConsoleDepth", "10000", Patterns::Integer());        
    }
    param.leave_subsection();
    
    
    data->declare_parameters(param);
}


//---------------------------------------------------------------------------
template<int dim>
void SimulatorBuilder<dim>::initialize(ParameterHandler& param)
{
    
    sim_selector->initialize(param);
            
    param.enter_subsection("Simulator");
    {
        simulator_parameter_file_name = param.get("simulator parameter file name");
        analysis_type = param.get("Analysis type");
        
        param.enter_subsection("Optimization");
        {
            optimization_parameter_file_name = param.get("optimization parameter file name");
            dakota_direct = param.get_bool("Dakota direct");
        }
        param.leave_subsection();       
    }
    param.leave_subsection();
    
    param_study.initialize(param);
    
    curve.initialize(param);
    
    param.enter_subsection("Logfile");
    {
        this->open_logfile(param.get("Logfile name"));
        FcstUtilities::log.depth_file(param.get_integer("FileDepth"));
        FcstUtilities::log.depth_console(param.get_integer("ConsoleDepth"));
    }
    param.leave_subsection();
    
    data->initialize(param);
}


//---------------------------------------------------------------------------
template<int dim>
void SimulatorBuilder<dim>::open_logfile (const std::string& extname)
{
    // If a parameter file was specified, use its base name for
    // the logfile.
    std::string name = extname;

    // Still no name, use hostname, allowing parallel programs to
    // write into different logfiles.
    if (name == "")
    {
        char log_name[40];
        gethostname(log_name, 19);
        name = log_name;
    }

    std::string::size_type pos = name.rfind(".log");
    if (pos < name.size())
    {
        name.erase(pos);
    }
    name += ".log";

    log_file = boost::shared_ptr<std::ofstream> (new std::ofstream(name.c_str()));

    // Log into file
    FcstUtilities::log.attach(*log_file);
}


//---------------------------------------------------------------------------
template<int dim>
void SimulatorBuilder<dim>::run_optimization()
{
    

#ifdef OPENFCST_WITH_PETSC
    FcstUtilities::log << "//=======================================================//" << std::endl;
    FcstUtilities::log << "//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//" << std::endl;
    FcstUtilities::log << "//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//" << std::endl;    
    FcstUtilities::log << ""  << std::endl;
    FcstUtilities::log<<" OpenFCST with Dakota not implemented in parallel. "<<std::endl;
    FcstUtilities::log<<" The same program will run in each processor providing no advantage. "<<std::endl;
    FcstUtilities::log<<"  "<<std::endl;
    FcstUtilities::log<<" Please run code without using mpirun. "<<std::endl;
    FcstUtilities::log << ""  << std::endl;    
    FcstUtilities::log << "//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//" << std::endl;
    FcstUtilities::log << "//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//" << std::endl;    
    FcstUtilities::log << "//=======================================================//" << std::endl;
        /*
    exit(-1);
        */
#endif


    if (dakota_direct)
    {
#ifdef _WITH_DAKOTA
        // NOTE: Must declare these in order for parameter handler to not complain when reading the parameter file specified.
        //        Not exclusively required for dakota application to run.
        Dakota::ParallelLibrary parallel_lib;
        shared_ptr<Dakota::ProblemDescDB> problem_db(new Dakota::ProblemDescDB (parallel_lib));
        SIM::DakotaApplication optimization(problem_db, optimization_parameter_file_name);
        optimization.declare_parameters(param);
        optimization.manage_inputs(param);

        Dakota::DirectApplicInterface* optimization_interface;

        if (optimization.use_NLS())
        {
            FcstUtilities::log<<"Entering DakotaLeastSquaresInterface"<<std::endl;
            optimization_interface = new SIM::DakotaLeastSquaresInterface<dim> (optimization, problem_db, param, data, sim_selector, simulator_parameter_file_name);
        }
        else
        {
            FcstUtilities::log<<"Entering DakotaDirectInterface"<<std::endl;
            optimization_interface = new SIM::DakotaDirectInterface<dim > (optimization, problem_db, param, data, sim_selector, simulator_parameter_file_name);
        }

        optimization.assign_interface(optimization_interface);
        optimization.run();

        FcstUtilities::log << "Optimization completed" << std::endl;
        if (!save_transfer_files)
        {
            remove( "._mesh_to_be_transfered.msh" );
            remove( ".transfer_solution.FEVector" );
        }
        //Delete
        // Allows results to be retrieved  --> TODO: needs to be redefined for multi-start methods
        //   Dakota::RealVector vars = optimization.vars_results();
        //   Dakota::RealVector resp = optimization.resp_results();
        //   FcstUtilities::log << "Results retrieved" << std::endl;
#else
        FcstUtilities::log<<"Dakota not available"<<std::endl;
#endif
    }
    else if (dakota_use)
    {

        app_lin = sim_selector->select_application();
        newton = sim_selector->select_solver(app_lin.get());
        solver = sim_selector->select_solver_method(app_lin.get(), newton.get());

        SIM::DakotaInterface<dim> optimization (simulator_parameter_file_name,
                                                param,
                                                dakota_parameters,
                                                dakota_results,
                                                *app_lin,
                                                *newton);
        optimization.run();
        FcstUtilities::log << "Optimization completed" << std::endl;
        if (!save_transfer_files)
        {
            remove( "._mesh_to_be_transfered.msh" );
            remove( ".transfer_solution.FEVector" );
        }

    }
    else
    {
        abort();
    }
}

//---------------------------------------------------------------------------
template<int dim>
void SimulatorBuilder<dim>::run_test()
{
    FcstUtilities::log <<    "==================Running Unit Tests===================" << std::endl;

    if (FcstTestSuite::run_tests())
        FcstUtilities::log << "===========Unit Tests Successfully Completed===========" << std::endl;
    else
        FcstUtilities::log << "==========Unit Tests Completed with Errors!============" << std::endl;
}


//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// XML converter routines:
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
template<int dim>
void SimulatorBuilder<dim>::convert_file(std::string main_file)
{
    //Create simulation selector:
    sim_selector = shared_ptr<SimulationSelector<dim> > (new SimulationSelector<dim>(data));
    
    // Check if the file exists:
    if (not FcstUtilities::file_exists(main_file)){
        FcstUtilities::log << "ERROR: Main file not found!" << std::endl;
        exit(1);
    }
    
    // Check the file extension and select conversion type:
    FileConversionOption conversion_type = NONE;
    
    std::size_t pos = main_file.find(".");
    std::string ext = main_file.substr(pos);
    
    //filepath for pathing the location of the data file
    std::string filePath;
    filePath = "";
    
    if( ext == ".prm")
    {
        FcstUtilities::log << "Converting .prm main, data and opt files to .xml. \n";
        conversion_type = PRM2XML;
        
    }
    else if (ext == ".xml")
    {
        FcstUtilities::log << "Converting .xml main, data and opt files to .prm. \n";
        conversion_type = XML2PRM;
        // set filepath only by xml to prm conversion because otherwise the reading function is not 
        // able to handle the further includes inside the prm file during a prm to xml conversion
        if(main_file.find_last_of("/") != std::string::npos)
        {
            filePath = main_file.substr(0,main_file.find_last_of("/")) + "/";
        }
        else
        {
            filePath = "";
        }
    }
    else
    {
        FcstUtilities::log << "ERROR: convert_file requires either a .prm or a .xml main file. \n";
        exit(1);
    }

    // Read the main.prm / .xml file and initialize parameter in SimulationBuilder
    param.clear();
    this->declare_parameters(param);
    curve.declare_parameters(param);
    FcstUtilities::read_parameter_files(param, main_file);
    this->initialize(param);
        
    //-- Modify necessary lines
    // Store actual data file location in main file:
    std::string data_file_path;
    param.enter_subsection("Simulator");
    {
        data_file_path = filePath + param.get("simulator parameter file name");
    }
    param.leave_subsection();
    // and modify to be data.xml / data.prm
    param.enter_subsection("Simulator");
    {
        if (conversion_type == XML2PRM)
            param.set("simulator parameter file name","data.prm");
        else if (conversion_type == PRM2XML)
            param.set("simulator parameter file name","data.xml");
        else
        {
            FcstUtilities::log << "ERROR: convert_file requires either a .prm or a .xml main file. \n";
            exit(1);
        }
                
    }
    param.leave_subsection();
    
    //-- Start crating new files:
    // Output the main
    if (conversion_type == XML2PRM)
    {
        FcstUtilities::log << "Creating converted main file in PRM" << std::endl;
        FcstUtilities::print_parameter_file_PRM(param, "main.prm");
    }     
    else if (conversion_type == PRM2XML)
    {
        FcstUtilities::log << "Creating converted main file in XML" << std::endl;
        FcstUtilities::print_parameter_file_XML(param, "main.xml");            
    }
    else
    {
        FcstUtilities::log << "ERROR: convert_file requires either a .prm or a .xml main file. \n";
        exit(1);
    }
    
    // Check if data file exists, if it does, read it, convert it and create new file:  
    if (FcstUtilities::file_exists(data_file_path))
    {
        //Create application to generate parameters
        app_lin = sim_selector->select_application();
        newton  = sim_selector->select_solver( app_lin.get() );
        solver  = sim_selector->select_solver_method( app_lin.get(), newton.get() );
        
        //Output data.xml, use separate param object to keep parameters separate
        ParameterHandler data_param;
        solver->declare_parameters(data_param);
        FcstUtilities::read_parameter_files(data_param,data_file_path);
        
        //Output the data file
        if (conversion_type == XML2PRM)
        {
            FcstUtilities::log << "Creating converted data file in PRM" << std::endl;
            FcstUtilities::print_parameter_file_PRM(data_param, + "data.prm");
        }
        else if (conversion_type == PRM2XML)
        {
            FcstUtilities::log << "Creating converted data file in XML" << std::endl;
            FcstUtilities::print_parameter_file_XML(data_param, "data.xml");
        }
        else
        {
            FcstUtilities::log << "ERROR: convert_file requires either a .prm or a .xml main file. \n";
            exit(1);
        }
    }        
    else
    {
        FcstUtilities::log << "ERROR: Data file '" << data_file_path << "' not found! \n";
        exit(1);
    }
    
    //Create Dakota components to generate
    #ifdef _WITH_DAKOTA
    {
        std::string opt_file_path;
        param.enter_subsection("Simulator");
        {
            param.enter_subsection("Optimization");
            {
                opt_file_path = param.get("optimization parameter file name");
            }
            param.leave_subsection();
        }
        param.leave_subsection();
        
        if (FcstUtilities::file_exists(opt_file_path))
        {
            
            Dakota::ParallelLibrary parallel_lib;
            ParameterHandler opt_param;
            shared_ptr<Dakota::ProblemDescDB> problem_db(new Dakota::ProblemDescDB (parallel_lib));
            SIM::DakotaApplication optimization(problem_db, optimization_parameter_file_name);
            optimization.declare_parameters(opt_param);
            
            //Output opt.xml, use separate param object to keep parameters separate
            FcstUtilities::read_parameter_files(opt_param,
                                                opt_file_path);
            
            if (conversion_type == XML2PRM)
            {
                FcstUtilities::log << "Creating converted opt file in XML" << std::endl;
                FcstUtilities::print_parameter_file_XML(opt_param, "opt.xml");
            }
            else if (conversion_type == PRM2XML)
            {
                FcstUtilities::log << "Creating converted opt file in XML" << std::endl;
                FcstUtilities::print_parameter_file_PRM(opt_param, "opt.xml");
            }
            else
            {
                FcstUtilities::log << "ERROR: convert_file requires either a .prm or a .xml main file. \n";
                exit(1);
            }
        }
        else
        {
            FcstUtilities::log << "No opt file detected, skipping." << std::endl;
        }
    }
    #endif
    
    sim_selector.reset();    
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
template<int dim>
void SimulatorBuilder<dim>::output_default_main(){
    
    FcstUtilities::log << "Creating default parameter file main.xml\n";
    sim_selector = shared_ptr<SimulationSelector<dim> > (new SimulationSelector<dim>(data));
    this->declare_parameters(param);
    curve.declare_parameters(param);
    
    // Modify main_file to make sure that the data file that is read is data.xml
    param.enter_subsection("Simulator");
    {
        param.set("simulator parameter file name","data.xml");
    }
    param.leave_subsection();
    
    FcstUtilities::print_parameter_file_XML(param, "main.xml");
    
    sim_selector.reset();
    
}

//---------------------------------------------------------------------------
template<int dim>
void SimulatorBuilder<dim>::output_default_other(std::string main_file){

    // Check if the file exists:
    if (not FcstUtilities::file_exists(main_file)){
        FcstUtilities::log << "ERROR: Main file not found!" << std::endl;
        exit(1);
    }


    FcstUtilities::log << "Creating default parameter file data.xml \n";

    sim_selector = shared_ptr<SimulationSelector<dim> > (new SimulationSelector<dim>(data));

    //Read in main
    param.clear();
    this->declare_parameters(param);
    curve.declare_parameters(param);    
    FcstUtilities::read_parameter_files(param, main_file);

    //Initialize main
    this->initialize(param);

    //Create application to generate parameters
    app_lin = sim_selector->select_application();
    newton  = sim_selector->select_solver( app_lin.get() );
    solver  = sim_selector->select_solver_method( app_lin.get(), newton.get() );

    //Output data.xml, use separate param object to keep parameters separate
    ParameterHandler param2;
    solver->declare_parameters(param2);
    FcstUtilities::print_parameter_file_XML(param2, "data.xml");

    
    //Create Dakota components to generate
    #ifdef _WITH_DAKOTA
    {
        FcstUtilities::log << "Creating default parameter file opt.xml \n";
        Dakota::ParallelLibrary parallel_lib;
        ParameterHandler param3;
        shared_ptr<Dakota::ProblemDescDB> problem_db(new Dakota::ProblemDescDB (parallel_lib));
        SIM::DakotaApplication optimization(problem_db, optimization_parameter_file_name);
        optimization.declare_parameters(param3);
        //Output opt.xml, use separate param object to keep parameters separate
        FcstUtilities::print_parameter_file_XML(param3, "opt.xml");
    }
    #endif
    
    sim_selector.reset();
    
    
}

//---------------------------------------------------------------------------
template<int dim>
void SimulatorBuilder <dim>::print_logo() const
{

     FcstUtilities::log << "//=====================================================================//" << std::endl;
     FcstUtilities::log << "//                                                                     " << std::endl;
     FcstUtilities::log << "//  OpenFCST: The open-source Fuel Cell Simulation Toolbox             " << std::endl;
     FcstUtilities::log << "//                                                                     " << std::endl;
     FcstUtilities::log << "//  Version: 0.3, Nov. 2016                                            " << std::endl;    
     FcstUtilities::log << "//                                                                     " << std::endl;     
     FcstUtilities::log << "//=====================================================================//" << std::endl;    
     FcstUtilities::log << "//                                                                     " << std::endl;
     FcstUtilities::log << "//  Current Developers                                                 " << std::endl;
     FcstUtilities::log << "//  ------------------                                                 " << std::endl;
     FcstUtilities::log << "//                                                                     " << std::endl;     
     FcstUtilities::log << "//  Marc Secanell, 2006- (a): Framework, cathode, anode, PEM, (...)    " << std::endl;
     FcstUtilities::log << "//  Alex Jarauta, 2016 (a): Local matrix info, enhancements            " << std::endl;
     FcstUtilities::log << "//  Aslan Kosakian, 2014- (a): Transient simulations.                  " << std::endl;     
     FcstUtilities::log << "//  Andreas Putz, 2012- (b): Framework                                 " << std::endl;
     FcstUtilities::log << "//  Mayank Sabharwal (a), 2014-: Microscale simulations                " << std::endl;
     FcstUtilities::log << "//  Jie Zhou (a), 2013-: Two-phase flow simuations                     " << std::endl;
     FcstUtilities::log << "//                                                                     " << std::endl;
     FcstUtilities::log << "//                                                                     " << std::endl;
     FcstUtilities::log << "//  Past Developers                                                    " << std::endl;
     FcstUtilities::log << "//  ---------------                                                    " << std::endl;     
     FcstUtilities::log << "//                                                                     " << std::endl;
     FcstUtilities::log << "//  Chad Balen, 2014-16 (a): CMake, fluid flow and many enhancements   " << std::endl;
     FcstUtilities::log << "//  Madhur Bhaiya, 2012-14 (a): Non-isothermal PEMFC                   " << std::endl;    
     FcstUtilities::log << "//  Peter Dobson, 2009-11 (a): Least-square parameter estimation       " << std::endl;
     FcstUtilities::log << "//  Michael Moore, 2010-12 (a): Double-trap ORR model                  " << std::endl;
     FcstUtilities::log << "//  Philip Wardlaw, 2012-14 (a): Micro-scale models and GUI            " << std::endl;
     FcstUtilities::log << "//  Valentin N. Zingan, 2012-14 (a): Fluid flow and enhancements       " << std::endl;
     FcstUtilities::log << "//                                                                     " << std::endl;
     FcstUtilities::log << "//  Other Contributors                                                 " << std::endl;     
     FcstUtilities::log << "//  ------------------                                                 " << std::endl;          
     FcstUtilities::log << "//                                                                     " << std::endl;    
     FcstUtilities::log << "//  Kailyn Domican, 2012-14 (a): Optimization docs and enhancements.   " << std::endl;
     FcstUtilities::log << "//  Simon Matterns, 2016-16 (a): Many GUI enhancements.                " << std::endl;
     FcstUtilities::log << "//                                                                     " << std::endl;
     FcstUtilities::log << "//                                                                     " << std::endl;    
     FcstUtilities::log << "//  a) Energy Systems Design Lab, University of Alberta                " << std::endl;
     FcstUtilities::log << "//  b) Automotive Fuel Cell Cooperation (AFCC) Corp.                   " << std::endl; 
     FcstUtilities::log << "//                                                                     " << std::endl;
     FcstUtilities::log << "//                                                                     " << std::endl;          
     FcstUtilities::log << "// (c) Copyright 2007-2016 by Energy Systems Design                    " << std::endl;
     FcstUtilities::log << "//     Laboratory, University of Alberta                               " << std::endl;  
     FcstUtilities::log << "//                                                                     " << std::endl; 
     FcstUtilities::log << "// This software is distributed under the MIT License.                 " << std::endl;
     FcstUtilities::log << "// For more information, see the README file in /doc/LICENSE           " << std::endl;
     FcstUtilities::log << "//                                                                     " << std::endl;     
     FcstUtilities::log << "//=====================================================================//" << std::endl;    
     FcstUtilities::log << "//=====================================================================//" << std::endl;    

}

//---------------------------------------------------------------------------
template<int dim>
void SimulatorBuilder<dim>::print_timer_info() const
{
      FcstUtilities::log.push("MAIN");
      FcstUtilities::log << "The program was executed in: " << timer.wall_time() << " seconds" << std::endl;
      FcstUtilities::log << "=============== END ====================" << std::endl;
      FcstUtilities::log.pop();
}

//---------------------------------------------------------------------------
// Explicit Instantiations
template class SimulatorBuilder<deal_II_dimension>;
