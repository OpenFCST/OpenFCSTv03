//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2006-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: geometry.cc
//    - Description: Used to store and output information for general fuel cell geometries
//    - Developers: M. Secanell and Peter Dobson
//                  Salome parser dealii::GridIn<dim, spacedim>::read_unv(): Valentin N. Zingan, University of Alberta, (C) 2012.
//    - $Id: geometry.cc 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#include <grid/geometry.h>

namespace NAME = FuelCellShop::Geometry;

//---------------------------------------------------------------------------
template <int dim>
NAME::GridBase<dim>::GridBase()
:
grid_in(new GridIn<dim>())
{}

//---------------------------------------------------------------------------
template <int dim>
NAME::GridBase<dim>::~GridBase()
{
    grid_in.reset();
}


//---------------------------------------------------------------------------
template <int dim>
void
NAME::GridBase<dim>::declare_GridGenerator_parameters(ParameterHandler& param)
{
    param.enter_subsection("Grid generation");
    {
        param.declare_entry("Type of mesh",
                            "CathodeMPL",
                            Patterns::Selection("GridExternal | CathodeMPL | PemfcMPL | Agglomerate | GridTest | HyperCube"),
                            "Type of mesh that you would like to generate.");
        param.declare_entry("File name",
                            "test.unv",
                            Patterns::Anything(),
                            "If the Type of mesh is set to GridExternal, then specify here the name of the mesh file");
        param.declare_entry("File type",
                            "Default",
                            Patterns::Selection("Default | vtk | unv | ucd | dbmesh |  xda | msh | netcdf | tecplot"),
                            "If the Type of mesh is set to GridExternal, then specify here the type of file (see deal.ii supported mesh types)");
        param.declare_entry("Initial refinement",
                            "0",
                            Patterns::Integer(),
                            "Number of times we want to globally refine the original grid before starting to solve");
    
        param.enter_subsection("Internal mesh generator parameters");
        {
            
            param.enter_subsection("Dimensions");
            {
                param.declare_entry("Cathode current collector width [cm]",      // thickness of the rib of the cathode bipolar plates (BPP)
                                    "0.1",
                                    Patterns::Double(),
                                    "Thickness of the rib of the bipolar plates (BPP)");
                param.declare_entry("Cathode channel width [cm]",
                                    "0.1",
                                    Patterns::Double(),
                                    "Thickness of the channel on the BPP");
                param.declare_entry("Cathode channel length [cm]",
                                    "0.2",
                                    Patterns::Double(),
                                    "Length of the channel on the Cathode in through plane direction ontop of fuel cell");
                param.declare_entry("Cathode channel height [cm]",
                                    "0.08",
                                    Patterns::Double(),
                                    "Height of the channel on the Cathode");
                param.declare_entry("Cathode channel entrance length [cm]",
                                    "0.1",
                                    Patterns::Double(),
                                    "Entrance Length of the channel on Cathode before fuel cell");
                param.declare_entry("Cathode channel exit length [cm]",
                                    "0.1",
                                    Patterns::Double(),
                                    "Exit Length of the Channel on Cathode after fuel cell");
                param.declare_entry ("Cathode CL thickness [cm]",    // thickness of the cathode catalyst layer
                                     "1E-3",
                                      Patterns::List(Patterns::Double()),
                                     "Thickness of the cathode catalyst layer [cm]");
                param.declare_entry ("Cathode GDL thickness [cm]",    // thickness of the cathode gas diffusion layer
                                     "2E-2",
                                     Patterns::Double(),
                                     "Thickness of the cathode gas diffusion layer [cm]");
                param.declare_entry ("Cathode MPL thickness [cm]",    // thickness of the cathode microporous layer
                                     "1.0e20",
                                     Patterns::Double(),
                                     "Thickness of the cathode microporous layer [cm]");
                param.declare_entry ("Membrane thickness [cm]",      // thickness of the membrane
                                     "2E-2",
                                     Patterns::Double(),
                                     "Thickness of the membrane");              
                param.declare_entry ("Anode MPL thickness [cm]",    // thickness of the anode microporous layer
                                     "0.0",
                                     Patterns::Double(),
                                     "Thickness of the anode microporous layer [cm]");                
                param.declare_entry ("Anode CL thickness [cm]",    // thickness of the anode catalyst layer
                                     "1E-3",
                                     Patterns::Double(),
                                     "Thickness of the anode catalyst layer [cm]");  
                param.declare_entry ("Anode GDL thickness [cm]",    // thickness of the anode gas diffusion layer
                                     "2E-2",
                                     Patterns::Double(),
                                     "Thickness of the anode gas diffusion layer [cm]");                
                param.declare_entry("Anode current collector width [cm]",      // thickness of the rib of the cathode bipolar plates (BPP)
                                    "0.1",
                                    Patterns::Double(),
                                    "Thickness of the rib of the bipolar plates (BPP) [cm]");
                param.declare_entry("Anode channel width [cm]",
                                    "0.1",
                                    Patterns::Double(),
                                    "Thickness of the channel on the BPP");
                
                param.declare_entry("Agglomerate radius [cm]",
                                    "100e-7",//cm
                                    Patterns::Double(),
                                    "Parameter used to define agglomerate size for spherical agglomerate mesh");
                param.declare_entry("Thin film thickness [cm]",
                                    "10e-7", //cm
                                    Patterns::Double(),
                                    "Parameter used to define thin film thickness for spherical agglomerate mesh");
                
                param.declare_entry("HyperCube edge",
                                    "1", //cm
                                    Patterns::Double(),
                                    "The cube will co from - this value to this value");

            }
            param.leave_subsection();
            
            param.enter_subsection("Mesh refinement parameters");
            {
                param.declare_entry("Initial vertical cell count",
                                    "6",
                                    Patterns::Integer(),
                                    "Number of cells we want in the y-dir of the original grid before starting to solve");
                param.declare_entry("Horizontal division of cathode Channel",
                                    "4",
                                    Patterns::Integer(),
                                    "Number of cells we want horizontally of the original grid before starting to solve");
                param.declare_entry("Horizontal division of cathode GDL",
                                    "4",
                                    Patterns::Integer(),
                                    "Number of cells we want horizontally in the cathode GDL layer");
                param.declare_entry("Horizontal division of cathode CL",
                                    "1",
                                    Patterns::Integer(),
                                    "Number of cells we want horizontally in the cathode CL layer");
                param.declare_entry("Horizontal division of cathode MPL",
                                    "1",
                                    Patterns::Integer(),
                                    "Number of cells we want horizontally in the cathode MPL layer");
                param.declare_entry("Horizontal division of membrane",
                                    "2",
                                    Patterns::Integer(),
                                    "Number of cells we want horizontally in the membrane layer");
                
                param.declare_entry("Horizontal division of anode MPL",
                                    "1",
                                    Patterns::Integer(),
                                    "Number of cells we want horizontally in the anode MPL layer");
                param.declare_entry("Horizontal division of anode CL",
                                    "1",
                                    Patterns::Integer(),
                                    "Number of cells we want horizontally in the anode CL layer");
                param.declare_entry("Horizontal division of anode GDL",
                                    "4",
                                    Patterns::Integer(),
                                    "Number of cells we want horizontally in the anode GDL layer");   
            }
            param.leave_subsection();
            
            param.enter_subsection("Material ID");
            {
                param.declare_entry("Test",// Current collector material_id
                                    "1",
                                    Patterns::Integer(),
                                    "GridTest material id");
                
                param.declare_entry("Cathode current collector",// Current collector material_id
                                    "0",
                                    Patterns::Integer(),
                                    "Current collector material_id");
                param.declare_entry("Cathode gas channel",      // Cathode gas channel material_id
                                    "11",
                                    Patterns::Integer(),
                                    "Cathode gas channel material_id");
                param.declare_entry("Cathode GDL",              // Cathode gas diffusion layer material_id
                                    "2",
                                    Patterns::Integer(),
                                    "Cathode gas diffusion layer material_id");
                param.declare_entry("Cathode MPL",              // Cathode microporous layer  material_id
                                    "3",
                                    Patterns::Integer(),
                                    "Cathode microporous layer material_id");
                param.declare_entry("Cathode CL",               // Cathode catalyst layer  material_id
                                    "4",
                                    Patterns::List(Patterns::Double()),
                                    //Patterns::Integer(),
                                    "Cathode catalyst layer material_id");
                
                param.declare_entry("Membrane",                 // Membrane material_id
                                    "5",
                                    Patterns::Integer(),
                                    "Membrane material_id");
                
                param.declare_entry("Anode MPL",                // Anode microporous layer material_id
                                    "7",
                                    Patterns::Integer(),
                                    "Anode microporous layer material_id");                
                param.declare_entry("Anode CL",                 // Anode catalyst layer material_id
                                    "6",
                                    Patterns::Integer(),
                                    "Anode catalyst layer material_id");
                param.declare_entry("Anode GDL",                // Anode gas diffusion layer material_id
                                    "8",
                                    Patterns::Integer(),
                                    "Anode gas diffusion layer material_id");
                param.declare_entry("Anode gas channel",        // Anode gas channel material_id
                                    "9",
                                    Patterns::Integer(),
                                    "Anode gas channel material_id");
                param.declare_entry("Anode current collector",  // Anode current collector material_id
                                    "10",
                                    Patterns::Integer(),
                                    "Anode current collector material_id");
                
                param.declare_entry("Agglomerate Material ID",
                                    "1",
                                    Patterns::Integer());
                param.declare_entry("Thin Film Material ID",
                                    "2",
                                    Patterns::Integer());
                
            }
            param.leave_subsection();
            
            param.enter_subsection("Boundary ID");
            {
                param.declare_entry("c_Ch_Inlet",
                                    "100",
                                    Patterns::Integer(),
                                    "Cathode gas channel inlet boundary_id");
                param.declare_entry("c_Ch_Outlet",
                                    "100",
                                    Patterns::Integer(),
                                    "Cathode gas channel outlet boundary_id");
                param.declare_entry("c_Ch_Base",
                                    "100",
                                    Patterns::Integer(),
                                    "Cathode gas channel entrance/exit base boundary_id");
                param.declare_entry("c_Ch_Roof",
                                    "100",
                                    Patterns::Integer(),
                                    "Cathode gas channel entrance/exit roof boundary_id");
                param.declare_entry("c_Ch/GDL",
                                    "100",
                                    Patterns::Integer(),
                                    "Cathode gas channel and gas diffusion layer boundary_id");
                param.declare_entry("c_BPP/GDL",
                                    "200",
                                    Patterns::Integer(),
                                    "Cathode bipolar plates and gas diffusion layer boundary_id");
                param.declare_entry("c_GDL/CL",
                                    "255",
                                    Patterns::Integer(),
                                    "Cathode gas diffusion layer and catalyst layer boundary_id");
                param.declare_entry("c_CL/Membrane",
                                    "255",
                                    Patterns::Integer(),
                                    "Cathode catalyst layer and membrane boundary_id");
                param.declare_entry("c_GDL/MPL",
                                    "255",
                                    Patterns::Integer(),
                                    "Cathode gas diffusion layer and microporous layer boundary_id");
                param.declare_entry("c_MPL/CL",
                                    "255",
                                    Patterns::Integer(),
                                    "Cathode microporous layer and catalyst layer boundary_id");
                param.declare_entry("Membrane/a_CL",
                                    "255",
                                    Patterns::Integer(),
                                    "Membrane and anode catalyst layer boundary_id");
                param.declare_entry("a_CL/GDL",
                                    "255",
                                    Patterns::Integer(),
                                    "Anode catalyst layer and gas diffusion layer boundary_id");
                param.declare_entry("a_CL/MPL",
                                    "255",
                                    Patterns::Integer(),
                                    "Anode catalyst layer and microporous layer boundary_id");
                param.declare_entry("a_MPL/GDL",
                                    "255",
                                    Patterns::Integer(),
                                    "Anode microporous layer and gas diffusion layer boundary_id");
                
                param.declare_entry("a_GDL/BPP",
                                    "8",
                                    Patterns::Integer(),
                                    "Anode gas diffusion layer and bipolar plate boundary_id");
                param.declare_entry("a_GDL/Ch",
                                    "7",
                                    Patterns::Integer(),
                                    "Anode gas diffusion layer and channel boundary_id");

                param.declare_entry("a_Ch_Inlet",
                                    "100",
                                    Patterns::Integer(),
                                    "Anode gas channel inlet boundary_id");
                param.declare_entry("a_Ch_Outlet",
                                    "100",
                                    Patterns::Integer(),
                                    "Anode gas channel outlet boundary_id");
                param.declare_entry("a_Ch_Base",
                                    "100",
                                    Patterns::Integer(),
                                    "Anode gas channel entrance/exit base boundary_id");
                param.declare_entry("a_Ch_Roof",
                                    "100",
                                    Patterns::Integer(),
                                    "Anode gas channel entrance/exit roof boundary_id");
                
                param.declare_entry("Agglomerate-Thin Film Boundary ID",
                                    "255",
                                    Patterns::Integer());                
                param.declare_entry("Thin Film Boundary ID",
                                    "0",
                                    Patterns::Integer());
                
            }
            param.leave_subsection();
            
            
        }
        param.leave_subsection();
        
    }
    param.leave_subsection();
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::GridBase<dim>::initialize(ParameterHandler& param)
{
    
    param.enter_subsection("Grid generation");
    {
        mesh_type_name= param.get("Type of mesh");
        if (mesh_type_name.compare("GridExternal") != 0)
            read_from_file = false;
        else 
            read_from_file = true;
        mesh_name = param.get("File name");
        mesh_type = param.get("File type");
        num_refine = param.get_integer("Initial refinement");
        
        param.enter_subsection("Internal mesh generator parameters");
        {
            
            param.enter_subsection("Dimensions");
            {
                l_land_c = 0.5*param.get_double("Cathode current collector width [cm]");
                l_channel_c = 0.5*param.get_double("Cathode channel width [cm]");
                l_channelLength_c = param.get_double("Cathode channel length [cm]");
                l_channelHeight_c = param.get_double("Cathode channel height [cm]");
                l_channelEntrance_c = param.get_double("Cathode channel entrance length [cm]");
                l_channelExit_c = param.get_double("Cathode channel exit length [cm]");
                l_gdl_c = param.get_double("Cathode GDL thickness [cm]");
                l_cat_c = FcstUtilities::string_to_number<double>( Utilities::split_string_list( param.get("Cathode CL thickness [cm]")));
                l_mpl_c = param.get_double("Cathode MPL thickness [cm]");
                l_mem = param.get_double("Membrane thickness [cm]");      
                l_cat_a = param.get_double("Anode CL thickness [cm]");
                l_mpl_a = param.get_double("Anode MPL thickness [cm]");               
                l_gdl_a = param.get_double("Anode GDL thickness [cm]");              
                l_land_a = 0.5*param.get_double("Anode current collector width [cm]");
                l_channel_a = 0.5*param.get_double("Anode channel width [cm]");
                
                if (this->l_land_c + this->l_channel_c != this->l_land_a + this->l_channel_a)
                {
                    FcstUtilities::log<<"Error: The sum of the current collector and channel in anode "
                    <<"and cathode should be the same. Otherwise, symmetry cannot "
                    <<"be guaranteed. Changing the value of the anode channel"
                    <<"width"<<std::endl;
                    this->l_channel_a = this->l_land_c + this->l_channel_c - this->l_land_a;
                }
                
                r_agg = param.get_double("Agglomerate radius [cm]");
                delta_agg = param.get_double("Thin film thickness [cm]");
                
                l_cube = param.get_double("HyperCube edge");
            }
            param.leave_subsection();
            
            param.enter_subsection("Mesh refinement parameters");
            {
                num_vert = param.get_integer("Initial vertical cell count");
                num_c_Channel = param.get_integer("Horizontal division of cathode Channel");
                num_c_GDL = param.get_integer("Horizontal division of cathode GDL");
                num_c_CL = param.get_integer("Horizontal division of cathode CL");
                num_c_MPL = param.get_integer("Horizontal division of cathode MPL");
                num_membrane = param.get_integer("Horizontal division of membrane");
                num_a_MPL = param.get_integer("Horizontal division of anode MPL");
                num_a_CL = param.get_integer("Horizontal division of anode CL");
                num_a_GDL = param.get_integer("Horizontal division of anode GDL");
            }
            param.leave_subsection();
            
            param.enter_subsection("Material ID");
            {
                test_mid = param.get_integer("Test");
                
                a_CC_mid = param.get_integer("Cathode current collector");
                c_GC_mid  = param.get_integer("Cathode gas channel");
                c_GDL_mid = param.get_integer("Cathode GDL");
                c_MPL_mid = param.get_integer("Cathode MPL");
                c_CL_mid  = FcstUtilities::string_to_number<unsigned int>( Utilities::split_string_list( param.get("Cathode CL")));                
                
                membrane_mid = param.get_integer("Membrane");
                
                a_CL_mid = param.get_integer("Anode CL");
                a_MPL_mid = param.get_integer("Anode MPL");               
                a_GDL_mid = param.get_integer("Anode GDL");
                a_GC_mid  = param.get_integer("Anode gas channel");
                a_CC_mid = param.get_integer("Anode current collector");
                
                r_agg_mid = param.get_integer("Agglomerate Material ID");
                delta_agg_mid = param.get_integer("Thin Film Material ID");
            }
            param.leave_subsection();
            
            param.enter_subsection("Boundary ID");
            {
                c_Ch_Inlet_bid = param.get_integer("c_Ch_Inlet");
                c_Ch_Outlet_bid = param.get_integer("c_Ch_Outlet");
                c_Ch_base_bid = param.get_integer("c_Ch_Base");
                c_Ch_roof_bid = param.get_integer("c_Ch_Roof");
                
                c_Ch_GDL_bid = param.get_integer("c_Ch/GDL");
                c_BPP_GDL_bid = param.get_integer("c_BPP/GDL");
                c_GDL_CL_bid = param.get_integer("c_GDL/CL");
                c_CL_Membrane_bid = param.get_integer("c_CL/Membrane");
                c_GDL_MPL_bid = param.get_integer("c_GDL/MPL");
                c_MPL_CL_bid = param.get_integer("c_MPL/CL");              
                c_CL_Membrane_bid = param.get_integer("c_CL/Membrane");
                
                a_Membrane_CL_bid = param.get_integer("Membrane/a_CL");
                a_CL_MPL_bid = param.get_integer("a_CL/MPL");
                a_MPL_GDL_bid = param.get_integer("a_MPL/GDL");               
                a_CL_GDL_bid = param.get_integer("a_CL/GDL");
                a_GDL_BPP_bid = param.get_integer("a_GDL/BPP");
                a_GDL_Ch_bid = param.get_integer("a_GDL/Ch");

                a_Ch_Inlet_bid = param.get_integer("a_Ch_Inlet");
                a_Ch_Outlet_bid = param.get_integer("a_Ch_Outlet");
                a_Ch_base_bid = param.get_integer("a_Ch_Base");
                a_Ch_roof_bid = param.get_integer("a_Ch_Roof");
                
                r_delta_bid = param.get_integer("Agglomerate-Thin Film Boundary ID");
                delta_bid = param.get_integer("Thin Film Boundary ID");
            }
            param.leave_subsection();
        }
        param.leave_subsection();
    }
    param.leave_subsection();
}

/*
//---------------------------------------------------------------------------
template<int dim>
void
NAME::GridBase<dim>::set_grid(Triangulation<dim>& triangulation,
                              bool& read_in_initial_solution)
{
    
  // If not reading from a grid from a previous solution, 
  //    check whether we can use one from external source
  if (read_in_initial_solution || read_from_file)
  {
    // TODO: Read grid_file name and format from parameter handler
    std::ifstream grid_file("._mesh_to_be_transfered.msh");
    if(!grid_file.fail())
    {
        FcstUtilities::log << "Reading grid from file" << std::endl;
        GridIn<dim> grid_in;
        grid_in.attach_triangulation(triangulation);
        grid_in.read_msh(grid_file);
        grid_file.close();
        read_from_file = true;
    }
    else
    {
        // if grid fails to read, generate standard grid, don't read initial solution
        FcstUtilities::log << "Grid file not available for reading" << std::endl;
        read_from_file = false;
        read_in_initial_solution = false;
    }
  }
  
  //if(!read_grid_from_file)
 // {
  //  FcstUtilities::log << "Generating new grid" << std::endl;  
  //  grid.generate_grid(*this->tr);
 // }
  
}
*/
//---------------------------------------------------------------------------
// --- Zingan ---
/*
template<int dim>
void
GridBase<dim>::generate_salome_grid(Triangulation<dim>&                triangulation,
                                    const types::boundary_id&          bdry_id,
                                    boost::shared_ptr< Boundary<dim> > boundary) const
{
  FcstUtilities::log << "Reading Salome grid from a UNV file" << std::endl;

  GridIn<dim> grid_in;
  grid_in.attach_triangulation(triangulation);
  std::ifstream input_file(salome_mesh_name.c_str());
  grid_in.read_unv(input_file);
  input_file.close();

  if( boundary )
      triangulation.set_boundary(bdry_id,
                                *boundary);

  triangulation.refine_global(this->num_refine);
}
*/


//---------------------------------------------------------------------------
template <int dim>
std::vector<unsigned int>
NAME::GridBase<dim>::get_material_id (const std::string material) const
{
    
    std::vector<unsigned int> material_IDs;
    
    if (material.compare("Cathode gas channel") == 0)
    {
        material_IDs.push_back(c_GC_mid);
    }    
    else if (material.compare("Cathode GDL") == 0)
    {
        material_IDs.push_back(c_GDL_mid);
    }
    else if (material.compare("Cathode MPL") == 0)
    {
        material_IDs.push_back(c_MPL_mid);
    }
    else if (material.compare("Cathode CL") == 0)
    {   
        // cathode is the only composite layer
        material_IDs = c_CL_mid;
    }                
    else if (material.compare("Membrane") == 0)
    {   
        material_IDs.push_back(membrane_mid);
    }
    else if (material.compare("Anode CL") == 0)
    {
        material_IDs.push_back(a_CL_mid);
    }
    else if (material.compare("Anode MPL") == 0)
    {
        material_IDs.push_back(a_MPL_mid);
    }
    else if (material.compare("Anode GDL") == 0)
    {   
        material_IDs.push_back(a_GDL_mid);
    }
    else if (material.compare("Anode gas channel") == 0)
    {
        material_IDs.push_back(a_GC_mid);
    }
    else if (material.compare("Agglomerate Material ID") == 0)
    {   
        material_IDs.push_back(r_agg_mid);
    }
    else if (material.compare("Thin Film Material ID") == 0)
    {   
        material_IDs.push_back(delta_agg_mid);
    }
    else
    {
        AssertThrow(true, ExcMessage("The material_id requested in GridBase<dim>::get_material_id is not defined"));
    }
    
    return material_IDs;
    
}


//---------------------------------------------------------------------------
template <int dim>
unsigned int
NAME::GridBase<dim>::get_boundary_id (const std::string boundary) const
{
    if (boundary == "c_Ch_Inlet")
        return this-> c_Ch_Inlet_bid;
    else if (boundary == "c_Ch_Outlet")
        return this-> c_Ch_Outlet_bid;
    else if (boundary == "c_Ch_Base")
        return this-> c_Ch_base_bid;
    else if (boundary == "c_Ch_Roof")
        return this-> c_Ch_roof_bid;
    else if (boundary == "c_Ch/GDL")
        return this->c_Ch_GDL_bid;
    else if (boundary == "c_BPP/GDL")
        return this->c_BPP_GDL_bid;
    else if (boundary == "c_GDL/CL")
        return this->c_GDL_CL_bid;
    else if (boundary == "c_GDL/MPL")
        return this->c_GDL_MPL_bid;
    else if (boundary == "c_MPL/CL")
        return this->c_MPL_CL_bid;
    else if (boundary == "c_CL/Membrane")
        return this->c_CL_Membrane_bid;
    else if (boundary == "Membrane/a_CL")
        return this->a_Membrane_CL_bid;
    else if (boundary == "a_CL/MPL")
        return this->a_CL_MPL_bid;
    else if (boundary == "a_MPL/GDL")
        return this->a_MPL_GDL_bid;
    else if (boundary == "a_CL/GDL")
        return this->a_CL_GDL_bid; 
    else if (boundary == "a_GDL/BPP")
        return this->a_GDL_BPP_bid;
    else if (boundary == "a_GDL/Ch")
        return this->a_GDL_Ch_bid;
    else if (boundary == "a_Ch_Inlet")
        return this->a_Ch_Inlet_bid;
    else if (boundary == "a_Ch_Outlet")
        return this->a_Ch_Outlet_bid;
    else if (boundary == "a_Ch_Base")
        return this->a_Ch_base_bid;
    else if (boundary == "a_Ch_Roof")
        return this->a_Ch_roof_bid;
    
    else if (boundary == "Agglomerate-Thin Film Boundary ID")
        return this->r_delta_bid; 
    else if (boundary == "Thin Film Boundary ID")
        return this->delta_bid;
    else
    {
        FcstUtilities::log<<"The boundary_id requested is not defined"<<std::endl;
        abort();
    }
}

//------------------------------------------------------------
// --- Zingan ---
template <int dim>
unsigned int
NAME::GridBase<dim>::get_boundary_id (ParameterHandler& param,
                                      const std::string boundary) const
{
  unsigned int id;
  
  param.enter_subsection("Mesh generation");
    {
        param.enter_subsection("Internal mesh generator parameters");
        {
            param.enter_subsection("Boundary ID");
            {
                id = param.get_integer(boundary);
                
            }
            param.leave_subsection();
            
        }
        param.leave_subsection();
  }
  param.leave_subsection();

  return id;
}

//-------------------------------------------------------------
template <int dim>
void 
NAME::GridBase<dim>::refine_area (Triangulation<dim> &triangulation,
                                  const unsigned int material_id)
{
    // Simply loop the cells looking for the material_id given and set refine flags
    typename Triangulation<dim>::active_cell_iterator
    cell = triangulation.begin_active(),
    endc = triangulation.end();
    for (; cell!=endc; ++cell)
    {
        if (cell->material_id() == material_id)
            cell->set_refine_flag();
    }
    triangulation.execute_coarsening_and_refinement();
}


//-------------------------------------------------------------
template <int dim>
void 
NAME::GridBase<dim>::output_grid (const Triangulation<dim> &triangulation,
                                 const std::string filename) const
{
    std::ofstream out;
    out.open(filename.c_str());
    GridOut grid_out;
    grid_out.write_eps (triangulation, out);
}

//------------------------------------------------------------------
template <int dim>
void
NAME::GridBase<dim>::print_material_id_and_boundary_id(const Triangulation<dim> &triangulation) const
{  
    #ifdef _1D_
    FcstUtilities::log << "Function " << __FUNCTION__ << " not defined in " << deal_II_dimension << "d \n";
    abort();
    #endif
    
    #ifdef _2D_
    //Used to test that the triangularization has the different boundary flags assigned.
    typename Triangulation<dim>::active_cell_iterator
    cellch = triangulation.begin_active(),
    endc = triangulation.end();
    int count = 0;
    for (; cellch!=endc; ++cellch)
    {
        count++;
        FcstUtilities::log << "Cell: " << count << "  M_id: " << int(cellch->material_id())
        << "  B1_id " << int((cellch->line(0))->boundary_indicator())
        << "  B2_id " << int((cellch->line(1))->boundary_indicator())
        << "  B3_id " << int((cellch->line(2))->boundary_indicator())
        << "  B4_id " << int((cellch->line(3))->boundary_indicator()) << "\n";
    }
    #endif
    
    #ifdef _3D_
    FcstUtilities::log << "Function " << __FUNCTION__ << " not defined in " << deal_II_dimension << "d \n";
    abort();
    #endif
}

//-------------------------------------------------------------
// Explicit instantiations.
template class NAME::GridBase<deal_II_dimension>;
