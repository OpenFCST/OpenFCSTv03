//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2006-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: geometry.h
//    - Description: Used to store and output information for general fuel cell geometries
//    - Developers: M. Secanell and Peter Dobson
//                  Salome parser dealii::GridIn<dim, spacedim>::read_unv(): Valentin N. Zingan, University of Alberta, (C) 2012.
//    - $Id: geometry.h 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#ifndef _FUELCELLSHOP__GEOMETRY_H
#define _FUELCELLSHOP__GEOMETRY_H

//------------------------------
// STANDARD LIBRARY DECLARATIONS
//------------------------------
#include <fstream>
#include <iostream>

//------------------------------
// DEAL.II DECLARATIONS
//------------------------------
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/exceptions.h>

#include <utils/fcst_utilities.h>

using namespace dealii;

namespace FuelCellShop
{
    namespace Geometry
    {

        /**
         * FuelCell Geometry information class
         *
         * \brief Class used to store the geometry and grid information for an application
         *
         * Derived classes should contain the specific implementation of
         * mesh generation and manipulation of the triangulation object.
         *
         * \note Concerning boundary id's:
         * 		All internal boundary id's must be set to 255 for deal.ii.
         * 		When implementing continuty over the domain, this is the
         * 		default parameter value deal.ii checks.  See the deal.ii
         * 		documentation for more information.
         *
         * <h3>Usage Details:</h3>
         *
         * GridBase is a virtual class that contains a map to all children applications. Using
         * the map, GridBase can pass the user a pointer to the desired class as specified in the input file.
         *
         * All children use the parameters in declare_GridGenerator_parameters; however, they have a re-implementation
         * of the function #generate_grid and other functions.
         *
         * To use the class, first create a pointer to the base object using:
         * @code
           boost::shared_ptr< FuelCellShop::Geometry::GridBase<dim> > grid;
         * @endcode
         *
         * Then, declare the data in the input file using:
         * @code
           FuelCellShop::Geometry::GridBase<dim>::declare_GridGenerator_parameters(param);
         * @endcode
         *
         * Finally, read the ParameterHandler object to initalize the grid. Then, use the objec to generate the desired grid:
         * @code
           grid = FuelCellShop::Geometry::GridBase<dim>::create_GridGenerator (param);
           grid->generate_grid(*this->tr);
         * @endcode
         *
         * \author M. Secanell, 2009-13
         * \author L. Brikett, 2009
         * \author P. Dobson, 2009-11
         */
        template <int dim>
        class GridBase
        {
        public:
            ///@name Instance Delivery (Public functions)
            /**
             * Function used to declare all the data necessary in the parameter files for
             * all GridBase children.
             *
             * @code
             *  section Grid generation
             *    +++++++ OPTIONS FROM APPFRAME +++++++++++
             *    set Sort by component = true
             *    set Sort Cuthill-McKee = false
             *    set Initial refinement = 1  # Number of times we want to globally refine the original grid before starting to solve
             *    set Refinement = global   # OPTIONS: global|adaptive|single|steps|circle
             *    set Refinement threshold = 0.3
             *    set Coarsening threshold = 0
             *    set Distort = 0
             *    +++++++ END OPTIONS FROM APPFRAME +++++++++++
             *    set Type of mesh = Cathode # OPTIONS: GridExternal | Cathode | Anode | CathodeMPL | Pemfc | PemfcMPL | Agglomerate -- Type of mesh that you would like to generate.")
             *    set File name = test.unv   # Name of the mesh file
             *    set File type = UNV        # Meshtype (see deal.ii supported mesh types
             *    ==== Section related to internal mesh generator: =======
             *    subsection Internal mesh generator parameters
             *      subsection Dimensions
             *       set Cathode current collector width [cm] = 0.1 # Thickness of the rib of the bipolar plates (BPP)
             *       set Cathode channel width [cm] = 0.1           # Thickness of the channel on the BPP
             *       set Cathode CL thickness [cm] = 1E-3           # Thickness of the cathode catalyst layer [cm]
             *       set Cathode GDL thickness [cm]  = 2E-2         # Thickness of the cathode gas diffusion layer [cm]
             *       set Cathode MPL thickness [cm] = 2E-3          # Thickness of the cathode microporous layer [cm]
             *       set Membrane thickness [cm] = 2E-2             # Thickness of the membrane [cm]
             *       set Anode MPL thickness [cm] = 2E-3            # Thickness of the anode microporous layer [cm]
             *       set Anode CL thickness [cm] = 1E-3             # Thickness of the anode catalyst layer [cm]
             *       set Anode GDL thickness [cm] = 2E-2            # Thickness of the anode gas diffusion layer [cm]
             *       set Anode current collector width [cm] = 0.1   # Thickness of the rib of the bipolar plates (BPP) [cm]
             *       set Anode channel width [cm] = 0.1             # Thickness of the channel on the BPP
             *       set Agglomerate radius [cm] = 100e-7           # Aggomerate radius (Used only for Agglomerata Type mesh)
             *       set Thin film thickness [cm] = 10e-7           # Agglomerate thin film (Used only for Agglomerata Type mesh)
             *     end
             *
             *     subsection Mesh refinement parameters
             *       set Initial vertical cell count = 6            # Number of cells we want in the y-dir of the original grid before starting to solve"
             *       set Horizontal division of cathode GDL = 4     # Number of cells we want horizontally in the cathode GDL layer
             *       set Horizontal division of cathode CL = 1      # Number of cells we want horizontally in the cathode CL layer
             *       set Horizontal division of cathode MPL = 1     # Number of cells we want horizontally in the cathode MPL layer
             *       set Horizontal division of membrane = 2        # Number of cells we want horizontally in the membrane layer
             *       set Horizontal division of anode MPL = 1       # Number of cells we want horizontally in the anode MPL layer
             *       set Horizontal division of anode CL = 1        # Number of cells we want horizontally in the anode CL layer
             *       set Horizontal division of anode GDL = 4       # Number of cells we want horizontally in the anode GDL layer
             *     end
             *
             *     subsection Material ID
             *       set Test = 1                                    # Material id for GridTest
             *       set Cathode current collector = 0               # Current collector material_id
             *       set Cathode gas channel = 1                     # Cathode gas channel material_id
             *       set Cathode GDL = 2                             # Cathode gas diffusion layer material_id
             *       set Cathode MPL = 3                             # Cathode microporous layer material_id
             *       set Cathode CL = 4                              # Cathode catalyst layer material_id
             *       set Membrane = 5                                # Membrane material_id
             *       set Anode MPL = 7                               # Anode microporous layer material_id
             *       set Anode CL = 6                                # Anode catalyst layer material_id
             *       set Anode GDL = 8                               # Anode gas diffusion layer material_id
             *       set Anode gas channel = 9                       # Anode gas channel material_id
             *       set Anode current collector = 10                # Anode current collector material_id
             *       set Agglomerate Material ID = 1
             *       set Thin Film Material ID = 2
             *     end
             *     subsection Boundary ID
             *       # NOTE: 255 defines an interior boundary conditions in deal.II. All internal boundaries
             *       # MUST have a 255 boundary id.
             *       set c_Ch/GDL = 100                              # Cathode gas channel and gas diffusion layer boundary_id
             *       set c_BPP/GDL = 200                             # Cathode bipolar plates and gas diffusion layer boundary_id
             *       set c_GDL/CL= 255                               # Cathode gas diffusion layer and catalyst layer boundary_id
             *       set c_CL/Membrane = 400                         # Cathode catalyst layer and membrane boundary_id
             *       set c_GDL/MPL = 255                             # Cathode gas diffusion layer and microporous layer boundary_id
             *       set c_MPL/CL = 255                              # Cathode microporous layer and catalyst layer boundary_id
             *       set c_CL/Membrane = 255                         # Cathode catalyst layer and membrane boundary_id
             *       set Membrane/a_CL = 255                         # Membrane and anode catalyst layer boundary_id
             *       set a_CL/GDL = 255                              # Anode catalyst layer and gas diffusion layer boundary_id
             *       set a_CL/MPL = 255                              # Anode catalyst layer and microporous layer boundary_id
             *       set a_MPL/GDL = 255                             # Anode microporous layer and gas diffusion layer boundary_id
             *       set a_GDL/BPP = 3                               # Anode gas diffusion layer and bipolar plate boundary_id
             *       set a_GDL/Ch = 4                                # Anode gas diffusion layer and channel boundary_id
             *       set Agglomerate-Thin Film Boundary ID = 255
             *       set Thin Film Boundary ID = 0
             *     end
             *   end
             * end
             * @endcode
             *
             */
            static void declare_GridGenerator_parameters (ParameterHandler &param);

            /**
             * Generate the appropriate mesh generator object based on the parameters in the input file
             */
            static boost::shared_ptr<FuelCellShop::Geometry::GridBase<dim> > create_GridGenerator (ParameterHandler &param)
            {
                boost::shared_ptr<FuelCellShop::Geometry::GridBase<dim> > pointer;
                std::string concrete_name;
                
                param.enter_subsection("Grid generation");
                {
                    concrete_name = param.get("Type of mesh");
                    FcstUtilities::log <<" "<<std::endl;
                    FcstUtilities::log << "name: "<<concrete_name.c_str()<<std::endl;
                }
                param.leave_subsection();
                //FcstUtilities::log <<FuelCellShop::Geometry::GridBase<dim>::get_mapFactory()->empty()<<std::endl;

                typename FuelCellShop::Geometry::GridBase<dim>::_mapFactory::iterator iterator = FuelCellShop::Geometry::GridBase<dim>::get_mapFactory()->find(concrete_name);
                if (iterator != FuelCellShop::Geometry::GridBase<dim>::get_mapFactory()->end())
                {
                    if (iterator->second)
                    {
                        pointer = iterator->second->create_replica(concrete_name);
                    }
                    else
                    {
                        FcstUtilities::log<<"Pointer not initialized"<<std::endl;
                        abort();
                    }
                }
                else
                {
                    FcstUtilities::log<<"Concrete name does not exist"<<std::endl;
                    abort();
                }
                
                pointer->initialize(param);
                
                return pointer;
            }
            //@}
            
            /**
             * Function is empty and must be reimplemented in derived classes.
             */
            virtual void generate_grid (Triangulation<dim>& ) 
            {
                const std::type_info& info = typeid ( *this );
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name()  << std::endl;
            }

            /**
             * This function is like the previous one \p generate_grid() but allows to
             * assign a curved boundary \p boundary to that with the boundary id
             * \p bdry_id.
             */
            void generate_grid_with_curved_boundaries(Triangulation<dim>&                triangulation,
                                                      const types::boundary_id&          bdry_id,
                                                      boost::shared_ptr< Boundary<dim> > boundary) const
            {
                   FcstUtilities::log << "Read mesh from file: " << std::endl;
                   FcstUtilities::log << "\t mesh name: " << this->mesh_name << std::endl;
                   FcstUtilities::log << "\t mesh type: " << this->mesh_type << std::endl;

                   GridIn<dim> grid_in;
                   grid_in.attach_triangulation(triangulation);
                   std::ifstream input_file( this->mesh_name.c_str() );
                   grid_in.read( input_file , grid_in.parse_format( this->mesh_type ) );
                   input_file.close();

                   triangulation.set_boundary(bdry_id,
                                             *boundary);

                   triangulation.refine_global(this->num_refine);
            }

            /** 
             * Return a vector with all material ids associated with string.
             * 
             * @note For composite layers this will be more than one value.
             */
            std::vector<unsigned int> get_material_id ( const std::string ) const;

            /**Return boundary id */
            unsigned int get_boundary_id ( ParameterHandler& param,
                                           const std::string boundary ) const;

            /**Return boundary id */
            unsigned int get_boundary_id ( const std::string ) const;

            /** Output the grid */
            void output_grid ( const Triangulation<dim> &triangulation, const std::string filename ) const;

            /**
             * Refine a certain part of the grid depending on material_ID given
             */
            void refine_area ( Triangulation<dim> &triangulation, const unsigned int material_id );

            /**
             * Return the type of mesh that is being used
             *
             * Options: GridExternal | CathodeMPL | PemfcMPL
             */
            std::string get_mesh_type()
            {
                return mesh_type_name;
            }
            ///@name Geometric data
            //@{
            /**
             * Return the Width of the cathode gas channel
             */
            inline double L_channel_c(){return 2.0*l_channel_c;}
            /**
             * Return the Length of the cathode gas channel
             */
            inline double L_channelLength_c(){return l_channelLength_c;}
            /**
             * Return the Height of the cathode gas channel
             */
            inline double L_channelHeight_c(){return l_channelHeight_c;}
            /**
             * Return the Entrance Length of the cathode gas channel
             */
            inline double L_channelEntrance_c(){return l_channelEntrance_c;}
            /**
             * Return the Exit Length of the cathode gas channel
             */
            inline double L_channelExit_c(){return l_channelExit_c;}
            /**
             * Return the Width of the cathode current collector
             */
            inline double L_land_c(){return 2.0*l_land_c;}
            /**
             * Return the Thickness of the cathode gas diffusion layer
             */
            inline double L_gdl_c(){return l_gdl_c;}
            /**
             * Return the thickness of the cathode catalyst layer
             */
            inline  std::vector<double> L_cat_c(){return l_cat_c;}

            /**
             * Return the Thickness of the cathode microporous layer
             */
            inline double L_mpl_c(){return l_mpl_c;}


            /**
             * Return the Thickness of the membrane
             */
            inline double L_mem(){return l_mem;}

            /**
             * Return the Thickness of the anode catalyst layer
             */
            inline double L_cat_a(){return l_cat_a;}
            /**
             * Return the Thickness of the anode microporous layer
             */
            inline double L_mpl_a(){return l_mpl_a;}
            /**
             * Return the Thickness of the anode gas diffusion layer
             */
            inline double L_gdl_a(){return l_gdl_a;}
            /**
             * Return the Width of the anode gas channel
             */
            inline double L_channel_a(){return 2.0*l_channel_a;}
            /**
             * Return the Width of the anode current collector
             */
            inline double L_land_a(){return 2.0*l_land_a;}
            //@}

            bool read_from_file;
            
            /**
             *  GridIn class object of dealii which allows to read various types of mesh formats.
             */
            boost::shared_ptr<GridIn<dim>> grid_in;
            

        protected:
             ///@name Constructors, destructor, and initalization
            //@{
            /** Constructor
             *
             */
            GridBase();

            /**
             * Destructor
             *
             */
            virtual ~GridBase();

            /**
             *
             */
            void initialize(ParameterHandler& param);

            //@}


            ///@name Instance Delivery (Types)
            //@{
            /**
             * This object is used to store all objects of type GasDiffusionLayer.
             */
            typedef std::map< std::string, GridBase<dim>* > _mapFactory;
            //@}

            ///@name Instance Delivery (Private and static)
            //@{
            /**
             *
             */
            static _mapFactory * get_mapFactory()
            {
                static _mapFactory mapFactory;
                return &mapFactory;
            }
            //@}
            ///@name Instance Delivery (Private functions)
            //@{
            /**
             * This member function is used to create an object of type gas diffusion layer
             *
             * \warning This class MUST be redeclared in every child.
             */
            virtual boost::shared_ptr<FuelCellShop::Geometry::GridBase<dim> > create_replica (const std::string &name)
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            }
            //@}
            ///@name Internal functions
            //@{
            /** */
            void print_material_id_and_boundary_id ( const Triangulation<dim> &triangulation ) const;
            //@}

            ///@name
            //@{
            /**
             * Specify if you would like to load a mesh from file or if you would like a module from the mesh generator
             * Options: File | Cathode | Anode | CathodeMPL | Pemfc | PemfcMPL
             */
            std::string mesh_type_name;
            /** Name of the mesh file */
            std::string mesh_name;
            /** Specify if it is a UNV file, MSH file, etc. */
            std::string mesh_type;
            /** Initial number of refinements */
            unsigned int num_refine;
            //@}

            ///@name  PEMFC Data (Dimensions)
            //@{

            /**
             * Width of the cathode gas channel
             */
            double l_channel_c;
            /**
             * Length of the cathode gas channel
             */
            double l_channelLength_c;
            /**
             * Height of the cathode gas channel
             */
            double l_channelHeight_c;
            /**
             * Entrance length of the cathode gas channel
             */
            double l_channelEntrance_c;
            /**
             * Exit Length of the cathode gas channel
             */
            double l_channelExit_c;
            /**
             * Width of the cathode current collector
             */
            double l_land_c;
            /**
             * Thickness of the cathode gas diffusion layer
             */
            double l_gdl_c;
            /**
             * Thickness of the cathode microporous layer
             */
            double l_mpl_c;
            /**
             * Thickness of the cathode catalyst layer
             */
            //double l_cat_c;								//##
            std::vector<double> l_cat_c;
            /**
             * Thickness of the membrane
             */
            double l_mem;
            /**
             * Thickness of the anode catalyst layer
             */
            double l_cat_a;
            /**
             * Thickness of the anode microporous layer
             */
            double l_mpl_a;
            /**
             * Thickness of the anode gas diffusion layer
             */
            double l_gdl_a;
            /**
             * Width of the anode gas channel
             */
            double l_channel_a;
            /**
             * Width of the anode current collector
             */
            double l_land_a;
            
            /** Cube edge for HyperCube mesh. Cube will go from -l_cube to l_cube */
            double l_cube;
             //@}

            ///@name Internal Mesh Generator Data
            //@{

            /**
             * Number of cells tall the initial grid
             */
            unsigned int num_vert;

            /**
             * Number of cells wide cathode channel
             */
            unsigned int num_c_Channel;
            /**
             * Number of cells wide cathode gas diffusion layer
             */
            unsigned int num_c_GDL;
            /**
             * Number of cells wide cathode microporous layer
             */
            unsigned int num_c_MPL;
            /**
             * Number of cells wide cathode catalyst layer
             */
            unsigned int num_c_CL;

            /**
             * Number of cells wide membrane layer
             */
            unsigned int num_membrane;
                        /**
             * Number of cells wide anode catalyst layer
             */
            unsigned int num_a_CL;
/**
             * Number of cells wide anode microporous layer
             */
            unsigned int num_a_MPL;
            /**
             * Number of cells wide anode gas diffusion layer
             */
            unsigned int num_a_GDL;

             //@}


            ///@name Material ID Data
            //@{
            /**
             * Material id for test cell (GridTest app)
             */
            unsigned int test_mid;
             /**
             * Material id cathode current collector
             */
            unsigned int c_CC_mid;
            /**
             * Material id cathode gas channel
             */
            unsigned int c_GC_mid;
                        /**
             * Material id cathode GDL
             */
            unsigned int c_GDL_mid;
            /**
             * Material id cathode MPL
             */
            unsigned int c_MPL_mid;
            /**
             * Material id cathode catalyst layer
             */
            std::vector<unsigned int> c_CL_mid;
            /**
             * Material id membrane
             */
            unsigned int membrane_mid;
            /**
             * Material id anode catalyst layer
             */
            unsigned int a_CL_mid;
            /**
             * Material id anode MPL
             */
            unsigned int a_MPL_mid;
            /**
             * Material id anode GDL
             */
            unsigned int a_GDL_mid;
            /**
             * Material id anode current collector
             */
            unsigned int a_CC_mid;
            /**
             * Material id anode gas channel
             */
            unsigned int a_GC_mid;


             //@}

            ///@name Boundary ID Data
            //@{
            /**
             * Boundary id cathode channel inlet
             */
            unsigned int c_Ch_Inlet_bid;    
            /**
             * Boundary id cathode channel outlet
             */
            unsigned int c_Ch_Outlet_bid;
            /**
             * Boundary id cathode channel entrance/exit base
             */
            unsigned int c_Ch_base_bid;
            /**
             * Boundary id cathode channel entrance/exit roof
             */
            unsigned int c_Ch_roof_bid;
            /**
             * Boundary id cathode channel and GDL
             */
            unsigned int c_Ch_GDL_bid;
            /**
             * Boundary id cathode BPP and GDL
             */
            unsigned int c_BPP_GDL_bid;
            /**
             * Boundary id cathode GDL and CL
             */
            unsigned int c_GDL_CL_bid;
            /**
             * Boundary id cathode CL and membrane
             */
            unsigned int c_CL_Membrane_bid;

            /**
             * Boundary id cathode GDL and MPL
             */
            unsigned int c_GDL_MPL_bid;
            /**
             * Boundary id cathode MPL and CL
             */
            unsigned int c_MPL_CL_bid;

            /**
             * Boundary id anode membrane and CL
             */
            unsigned int a_Membrane_CL_bid;
            /**
             * Boundary id anode CL and GDL
             */
            unsigned int a_CL_GDL_bid;
            /**
             * Boundary id anode GDL and BPP
             */
            unsigned int a_GDL_BPP_bid;
            /**
             * Boundary id anode GDL and channel
             */
            unsigned int a_GDL_Ch_bid;

            /**
             * Boundary id anode GDL and MPL
             */
            unsigned int a_MPL_GDL_bid;
            /**
             * Boundary id anode MPL and CL
             */
            unsigned int a_CL_MPL_bid;
            /**
             * Boundary id cathode channel inlet
             */
            unsigned int a_Ch_Inlet_bid;    
            /**
             * Boundary id cathode channel outlet
             */
            unsigned int a_Ch_Outlet_bid;
            /**
             * Boundary id cathode channel entrance/exit base
             */
            unsigned int a_Ch_base_bid;
            /**
             * Boundary id cathode channel entrance/exit roof
             */
            unsigned int a_Ch_roof_bid;
            //@}

            ///@name Agglomerate Data
            //@{
            /** Agglomerate radius */
            double r_agg;

            /** Electrolyte thin film thickness */
            double delta_agg;

            /** Porous Agglomerate Material ID*/
            unsigned int r_agg_mid;

            /** Electrolyte Thin Film Material ID */
            unsigned int delta_agg_mid;

            /** Boundary ID for agglomerate/thin film boundary */
            unsigned int r_delta_bid;

            /** Boundary ID for thin film external boundary */
            unsigned int delta_bid;

            /** Centre point of the circular/spherical domain */
            Point<dim> center;
            //@}
        };

    }//namespace Geometry

}//namespace FuelCell

#endif
