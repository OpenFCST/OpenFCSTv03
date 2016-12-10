//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2009-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: geometries.cc
//    - Description: Geometry definition for several fuel cell elements.
//    - Developers: L. Birkett, P. Dobson, M. Secanell and V. Zingan
//    - $Id: geometries.cc 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#include <grid/geometries.h>

namespace NAME = FuelCellShop::Geometry;

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//----------------------------- External mesh -------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
template <int dim>
const std::string NAME::GridExternal<dim>::concrete_name ("GridExternal");

template <int dim>
NAME::GridExternal<dim> const* NAME::GridExternal<dim>::PROTOTYPE = new NAME::GridExternal<dim>();

//---------------------------------------------------------------------------
template <int dim>
void
NAME::GridExternal<dim>::generate_grid( Triangulation<dim>& triangulation) 
{
    FcstUtilities::log << "Read mesh from file: " << std::endl; 
    FcstUtilities::log << "\t mesh name: " << this->mesh_name << std::endl;
    FcstUtilities::log << "\t mesh type: " << this->mesh_type << std::endl;
    
    
    //GridIn<dim> grid_in;
    this->grid_in->attach_triangulation(triangulation);
    
    // File reading version. Does not work yet
    // grid_in.read(  mesh_name, mesh_type );
    // Read from stream
    std::ifstream input_file( this->mesh_name.c_str() );
    this->grid_in->read( input_file ,  this->grid_in->parse_format( this->mesh_type ) );
    input_file.close();
    
    triangulation.refine_global(this->num_refine);
 }

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//----------------------------- Cube -------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
template <int dim>
const std::string NAME::HyperCube<dim>::concrete_name ("HyperCube");

template <int dim>
NAME::HyperCube<dim> const* NAME::HyperCube<dim>::PROTOTYPE = new NAME::HyperCube<dim>();

//---------------------------------------------------------------------------
template <int dim>
void
NAME::HyperCube<dim>::generate_grid( Triangulation<dim>& triangulation) 
{
    GridGenerator::hyper_cube (triangulation, -1.0*this->l_cube, this->l_cube);
    triangulation.refine_global(this->num_refine);
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//----------------------------- Cathdode w/ MPL -----------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
template <int dim>
const std::string NAME::CathodeMPL<dim>::concrete_name ("CathodeMPL");

template <int dim>
NAME::CathodeMPL<dim> const* NAME::CathodeMPL<dim>::PROTOTYPE = new NAME::CathodeMPL<dim>();

//---------------------------------------------------------------------------
template <int dim>
void
NAME::CathodeMPL<dim>::generate_grid(Triangulation<dim> &triangulation) 
{
    
    #ifdef _1D_
    FcstUtilities::log << "Function " << __FUNCTION__ << " not defined in " << deal_II_dimension << "d \n";
    abort();
    #endif
    
    #ifdef _3D_
    FcstUtilities::log << "Function " << __FUNCTION__ << " not defined in " << deal_II_dimension << "d \n";
    abort();
    #endif
    
    FcstUtilities::log << "Using FCST grid generator." << std::endl;
    
    #ifdef _2D_
    // First need to create a vector containing cells in x and y directions
    std::vector< std::vector <double> > mesh;
    std::vector<double> part;
    for(unsigned int j = 0; j < this->l_cat_c.size(); j++)
    {
        for(unsigned int i = 0; i < this->num_c_CL; i++)
        {
            if(this->l_cat_c.at(j) > 0)
                part.push_back((this->l_cat_c.at(j))/this->num_c_CL);
        }
    }
    for(unsigned int i = 0; i < this->num_c_MPL; i++)
    {
        if(this->l_mpl_c > 0)
            part.push_back((this->l_mpl_c)/this->num_c_MPL);
    }
    for(unsigned int i = 0; i < this->num_c_GDL; i++)
    {
        if(this->l_gdl_c > 0)
            part.push_back((this->l_gdl_c)/this->num_c_GDL);
    }
    mesh.push_back(part);
    part.clear();
    for(unsigned int i = 0; i < this->num_vert; i++)
    {
        part.push_back((this->l_land_c + this->l_channel_c)/this->num_vert);
    }
    mesh.push_back(part);
    
    // Now create the mesh using the lengths to make the upper corner point
    GridGenerator::subdivided_hyper_rectangle(triangulation, mesh, Point<2> (0,0), Point<2> (this->l_gdl_c + std::accumulate(this->l_cat_c.begin(),this->l_cat_c.end(),0.0) + this->l_mpl_c, this->l_land_c + this->l_channel_c), false);
    
    // Loop over the cells checking to see if the cell vertices are close to the FC BPP/Ch boundary
    // Move these vertices to the FC boundary
    typename Triangulation<dim>::active_cell_iterator
    cell = triangulation.begin_active(),
    endc = triangulation.end();
    typename Triangulation<dim>::active_line_iterator
    line;
    //point variable to store cell vertices as we loop over them
    Point<dim> vert;
    //calculate current cell height
    const double initial_cell_height = (this->l_land_c + this->l_channel_c)/this->num_vert;
    
    for (; cell!=endc; ++cell)
    {
        for (unsigned int i = 0; i < (GeometryInfo<dim>::vertices_per_cell); i++)
        {
            vert = cell->vertex(i);
            //now check to see if current vertex is close to current collector/gas channel boundary
            if (std::abs(vert(1) - this->l_land_c - 0.00001) <= initial_cell_height/2)
            {
                vert(1) = this->l_land_c;
                cell->vertex(i) = vert;
            }
        }
        // Set material id
        vert = cell->center();
        if(vert(0) < std::accumulate(this->l_cat_c.begin(),this->l_cat_c.end(), 0.0))
        {
            for(unsigned int j = 0; j < this->l_cat_c.size(); j++)
            {
                if (vert(0) > std::accumulate(this->l_cat_c.begin(),this->l_cat_c.begin() + j, 0.0) && 
                    vert(0) <= std::accumulate(this->l_cat_c.begin(),this->l_cat_c.begin() + j + 1, 0.0))
                    cell->set_material_id(this->c_CL_mid.at(j));
            }
        }
        else if ((std::accumulate(this->l_cat_c.begin(),this->l_cat_c.end(),0.0) + this->l_mpl_c)>=vert(0))
            cell->set_material_id(this->c_MPL_mid);
        else
            cell->set_material_id(this->c_GDL_mid);
        
        //loop through lines of the cell and change boundary_ids as necessary
        
        for (unsigned int j=0; j < GeometryInfo<dim>::lines_per_cell; j++)
        {
            line = cell->line(j);
            vert = line->vertex(0);
            if (vert(0) == 0.0)
            {
                vert = line->vertex(1);
                if (vert(0) == 0.0)
                {
                    line->set_boundary_indicator(this->c_CL_Membrane_bid);
                    continue;
                }
                vert = line->vertex(0);
            }
            for(unsigned int j = 0; j < this->l_cat_c.size(); j++)
                {
                if (vert(0) == std::accumulate(this->l_cat_c.begin(),this->l_cat_c.begin() + j,0.0))
                {
                    vert = line->vertex(1);
                    if (vert(0) == std::accumulate(this->l_cat_c.begin(),this->l_cat_c.begin() + j,0.0))
                    {
                        line->set_boundary_indicator(this->c_MPL_CL_bid);                                                         // CCL|CCL or c_MPL_CL_bid
                        continue;
                    }
                    vert = line->vertex(0);
                }
                }
            if (std::abs(vert(0) - std::accumulate(this->l_cat_c.begin(),this->l_cat_c.end(),0.0) - this->l_mpl_c) < 0.00001)
            {
                vert = line->vertex(1);
                if (std::abs(vert(0) - std::accumulate(this->l_cat_c.begin(),this->l_cat_c.end(),0.0) - this->l_mpl_c) < 0.00001)
                {
                    line->set_boundary_indicator(this->c_GDL_MPL_bid);
                    continue;
                }
                vert = line->vertex(0);
            }
            if (std::abs(vert(0) - std::accumulate(this->l_cat_c.begin(),this->l_cat_c.end(),0.0) - this->l_mpl_c - this->l_gdl_c) < 0.00001)
            {
                if (std::abs(vert(0) - std::accumulate(this->l_cat_c.begin(),this->l_cat_c.end(),0.0) - this->l_mpl_c - this->l_gdl_c) < 0.00001)
                {
                    if (vert(1) < this->l_land_c)
                        line->set_boundary_indicator(this->c_BPP_GDL_bid);
                    else
                        line->set_boundary_indicator(this->c_Ch_GDL_bid);
                    continue;
                }
            }
        }
    }
    #endif
    
    
    triangulation.refine_global(this->num_refine);
}


//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//----------------------------- PEMFC w/ MPL --------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
template <int dim>
const std::string NAME::PemfcMPL<dim>::concrete_name ("PemfcMPL");

template <int dim>
NAME::PemfcMPL<dim> const* NAME::PemfcMPL<dim>::PROTOTYPE = new NAME::PemfcMPL<dim>();

//---------------------------------------------------------------------------
template <int dim>
void
NAME::PemfcMPL<dim>::generate_grid(Triangulation<dim> &triangulation) 
{
    
    #ifdef _1D_
    FcstUtilities::log << "Function " << __FUNCTION__ << " not defined in " << deal_II_dimension << "d \n";
    abort();
    #endif
    
    #ifdef _3D_
    FcstUtilities::log << "Function " << __FUNCTION__ << " not defined in " << deal_II_dimension << "d \n";
    abort();
    #endif
    
    FcstUtilities::log << "Using OpenFCST grid generator::PemfcMPL." << std::endl;
    
    #ifdef _2D_
    // First need to create a vector containing cells in x and y directions
    std::vector< std::vector <double> > mesh;
    std::vector<double> part;
    for(unsigned int i = 0; i < this->num_a_GDL; i++)
    {
        if(this->l_gdl_a > 0)
            part.push_back((this->l_gdl_a)/this->num_a_GDL);
    }
    for(unsigned int i = 0; i < this->num_a_MPL; i++)
    {
        if(this->l_mpl_a > 0)
            part.push_back((this->l_mpl_a)/this->num_a_MPL);
    }
    for(unsigned int i = 0; i < this->num_a_CL; i++)
    {
        if(this->l_cat_a > 0)
            part.push_back((this->l_cat_a)/this->num_a_CL);
    }
    for(unsigned int i = 0; i < this->num_membrane; i++)
    {
        part.push_back((this->l_mem)/this->num_membrane);
    }
    for(unsigned int j = 0; j < this->l_cat_c.size(); j++)
    {
        for(unsigned int i = 0; i < this->num_c_CL; i++)
        {
            if(this->l_cat_c.at(j) > 0)
                part.push_back((this->l_cat_c.at(j))/this->num_c_CL);
        }
    }
    for(unsigned int i = 0; i < this->num_c_MPL; i++)
    {
        if(this->l_mpl_c > 0)
            part.push_back((this->l_mpl_c)/this->num_c_MPL);
    }
    for(unsigned int i = 0; i < this->num_c_GDL; i++)
    {
        if(this->l_gdl_c > 0)
            part.push_back((this->l_gdl_c)/this->num_c_GDL);
    }
    mesh.push_back(part);
    part.clear();
    for(unsigned int i = 0; i <this-> num_vert; i++)
    {
        part.push_back((this->l_land_c + this->l_channel_c)/this->num_vert);
    }
    mesh.push_back(part);
    
    // Now create the mesh using the lengths to make the upper corner point
    GridGenerator::subdivided_hyper_rectangle(triangulation, mesh, Point<2> (0,0),
                                              Point<2> (this->l_gdl_a + this->l_mpl_a + this->l_cat_a + this->l_mem + std::accumulate(this->l_cat_c.begin(),this->l_cat_c.end(),0.0) + this->l_mpl_c + this->l_gdl_c, this->l_land_c + this->l_channel_c), false);
    
    // Loop over the cells checking to see if the cell vertices are close to the FC BPP/Ch boundary
    // Move these vertices to the FC boundary
    typename Triangulation<dim>::active_cell_iterator
    cell = triangulation.begin_active(),
    endc = triangulation.end();
    typename Triangulation<dim>::active_line_iterator line;
    //point variable to store cell vertices as we loop over them
    Point<dim> vert;
    //calculate current cell width
    const double initial_cell_height = (this->l_land_c + this->l_channel_c)/this->num_vert;
    
    for (; cell!=endc; ++cell)
    {
        for (unsigned int i = 0; i < (GeometryInfo<dim>::vertices_per_cell); i++)
        {
            //if plate boundaries are far apart move cells to match plates
            if (std::abs(this->l_land_c - this->l_land_a) > initial_cell_height)
            {
                vert = cell->vertex(i);
                //now check to see if current vertex is close to cathode current collector/gas channel boundary
                if (std::abs(vert(1) - this->l_land_c - 0.0001) <= initial_cell_height/2)
                {
                    vert(1) = this->l_land_c;
                    cell->vertex(i) = vert;
                }
                //now check if currrent vertex is close to anode current collector/gas channel boundary
                else if (std::abs(vert(1) - this->l_land_a - 0.0001) <= initial_cell_height/2)
                {
                    vert(1) = this->l_land_a;
                    cell->vertex(i) = vert;
                }
            }
            //if boundaries are close keep boundary on same cell division
            else if (std::abs(this->l_land_c - this->l_land_a) < initial_cell_height/2)
            {
                vert = cell->vertex(i);
                //now check to see if current vertex is close to cathode current collector/gas channel boundary
                if (std::abs(vert(1) - this->l_land_c - 0.0001) <= initial_cell_height/2)
                {
                    if ((this->l_gdl_c + this->l_mpl_c + std::accumulate(this->l_cat_c.begin(),this->l_cat_c.end(),0.0)) > vert(0))                                           // go over statement.
                        vert(1) = this->l_land_c;
                    else if ((this->l_gdl_c + this->l_mpl_c + std::accumulate(this->l_cat_c.begin(),this->l_cat_c.end(),0.0) + this->l_mem) < vert(0))
                        vert(1) = this->l_land_a;
                    else
                        vert(1) = this->l_land_c + (vert(0) - (this->l_gdl_c + this->l_mpl_c + std::accumulate(this->l_cat_c.begin(),this->l_cat_c.end(),0.0)))*(this->l_land_a - this->l_land_c)/this->l_mem;
                    cell->vertex(i) = vert;
                }
                
            }
            //if boundaries are in between those two criterion we need to refine cells first and then go back to move boundaries
            else
            {
                vert = cell->vertex(i);
                //now check to see if current vertex is close to cathode current collector/gas channel boundary
                if ((std::abs(vert(1) - this->l_land_c - 0.0001) <= initial_cell_height/2) &&
                    (std::abs(vert(1) - this->l_land_a - 0.0001) > initial_cell_height/2) &&
                    (vert(1) != this->l_land_c) && (i < GeometryInfo<dim>::vertices_per_cell/2))
                {
                    vert(1) = this->l_land_c;
                    cell->vertex(i) = vert;
                    vert = cell->vertex(i+2);
                    vert(1) = this->l_land_a;
                    cell->vertex(i+2) = vert;
                }
                else if ((std::abs(vert(1) - this->l_land_a - 0.0001) <= initial_cell_height/2) &&
                         (std::abs(vert(1) - this->l_land_c - 0.0001) > initial_cell_height/2) &&
                    (vert(1) != this->l_land_a) && (i < GeometryInfo<dim>::vertices_per_cell/2))
                {
                    vert(1) = this->l_land_a;
                    cell->vertex(i) = vert;
                    vert = cell->vertex(i+2);
                    vert(1) = this->l_land_c;
                    cell->vertex(i+2) = vert;
                }
                else if ((std::abs(vert(1) - this->l_land_c - 0.0001) <= initial_cell_height/2) &&
                         (std::abs(vert(1) - this->l_land_a - 0.0001) <= initial_cell_height/2) &&
                    (vert(1) != this->l_land_c) && (vert(1) != this->l_land_a) && (i < GeometryInfo<dim>::vertices_per_cell/2))
                {
                    if (this->l_land_c > this->l_land_a)
                    {
                        vert(1) = this->l_land_a;
                        cell->vertex(i) = vert;
                        vert = cell->vertex(i+2);
                        vert(1) = this->l_land_c;
                        cell->vertex(i+2) = vert;
                    }
                    else
                    {
                        vert(1) = this->l_land_c;
                        cell->vertex(i) = vert;
                        vert = cell->vertex(i+2);
                        vert(1) = this->l_land_a;
                        cell->vertex(i+2) = vert;
                    }
                }
            }
        }
        
        // Set material id based on the center of the cell
        vert = cell->center();
        if (vert(0) >= 0 &&
            vert(0) <= this->l_gdl_a)
            cell->set_material_id(this->a_GDL_mid);
        else if (vert(0) > this->l_gdl_a &&
                 vert(0) <= this->l_gdl_a + this->l_mpl_a)
            cell->set_material_id(this->a_MPL_mid);
        else if (vert(0) > this->l_gdl_a + this->l_mpl_a &&
                 vert(0) <= (this->l_gdl_a + this->l_mpl_a + this->l_cat_a))
            cell->set_material_id(this->a_CL_mid);
        else if (vert(0) > (this->l_gdl_a + this->l_mpl_a + this->l_cat_a) &&
                 vert(0) <= (this->l_gdl_a + this->l_mpl_a + this->l_cat_a + this->l_mem ))
            cell->set_material_id(this->membrane_mid);
        else if(vert(0) > this->l_gdl_a + this->l_mpl_a + this->l_cat_a + this->l_mem  &&
                vert(0) <= (this->l_gdl_a + this->l_mpl_a + this->l_cat_a + this->l_mem + std::accumulate(this->l_cat_c.begin(),this->l_cat_c.end(),0.0)))
        {
            for(unsigned int j = 0; j < this->l_cat_c.size(); j++)
            {
                if (vert(0) > this->l_gdl_a + this->l_mpl_a + this->l_cat_a + this->l_mem + std::accumulate(this->l_cat_c.begin(),this->l_cat_c.begin() + j, 0.0) && 
                    vert(0) <= this->l_gdl_a + this->l_mpl_a + this->l_cat_a + this->l_mem + std::accumulate(this->l_cat_c.begin(),this->l_cat_c.begin() + j + 1, 0.0))
                    cell->set_material_id(this->c_CL_mid.at(j));
            }
        }
        else if (vert(0) > this->l_gdl_a + this->l_mpl_a + this->l_cat_a + this->l_mem  + std::accumulate(this->l_cat_c.begin(),this->l_cat_c.end(),0.0) &&
                 vert(0) <= this->l_gdl_a + this->l_mpl_a + this->l_cat_a + this->l_mem + std::accumulate(this->l_cat_c.begin(),this->l_cat_c.end(),0.0) + this->l_mpl_c)
            cell->set_material_id(this->c_MPL_mid);
        else
            cell->set_material_id(this->c_GDL_mid);
        
        //loop through lines of the cell and change boundary_ids as necessary
        
        for (unsigned int j=0; j < GeometryInfo<dim>::lines_per_cell; j++)
        {
            line = cell->line(j);
            vert = line->vertex(0);
            if (line->vertex(0)(0) == 0.0 &&
                line->vertex(1)(0) == 0.0)
            {
                if (line->vertex(0)(1) <= this->l_land_a &&
                    line->vertex(1)(1) <= this->l_land_a)
                    line->set_boundary_indicator(this->a_GDL_BPP_bid);
                else
                    line->set_boundary_indicator(this->a_GDL_Ch_bid);
                continue;
            }
            else if (line->vertex(0)(0) == this->l_gdl_a &&
                     line->vertex(1)(0) == this->l_gdl_a)
            {
                line->set_boundary_indicator(this->a_MPL_GDL_bid);
                continue;
            }
            else if (std::abs(line->vertex(0)(0) - this->l_gdl_a - this->l_mpl_a) < 0.00001 &&
                     std::abs(line->vertex(1)(0) - this->l_gdl_a - this->l_mpl_a) < 0.00001)
            {
                line->set_boundary_indicator(this->a_CL_MPL_bid);
                continue;
            }
            else if (std::abs(line->vertex(0)(0) - this->l_gdl_a - this->l_mpl_a - this->l_cat_a) < 0.00001 &&
                     std::abs(line->vertex(1)(0) - this->l_gdl_a - this->l_mpl_a - this->l_cat_a) < 0.00001)
            {
                line->set_boundary_indicator(this->a_Membrane_CL_bid);
                continue;
            }
            else if (std::abs(line->vertex(0)(0) - this->l_gdl_a - this->l_mpl_a - this->l_cat_a - this->l_mem) < 0.00001 &&
                     std::abs(line->vertex(1)(0) - this->l_gdl_a - this->l_mpl_a - this->l_cat_a - this->l_mem) < 0.00001)
            {
                line->set_boundary_indicator(this->c_CL_Membrane_bid);
                continue;
            }
            else if (std::abs(line->vertex(0)(0) - this->l_gdl_a - this->l_mpl_a - this->l_cat_a - this->l_mem - std::accumulate(this->l_cat_c.begin(),this->l_cat_c.end(),0.0)) < 0.00001 &&
                     std::abs(line->vertex(1)(0) - this->l_gdl_a - this->l_mpl_a - this->l_cat_a - this->l_mem - std::accumulate(this->l_cat_c.begin(),this->l_cat_c.end(),0.0)) < 0.00001)
            {
                for(unsigned int j = 0; j < this->l_cat_c.size(); j++)
                { 
                    if (std::abs(line->vertex(0)(0) - this->l_gdl_a - this->l_mpl_a - this->l_cat_a - this->l_mem - std::accumulate(this->l_cat_c.begin(),this->l_cat_c.begin() + j,0.0)) < 0.00001 &&
                        std::abs(line->vertex(1)(0) - this->l_gdl_a - this->l_mpl_a - this->l_cat_a - this->l_mem - std::accumulate(this->l_cat_c.begin(),this->l_cat_c.begin() + j,0.0)) < 0.00001)
                    {
                        line->set_boundary_indicator(this->c_MPL_CL_bid);                           // CCL|CCL or c_MPL_CL_bid
                        continue;  
                    }
                }   
            }
            else if (std::abs(line->vertex(0)(0) - this->l_gdl_a - this->l_mpl_a - this->l_cat_a - this->l_mem - std::accumulate(this->l_cat_c.begin(),this->l_cat_c.end(),0.0) - this->l_mpl_c) < 0.00001 &&
                     std::abs(line->vertex(1)(0) - this->l_gdl_a - this->l_mpl_a - this->l_cat_a - this->l_mem - std::accumulate(this->l_cat_c.begin(),this->l_cat_c.end(),0.0) - this->l_mpl_c) < 0.00001)
            {
                line->set_boundary_indicator(this->c_GDL_MPL_bid);
                continue;
            }
            else if (std::abs(line->vertex(0)(0) - this->l_gdl_a - this->l_mpl_a - this->l_cat_a - this->l_mem - std::accumulate(this->l_cat_c.begin(),this->l_cat_c.end(),0.0) - this->l_mpl_c - this->l_gdl_c) < 0.00001 &&
                     std::abs(line->vertex(1)(0) - this->l_gdl_a - this->l_mpl_a - this->l_cat_a - this->l_mem - std::accumulate(this->l_cat_c.begin(),this->l_cat_c.end(),0.0) - this->l_mpl_c - this->l_gdl_c) < 0.00001)
            {
                if (line->vertex(0)(1) <= this->l_land_c &&
                    line->vertex(1)(1) <= this->l_land_c)
                    line->set_boundary_indicator(this->c_BPP_GDL_bid);
                else
                    line->set_boundary_indicator(this->c_Ch_GDL_bid);
                continue;
            }
        }
    }
    #endif
    
    
    triangulation.refine_global(this->num_refine);
}


//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//----------------------------- AGGLOMERATE ---------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
template <int dim>
const std::string NAME::Agglomerate<dim>::concrete_name ("Agglomerate");

template <int dim>
NAME::Agglomerate<dim> const* NAME::Agglomerate<dim>::PROTOTYPE = new NAME::Agglomerate<dim>();

//---------------------------------------------------------------------------
template <int dim>
void
NAME::Agglomerate<dim>::generate_grid(Triangulation<dim> &triangulation) 
{
    
    #ifdef _1D_
    FcstUtilities::log << "Function " << __FUNCTION__ << " not defined in " << deal_II_dimension << "d \n";
    abort();
    #endif
    FcstUtilities::log << "Using FCST grid generator." << std::endl;
    this->center(true); // Set the center to the origin
    double domain_radius = (this->r_agg + this->delta_agg)/this->r_agg; // radius - normalized
    
    // Create a circle at 'center' with 'radius'
    // By default creates 5 cells in 2D - 1 sq at the center + 4 trapezoid
    GridGenerator::hyper_ball (triangulation,this->center,domain_radius);
    
    // Define a circular boundary at the radius.
    // Otherwise, refinements will be done assuming the outer boundary is a straight line
    static const HyperBallBoundary<dim> agglomerate_radius(this->center, domain_radius);
    triangulation.set_boundary (0, agglomerate_radius);
    
    reset_material_ids(triangulation);
    
    
    triangulation.refine_global(this->num_refine);
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::Agglomerate<dim>::reset_material_ids(Triangulation<dim> &triangulation) const
{
    // Set the material ID's based on the agglomerate radius and thin film thickness.
    typename Triangulation<dim>::active_cell_iterator
    cell = triangulation.begin_active(),
    endc = triangulation.end();
    // points to store the coordinates of the vertex
    Point<dim> vert;
    
    // Set the material_id's and the boundary_id's for all the cells
    for (; cell!=endc; ++cell) // Loop over all cells and vertices
    {
        double max_dist = 0.; // Variable to store the maximum distance of the vertex from the center
        for (unsigned int i=0; i < GeometryInfo<dim>::vertices_per_cell; ++i)
        {
            vert = cell->vertex(i);
            double vert_dist = vert.distance(this->center); //calculate the distance between the vertex and the center
            if (vert_dist > max_dist)
                max_dist = vert_dist;
        }
        // If any of the vertices of the cell are larger than the agglomerate radius, it belongs to the thin film
        if (max_dist > 1.0)
        {
            // Set the material_id for the thin film
            cell->set_material_id(this->delta_agg_mid);
        }
        else
        {
            // Set all other cells to the porous agglomerate initially
            cell->set_material_id(this->r_agg_mid);
        }
    }
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//----------------------------- External mesh -------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
template <int dim>
const std::string NAME::GridTest<dim>::concrete_name ("GridTest");

template <int dim>
NAME::GridTest<dim> const* NAME::GridTest<dim>::PROTOTYPE = new NAME::GridTest<dim>();

//---------------------------------------------------------------------------
template <int dim>
void
NAME::GridTest<dim>::generate_grid( Triangulation<dim>& triangulation) 
{
    // Iterator to loop over lines:
    typename Triangulation<dim>::active_line_iterator line;
    // Point variable to store cell vertices as we loop over them
    Point<dim> vert;
    
    // deal.II grid generator:
    GridGenerator::hyper_cube(triangulation);
    
    typename Triangulation<dim>::active_cell_iterator
    cell = triangulation.begin_active(),
    endc = triangulation.end();
    
    for (; cell!=endc; ++cell)
    {
        cell->set_material_id(this->test_mid);

        //loop through lines of the cell and change boundary_ids as necessary
        for (unsigned int j=0; j < GeometryInfo<dim>::lines_per_cell; j++)
        {
            line = cell->line(j);
            
            vert = line->vertex(0);
            if (vert(0) == 0.0)
            {
                vert = line->vertex(1);
                if (vert(0) == 0.0)
                {
                    line->set_boundary_indicator(1);
                    continue;
                }
                vert = line->vertex(0);
            }
            if (vert(0) == 1.0)
            {
                vert = line->vertex(1);
                if (vert(0) == 1.0)
                {
                    line->set_boundary_indicator(2);
                    continue;
                }
                vert = line->vertex(0);
            }            
        }
    }
    //triangulation.refine_global(this->num_refine);
}

//-------------------------------------------------------------
// Explicit instantiations.
template class NAME::GridExternal<deal_II_dimension>;
template class NAME::HyperCube<deal_II_dimension>;

template class NAME::CathodeMPL<deal_II_dimension>;
template class NAME::PemfcMPL<deal_II_dimension>;
template class NAME::Agglomerate<deal_II_dimension>;

template class NAME::GridTest<deal_II_dimension>;