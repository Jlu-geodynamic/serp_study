/*
  Copyright (C) 2011 - 2017 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include "tag_composition.h"

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

namespace aspect
{
  namespace MeshRefinement
  {
    template <int dim>
    void
    TagComposition<dim>::tag_additional_cells() const
    {
      if (this->get_dof_handler().n_dofs() != 0)
        {

          AssertThrow (this->n_compositional_fields() >= 1,
                       ExcMessage ("This refinement criterion can not be used when no "
                                   "compositional fields are active!"));

          QGauss<dim> quadrature (this->get_fe().base_element(this->introspection().base_elements.compositional_fields).degree+1);

          FEValues<dim> fe_values (this->get_mapping(),
                                   this->get_fe(),
                                   quadrature,
                                   update_quadrature_points | update_values);

          // the values of the compositional fields are stored as blockvectors for each field
          // we have to extract them in this structure
          std::vector<std::vector<double> > prelim_composition_values
          (this->n_compositional_fields(),
           std::vector<double> (quadrature.size()));

          typename DoFHandler<dim>::active_cell_iterator
          cell = this->get_dof_handler().begin_active(),
          endc = this->get_dof_handler().end();
          for (; cell!=endc; ++cell)
            if (cell->is_locally_owned())
              {
                bool coarsen = false;
                bool refine = false;
                bool clear_refine = false;
                bool clear_coarsen = false;

                bool sed1_present = false;
                bool ly_1_present = false;
                bool ly_2_present = false;
                bool ly_3_present = false;
                bool ly_4_present = false;
				bool ly_5_present = false;
				bool melt_present = false;
				bool new_crust_present = false;

                fe_values.reinit(cell);

                for (unsigned int c = 0; c<this->n_compositional_fields(); c++)
                  {
                    fe_values[this->introspection().extractors.compositional_fields[c]].get_function_values (this->get_solution(),
                        prelim_composition_values[c]);
                  }


                for (unsigned int p=0; p<quadrature.size(); ++p)
                  {
                    if (this->introspection().compositional_name_exists("sediment_1") && prelim_composition_values[sed1_refinement[0]][p] > 0.2)
                      {
                        sed1_present = true;
                        //Crust will have smallest res, so not interested in other fields
                        break;
                      }
                    if ((this->introspection().compositional_name_exists("sediment_init") || this->introspection().compositional_name_exists("upper")) && prelim_composition_values[ly_1_refinement[0]][p] > 0.2)
					  {
						ly_1_present = true;
						break;
					  }
                    if ((this->introspection().compositional_name_exists("salt") || this->introspection().compositional_name_exists("lower")) && prelim_composition_values[ly_2_refinement[0]][p] > 0.2)
					  {
						ly_2_present = true;
						break;
					  }
					
                    if ((this->introspection().compositional_name_exists("crust") || this->introspection().compositional_name_exists("mantle_upper")) && prelim_composition_values[ly_3_refinement[0]][p] > 0.2)
					  {
						ly_3_present = true;
						break;
					  }
					if (this->introspection().compositional_name_exists("new_crust") && prelim_composition_values[new_crust_refinement[0]][p] > 0.2)
					  {
						new_crust_present = true;
						break;
					  }
					if (this->introspection().compositional_name_exists("melt_fraction") && prelim_composition_values[melt_refinement[0]][p] > 5e-6)
                        melt_present = true;
                    if ((this->introspection().compositional_name_exists("mantle_L") || this->introspection().compositional_name_exists("mantle_middle")) && prelim_composition_values[ly_4_refinement[0]][p] > 0.2)
                        ly_4_present = true;
					if ((this->introspection().compositional_name_exists("mantle_D") || this->introspection().compositional_name_exists("mantle_lower")) && prelim_composition_values[ly_5_refinement[0]][p] > 0.2)
						ly_5_present = true;
					
                  }


                //Only continue if at least one is true

                    int maximum_refinement_level = max_level;
                    int minimum_refinement_level = min_level;

                    if (sed1_present)
                      {
                        minimum_refinement_level = sed1_refinement[1];
                        maximum_refinement_level = sed1_refinement[2];
                      }
                    else if (ly_1_present)
                      {
                        minimum_refinement_level = ly_1_refinement[1];
                        maximum_refinement_level = ly_1_refinement[2];
                      }
                    else if (ly_2_present)
                      {
                        minimum_refinement_level = ly_2_refinement[1];
                        maximum_refinement_level = ly_2_refinement[2];
                      }
                    else if (ly_3_present)
                      {
                        minimum_refinement_level = ly_3_refinement[1];
                        maximum_refinement_level = ly_3_refinement[2];
                      }
                    else if (ly_4_present)
                      {
                        minimum_refinement_level = ly_4_refinement[1];
                        maximum_refinement_level = ly_4_refinement[2];
                      }
					else if (ly_5_present)
                      {
                        minimum_refinement_level = ly_5_refinement[1];
                        maximum_refinement_level = ly_5_refinement[2];
                      }
					else if (melt_present)
					{
						minimum_refinement_level = melt_refinement[1];
                        maximum_refinement_level = melt_refinement[2];
					}
					else if (new_crust_present)
					{
						minimum_refinement_level = new_crust_refinement[1];
                        maximum_refinement_level = new_crust_refinement[2];
					}

                    const int cell_level = cell->level();
                    if (cell_level >= maximum_refinement_level)
                      {
                        clear_refine = true;
                      }
                    if (cell_level >  maximum_refinement_level)
                      {
                        coarsen = true;
                      }
                    if (cell_level <= minimum_refinement_level)
                      {
                        clear_coarsen = true;
                      }
                    if (cell_level < minimum_refinement_level)
                      {
                        refine = true;
                      }

                    if (clear_refine)
                      cell->clear_refine_flag ();
                    if (clear_coarsen)
                      cell->clear_coarsen_flag ();
                    if (refine)
                      cell->set_refine_flag ();
                    if (coarsen)
                      cell->set_coarsen_flag ();

              }
        }
    }

    template <int dim>
    void
    TagComposition<dim>::
    declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Mesh refinement");
      {
        prm.enter_subsection("Composition");
        {
          prm.declare_entry("Sediment refinement","",
                            Patterns::List (Patterns::Integer(0)),
                            "The compositional field number of the crust, its minimum refinement level and "
                            "its maximum refinement level.");
          prm.declare_entry("Layer 1 refinement","",
                            Patterns::List (Patterns::Integer(0)),
                            "The compositional field number of the slab mantle, its minimum refinement level and "
                            "its maximum refinement level.");
          prm.declare_entry("Layer 2 refinement","",
                            Patterns::List (Patterns::Integer(0)),
                            "The compositional field number of the slab mantle, its minimum refinement level and "
                            "its maximum refinement level.");
          prm.declare_entry("Layer 3 refinement","",
                            Patterns::List (Patterns::Integer(0)),
                            "The compositional field number of the slab mantle, its minimum refinement level and "
                            "its maximum refinement level.");
          prm.declare_entry("Layer 4 refinement","",
                            Patterns::List (Patterns::Integer(0)),
                            "The compositional field number of the slab mantle, its minimum refinement level and "
                            "its maximum refinement level.");
		  prm.declare_entry("Layer 5 refinement","",
                            Patterns::List (Patterns::Integer(0)),
                            "The compositional field number of the slab mantle, its minimum refinement level and "
                            "its maximum refinement level.");
		  prm.declare_entry("Melt refinement","",
                            Patterns::List (Patterns::Integer(0)),
                            "The compositional field number of the slab mantle, its minimum refinement level and "
                            "its maximum refinement level.");
		  prm.declare_entry("New crust refinement","",
                            Patterns::List (Patterns::Integer(0)),
                            "The compositional field number of the slab mantle, its minimum refinement level and "
                            "its maximum refinement level.");

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    TagComposition<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Mesh refinement");
      {
        min_level = prm.get_integer("Minimum refinement level");
        max_level = prm.get_integer("Initial adaptive refinement") + prm.get_integer("Initial global refinement");
        prm.enter_subsection("Composition");
        {

          const std::vector<int> sed1
            = Utilities::string_to_int(
                Utilities::split_string_list(prm.get("Sediment refinement")));

          sed1_refinement = std::vector<unsigned int> (sed1.begin(),sed1.end());

          AssertThrow (sed1_refinement.size() == 3,
                       ExcMessage ("The number of refinement data given here must be "
                                   "equal to 3 (field number + min level + max level). "));

          AssertThrow (sed1_refinement[0] < this->n_compositional_fields(),
                       ExcMessage ("The number of compositional field to refine (starting "
                                   "from 0) should be smaller than the number of fields. "));

          AssertThrow (sed1_refinement[1] >= min_level,
                       ExcMessage ("The minimum refinement for the crust cannot be "
                                   "smaller than the minimum level of the whole model. "));

          AssertThrow (sed1_refinement[2] <= max_level,
                       ExcMessage ("The maximum refinement for the crust cannot be "
                                   "greater than the maximum level of the whole model. "));

          const std::vector<int> ly_1
            = Utilities::string_to_int(
                Utilities::split_string_list(prm.get("Layer 1 refinement")));

          ly_1_refinement = std::vector<unsigned int> (ly_1.begin(),
                                                              ly_1.end());

          AssertThrow (ly_1_refinement.size() == 3,
                       ExcMessage ("The number of refinement data given here must be "
                                   "equal to 3 (field number + min level + max level). "));

          AssertThrow (ly_1_refinement[0] < this->n_compositional_fields(),
                       ExcMessage ("The number of compositional field to refine (starting "
                                   "from 0) should be smaller than the number of fields. "));

          AssertThrow (ly_1_refinement[1] >= min_level,
                       ExcMessage ("The minimum refinement for the slab mantle cannot be "
                                   "smaller than the minimum level of the whole model. "));

          AssertThrow (ly_1_refinement[2] <= max_level,
                       ExcMessage ("The maximum refinement for the slab mantle cannot be "
                                   "greater than the maximum level of the whole model. "));
        

          const std::vector<int> ly_2
            = Utilities::string_to_int(
                Utilities::split_string_list(prm.get("Layer 2 refinement")));

          ly_2_refinement = std::vector<unsigned int> (ly_2.begin(),
                                                              ly_2.end());

          AssertThrow (ly_2_refinement.size() == 3,
                       ExcMessage ("The number of refinement data given here must be "
                                   "equal to 3 (field number + min level + max level). "));

          AssertThrow (ly_2_refinement[0] < this->n_compositional_fields(),
                       ExcMessage ("The number of compositional field to refine (starting "
                                   "from 0) should be smaller than the number of fields. "));

          AssertThrow (ly_2_refinement[1] >= min_level,
                       ExcMessage ("The minimum refinement for the slab mantle cannot be "
                                   "smaller than the minimum level of the whole model. "));

          AssertThrow (ly_2_refinement[2] <= max_level,
                       ExcMessage ("The maximum refinement for the slab mantle cannot be "
                                   "greater than the maximum level of the whole model. "));
        


          const std::vector<int> ly_3
            = Utilities::string_to_int(
                Utilities::split_string_list(prm.get("Layer 3 refinement")));

          ly_3_refinement = std::vector<unsigned int> (ly_3.begin(),
                                                              ly_3.end());

          AssertThrow (ly_3_refinement.size() == 3,
                       ExcMessage ("The number of refinement data given here must be "
                                   "equal to 3 (field number + min level + max level). "));

          AssertThrow (ly_3_refinement[0] < this->n_compositional_fields(),
                       ExcMessage ("The number of compositional field to refine (starting "
                                   "from 0) should be smaller than the number of fields. "));

          AssertThrow (ly_3_refinement[1] >= min_level,
                       ExcMessage ("The minimum refinement for the slab mantle cannot be "
                                   "smaller than the minimum level of the whole model. "));

          AssertThrow (ly_3_refinement[2] <= max_level,
                       ExcMessage ("The maximum refinement for the slab mantle cannot be "
                                   "greater than the maximum level of the whole model. "));
        


          const std::vector<int> ly_4
            = Utilities::string_to_int(
                Utilities::split_string_list(prm.get("Layer 4 refinement")));

          ly_4_refinement = std::vector<unsigned int> (ly_4.begin(),
                                                              ly_4.end());

          AssertThrow (ly_4_refinement.size() == 3,
                       ExcMessage ("The number of refinement data given here must be "
                                   "equal to 3 (field number + min level + max level). "));

          AssertThrow (ly_4_refinement[0] < this->n_compositional_fields(),
                       ExcMessage ("The number of compositional field to refine (starting "
                                   "from 0) should be smaller than the number of fields. "));

          AssertThrow (ly_4_refinement[1] >= min_level,
                       ExcMessage ("The minimum refinement for the slab mantle cannot be "
                                   "smaller than the minimum level of the whole model. "));

          AssertThrow (ly_4_refinement[2] <= max_level,
                       ExcMessage ("The maximum refinement for the slab mantle cannot be "
                                   "greater than the maximum level of the whole model. "));
								   
		
		  const std::vector<int> ly_5
            = Utilities::string_to_int(
                Utilities::split_string_list(prm.get("Layer 5 refinement")));

          ly_5_refinement = std::vector<unsigned int> (ly_5.begin(),
                                                              ly_5.end());

          AssertThrow (ly_5_refinement.size() == 3,
                       ExcMessage ("The number of refinement data given here must be "
                                   "equal to 3 (field number + min level + max level). "));

          AssertThrow (ly_5_refinement[0] < this->n_compositional_fields(),
                       ExcMessage ("The number of compositional field to refine (starting "
                                   "from 0) should be smaller than the number of fields. "));

          AssertThrow (ly_5_refinement[1] >= min_level,
                       ExcMessage ("The minimum refinement for the slab mantle cannot be "
                                   "smaller than the minimum level of the whole model. "));

          AssertThrow (ly_5_refinement[2] <= max_level,
                       ExcMessage ("The maximum refinement for the slab mantle cannot be "
                                   "greater than the maximum level of the whole model. "));	  
		  
		  
		  const std::vector<int> melt
            = Utilities::string_to_int(
                Utilities::split_string_list(prm.get("Melt refinement")));

          melt_refinement = std::vector<unsigned int> (melt.begin(),
                                                              melt.end());

          AssertThrow (melt_refinement.size() == 3,
                       ExcMessage ("The number of refinement data given here must be "
                                   "equal to 3 (field number + min level + max level). "));

          AssertThrow (melt_refinement[0] < this->n_compositional_fields(),
                       ExcMessage ("The number of compositional field to refine (starting "
                                   "from 0) should be smaller than the number of fields. "));

          AssertThrow (melt_refinement[1] >= min_level,
                       ExcMessage ("The minimum refinement for the slab mantle cannot be "
                                   "smaller than the minimum level of the whole model. "));

          AssertThrow (melt_refinement[2] <= max_level,
                       ExcMessage ("The maximum refinement for the slab mantle cannot be "
                                   "greater than the maximum level of the whole model. "));
			
          const std::vector<int> new_crust
            = Utilities::string_to_int(
                Utilities::split_string_list(prm.get("New crust refinement")));
				
		  new_crust_refinement = std::vector<unsigned int> (new_crust.begin(),
                                                              new_crust.end());
		  AssertThrow (new_crust_refinement.size() == 3,
                       ExcMessage ("The number of refinement data given here must be "
                                   "equal to 3 (field number + min level + max level). "));

          AssertThrow (new_crust_refinement[0] < this->n_compositional_fields(),
                       ExcMessage ("The number of compositional field to refine (starting "
                                   "from 0) should be smaller than the number of fields. "));

          AssertThrow (new_crust_refinement[1] >= min_level,
                       ExcMessage ("The minimum refinement for the slab mantle cannot be "
                                   "smaller than the minimum level of the whole model. "));

          AssertThrow (new_crust_refinement[2] <= max_level,
                       ExcMessage ("The maximum refinement for the slab mantle cannot be "
                                   "greater than the maximum level of the whole model. "));
        }

        prm.leave_subsection();
      }
      prm.leave_subsection();
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MeshRefinement
  {
    ASPECT_REGISTER_MESH_REFINEMENT_CRITERION(TagComposition,
                                              "tag composition",
                                              "A mesh refinement criterion that "
                                              "(de)flags cells for refinement and coarsening "
                                              "based on what composition is present. Different max "
                                              "and min refinement levels are set for the mantle, the crustal "
                                              "field, the slab mantle and the overriding plate. ")
  }
}
