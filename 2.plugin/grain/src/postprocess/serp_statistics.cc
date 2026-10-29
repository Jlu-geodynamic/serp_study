/*
  Copyright (C) 2011 - 2019 by the authors of the ASPECT code.

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
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include </fs2/home/liuzhonglan/wy/lib_extra/melt20260814/grain/include/postprocess/serp_statistics.h>

#include <aspect/geometry_model/interface.h>
#include <aspect/geometry_model/box.h>
#include <aspect/utilities.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <iostream>
#include <dirent.h>
#include <sys/stat.h>

namespace aspect
{
  namespace Postprocess
  {

    template <int dim>
    std::pair<std::string,std::string>
    SerpStatistics<dim>::execute (TableHandler &statistics)
    {
      if (this->n_compositional_fields() == 0)
        return std::pair<std::string,std::string>();

      // create a quadrature formula based on the compositional element alone.
      // be defensive about determining that a compositional field actually exists
      AssertThrow (this->introspection().base_elements.compositional_fields
                   != numbers::invalid_unsigned_int,
                   ExcMessage("This postprocessor cannot be used without compositional fields."));

      const QGauss<dim> quadrature_formula (this->get_fe().base_element(this->introspection().base_elements.compositional_fields).degree+1);
      const unsigned int n_q_points = quadrature_formula.size();

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values   |
                               update_quadrature_points |
                               update_JxW_values);

	  MaterialModel::MaterialModelInputs<dim> in(fe_values.n_quadrature_points, this->n_compositional_fields());
	  double local_M_t = 0.;
      // compute the integral quantities by quadrature
      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          {
			//填充输入
            fe_values.reinit (cell);
			in.reinit(fe_values, cell, this->introspection(), this->get_solution());
			std::vector<double> d_serp_values(n_q_points, 0.);
			const unsigned int d_serp_id = this->introspection().compositional_index_for_name("d_serp");
		    fe_values[this->introspection().extractors.compositional_fields[d_serp_id]].get_function_values (this->get_solution(), d_serp_values);
			//20260806
			//获取洋壳数据，含洋壳的区域认为不发生蛇纹岩化
			std::vector<double> new_crust_values(n_q_points, 0.);
			const unsigned int new_crust_id = this->introspection().compositional_index_for_name("new_crust");
			fe_values[this->introspection().extractors.compositional_fields[new_crust_id]].get_function_values (this->get_solution(), new_crust_values);
			for (unsigned int q=0; q<n_q_points; ++q)
			{
				const double posx = in.position[q][0];
				//20260126 考虑到挤压背景下的边界问题，不统计左右边界区域的蛇纹岩化
				//20260806
				//获取洋壳数据，含洋壳的区域认为不发生蛇纹岩化
				if(posx > 20e3 && posx < max_x - 20e3 && d_serp_values[q] > 0 && new_crust_values[q] <= 0.1)
				{
					//所有地幔相放在一起计算
					//根据温度和参考密度估算当前区域的密度
					const double temperature = in.temperature[q];
					const double mantle_density_ref = densities[0];
					const double mantle_thermal_expansivity = thermal_expansivities[0];
					//默认温度小于Tp
					const double mantle_density_with_T = mantle_density_ref * (1 - (temperature - Tp) * mantle_thermal_expansivity);
					//面积 fe_values.JxW(q) / 1e6
					//20260127 将输出结果的单位修改为 kg/yr*m
					//local_M_t += d_serp_values[q] * mantle_density_with_T * fe_values.JxW(q) / 1e6;
					local_M_t += d_serp_values[q] * mantle_density_with_T * fe_values.JxW(q);
				}
			}
				
          }
		  
	  const double global_M_t = Utilities::MPI::sum (local_M_t, this->get_mpi_communicator());

	  
      statistics.add_value ("Serp per year(kg/(yr*km))", global_M_t);
	  {
        const char *columns[] = {"Serp per year(kg/(yr*km))"};
        for (unsigned int i=0; i<sizeof(columns)/sizeof(columns[0]); ++i)
          {
            statistics.set_precision (columns[i], 8);
            statistics.set_scientific (columns[i], true);
          }
      }
      std::ostringstream output;
      output.precision(4);
	  output << global_M_t << "   kg/(yr*km). ";
      return std::pair<std::string, std::string> ("Serp per year:",
                                                  output.str());
    }
	
	template <int dim>
    void
    SerpStatistics<dim>::parse_parameters (ParameterHandler &prm)
	{
		Tp = prm.get_double ("Adiabatic surface temperature");
		
		prm.enter_subsection("Material model");
		{
			prm.enter_subsection ("Visco Plastic Grain Size");
			{
				const unsigned int n_fields = this->n_compositional_fields() + 1;
				densities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Densities"))),
                                                     n_fields,
                                                     "Densities");
				thermal_expansivities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Thermal expansivities"))),
                                                                           n_fields,
                                                                           "Thermal expansivities");
			}
			prm.leave_subsection();
		}
		prm.leave_subsection();
		
		prm.enter_subsection("Geometry model");
		{
			prm.enter_subsection("Box");
			{
				max_x = prm.get_double ("X extent");
			}
			prm.leave_subsection ();
		}
		prm.leave_subsection ();
		
	}
	
	
	
	
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(SerpStatistics,
                                  "serp statistics",
                                  "A postprocessor that computes some statistics about "
                                  "the compositional fields, if present in this simulation. "
                                  "In particular, it computes maximal and minimal values of "
                                  "each field, as well as the total mass contained in this "
                                  "field as defined by the integral "
                                  "$m_i(t) = \\int_\\Omega c_i(\\mathbf x,t) \\; \\text{d}x$.")
  }
}
