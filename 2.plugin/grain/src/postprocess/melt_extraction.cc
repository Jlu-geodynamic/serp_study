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


#include </fs2/home/liuzhonglan/wy/lib_extra/melt20260814/grain/include/postprocess/melt_extraction.h>
#include </fs2/home/liuzhonglan/wy/lib_extra/melt20260814/grain/include/postprocess/pre_melt_extraction.h>

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

//#include <mpi.h>

namespace aspect
{
  namespace Postprocess
  {
	  
	  
	  
	template <int dim>
	double
	MeltExtraction<dim>::
	read_from_txt(const std::string path,
	              const unsigned int mpi_communicator)const
	{
		double temp_input = 0;
		if(this->simulator_is_past_initialization() && this->get_timestep_number() > 0 && Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == mpi_communicator)
		{
			std::ifstream data_from_txt(path.c_str());
			std::string line;
			getline(data_from_txt, line);
			while(line == "")
			{
				data_from_txt.close();
				sleep(1.);
				std::ifstream data_from_txt(path.c_str());
				getline(data_from_txt, line);
			}
			std::stringstream ss(line);
			ss >> temp_input;
			data_from_txt.close();
			return temp_input;
		} 
		return temp_input;
	}	

    template <int dim>
    std::pair<std::string,std::string>
    MeltExtraction<dim>::execute (TableHandler &statistics)
    {
      if (this->n_compositional_fields() == 0)
        return std::pair<std::string,std::string>();

      // create a quadrature formula based on the compositional element alone.
      // be defensive about determining that a compositional field actually exists
      AssertThrow (this->introspection().base_elements.compositional_fields
                   != numbers::invalid_unsigned_int,
                   ExcMessage("This postprocessor cannot be used without compositional fields."));
	  AssertThrow (!dike_mode || ((dike_mode && include_mantle_melting) && !(dike_mode && include_crustal_melting) && !(dike_mode && include_carbonate_melting)),
                   ExcMessage("Cannot use dike mode with crustal melting or carbonate melting together."));
				   
	  
      const QGauss<dim> quadrature_formula (this->get_fe().base_element(this->introspection().base_elements.compositional_fields).degree+1);
      const unsigned int n_q_points = quadrature_formula.size();

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values   |
                               update_quadrature_points |
                               update_JxW_values);

	  MaterialModel::MaterialModelInputs<dim> in(fe_values.n_quadrature_points, this->n_compositional_fields());
	  unsigned int mantle_local_min_x_crood = 9999;
	  unsigned int mantle_local_max_x_crood = 0;
	  unsigned int crustal_local_min_x_crood = 9999;
	  unsigned int crustal_local_max_x_crood = 0;
	  //20251104
	  //20251113 新增岩墙位置判定
	  double dike_emplace_location_xmax_max = 0;
	  double dike_emplace_location_xmax_min = 99999;
	  double dike_emplace_location_ymin = 0;
	  //20251201 新增岩墙顶部y坐标
	  //岩墙长度不固定时启用
	  double dike_emplace_location_ymax = 0;

      // compute the integral quantities by quadrature
      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          {
			//获取mpi通道编号
			unsigned int mpi_comm_id = 0;
			unsigned int out_comm_id = 0;
		    while (out_comm_id < 5000)
		    {
			    if (Utilities::MPI::this_mpi_process (this->get_mpi_communicator()) == out_comm_id)
			    {
				    mpi_comm_id = out_comm_id;
				    break;
			    } 
			    out_comm_id ++;
		    } 
			//填充输入
            fe_values.reinit (cell);
			in.reinit(fe_values, cell, this->introspection(), this->get_solution());
			const unsigned int mid_point = n_q_points / 2;
			const double posx = in.position[mid_point][0];
	        const unsigned int crood_x = std::floor( posx / ((max_x / repetitions_x) / std::pow(2, melt_fraction_refinement)));
			//20240703
			if(!vertical_trans && !dike_mode)
			{
				if(include_mantle_melting)
				{
					const double mantle_emplace_location_y = read_from_txt(melt_info_path + "mantle_emplace_location_y.txt", mpi_comm_id);
					for (unsigned int q=0; q<n_q_points; ++q)
					{
						const double posy = in.position[q][1];
						//20260314 宽度降低
						  if (posy / 1e3 >= mantle_emplace_location_y - 6 && posy / 1e3 < mantle_emplace_location_y  && in.temperature[q] >= 1473)
					      //if (posy / 1e3 >= mantle_emplace_location_y - 2 && posy / 1e3 < mantle_emplace_location_y  && in.temperature[q] >= 1473)
						{
							mantle_local_min_x_crood = crood_x;
							mantle_local_max_x_crood = crood_x;
						}
					}
				}
				if(include_crustal_melting)
				{
					const double crustal_emplace_location_y = read_from_txt(melt_info_path + "crust_emplace_location_y.txt", mpi_comm_id);
					for (unsigned int q=0; q<n_q_points; ++q)
					{
						const double posy = in.position[q][1];
						  if (posy / 1e3 >= crustal_emplace_location_y - 6 && posy / 1e3 < crustal_emplace_location_y  && in.temperature[q] >= 873)
						{
							crustal_local_min_x_crood = crood_x;
							crustal_local_max_x_crood = crood_x;
						}
					}
				}
			}
			else if(dike_mode)
			{
				const double mantle_emplace_location_y = read_from_txt(melt_info_path + "mantle_emplace_location_y.txt", mpi_comm_id);
				for (unsigned int q=0; q<n_q_points; ++q)
					{
						const double posx_in_km = in.position[q][0] / 1e3;
						const double posy_in_km = in.position[q][1] / 1e3;
						//20251112 修改岩墙位置判定
						  //if (std::abs(posy_in_km - mantle_emplace_location_y) <= 2 * max_refinement_length_y / 1e3 && in.temperature[q] >= 1473 && in.temperature[q] <= 1483) 
						  if (std::abs(posy_in_km - mantle_emplace_location_y) <= 10 && in.temperature[q] >= 1453) 
						    {
								//20251113 修改岩墙位置判定
							    dike_emplace_location_xmax_max = posx_in_km;
								dike_emplace_location_xmax_min = posx_in_km;
						    }
					}
			}
			else
			{
				mantle_local_min_x_crood = 0;
				mantle_local_max_x_crood = 0;
				crustal_local_min_x_crood = 0.;
				crustal_local_max_x_crood = 0.;
			}
			
          }
		  
	  const unsigned int global_mantle_min_x_crood = Utilities::MPI::min (mantle_local_min_x_crood, this->get_mpi_communicator());
	  const unsigned int global_mantle_max_x_crood = Utilities::MPI::max (mantle_local_max_x_crood, this->get_mpi_communicator());
	  
	  const unsigned int global_crustal_min_x_crood = Utilities::MPI::min (crustal_local_min_x_crood, this->get_mpi_communicator());
	  const unsigned int global_crustal_max_x_crood = Utilities::MPI::max (crustal_local_max_x_crood, this->get_mpi_communicator());
	  const double global_dike_emplace_location_xmax_max = Utilities::MPI::max (dike_emplace_location_xmax_max, this->get_mpi_communicator());
	  const double global_dike_emplace_location_xmax_min = Utilities::MPI::min (dike_emplace_location_xmax_min, this->get_mpi_communicator());

	  //根据具体情况将熔体传送至地壳增生窗口，模型最大熔融程度小于0.2%时无法开启传输，期间产生的熔体会流失
	  //这部分熔体量非常少，对地壳增生的影响可以忽略不计。
	  
	  //当最大值大于0时，代表存在地壳增生窗口
	  const double min_mantle_emplace_window = global_mantle_max_x_crood > 0 ? global_mantle_min_x_crood : 0;
	  const double max_mantle_emplace_window = global_mantle_max_x_crood > 0 ? global_mantle_max_x_crood : 0;
	  const double min_mantle_emplace_window_in_km = min_mantle_emplace_window * ((max_x / repetitions_x) / std::pow(2, melt_fraction_refinement)) / 1e3;
	  const double max_mantle_emplace_window_in_km = max_mantle_emplace_window * ((max_x / repetitions_x) / std::pow(2, melt_fraction_refinement)) / 1e3;
	  
	  const double min_crustal_emplace_window = global_crustal_max_x_crood > 0 ? global_crustal_min_x_crood : 0;
	  const double max_crustal_emplace_window = global_crustal_max_x_crood > 0 ? global_crustal_max_x_crood : 0;
	  const double min_crustal_emplace_window_in_km = min_crustal_emplace_window * ((max_x / repetitions_x) / std::pow(2, melt_fraction_refinement)) / 1e3;
	  const double max_crustal_emplace_window_in_km = max_crustal_emplace_window * ((max_x / repetitions_x) / std::pow(2, melt_fraction_refinement)) / 1e3;
	  
	  if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
	  {
		  if(include_mantle_melting)
		  {
			  if(!dike_mode)
			  {
				  std::ofstream mantle_write_min_win((melt_info_path + "mantle_min_emplace_window.txt").c_str());
			      if(mantle_write_min_win)
				      mantle_write_min_win << min_mantle_emplace_window;
			      mantle_write_min_win.close();
			  
			      std::ofstream mantle_write_max_win((melt_info_path + "mantle_max_emplace_window.txt").c_str());
			      if(mantle_write_max_win)
				      mantle_write_max_win << max_mantle_emplace_window;
			      mantle_write_max_win.close();
			  }
			  else
			  {
				  //20251113 修改岩墙位置判定
				  //20251118 岩墙X位置仅在岩墙形成前更新，岩墙形成后X坐标固定
				  //岩浆判定条件临时使用10kg/m2
				  //20251124 功能：可控制岩墙x轴是否更新，用vertical_trans临时代替
				  const double temp_mantle_total_melt_extraction = read_from_txt(melt_info_path + "mantle_total_melt_extraction.txt", 0);
				  const double out_dike_emplace_location_xmax = 0.5 * (global_dike_emplace_location_xmax_max + global_dike_emplace_location_xmax_min);
				  //20251201
				  //岩墙x坐标固定时，当熔融提取面积到达岩墙面积90%则停止x轴位置更新
				  //岩墙长度不固定时，使用上一步的y坐标计算面积。
				  //y坐标信息为0时（初始步或没有明显软流圈上涌）或者使用上一步y坐标信息，对最终位置影响不大，因此不针对此种情况进行判定。
				  const double dike_ymax_km = constant_dike_length ? 0. : read_from_txt(melt_info_path + "dike_emplace_location_ymax.txt", 0);
				  const double dike_ymin_km = constant_dike_length ? 0. : read_from_txt(melt_info_path + "dike_emplace_location_ymin.txt", 0);
				  const double dike_length_value = constant_dike_length ? dike_length : (dike_ymax_km - dike_ymin_km) * 1e3;
				  if(!constant_dike_x || (constant_dike_x && temp_mantle_total_melt_extraction <= 0.9 * dike_width * dike_length_value / 1e6))
				  {
					  std::ofstream dike_write_location_xmax((melt_info_path + "dike_emplace_location_xmax.txt").c_str());
				      if(dike_write_location_xmax)
					      dike_write_location_xmax << out_dike_emplace_location_xmax;
				      dike_write_location_xmax.close();
				  }
				  
			  }
			  
		  }
		  
		  if(include_crustal_melting)
		  {
			  std::ofstream crustal_write_min_win((melt_info_path + "crust_min_emplace_window.txt").c_str());
			  if(crustal_write_min_win)
				  crustal_write_min_win << min_crustal_emplace_window;
			  crustal_write_min_win.close();
			  
			  std::ofstream crustal_write_max_win((melt_info_path + "crust_max_emplace_window.txt").c_str());
			  if(crustal_write_max_win)
				  crustal_write_max_win << max_crustal_emplace_window;
			  crustal_write_max_win.close();
		  }
		  
	  }
	  //20251104
	  //岩墙启用时，再次遍历网格，确定岩墙底部的y坐标
	  if(dike_mode)
	  {
		  for (const auto &cell : this->get_dof_handler().active_cell_iterators())
			  if (cell->is_locally_owned())
			  {
				  //获取mpi通道编号
			      unsigned int dike_mpi_comm_id = 0;
			      unsigned int dike_out_comm_id = 0;
		          while (dike_out_comm_id < 5000)
		          {
			          if (Utilities::MPI::this_mpi_process (this->get_mpi_communicator()) == dike_out_comm_id)
			          {
				          dike_mpi_comm_id = dike_out_comm_id;
				          break;
			          } 
			          dike_out_comm_id ++;
		          } 
			      //填充输入
                  fe_values.reinit (cell);
			      in.reinit(fe_values, cell, this->introspection(), this->get_solution());
				  const double dike_emplace_location_x = read_from_txt(melt_info_path + "dike_emplace_location_xmax.txt", dike_mpi_comm_id);
				  for (unsigned int q=0; q<n_q_points; ++q)
					{
						const double dike_posx_in_km = in.position[q][0] / 1e3;
						const double dike_posy_in_km = in.position[q][1] / 1e3;
						//20251112 修改位置判定，x轴范围缩减，y轴只要大于设置的岩墙底部值便置入，随后取最大值
						//if (std::abs(dike_posx_in_km - dike_emplace_location_x) <= 2 * max_refinement_length_x / 1e3 && in.temperature[q] >= dike_bottom_temperature && in.temperature[q] <= dike_bottom_temperature + 10)
						if (std::abs(dike_posx_in_km - dike_emplace_location_x) <= max_refinement_length_x / 1e3 && in.temperature[q] >= dike_bottom_temperature)
							dike_emplace_location_ymin = dike_posy_in_km;
						//20251201 岩墙最大y坐标
						if (!constant_dike_length && std::abs(dike_posx_in_km - dike_emplace_location_x) <= max_refinement_length_x / 1e3 && in.temperature[q] >= dike_top_temperature)
							dike_emplace_location_ymax = dike_posy_in_km;
					}
			  }
	  }
	  const double global_dike_emplace_location_ymin = Utilities::MPI::max (dike_emplace_location_ymin, this->get_mpi_communicator());
	  //20251201
	  const double global_dike_emplace_location_ymax = Utilities::MPI::max (dike_emplace_location_ymax, this->get_mpi_communicator());
	  if (dike_mode && include_mantle_melting && Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
	  {
		  std::ofstream dike_write_location_ymin((melt_info_path + "dike_emplace_location_ymin.txt").c_str());
		  if(dike_write_location_ymin)
			  dike_write_location_ymin << global_dike_emplace_location_ymin;
		  dike_write_location_ymin.close();
		  //20251201
		  if(!constant_dike_length)
		  {
			  std::ofstream dike_write_location_ymax((melt_info_path + "dike_emplace_location_ymax.txt").c_str());
			  if(dike_write_location_ymax)
			      dike_write_location_ymax << global_dike_emplace_location_ymax;
		      dike_write_location_ymax.close();
		  }
	  }
	  //20251202 添加岩墙相关输出
	  const double dike_length_dynamic = constant_dike_length ? dike_length : (global_dike_emplace_location_ymax - global_dike_emplace_location_ymin) * 1e3;
	  
      statistics.add_value ("Min emplacing window(mantle,km)", min_mantle_emplace_window_in_km);
	  statistics.add_value ("Max emplacing window(mantle,km)", max_mantle_emplace_window_in_km);
	  statistics.add_value ("Min emplacing window(crustal,km)", min_crustal_emplace_window_in_km);
	  statistics.add_value ("Max emplacing window(crustal,km)", max_crustal_emplace_window_in_km);
	  statistics.add_value ("Dike length(m)", dike_length_dynamic);
	  {
        const char *columns[] = {"Min emplacing window(mantle,km)", "Max emplacing window(mantle,km)", "Min emplacing window(crustal,km)", "Max emplacing window(crustal,km)"};
        for (unsigned int i=0; i<sizeof(columns)/sizeof(columns[0]); ++i)
          {
            statistics.set_precision (columns[i], 8);
            statistics.set_scientific (columns[i], true);
          }
      }
      std::ostringstream output;
      output.precision(4);
	  output << min_mantle_emplace_window_in_km << "km to "
	         << max_mantle_emplace_window_in_km << "km /"
			 << min_crustal_emplace_window_in_km << "km to"
			 << max_crustal_emplace_window_in_km << "km, dike x range "
			 << global_dike_emplace_location_xmax_min << "to"
			 << global_dike_emplace_location_xmax_max <<"km, dike length "
			 << dike_length_dynamic <<"m.";
      return std::pair<std::string, std::string> ("Extracting melt, emplacing window mantle/crustal:",
                                                  output.str());
    }
	
	template <int dim>
    void
    MeltExtraction<dim>::parse_parameters (ParameterHandler &prm)
	{
		//20240703
		 prm.enter_subsection("Postprocess");
		{
			prm.enter_subsection("Melt extraction");
			{
				vertical_trans = prm.get_bool("Vertical transportation mode");
				//20251104
				dike_mode = prm.get_bool("Use dike mode");
				dike_bottom_temperature = prm.get_double("Temperature of the bottom of dike(K)");
				//20251201
				constant_dike_length = prm.get_bool("Constant dike length");
			    constant_dike_x = prm.get_bool("Constant dike X position");
				dike_top_temperature = prm.get_double("Temperature of the top of dike(K)");
				dike_width = prm.get_double ("Width of dike(m)");
			    dike_length = prm.get_double ("Length of dike(m)");
			}
			prm.leave_subsection ();
		}
		prm.leave_subsection (); 
		
		
		//20240614
		prm.enter_subsection("Mesh refinement");
		{
			melt_fraction_refinement = prm.get_integer("Initial adaptive refinement") + prm.get_integer("Initial global refinement");
			
		}
		prm.leave_subsection ();
		
		prm.enter_subsection("Geometry model");
		{
			prm.enter_subsection("Box");
			{
				max_x = prm.get_double ("X extent");
				repetitions_x = prm.get_integer ("X repetitions");
				max_y = prm.get_double ("Y extent");
				repetitions_y = prm.get_integer ("Y repetitions");
				melt_info_path = this->get_output_directory() + "melt_info/";
				max_refinement_length_x = (max_x / repetitions_x) / std::pow(2, melt_fraction_refinement);
				max_refinement_length_y = (max_y / repetitions_y) / std::pow(2, melt_fraction_refinement);
			}
			prm.leave_subsection ();
		}
		prm.leave_subsection ();
		
		prm.enter_subsection("Material model");
		{
			prm.enter_subsection ("Visco Plastic Grain Size");
			{
				include_mantle_melting = prm.get_bool ("Include mantle melting");
				include_crustal_melting = prm.get_bool("Include crustal melting");
				//20251104
				include_carbonate_melting = prm.get_bool("Include carbonate melting");
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
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(MeltExtraction,
                                  "melt extraction",
                                  "A postprocessor that computes some statistics about "
                                  "the compositional fields, if present in this simulation. "
                                  "In particular, it computes maximal and minimal values of "
                                  "each field, as well as the total mass contained in this "
                                  "field as defined by the integral "
                                  "$m_i(t) = \\int_\\Omega c_i(\\mathbf x,t) \\; \\text{d}x$.")
  }
}
