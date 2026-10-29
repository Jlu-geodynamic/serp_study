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


#include </fs2/home/liuzhonglan/wy/lib_extra/melt20260814/grain/include/postprocess/melt_emplacement.h>
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
	void
	MeltEmplacement<dim>::
	write_to_container(const std::string path,
	                  const unsigned int crood,
	                  const unsigned int mpi_communicator,
				      const double write_information)const
	{
		if(this->simulator_is_past_initialization() && this->get_timestep_number() > 0 && Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == mpi_communicator)
		{
			std::stringstream get_name;
		    get_name << crood;
		    std::string temp_txt_name = get_name.str();
		    std::ofstream write((path + temp_txt_name + ".txt").c_str());
			while(!write)
			{
				write.close();
				sleep(1.);
				std::ofstream write((path + temp_txt_name + ".txt").c_str());
			}
		    write << write_information;
		    write.close();
		}
	}
	
	template <int dim>
	double
	MeltEmplacement<dim>::
	read_from_container(const std::string path,
	                   const unsigned int crood,
	                   const unsigned int mpi_communicator)const
	{
		double temp_input = 0;
		if(this->simulator_is_past_initialization() && this->get_timestep_number() > 0 && Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == mpi_communicator)
		{
			std::stringstream get_name;
		    get_name << crood;
		    std::string temp_txt_name = get_name.str();
			std::ifstream data_from_txt((path + temp_txt_name + ".txt").c_str());
			std::string line;
			getline(data_from_txt, line);
			while(line == "")
			{
				data_from_txt.close();
				sleep(1.);
				std::ifstream data_from_txt((path + temp_txt_name + ".txt").c_str());
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
	double
	MeltEmplacement<dim>::
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
	void
	MeltEmplacement<dim>::
    convey_flag_update(const std::string path,
	                   const unsigned int crood_x,
					   const unsigned int crood_y_in_flag_container,
					   const unsigned int mpi_communicator,
					   const unsigned int flag)const
	{
		if(this->simulator_is_past_initialization() && this->get_timestep_number() > 0 && Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == mpi_communicator)
		{
			const int flag_place = 2 * crood_y_in_flag_container;
			const double write_information = flag;
			std::stringstream get_name;
			get_name << crood_x;
			std::string temp_txt_name = get_name.str();
			std::fstream flag_container((path + temp_txt_name + ".txt").c_str(), std::ios::in | std::ios::out);
			while(!flag_container)
			{
				flag_container.close();
				sleep(1.);
				std::fstream flag_container((path + temp_txt_name + ".txt").c_str(), std::ios::in | std::ios::out);
			}
			std::stringstream ss;
			ss << write_information;
			flag_container.seekp(flag_place, std::ios::beg);
			flag_container << ss.str();
			flag_container.close(); 
		}
	}
    	
    template <int dim>
    std::pair<std::string,std::string>
    MeltEmplacement<dim>::execute (TableHandler &statistics)
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
			const double posy = in.position[mid_point][1];
	        const unsigned int crood_x = std::floor( posx / ((max_x / repetitions_x) / std::pow(2, melt_fraction_refinement)));
			const unsigned int crood_y = std::floor( posy / ((max_y / repetitions_y) / std::pow(2, melt_fraction_refinement)));
			const int crood_y_in_flag_container = crood_y - min_container_amount_y;
			
			//清空搬运标记
			if (crood_y_in_flag_container >= 0)
			{
				if(include_mantle_melting)
					convey_flag_update(convey_flag_mantle_melting_path, crood_x, crood_y_in_flag_container, mpi_comm_id, 0);
				if(include_crustal_melting)
					convey_flag_update(convey_flag_crustal_melting_path, crood_x, crood_y_in_flag_container, mpi_comm_id, 0);
			}
			
			//仅在混合熔融中需要判定成分
			//以中间位置的占比为判定依据
			double crust_fraction = 0.;
			double mantle_fraction = 0.;
			if(include_mantle_melting && include_crustal_melting)
			{
				//20250701
		        //调整成分场调用逻辑，可在prm文件中删除不需要的成分场
				//20260127 膏盐模式与地壳熔融无法同时使用，因此膏盐模式中新增的相不在此添加
				//vector
				std::vector<double> sd1_values(n_q_points, 0.);
				std::vector<double> sd2_values(n_q_points, 0.);
				std::vector<double> upper_values(n_q_points, 0.);
				std::vector<double> lower_values(n_q_points, 0.);
				std::vector<double> mantle_upper_values(n_q_points, 0.);
				//std::vector<double> mantle_middle_values(n_q_points, 0.);
				//std::vector<double> mantle_lower_values(n_q_points, 0.);
				std::vector<double> target_mantle_melting_values(n_q_points, 0.);
				std::vector<double> target_crustal_melting_values(n_q_points, 0.);
				//mid value
				double sd1_mid = 0.;
				double sd2_mid = 0.;
				double upper_mid = 0.;
				double lower_mid = 0.;
				double mantle_upper_mid = 0.;
				//double mantle_middle_mid = 0.;
				//double mantle_lower_mid = 0.;
				double target_mantle_melting_mid = 0.;
				double target_crustal_melting_mid = 0.;
				if (this->introspection().compositional_name_exists("sediment_1"))
				{
					const unsigned int sd1_id = this->introspection().compositional_index_for_name("sediment_1");
				    fe_values[this->introspection().extractors.compositional_fields[sd1_id]].get_function_values (this->get_solution(), sd1_values);
					sd1_mid = sd1_values[mid_point];
				}
				if(this->introspection().compositional_name_exists("sediment_2"))
				{
					const unsigned int sd2_id = this->introspection().compositional_index_for_name("sediment_2");
				    fe_values[this->introspection().extractors.compositional_fields[sd2_id]].get_function_values (this->get_solution(), sd2_values); 
					sd2_mid = sd2_values[mid_point];
				}
				if(this->introspection().compositional_name_exists("upper"))
				{
					const unsigned int upper_id = this->introspection().compositional_index_for_name("upper");
				    fe_values[this->introspection().extractors.compositional_fields[upper_id]].get_function_values (this->get_solution(), upper_values);
					upper_mid = upper_values[mid_point];
				}
			    if(this->introspection().compositional_name_exists("lower"))
				{
					const unsigned int lower_id = this->introspection().compositional_index_for_name("lower");
				    fe_values[this->introspection().extractors.compositional_fields[lower_id]].get_function_values (this->get_solution(), lower_values);
					lower_mid = lower_values[mid_point];
				}
				if(this->introspection().compositional_name_exists("mantle_upper"))
				{
					const unsigned int mantle_upper_id = this->introspection().compositional_index_for_name("mantle_upper");
				    fe_values[this->introspection().extractors.compositional_fields[mantle_upper_id]].get_function_values (this->get_solution(), mantle_upper_values); 
					mantle_upper_mid = mantle_upper_values[mid_point];
				}
				/* if(this->introspection().compositional_name_exists("mantle_middle"))
				{
					const unsigned int mantle_middle_id = this->introspection().compositional_index_for_name("mantle_middle");
				    fe_values[this->introspection().extractors.compositional_fields[mantle_middle_id]].get_function_values (this->get_solution(), mantle_middle_values); 
					mantle_middle_mid = mantle_middle_values[mid_point];
				}
				if(this->introspection().compositional_name_exists("mantle_lower"))
				{
					const unsigned int mantle_lower_id = this->introspection().compositional_index_for_name("mantle_lower");
				    fe_values[this->introspection().extractors.compositional_fields[mantle_lower_id]].get_function_values (this->get_solution(), mantle_lower_values); 
					mantle_lower_mid = mantle_lower_values[mid_point];
				} */
			    if(this->introspection().compositional_name_exists("new_crust"))
				{
					const unsigned int target_mantle_melting_id = this->introspection().compositional_index_for_name("new_crust");
				    fe_values[this->introspection().extractors.compositional_fields[target_mantle_melting_id]].get_function_values (this->get_solution(), target_mantle_melting_values); 
					target_mantle_melting_mid = target_mantle_melting_values[mid_point];
				}					
				if(this->introspection().compositional_name_exists("core_complex"))
				{
					const unsigned int target_crustal_melting_id = this->introspection().compositional_index_for_name("core_complex");
				    fe_values[this->introspection().extractors.compositional_fields[target_crustal_melting_id]].get_function_values (this->get_solution(), target_crustal_melting_values); 
					target_crustal_melting_mid = target_crustal_melting_values[mid_point];
				}			
				crust_fraction = upper_mid + lower_mid;
				mantle_fraction = 1. - sd1_mid - sd2_mid - upper_mid - lower_mid - mantle_upper_mid - target_mantle_melting_mid - target_crustal_melting_mid;
			}
			

			//1.熔体抽取
			//20240703
			//模式1：集中抽取
			//将当前步提取的熔体平均输入到窗口中
			if (!vertical_trans)
			{
				if(include_mantle_melting)
				{
					//岩浆侵位窗口
					const double temp_mantle_min_window = read_from_txt(melt_info_path + "mantle_min_emplace_window.txt", mpi_comm_id);
					const unsigned int mantle_min_window = temp_mantle_min_window;
					const double temp_mantle_max_window = read_from_txt(melt_info_path + "mantle_max_emplace_window.txt", mpi_comm_id);
					const unsigned int mantle_max_window = temp_mantle_max_window;
					//每一列提取一次即可，此处取熔融传输最深的位置之下(min_container_amount_y - 1)作为参考
					const double mantle_emplace_range = mantle_max_window - mantle_min_window + 1;
					if (mantle_max_window > 0 && crood_y == min_container_amount_y - 1 && crood_x >= mantle_min_window && crood_x <= mantle_max_window)
					{
						const double current_mantle_melt_extraction = read_from_txt(melt_info_path + "mantle_current_melt_extraction.txt", mpi_comm_id);
						if (current_mantle_melt_extraction > 0.)
						{
							const double container_mantle_melt_fraction = read_from_container(container_mantle_melting_path, crood_x, mpi_comm_id);
							const double mantle_write_information = container_mantle_melt_fraction + current_mantle_melt_extraction / mantle_emplace_range;
							write_to_container(container_mantle_melting_path, crood_x, mpi_comm_id, mantle_write_information);
						}
					}
				}
				if(include_crustal_melting)
				{
					const double temp_crustal_min_window = read_from_txt(melt_info_path + "crust_min_emplace_window.txt", mpi_comm_id);
					const unsigned int crustal_min_window = temp_crustal_min_window;
					const double temp_crustal_max_window = read_from_txt(melt_info_path + "crust_max_emplace_window.txt", mpi_comm_id);
					const unsigned int crustal_max_window = temp_crustal_max_window;

					const double crustal_emplace_range = crustal_max_window - crustal_min_window + 1;
					if (crustal_max_window > 0 && crood_y == min_container_amount_y - 1 && crood_x >= crustal_min_window && crood_x <= crustal_max_window)
					{
						const double current_crustal_melt_extraction = read_from_txt(melt_info_path + "crust_current_melt_extraction.txt", mpi_comm_id);
						if (current_crustal_melt_extraction > 0.)
						{
							const double container_crustal_melt_fraction = read_from_container(container_crustal_melting_path, crood_x, mpi_comm_id);
							const double crustal_write_information = container_crustal_melt_fraction + current_crustal_melt_extraction / crustal_emplace_range;
							write_to_container(container_crustal_melting_path, crood_x, mpi_comm_id, crustal_write_information);
						}
					}
				}
			}
			//模式2：垂直抽取
			//当前步提取的熔体将在其产生的位置垂直上升到地表。
			else if (vertical_trans)
			{
				std::vector<double> melt_record_values(n_q_points, 0.);
				std::vector<double> melt_fraction_values(n_q_points, 0.);
				std::vector<double> max_melt_values(n_q_points, 0.);
				std::vector<double> mantle_lower_values(n_q_points, 0.);
				if(this->introspection().compositional_name_exists("melt_record"))
				{
					unsigned int melt_record_id = this->introspection().compositional_index_for_name("melt_record");
			        fe_values[this->introspection().extractors.compositional_fields[melt_record_id]].get_function_values (this->get_solution(), melt_record_values);
				}
				if(this->introspection().compositional_name_exists("melt_fraction"))
				{
					unsigned int melt_fraction_id = this->introspection().compositional_index_for_name("melt_fraction");
			        fe_values[this->introspection().extractors.compositional_fields[melt_fraction_id]].get_function_values (this->get_solution(), melt_fraction_values);
				}
				//20250428 含碳熔融相关
			    if(this->introspection().compositional_name_exists("max_melt"))
				{
					unsigned int max_melt_id = this->introspection().compositional_index_for_name("max_melt");
			        fe_values[this->introspection().extractors.compositional_fields[max_melt_id]].get_function_values (this->get_solution(), max_melt_values);
				}
				//默认含碳熔融在‘mantle_lower’相中产生。
				//如有需求可酌情修改。
				if(this->introspection().compositional_name_exists("mantle_lower"))
			    {
					unsigned int mantle_lower_id = this->introspection().compositional_index_for_name("mantle_lower");
			        fe_values[this->introspection().extractors.compositional_fields[mantle_lower_id]].get_function_values (this->get_solution(), mantle_lower_values);
				}
			    double mantle_melt_add = 0;
				double crustal_melt_add = 0;
				//20250428
				double carbon_melt_add = 0;
			    for (unsigned int q=0; q<n_q_points; ++q)
		        {
					const double melt_record_in_point = melt_record_values[q];
				    const double melt_fraction_in_point = melt_fraction_values[q];
					//20250428
					const double max_melt_fraction_in_point = max_melt_values[q];
					const double mantle_lower_in_point = mantle_lower_values[q];
				    if (melt_fraction_in_point > melt_limit && melt_record_in_point > 0)
					{
						if((include_mantle_melting && !include_crustal_melting && !include_carbonate_melting) 
						|| (include_mantle_melting && include_crustal_melting && !include_carbonate_melting && mantle_fraction >= crust_fraction)
					    || (include_mantle_melting && include_carbonate_melting && (mantle_lower_in_point < 0.3 || mantle_lower_in_point >= 0.3 && max_melt_fraction_in_point >= 0.02)))
							mantle_melt_add += std::min(melt_record_in_point * melt_factor * fe_values.JxW(q) / 1e6, 5e-2 * fe_values.JxW(q) / 1e6);
						else if((!include_mantle_melting && include_crustal_melting) || (include_mantle_melting && include_crustal_melting && mantle_fraction < crust_fraction))
							crustal_melt_add += std::min(melt_record_in_point * melt_factor * fe_values.JxW(q) / 1e6, 5e-2 * fe_values.JxW(q) / 1e6);
						else if(include_mantle_melting && include_carbonate_melting && mantle_lower_in_point >= 0.3 && max_melt_fraction_in_point < 0.03)
							carbon_melt_add += std::min(melt_record_in_point * melt_factor * fe_values.JxW(q) / 1e6, 5e-2 * fe_values.JxW(q) / 1e6);
						else
						{
							mantle_melt_add = 0;
							crustal_melt_add = 0;
							carbon_melt_add = 0;
						}
					}
		        }
				
				if (mantle_melt_add > 0 && this->get_timestep_number() > 0 && mpi_comm_id > 0)
			    {
				    const double container_mantle_melt_fraction = read_from_container(container_mantle_melting_path, crood_x, mpi_comm_id);
				    const double mantle_write_information = container_mantle_melt_fraction + mantle_melt_add;
				    write_to_container(container_mantle_melting_path, crood_x, mpi_comm_id, mantle_write_information);
			    }
				if (crustal_melt_add > 0 && this->get_timestep_number() > 0 && mpi_comm_id > 0)
				{
					const double container_crustal_melt_fraction = read_from_container(container_crustal_melting_path, crood_x, mpi_comm_id);
				    const double crustal_write_information = container_crustal_melt_fraction + crustal_melt_add;
				    write_to_container(container_crustal_melting_path, crood_x, mpi_comm_id, crustal_write_information);
				}
				if (carbon_melt_add > 0 && this->get_timestep_number() > 0 && mpi_comm_id > 0)
			    {
				    const double carbon_container_melt_fraction = read_from_container(carbon_container_path, crood_x, mpi_comm_id);
				    const double carbon_write_information = carbon_container_melt_fraction + carbon_melt_add;
				    write_to_container(carbon_container_path, crood_x, mpi_comm_id, carbon_write_information);
			    }
			}
          }
	  //预留位置，如有需要可视化的参数再进行添加
      statistics.add_value ("Melt extraction", 0);
	  {
        const char *columns[] = {"Melt extraction"};
        for (unsigned int i=0; i<sizeof(columns)/sizeof(columns[0]); ++i)
          {
            statistics.set_precision (columns[i], 8);
            statistics.set_scientific (columns[i], true);
          }
      }
      std::ostringstream output;
      output.precision(4);
	  output << "complete";
      return std::pair<std::string, std::string> ("Processing of melt emplacement,",
                                                  output.str());
    }
	
	template <int dim>
    void
    MeltEmplacement<dim>::parse_parameters (ParameterHandler &prm)
	{
		prm.enter_subsection("Postprocess");
		{
			prm.enter_subsection("Melt extraction");
			{
				depth_max = prm.get_double ("Max depth of melt convey function");
				//20240703
				vertical_trans = prm.get_bool("Vertical transportation mode");
				melt_limit = prm.get_double ("Limit melt fraction of extracting");
				melt_factor = prm.get_double ("Coefficient of melt extraction");
				//20240614
				//melt_fraction_refinement = prm.get_integer ("Mesh refinement of melt fraction");
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
				
				min_container_amount_y = std::ceil((max_y - depth_max - 5e3) / ((max_y / repetitions_y) / std::pow(2, melt_fraction_refinement)));
				convey_flag_mantle_melting_path = this->get_output_directory() + "convey_flag_mantle_melting/";
				container_mantle_melting_path = this->get_output_directory() + "container_mantle_melting/";
				convey_flag_crustal_melting_path = this->get_output_directory() + "convey_flag_crustal_melting/";
				container_crustal_melting_path = this->get_output_directory() + "container_crustal_melting/";
				melt_info_path = this->get_output_directory() + "melt_info/";
				//20250428
				carbon_container_path = this->get_output_directory() + "container_carbonate_melting/";
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
    ASPECT_REGISTER_POSTPROCESSOR(MeltEmplacement,
                                  "melt emplacement",
                                  "A postprocessor that computes some statistics about "
                                  "the compositional fields, if present in this simulation. "
                                  "In particular, it computes maximal and minimal values of "
                                  "each field, as well as the total mass contained in this "
                                  "field as defined by the integral "
                                  "$m_i(t) = \\int_\\Omega c_i(\\mathbf x,t) \\; \\text{d}x$.")
  }
}
