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
	PreMeltExtraction<dim>::
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
    PreMeltExtraction<dim>::execute (TableHandler &statistics)
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
				   
      //初始化阶段创建相应txt容器
	  if (this->get_timestep_number() == 0 && Utilities::MPI::this_mpi_process (this->get_mpi_communicator()) == 0)
	  {
		  //20251104
		  //dike模式不需要flag和容器
		  if(!dike_mode)
		  {
			  DIR *convey_flag_mantle_melting_directory = opendir(convey_flag_mantle_melting_path.c_str());
		      if (include_mantle_melting && convey_flag_mantle_melting_directory == nullptr)
			 {
				mkdir (convey_flag_mantle_melting_path.c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
				for (unsigned int m = 0; m <= melt_container_amount_x; m++)
					{
						std::stringstream init_get_name;
						init_get_name << m;
						std::string temp_txt_name = init_get_name.str();
						std::ofstream init_write((convey_flag_mantle_melting_path + temp_txt_name + ".txt").c_str());
						for (unsigned int a = 0; a <= flag_amount; a++) 
							init_write << "0" << " " ;

						init_write.close();
					}
			 } 
			 closedir(convey_flag_mantle_melting_directory);
			 
			 DIR *container_mantle_directory = opendir(container_mantle_melting_path.c_str());
		     if (include_mantle_melting && container_mantle_directory == nullptr)
			 {
				 mkdir (container_mantle_melting_path.c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
				 std::ifstream exist(container_mantle_melting_path + "0" + ".txt");
				if (!exist.good())
				{
					for (unsigned int p = 0; p < melt_container_amount_x; p++)
					 {
						 std::stringstream init_get_name;
						 init_get_name << p;
						 std::string temp_txt_name = init_get_name.str();
						 std::ofstream init_write((container_mantle_melting_path + temp_txt_name + ".txt").c_str());
						 init_write << "0" ;
						 init_write.close();
					 }
				}
			 } 
			closedir(container_mantle_directory);

			DIR *container_carbon_directory = opendir(container_carbon_melting_path.c_str());
		     if (include_carbonate_melting && container_carbon_directory == nullptr)
			 {
				 mkdir (container_carbon_melting_path.c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
				 std::ifstream exist(container_carbon_melting_path + "0" + ".txt");
				if (!exist.good())
				{
					for (unsigned int p = 0; p < melt_container_amount_x; p++)
					 {
						 std::stringstream init_get_name;
						 init_get_name << p;
						 std::string temp_txt_name = init_get_name.str();
						 std::ofstream init_write((container_carbon_melting_path + temp_txt_name + ".txt").c_str());
						 init_write << "0" ;
						 init_write.close();
					 }
				}
			 } 
			closedir(container_carbon_directory); 
			DIR *convey_flag_crustal_melting_directory = opendir(convey_flag_crustal_melting_path.c_str());
		    if (include_crustal_melting && convey_flag_crustal_melting_directory == nullptr)
			   {
				  mkdir (convey_flag_crustal_melting_path.c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
				  for (unsigned int m = 0; m <= melt_container_amount_x; m++)
					  {
						  std::stringstream init_get_name;
						  init_get_name << m;
						  std::string temp_txt_name = init_get_name.str();
						  std::ofstream init_write((convey_flag_crustal_melting_path + temp_txt_name + ".txt").c_str());
						  for (unsigned int a = 0; a <= flag_amount; a++) 
							  init_write << "0" << " " ;
						  
						  init_write.close();
					  }
			   } 
			 closedir(convey_flag_crustal_melting_directory);
			 
			 DIR *container_crustal_directory = opendir(container_crustal_melting_path.c_str());
		     if (include_crustal_melting && container_crustal_directory == nullptr)
			 {
				 mkdir (container_crustal_melting_path.c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
				 std::ifstream exist(container_crustal_melting_path + "0" + ".txt");
				if (!exist.good())
				{
					for (unsigned int p = 0; p < melt_container_amount_x; p++)
					 {
						 std::stringstream init_get_name;
						 init_get_name << p;
						 std::string temp_txt_name = init_get_name.str();
						 std::ofstream init_write((container_crustal_melting_path + temp_txt_name + ".txt").c_str());
						 init_write << "0" ;
						 init_write.close();
					 }
				}
			 } 
			closedir(container_crustal_directory);
		  }
		  
			
			DIR *melt_info_directory = opendir(melt_info_path.c_str());
			if (melt_info_directory == nullptr)
			{
				mkdir (melt_info_path.c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
				std::ifstream exist_mantle(melt_info_path + "mantle_max_melt_fraction.txt");
				if (!exist_mantle.good() && include_mantle_melting)
				{
					//当前步最大熔融程度
				    std::ofstream init_max_fraction_write((melt_info_path + "mantle_max_melt_fraction.txt").c_str());
				    init_max_fraction_write << "0" ;
				    init_max_fraction_write.close();
				    
				    //当前步提取熔体量
				    std::ofstream init_cur_ext_write((melt_info_path + "mantle_current_melt_extraction.txt").c_str());
				    init_cur_ext_write << "0" ;
				    init_cur_ext_write.close();
					
					//等温线1200°C最高点
					std::ofstream init_emplace_location_y_write((melt_info_path + "mantle_emplace_location_y.txt").c_str());
				    init_emplace_location_y_write << "0" ;
				    init_emplace_location_y_write.close();
					if(!dike_mode)
					{
						//最小地壳增生窗口
				        std::ofstream init_min_win_write((melt_info_path + "mantle_min_emplace_window.txt").c_str());
				        init_min_win_write << "0" ;
				        init_min_win_write.close();
					
					    //最大地壳增生窗口
				        std::ofstream init_max_win_write((melt_info_path + "mantle_max_emplace_window.txt").c_str());
				        init_max_win_write << "0" ;
				        init_max_win_write.close();
					}
					else
					{
						//总熔融提取量，用于确定岩墙打开时间（暂定）
						std::ofstream init_min_win_write((melt_info_path + "mantle_total_melt_extraction.txt").c_str());
				        init_min_win_write << "0" ;
				        init_min_win_write.close();
						//岩墙最大x坐标，即等温线1200°C最高点x坐标
						std::ofstream init_emplace_location_x_write((melt_info_path + "dike_emplace_location_xmax.txt").c_str());
				        init_emplace_location_x_write << "0" ;
				        init_emplace_location_x_write.close();
						//岩墙最小y坐标，即根据等温线以及x坐标确定的y坐标
						std::ofstream init_emplace_location_ymin_write((melt_info_path + "dike_emplace_location_ymin.txt").c_str());
				        init_emplace_location_ymin_write << "0" ;
				        init_emplace_location_ymin_write.close();
						//20251201 
						//岩墙最大y坐标，根据岩墙顶部等温线以及x坐标确定
						//岩墙长度不固定时启用
						std::ofstream init_emplace_location_ymax_write((melt_info_path + "dike_emplace_location_ymax.txt").c_str());
				        init_emplace_location_ymax_write << "0" ;
				        init_emplace_location_ymax_write.close();
					}
				}
				std::ifstream exist_crust(melt_info_path + "crust_max_melt_fraction.txt");
				if (!exist_crust.good() && include_crustal_melting)
				{
					//当前步最大熔融程度
				    std::ofstream init_max_fraction_write((melt_info_path + "crust_max_melt_fraction.txt").c_str());
				    init_max_fraction_write << "0" ;
				    init_max_fraction_write.close();
				    

				    //当前步提取熔体量
				    std::ofstream init_cur_ext_write((melt_info_path + "crust_current_melt_extraction.txt").c_str());
				    init_cur_ext_write << "0" ;
				    init_cur_ext_write.close();
					
					//最小地壳增生窗口
				    std::ofstream init_min_win_write((melt_info_path + "crust_min_emplace_window.txt").c_str());
				    init_min_win_write << "0" ;
				    init_min_win_write.close();
					
					//最大地壳增生窗口
				    std::ofstream init_max_win_write((melt_info_path + "crust_max_emplace_window.txt").c_str());
				    init_max_win_write << "0" ;
				    init_max_win_write.close();
					
					//等温线°C最高点（暂定600）
					std::ofstream init_emplace_location_y_write((melt_info_path + "crust_emplace_location_y.txt").c_str());
				    init_emplace_location_y_write << "0" ;
				    init_emplace_location_y_write.close();
				}
			}
			closedir(melt_info_directory);
	  }
	  
      const QGauss<dim> quadrature_formula (this->get_fe().base_element(this->introspection().base_elements.compositional_fields).degree+1);
      const unsigned int n_q_points = quadrature_formula.size();

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values   |
                               update_quadrature_points |
                               update_JxW_values);
	  MaterialModel::MaterialModelInputs<dim> in(fe_values.n_quadrature_points, this->n_compositional_fields());
	  
	  double local_mantle_max_melt_fraction = -std::numeric_limits<double>::max();
	  double local_mantle_add_melt_fraction = 0.;
	  double local_mantle_current_melt_extraction = 0.;
	  double local_mantle_emplace_location_y = 0.;
	  
	  double local_crustal_max_melt_fraction = -std::numeric_limits<double>::max();
	  double local_crustal_add_melt_fraction = 0.;
	  double local_crustal_current_melt_extraction = 0.;
	  double local_crustal_emplace_location_y = 0.;

      // compute the integral quantities by quadrature
      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          {
			//填充输入
            fe_values.reinit (cell);
			//20240620
			in.reinit(fe_values, cell, this->introspection(), this->get_solution());
			const unsigned int mid_point = n_q_points / 2;
			const double posy = in.position[mid_point][1];
			
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
            std::vector<double> melt_record_values(n_q_points, 0.);
			std::vector<double> melt_fraction_values(n_q_points, 0.);
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
			double melt_record_in_point = 0.;
			double melt_fraction_in_point = 0.;
			for (unsigned int q=0; q<n_q_points; ++q)
		    {
				//仅有地幔熔融
				if(include_mantle_melting && !include_crustal_melting)
				{
					local_crustal_max_melt_fraction = 0.;
	                local_crustal_add_melt_fraction = 0.;
	                local_crustal_current_melt_extraction = 0.;
	                local_crustal_emplace_location_y = 0.;
					//20240620
				    //定位熔融传输位置的y坐标
				    if (in.temperature[q] >= 1473 && in.temperature[q] <= 1483)
					    local_mantle_emplace_location_y = posy;
				    local_mantle_max_melt_fraction = std::max(local_mantle_max_melt_fraction, melt_fraction_values[q]);
				    local_mantle_add_melt_fraction += std::max(0., melt_record_values[q] * fe_values.JxW(q) / 1e6);
				    melt_record_in_point = melt_record_values[q];
				    melt_fraction_in_point = melt_fraction_values[q];
				    //单次最大提取量保持5%
			        if (melt_fraction_in_point > melt_limit && melt_record_in_point > 0)
					    local_mantle_current_melt_extraction += std::min(melt_record_in_point * melt_factor * fe_values.JxW(q) / 1e6, 5e-2 * fe_values.JxW(q) / 1e6);
				}
			    //仅有地壳熔融
				else if(!include_mantle_melting && include_crustal_melting)
				{
					local_mantle_max_melt_fraction = 0.;
	                local_mantle_add_melt_fraction = 0.;
	                local_mantle_current_melt_extraction = 0.;
	                local_mantle_emplace_location_y = 0.;
					//暂定600°C，根据熔融情况调整
					if (in.temperature[q] >= 873 && in.temperature[q] <= 883)
						local_crustal_emplace_location_y = posy;
					local_crustal_max_melt_fraction = std::max(local_crustal_max_melt_fraction, melt_fraction_values[q]);
				    local_crustal_add_melt_fraction += std::max(0., melt_record_values[q] * fe_values.JxW(q) / 1e6);
				    melt_record_in_point = melt_record_values[q];
				    melt_fraction_in_point = melt_fraction_values[q];
					if (melt_fraction_in_point > melt_limit && melt_record_in_point > 0)
						local_crustal_current_melt_extraction += std::min(melt_record_in_point * melt_factor * fe_values.JxW(q) / 1e6, 5e-2 * fe_values.JxW(q) / 1e6);
				}
				//混合熔融
				else if(!dike_mode && include_mantle_melting && include_crustal_melting)
				{
					//始终统计地幔熔融与地壳熔融的最高点
					if (in.temperature[q] >= 1473 && in.temperature[q] <= 1483)
					    local_mantle_emplace_location_y = posy;
					if (in.temperature[q] >= 873 && in.temperature[q] <= 883)
						local_crustal_emplace_location_y = posy;
					melt_record_in_point = melt_record_values[q];
				    melt_fraction_in_point = melt_fraction_values[q];
					
					//根据含量确定熔融类型
					if(mantle_fraction >= crust_fraction)
					{
						local_mantle_max_melt_fraction = std::max(local_mantle_max_melt_fraction, melt_fraction_values[q]);
				        local_mantle_add_melt_fraction += std::max(0., melt_record_values[q] * fe_values.JxW(q) / 1e6);
				        //单次最大提取量保持5%
			            if (melt_fraction_in_point > melt_limit && melt_record_in_point > 0)
					        local_mantle_current_melt_extraction += std::min(melt_record_in_point * melt_factor * fe_values.JxW(q) / 1e6, 5e-2 * fe_values.JxW(q) / 1e6);
					}
					else
					{
						local_crustal_max_melt_fraction = std::max(local_crustal_max_melt_fraction, melt_fraction_values[q]);
				        local_crustal_add_melt_fraction += std::max(0., melt_record_values[q] * fe_values.JxW(q) / 1e6);
						if (melt_fraction_in_point > melt_limit && melt_record_in_point > 0)
						    local_crustal_current_melt_extraction += std::min(melt_record_in_point * melt_factor * fe_values.JxW(q) / 1e6, 5e-2 * fe_values.JxW(q) / 1e6);
					}
				}
				else
				{
					local_mantle_max_melt_fraction = 0.;
	                local_mantle_add_melt_fraction = 0.;
	                local_mantle_current_melt_extraction = 0.;
	                local_mantle_emplace_location_y = 0.;
	                local_crustal_max_melt_fraction = 0.;
	                local_crustal_add_melt_fraction = 0.;
	                local_crustal_current_melt_extraction = 0.;
	                local_crustal_emplace_location_y = 0.;
				}
		    }
          }
		  
	  const double global_mantle_current_melt_extraction = Utilities::MPI::sum (local_mantle_current_melt_extraction, this->get_mpi_communicator());
	  const double global_mantle_add_melt_fraction = Utilities::MPI::sum (local_mantle_add_melt_fraction, this->get_mpi_communicator());
	  const double global_mantle_max_melt_fraction = Utilities::MPI::max (local_mantle_max_melt_fraction, this->get_mpi_communicator());
	  const double global_mantle_max_emplace_location_y = Utilities::MPI::max (local_mantle_emplace_location_y, this->get_mpi_communicator());
	  
	  const double global_crustal_current_melt_extraction = Utilities::MPI::sum (local_crustal_current_melt_extraction, this->get_mpi_communicator());
	  const double global_crustal_add_melt_fraction = Utilities::MPI::sum (local_crustal_add_melt_fraction, this->get_mpi_communicator());
	  const double global_crustal_max_melt_fraction = Utilities::MPI::max (local_crustal_max_melt_fraction, this->get_mpi_communicator());
	  const double global_crustal_max_emplace_location_y = Utilities::MPI::max (local_crustal_emplace_location_y, this->get_mpi_communicator());
	  
	  if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
	  {
		  if(include_mantle_melting)
		  {
			  std::ofstream mantle_write_max_melt((melt_info_path + "mantle_max_melt_fraction.txt").c_str());
			  if(mantle_write_max_melt)
				  mantle_write_max_melt << global_mantle_max_melt_fraction;
			  mantle_write_max_melt.close();
			  
			  std::ofstream mantle_write_curr_extract((melt_info_path + "mantle_current_melt_extraction.txt").c_str());
			  if(mantle_write_curr_extract)
				  mantle_write_curr_extract << global_mantle_current_melt_extraction;
			  mantle_write_curr_extract.close();
			  
			  const double mantle_location_y_in_km = global_mantle_max_emplace_location_y / 1e3;
			  std::ofstream mantle_write_location_y((melt_info_path + "mantle_emplace_location_y.txt").c_str());
			  if(mantle_write_location_y)
				  mantle_write_location_y << mantle_location_y_in_km;
			  mantle_write_location_y.close();
			  //20251104
			  //固定M值不需要记录最大熔融，暂时默认所有岩墙模式都记录熔融情况
			  if(dike_mode)
			  {
				  const double temp_mantle_total_melt_extraction = read_from_txt(melt_info_path + "mantle_total_melt_extraction.txt", 0);
				  const double out_mantle_total_melt_extraction = temp_mantle_total_melt_extraction + global_mantle_current_melt_extraction;
				  std::ofstream mantle_write_total_melt_extract((melt_info_path + "mantle_total_melt_extraction.txt").c_str());
			      if(mantle_write_total_melt_extract)
				      mantle_write_total_melt_extract << out_mantle_total_melt_extraction;
			      mantle_write_total_melt_extract.close();
			  }
		  }
		  
		  if(include_crustal_melting)
		  {
			  std::ofstream crustal_write_max_melt((melt_info_path + "crust_max_melt_fraction.txt").c_str());
			  if(crustal_write_max_melt)
				  crustal_write_max_melt << global_crustal_max_melt_fraction;
			  crustal_write_max_melt.close();
			  
			  std::ofstream crustal_write_curr_extract((melt_info_path + "crust_current_melt_extraction.txt").c_str());
			  if(crustal_write_curr_extract)
				  crustal_write_curr_extract << global_crustal_current_melt_extraction;
			  crustal_write_curr_extract.close();
			  
			  const double crustal_location_y_in_km = global_crustal_max_emplace_location_y / 1e3;
			  std::ofstream crustal_write_location_y((melt_info_path + "crust_emplace_location_y.txt").c_str());
			  if(crustal_write_location_y)
				  crustal_write_location_y << crustal_location_y_in_km;
			  crustal_write_location_y.close();
		  }
	  }
	  
		  
      statistics.add_value ("Max_melt_fraction(mantle)", global_mantle_max_melt_fraction);
	  statistics.add_value ("Current_melt_extraction(mantle, km^2)", global_mantle_current_melt_extraction);
	  statistics.add_value ("Current_delta_melt_fraction(mantle, km^2)", global_mantle_add_melt_fraction);
	  statistics.add_value ("Max_melt_fraction(crustal)", global_mantle_max_melt_fraction);
	  statistics.add_value ("Current_melt_extraction(crustal, km^2)", global_mantle_current_melt_extraction);
	  statistics.add_value ("Current_delta_melt_fraction(crustal, km^2)", global_mantle_add_melt_fraction);
	  {
        const char *columns[] = {"Max_melt_fraction(mantle)", "Current_melt_extraction(mantle, km^2)", "Current_delta_melt_fraction(mantle, km^2)", "Max_melt_fraction(crustal)", "Current_melt_extraction(crustal, km^2)", "Current_delta_melt_fraction(crustal, km^2)"};
        for (unsigned int i=0; i<sizeof(columns)/sizeof(columns[0]); ++i)
          {
            statistics.set_precision (columns[i], 8);
            statistics.set_scientific (columns[i], true);
          }
      }
      std::ostringstream output;
      output.precision(4);
	  output << global_mantle_max_melt_fraction << "/"
	         << global_mantle_current_melt_extraction << "/"
	         << global_mantle_add_melt_fraction << "/"
			 << global_crustal_max_melt_fraction << "/"
			 << global_crustal_current_melt_extraction << "/"
			 << global_crustal_add_melt_fraction << "/";
      return std::pair<std::string, std::string> ("Pre-Processing of melt extraction, max melt fraction(mantle) / current melt extraction(mantle, km^2) / current_delta_melt_fraction(mantle, km^2) /max melt fraction(crustal) / current melt extraction(crustal, km^2) / current_delta_melt_fraction(crustal, km^2),",
                                                  output.str());
    }
	
	template <int dim>
    void
    PreMeltExtraction<dim>::declare_parameters (ParameterHandler &prm)
	{
		prm.enter_subsection("Postprocess");
		{
			prm.enter_subsection("Melt extraction");
			{
				prm.declare_entry ("Max depth of melt convey function", "20e3",
                                   Patterns::Double (0.),
                                   "Max depth of melt convey function"
                                   "Units: meters.");
				prm.declare_entry ("Max topography above the initial surface", "10e3",
                                   Patterns::Double (0.),
                                   "Max topography above the initial surface"
                                   "Units: meters.");
				prm.declare_entry ("Limit melt fraction of extracting", "1e-5",
                                   Patterns::Double (0.),
                                   "Limit melt fraction of extracting"
                                   "Units: none.");
				prm.declare_entry ("Coefficient of melt extraction", "0.7",
                                   Patterns::Double (0.),
                                   "Coefficient of melt extraction"
                                   "Units: none.");
				prm.declare_entry ("Vertical transportation mode", "false",
                                   Patterns::Bool (),
                                   "Whether or not to choose vertical transportation mode. "
                                   "Units: /.");
				prm.declare_entry ("Container limit of new crust generation", "0.",
                                   Patterns::Double (0.),
                                   "Container limit of melt emplacement. "
								   "E.g. If this value is set to 2, the function of new crust generation "
								   "will start until the melt volume is more than 2 km^3/km in the container. "
                                   "Units: km^3/km.");
			   //20251104
			   prm.declare_entry ("Use dike mode", "false",
                                   Patterns::Bool (),
                                   "Whether or not to use dike mode. "
                                   "Units: /.");
			   prm.declare_entry ("Constant M", "false",
                                   Patterns::Bool (),
                                   "Whether or not to use constant M. "
                                   "Units: /.");
			   prm.declare_entry ("M for dike", "0.5",
                                   Patterns::Double (0.),
                                   "M for dike. "
                                   "Units: /.");
			   prm.declare_entry ("Temperature of the bottom of dike(K)", "923",
                                   Patterns::Double (0.),
                                   "Temperature of the bottom of dike. "
                                   "Units: Kelvin");
			   prm.declare_entry ("Width of dike(m)", "1.25e3",
                                   Patterns::Double (0.),
                                   "Width of dike. "
                                   "Units: m");
			   prm.declare_entry ("Length of dike(m)", "5e3",
                                   Patterns::Double (0.),
                                   "Length of dike. "
                                   "Units: m");
			   prm.declare_entry ("Full spreading rate of model(m/yr)", "0.02",
                                   Patterns::Double (0.),
                                   "Full spreading rate of model, used in calculating phi. "
                                   "Units: m");
			  //20251201
			  prm.declare_entry ("Constant dike length", "false",
                                   Patterns::Bool (),
                                   "Whether or not to use constant dike length. "
                                   "Units: /.");
			  prm.declare_entry ("Constant dike X position", "false",
                                   Patterns::Bool (),
                                   "Whether or not to use constant dike X position. "
                                   "Units: /.");
			  prm.declare_entry ("Temperature of the top of dike(K)", "283",
                                   Patterns::Double (0.),
                                   "Temperature of the top of dike, used when 'Constant dike length' is set false. "
                                   "Units: Kelvin");
			}
			prm.leave_subsection ();
		}
		prm.leave_subsection ();
	}
	
	
	template <int dim>
    void
    PreMeltExtraction<dim>::parse_parameters (ParameterHandler &prm)
	{
		prm.enter_subsection("Postprocess");
		{
			prm.enter_subsection("Melt extraction");
			{
				depth_max = prm.get_double ("Max depth of melt convey function");
				y_add = prm.get_double ("Max topography above the initial surface");
				melt_limit = prm.get_double ("Limit melt fraction of extracting");
				melt_factor = prm.get_double ("Coefficient of melt extraction");
				//20251104
				dike_mode = prm.get_bool("Use dike mode");
				
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
				melt_container_amount_x = std::ceil(max_x / ((max_x / repetitions_x) / std::pow(2, melt_fraction_refinement)));
				melt_container_amount_y = std::ceil((max_y + y_add) / ((max_y / repetitions_y) / std::pow(2, melt_fraction_refinement)));
				flag_amount = melt_container_amount_y - min_container_amount_y;
				convey_flag_mantle_melting_path = this->get_output_directory() + "convey_flag_mantle_melting/";
				container_mantle_melting_path = this->get_output_directory() + "container_mantle_melting/";
				convey_flag_crustal_melting_path = this->get_output_directory() + "convey_flag_crustal_melting/";
				container_crustal_melting_path = this->get_output_directory() + "container_crustal_melting/";
				melt_info_path = this->get_output_directory() + "melt_info/";
			    //20250428
				container_carbon_melting_path = this->get_output_directory() + "container_carbonate_melting/";
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
			    //20250428
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
    ASPECT_REGISTER_POSTPROCESSOR(PreMeltExtraction,
                                  "pre melt extraction",
                                  "A postprocessor that computes some statistics about "
                                  "the compositional fields, if present in this simulation. "
                                  "In particular, it computes maximal and minimal values of "
                                  "each field, as well as the total mass contained in this "
                                  "field as defined by the integral "
                                  "$m_i(t) = \\int_\\Omega c_i(\\mathbf x,t) \\; \\text{d}x$.")
  }
}
