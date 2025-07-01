program only_analysis
  use analysis
  implicit none

  character(256) :: model, dir

  real(8) :: time_map
  logical :: old_format, format_2d

  !time_map = 2.926d-2
  !model = "SFHoTim276_13_14_0025_150mstg_B0_HLLC"
  !dir = "/sakura/ptmp/shofu/SFHoTim276_13_14_0025_150mstg_B0_HLLC/Analysis_ptr/data_std"

  ! time_map = 3.096E-2; old_format=.true.; format_2d=.false.
  ! model = "SFHoTim276_12_15_0025_150mstg_B0_HLLC_rerun"
  ! dir = "/sakura/ptmp/shofu/SFHoTim276_12_15_0025_150mstg_B0_HLLC_rerun/Analysis_ptr/data"

  ! time_map = 2.715d-2; old_format=.true.; format_2d=.false.
  ! model = "SFHoTim276_135_135_0025_150mstg_B0_HLLC_rerun"
  ! dir = "/sakura/ptmp/shofu/SFHoTim276_135_135_0025_150mstg_B0_HLLC_rerun/Analysis_ptr/data"

  ! time_map = 3.027d-2; old_format=.true.; format_2d=.false.
  ! model = "SFHoTim276_125_145_0025_150mstg_B0_HLLC_rerun"
  ! dir = "/sakura/ptmp/shofu/SFHoTim276_125_145_0025_150mstg_B0_HLLC_rerun/Analysis_ptr/data"
  
  ! time_map = 1.871d-02; old_format=.true.; format_2d=.false.
  ! model = "SFHoTim276_125_155_0025_150mstg_B0_HLLC_inspiral"
  ! dir = "/sakura/ptmp/shofu/SFHoTim276_125_155_0025_150mstg_B0_HLLC_inspiral/Analysis_ptr/data_std"

  !time_map = 1.871d-02; old_format=.true.; format_2d=.false.
  !model = "DD2Tim326_Q4_M135_a75_0056_270m_B5e16_Hon5png"
  !dir = "/sakura/ptmp/shofu/DD2Tim326_Q4_M135_a75_0056_270m_B5e16_Hon5png/Analysis_ptr/data_small_set/data_post"
  !dir = "/sakura/ptmp/shofu/DD2Tim326_Q4_M135_a75_0056_270m_B5e16_Hon5png/Analysis_ptr/data_small_set/data_dyn_unified"
  !dir = "/sakura/ptmp/shofu/DD2Tim326_Q4_M135_a75_0056_270m_B5e16_Hon5png/Analysis_ptr/data_3e8km"; time_map = 0d0 ; old_format=.false.

  !time_map = 1.871d-02; old_format=.false.; format_2d=.false.
  model = "DD2Tim326_135_135_0028_12.5mstg_B15.5_HLLD_CT_GS"
  dir = "/scratch/sfujibayashi/DD2Tim326_135_135_0028_12.5mstg_B15.5_HLLD_CT_GS/Analysis_ptr/data_3e8cm/"
  
  call tr_analysis(model,dir,old_format=old_format,format_2d=format_2d, time_map=time_map)
  stop
  

  format_2d = .true.; old_format=.true.

  ! !!! MHD+dynamo seriese
  ! model = "DD2_135H_DDMHD4p-5P75-4a3"
  ! model="DD2_135H_DDMHD4p-5P7-49s"
  ! model="DD2_135H_DDMHD4p-5P7-46s"
  
  !! model="DD2_135H_DDMHD4p-5P75-57s"
  ! model="DD2_135H_DDMHD4p-5P8-4s"
  ! model="DD2_135H_DDMHD4p-5P75-4s"
  ! model="DD2_135H_DDMHD4p-5P75-43s"
  model="DD2_135H_DDMHD4p-5P7-4s"
  ! model = "DD2_135H_DDMHD4p-5P7-43s"

  dir = "/sakura/ptmp/shofu/"//trim(model)//"/Analysis_ptr/data"

  ! model = "SFHo_13_14_150m_t10_070_01000_15_0.04_1.0e4_vcut13_RegCowl"
  ! dir = "/sakura/ptmp/shofu/SFHo_13_14_150m_t10_070_01000_15_0.04_1.0e4_vcut13_RegCowl/Analysis_ptr/data"

  !model = "SFHo_12_15_150m_t20_070_01000_15_0.04_1.0e4_vcut13_RegCowl"
  !dir = "/sakura/ptmp/shofu/SFHo_12_15_150m_t20_070_01000_15_0.04_1.0e4_vcut13_RegCowl/Analysis_ptr/data"

  ! model = "SFHo_135_135_150m_t9_070_01000_15_0.04_1.0e4_vcut12_RegCowl"
  ! dir = "/sakura/ptmp/shofu/SFHo_135_135_150m_t9_070_01000_15_0.04_1.0e4_vcut12_RegCowl/Analysis_ptr/data_std"

  ! model = "SFHo_125_145_150m_t8_070_01000_15_0.04_1.0e4_vcut13_RegCowl"
  ! dir = "/sakura/ptmp/shofu/SFHo_125_145_150m_t8_070_01000_15_0.04_1.0e4_vcut13_RegCowl/Analysis_ptr/data"

  ! model = "SFHo_125_155_150m_t9_070_01000_15_0.04_1.0e4_vcut13_RegCowl"
  ! dir = "/sakura/ptmp/shofu/SFHo_125_155_150m_t9_070_01000_15_0.04_1.0e4_vcut13_RegCowl/Analysis_ptr/data"
  
  call tr_analysis(model,dir,old_format=old_format,format_2d=format_2d, time_map=time_map)

end program only_analysis
