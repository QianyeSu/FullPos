# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

##################################################################
## Project dependencies

include( CMakeFindDependencyMacro )

set( fiat_VERSION_STR 1.0.0 )
set( fiat_HAVE_MPI    0 )
set( fiat_HAVE_OMP    True )
set( fiat_HAVE_FCKIT   )
set( fiat_SOURCE_FILENAMES      drhook.c;crc.c;addrdiff.c;countingsort.c;ecqsort.c;gnomesort.c;rsort32.c;rsort64.c;fiat_constructor.c;ec_mpi_atexit.c;ec_get_cycles.c;getcurheap.c;gethwm.c;getmaxrss.c;getpag.c;getrss.c;getstackusage.c;getstatm.c;getstk.c;fnec.c;get_tcmalloc_info.c;linux_bind.c;linuxtrbk.c;memory_hook.c;opfla_perfmon.c;pthread_attr_init.c;tabort.c;wrap_ftn.c;abor1_c.c;bytes_io.c;ec_args.c;ec_datetime.c;ec_endian.c;ec_env.c;ec_exit.c;ec_raise.c;ec_set_umask.c;ecmpi_version.c;ecmwf_transfer.c;ecomp_version.c;loc_addr.c;dr_hook_end.F90;dr_hook_init.F90;dr_hook_watch_mod.F90;cdrhookinit.F90;dr_hack_mod.F90;dr_hook_procinfo.F90;dr_hook_prt.F90;dr_hook_stackcheck_mod.F90;dr_hook_util.F90;dr_hook_util_multi.F90;drhook_run_omp_parallel.F90;yomhook.F90;ecsort_mix.F90;gstats.F90;gstats_barrier.F90;gstats_barrier2.F90;gstats_label.F90;gstats_print.F90;gstats_psut.F90;gstats_setup.F90;yomgstats.F90;ec_mpi_finalize.F90;mpi4to8.F90;mpi4to8_m.F90;mpi4to8_s.F90;mpl_abort_mod.F90;mpl_allgather_mod.F90;mpl_allgatherv_mod.F90;mpl_allreduce_mod.F90;mpl_alltoallv_mod.F90;mpl_arg_mod.F90;mpl_barrier_mod.F90;mpl_broadcast_mod.F90;mpl_buffer_method_mod.F90;mpl_bytes_mod.F90;mpl_close_mod.F90;mpl_comm_create_mod.F90;mpl_comm_free_mod.F90;mpl_comm_split_mod.F90;mpl_data_module.F90;mpl_end_mod.F90;mpl_gatherv_mod.F90;mpl_groups.F90;mpl_init_mod.F90;mpl_ioinit_mod.F90;mpl_locomm_create_mod.F90;mpl_message_mod.F90;mpl_mpif.F90;mpl_mygatherv_mod.F90;mpl_myrank_mod.F90;mpl_nproc_mod.F90;mpl_open_mod.F90;mpl_probe_mod.F90;mpl_read_mod.F90;mpl_recv_mod.F90;mpl_scatterv_mod.F90;mpl_send_mod.F90;mpl_setdflt_comm_mod.F90;mpl_stats_mod.F90;mpl_testsome_mod.F90;mpl_tour_table_mod.F90;mpl_wait_mod.F90;mpl_waitany_mod.F90;mpl_write_mod.F90;yommplstats.F90;mpl_bindc.F90;mpl_module.F90;oml_mod.F90;getheapstat.F90;getmemstat.F90;getmemvals.F90;gentrbk.F90;sdl_mod.F90;abor1.F90;bytes_io_mod.F90;ec_args_mod.F90;ec_datetime_mod.F90;ec_env_mod.F90;ec_flush.F90;ec_khz.F90;ec_lun.F90;ec_meminfo.F90;ec_parkind.F90;ec_pmon.F90;getopt.F90;cptime.F90;get_openmp.F90;qsortc.F90;strhandler_mod.F90;timef.F90;user_clock.F90;cxxdemangle.cc;parkind1.F90;parkind2.F90 )

if( fiat_HAVE_OMP AND NOT TARGET OpenMP::OpenMP_Fortran )
    if( NOT CMAKE_Fortran_COMPILER_LOADED )
        enable_language( Fortran )
    endif()
    find_dependency( OpenMP COMPONENTS Fortran )
endif()

if( fiat_HAVE_FCKIT AND NOT TARGET fckit )
  find_dependency( fckit HINTS ${CMAKE_CURRENT_LIST_DIR}/../fckit  )
endif()

if( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(_whole_archive "-Wl,-force_load")
  set(_no_whole_archive "")
else()
  set(_whole_archive "-Wl,--whole-archive")
  set(_no_whole_archive "-Wl,--no-whole-archive")
endif()
set(MPI_SERIAL_LIBRARIES ${_whole_archive} mpi_serial ${_no_whole_archive})
