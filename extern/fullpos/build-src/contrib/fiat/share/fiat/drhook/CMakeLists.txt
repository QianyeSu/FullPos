list( APPEND files
  drhook_merge_walltime.pl
  drhook_merge_walltime_max.pl
  drhook_merge_walltime_total_max.pl
)
set( destination share/fiat/drhook )

install(
    FILES       ${files}
    DESTINATION ${destination}
    PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)

file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/${destination})

foreach( file ${files} )
  if( WIN32 )
    execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "copy_if_different"
        "${CMAKE_CURRENT_SOURCE_DIR}/${file}"
        "${CMAKE_BINARY_DIR}/${destination}/${file}")
  else()
    execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink"
        "${CMAKE_CURRENT_SOURCE_DIR}/${file}"
        "${CMAKE_BINARY_DIR}/${destination}/${file}")
  endif()
endforeach()
