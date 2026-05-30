project(
  'fullpos',
  ['c', 'fortran'],
  version: '0.1.0a0',
  default_options: ['warning_level=2', 'buildtype=release', 'b_lto=false'],
)

py_mod = import('python')
fs = import('fs')
py = py_mod.find_installation(pure: false)
numpy_inc = run_command(
  py,
  ['-c', 'import numpy; print(numpy.get_include())'],
  check: true,
).stdout().strip()
fortran_compiler = meson.get_compiler('fortran')

native_prefix_opt = get_option('fullpos_native_prefix')
if native_prefix_opt == ''
  fullpos_native_prefix = meson.current_source_dir() / 'extern' / 'fullpos' / 'local'
  ectrans_inc = include_directories(
    numpy_inc,
    'extern/fullpos/local/include',
  )
else
  fullpos_native_prefix = native_prefix_opt
  ectrans_inc = include_directories(numpy_inc)
endif
native_libdir_name = 'lib'
if host_machine.system() != 'windows' and fs.exists(fullpos_native_prefix / 'lib64')
  native_libdir_name = 'lib64'
endif
ectrans_includedir = fullpos_native_prefix / 'include'
ectrans_libdir = fullpos_native_prefix / native_libdir_name

ectrans_c_args = [
  '-I' + ectrans_includedir,
]
ectrans_link_args = [
  '-L' + ectrans_libdir,
  '-ltransi_dp',
  '-ltrans_dp',
  '-lfiat',
  '-lparkind_dp',
  '-lparkind_sp',
  '-lopenblas',
]

if host_machine.system() == 'windows'
  ectrans_link_args += ['-lws2_32']
endif

if host_machine.system() == 'windows'
  native_runtime_source_dir = fullpos_native_prefix / 'bin'
  native_runtime_install_dir = py.get_install_dir() / 'fullpos' / '_native' / 'bin'
  native_install_rpath = ''
  native_runtime_libraries = [
    'libtransi_dp.dll',
    'libtrans_dp.dll',
    'libfiat.dll',
    'libparkind_dp.dll',
    'libparkind_sp.dll',
  ]
elif host_machine.system() == 'darwin'
  native_runtime_source_dir = fullpos_native_prefix / 'lib'
  native_runtime_install_dir = py.get_install_dir() / 'fullpos' / '_native' / 'lib'
  native_install_rpath = '@loader_path/_native/lib'
  native_runtime_libraries = [
    'libtransi_dp.dylib',
    'libtrans_dp.dylib',
    'libfiat.dylib',
    'libparkind_dp.dylib',
    'libparkind_sp.dylib',
  ]
else
  native_runtime_source_dir = fullpos_native_prefix / native_libdir_name
  native_runtime_install_dir = py.get_install_dir() / 'fullpos' / '_native' / 'lib'
  native_install_rpath = '$ORIGIN/_native/lib'
  native_runtime_libraries = [
    'libtransi_dp.so',
    'libtrans_dp.so',
    'libfiat.so',
    'libparkind_dp.so',
    'libparkind_sp.so',
  ]
endif

fullpos_fortran_args = [
  '-I' + fullpos_native_prefix / 'module' / 'parkind_dp',
  '-I' + fullpos_native_prefix / 'module' / 'fiat',
  '-I' + fullpos_native_prefix / 'include' / 'fiat',
]

if fortran_compiler.get_id() == 'gcc'
  fullpos_fortran_args += ['-ffree-line-length-none']
  fullpos_fortran_args += ['-fopenmp']
  ectrans_link_args += ['-fopenmp']
endif

py.extension_module(
  '_ectrans',
  [
    'src/fullpos/_ectrans.c',
    'src/fullpos/_fpint_wrappers.F90',
    'extern/fullpos/upstream/arpifs/interpol/fpint4.F90',
    'extern/fullpos/upstream/arpifs/interpol/fpint12.F90',
    'extern/fullpos/upstream/arpifs/interpol/fpavg.F90',
    'extern/fullpos/upstream/arpifs/interpol/fpnear.F90',
    'extern/fullpos/upstream/arpifs/interpol/suhow1.F90',
    'extern/fullpos/upstream/arpifs/interpol/suhow2.F90',
    'extern/fullpos/upstream/arpifs/interpol/suhox1.F90',
  ],
  include_directories: ectrans_inc,
  c_args: ectrans_c_args,
  fortran_args: fullpos_fortran_args,
  link_args: ectrans_link_args,
  build_rpath: ectrans_libdir,
  install_rpath: native_install_rpath,
  install_dir: py.get_install_dir() / 'fullpos',
  install: true,
)

py.extension_module(
  '_vertical_native',
  [
    'src/fullpos/_vertical_native.c',
    'src/fullpos/_vertical_support.F90',
    'src/fullpos/_vertical_native_wrappers.F90',
    'extern/fullpos/upstream/arpifs/pp_obs/ppinit.F90',
    'extern/fullpos/upstream/arpifs/pp_obs/ppflev.F90',
    'extern/fullpos/upstream/arpifs/pp_obs/ppintp.F90',
    'extern/fullpos/upstream/arpifs/pp_obs/ppitpq.F90',
    'extern/fullpos/upstream/arpifs/pp_obs/ppsta.F90',
    'extern/fullpos/upstream/arpifs/pp_obs/ppq.F90',
    'extern/fullpos/upstream/arpifs/pp_obs/ppuv.F90',
    'extern/fullpos/upstream/arpifs/pp_obs/ppt.F90',
    'extern/fullpos/upstream/arpifs/pp_obs/ppt_old.F90',
    'extern/fullpos/upstream/arpifs/pp_obs/apache.F90',
    'extern/fullpos/upstream/arpifs/fullpos/fpview.F90',
    'extern/fullpos/upstream/arpifs/pp_obs/pp2dint.F90',
    'extern/fullpos/upstream/arpifs/pp_obs/ppgeop.F90',
    'extern/fullpos/upstream/arpifs/pp_obs/ppleta.F90',
    'extern/fullpos/upstream/arpifs/pp_obs/pplteta.F90',
    'extern/fullpos/upstream/arpifs/pp_obs/ppltp.F90',
    'extern/fullpos/upstream/arpifs/pp_obs/ppltw.F90',
    'extern/fullpos/upstream/arpifs/adiab/gphpre.F90',
    'extern/fullpos/upstream/arpifs/adiab/gpgeo.F90',
    'extern/fullpos/upstream/arpifs/adiab/gptet.F90',
    'extern/fullpos/upstream/arpifs/adiab/gprcp.F90',
    'extern/fullpos/upstream/arpifs/adiab/gppvo.F90',
    'extern/fullpos/upstream/arpifs/fullpos/fpps.F90',
    'extern/fullpos/upstream/arpifs/pp_obs/pppmer.F90',
    'extern/fullpos/upstream/ifsaux/support/isrchflt.F',
  ],
  include_directories: include_directories(
    numpy_inc,
    'src/fullpos',
    'extern/fullpos/upstream/ifsaux/include',
  ),
  fortran_args: fullpos_fortran_args,
  install_rpath: native_install_rpath,
  install_dir: py.get_install_dir() / 'fullpos',
  install: true,
)

build_config = configuration_data()
build_config.set_quoted('NATIVE_PREFIX', fullpos_native_prefix)
configure_file(
  input: 'src/fullpos/_build_config.py.in',
  output: '_build_config.py',
  configuration: build_config,
  install_dir: py.get_install_dir() / 'fullpos',
)

py.install_sources(
  [
    'src/fullpos/__init__.py',
    'src/fullpos/api.py',
    'src/fullpos/config.py',
    'src/fullpos/diagnostics.py',
    'src/fullpos/errors.py',
    'src/fullpos/filtering.py',
    'src/fullpos/grib.py',
    'src/fullpos/grids.py',
    'src/fullpos/masks.py',
    'src/fullpos/metadata.py',
    'src/fullpos/native.py',
    'src/fullpos/reduced.py',
    'src/fullpos/regridder.py',
    'src/fullpos/spectral.py',
    'src/fullpos/surface.py',
    'src/fullpos/vertical.py',
    'src/fullpos/xarray.py',
  ],
  subdir: 'fullpos',
)

py.install_sources(
  [
    'src/fullpos/filters/__init__.py',
    'src/fullpos/filters/core.py',
    'src/fullpos/filters/io.py',
    'src/fullpos/filters/profiles.py',
    'src/fullpos/filters/spectral.py',
  ],
  subdir: 'fullpos/filters',
)

py.install_sources(
  [
    'src/fullpos/_vertical/__init__.py',
    'src/fullpos/_vertical/api.py',
    'src/fullpos/_vertical/common.py',
    'src/fullpos/_vertical/diagnostics.py',
    'src/fullpos/_vertical/eta.py',
    'src/fullpos/_vertical/height.py',
    'src/fullpos/_vertical/model_level.py',
    'src/fullpos/_vertical/pos.py',
    'src/fullpos/_vertical/potential_temperature.py',
    'src/fullpos/_vertical/potential_vorticity.py',
    'src/fullpos/_vertical/pressure.py',
    'src/fullpos/_vertical/temperature.py',
    'src/fullpos/_vertical/validation.py',
  ],
  subdir: 'fullpos/_vertical',
)

py.install_sources(
  [
    'src/fullpos/interpolation/__init__.py',
    'src/fullpos/interpolation/api.py',
    'src/fullpos/interpolation/horizontal.py',
    'src/fullpos/interpolation/kernels.py',
    'src/fullpos/interpolation/native.py',
    'src/fullpos/interpolation/regridder.py',
    'src/fullpos/interpolation/surface.py',
    'src/fullpos/interpolation/xarray.py',
  ],
  subdir: 'fullpos/interpolation',
)

if get_option('fullpos_bundle_native_runtime')
  foreach name : native_runtime_libraries
    candidate = native_runtime_source_dir / name
    if fs.exists(candidate)
      install_data(candidate, install_dir: native_runtime_install_dir)
    endif
  endforeach
endif
