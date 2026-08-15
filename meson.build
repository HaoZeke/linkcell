project(
  'linkcell',
  ['c', 'cpp'],
  version: '0.1.0',
  meson_version: '>=1.2.0',
  license: 'MIT',
  default_options: [
    'buildtype=debugoptimized',
    'warning_level=2',
    'c_std=c11',
    'cpp_std=c++17',
  ],
)

cc = meson.get_compiler('c')
cargo_prog = find_program('cargo', required: true)
py = find_program('python3', 'python')
_install = not meson.is_subproject()
_msproot = meson.project_source_root()

if host_machine.system() == 'windows'
  cargo_lib_name = 'linkcell.lib'
else
  cargo_lib_name = 'liblinkcell.a'
endif

linkcell_cargo = custom_target(
  'linkcell-cargo',
  input: files(
    'Cargo.toml',
    'Cargo.lock',
    'src/lib.rs',
    'src/knearest.rs',
    'src/cell.rs',
    'src/error.rs',
    'src/capi.rs',
  ),
  output: cargo_lib_name,
  command: [
    py, '-c',
    'import subprocess, shutil, sys; ' +
    'subprocess.check_call([sys.argv[1], "build", "--release", "--features", "capi", "--manifest-path", sys.argv[2], "--target-dir", sys.argv[3] + "/cargo-target"]); ' +
    'shutil.copy2(sys.argv[3] + "/cargo-target/release/" + sys.argv[4], sys.argv[5])',
    cargo_prog.full_path(), '@INPUT0@', '@PRIVATE_DIR@', cargo_lib_name, '@OUTPUT0@',
  ],
  build_by_default: true,
  install: _install,
  install_dir: get_option('libdir'),
)

rust_sys = [dependency('threads')]
if host_machine.system() == 'linux'
  rust_sys += [
    cc.find_library('dl', required: true),
    cc.find_library('m', required: true),
  ]
elif host_machine.system() == 'darwin'
  rust_sys += [cc.find_library('m', required: true)]
elif host_machine.system() == 'windows'
  rust_sys += [
    cc.find_library('ws2_32', required: true),
    cc.find_library('userenv', required: true),
    cc.find_library('bcrypt', required: true),
    cc.find_library('ntdll', required: true),
  ]
endif

# link_args rather than link_with: Meson promotes installed static
# consumers to link_whole, and link_whole rejects custom_targets.
linkcell_dep = declare_dependency(
  link_args: [linkcell_cargo.full_path()],
  sources: linkcell_cargo,
  include_directories: include_directories('include'),
  dependencies: rust_sys,
)

if _install
  install_headers('include/linkcell.h', 'include/linkcell.hpp')
  pkg = import('pkgconfig')
  pkg.generate(
    libraries: [linkcell_cargo] + rust_sys,
    name: 'linkcell',
    filebase: 'linkcell',
    description: 'Periodic linked-cell k-nearest neighbour search',
    version: meson.project_version(),
    url: 'https://github.com/HaoZeke/linkcell',
    subdirs: '.',
  )
endif

if get_option('with_tests')
  test(
    'cargo-tests',
    cargo_prog,
    args: ['test', '--manifest-path', _msproot / 'Cargo.toml'],
    timeout: 300,
  )
endif

if get_option('with_examples')
  subdir('examples')
endif
