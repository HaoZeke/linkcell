project(
  'linkcell',
  ['c'],
  version: '0.1.0',
  meson_version: '>=1.2.0',
  default_options: [
    'buildtype=debugoptimized',
    'warning_level=2',
    'c_std=c11',
  ],
)

cc = meson.get_compiler('c')
cargo_prog = find_program('cargo', required: true)
py = find_program('python3', 'python')

if host_machine.system() == 'windows'
  cargo_lib_name = 'linkcell.lib'
else
  cargo_lib_name = 'liblinkcell.a'
endif

linkcell_cargo = custom_target(
  'linkcell-cargo',
  input: files('Cargo.toml', 'src/lib.rs', 'src/knearest.rs', 'src/cell.rs', 'src/error.rs', 'src/capi.rs'),
  output: cargo_lib_name,
  command: [
    py, '-c',
    'import subprocess, shutil, sys; ' +
    'subprocess.check_call([sys.argv[1], "build", "--release", "--features", "capi", "--manifest-path", sys.argv[2], "--target-dir", sys.argv[3] + "/cargo-target"]); ' +
    'shutil.copy2(sys.argv[3] + "/cargo-target/release/" + sys.argv[4], sys.argv[5])',
    cargo_prog.full_path(), '@INPUT0@', '@PRIVATE_DIR@', cargo_lib_name, '@OUTPUT0@',
  ],
  build_by_default: true,
)

rust_sys = [dependency('threads')]
if host_machine.system() == 'linux'
  rust_sys += [
    cc.find_library('dl', required: true),
    cc.find_library('m', required: true),
  ]
endif

linkcell_dep = declare_dependency(
  sources: linkcell_cargo,
  include_directories: include_directories('include'),
  dependencies: rust_sys,
)
