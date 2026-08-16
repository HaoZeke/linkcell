project(
  'linkcell',
  ['c', 'cpp'],
  # Keep in lockstep with Cargo.toml [package].version (source of truth).
  version: '0.2.3',
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
_mbproot = meson.project_build_root()

if host_machine.system() == 'windows'
  cargo_shared = 'linkcell.dll'
  cargo_static = 'linkcell.lib'
else
  _so = host_machine.system() == 'darwin' ? '.dylib' : '.so'
  cargo_shared = 'liblinkcell' + _so
  cargo_static = 'liblinkcell.a'
endif

_profile = get_option('buildtype') == 'debug' ? 'debug' : 'release'

# cargo rustc emits both crate types; the target lists both so Meson
# does not pretend the static archive is the only artifact.
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
  output: [cargo_shared, cargo_static],
  command: [
    py, '-c',
    'import shutil, subprocess, sys; from pathlib import Path; ' +
    'cargo, src, target, profile, shared, archive, out_s, out_a = sys.argv[1:9]; ' +
    'cmd = [cargo, "rustc", "--lib", "--manifest-path", str(Path(src) / "Cargo.toml"), ' +
    '"--target-dir", target, "--features", "capi"]; ' +
    'Path(src, "Cargo.lock").is_file() and cmd.append("--locked"); ' +
    'profile == "release" and cmd.append("--release"); ' +
    'rustc = ["--crate-type=cdylib", "--crate-type=staticlib"]; ' +
    'sys.platform.startswith("linux") and rustc.append("-Clink-arg=-Wl,-soname," + shared); ' +
    'sys.platform == "darwin" and rustc.append("-Clink-arg=-Wl,-install_name," + chr(64) + "rpath/" + shared); ' +
    'cmd.extend(["--"] + rustc); subprocess.check_call(cmd); ' +
    'built = Path(target) / profile; ' +
    'shutil.copy2(built / shared, out_s); shutil.copy2(built / archive, out_a)',
    cargo_prog.full_path(),
    _msproot,
    _mbproot / 'cargo-target',
    _profile,
    cargo_shared,
    cargo_static,
    '@OUTPUT0@',
    '@OUTPUT1@',
  ],
  build_by_default: true,
  install: _install,
  install_dir: [get_option('libdir'), get_option('libdir')],
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
# Do not attach the cargo outputs as sources: a static consumer then
# passes liblinkcell.so to `ar` ("file format not recognized").
linkcell_built_h = custom_target(
  'linkcell-built-h',
  input: linkcell_cargo,
  output: 'linkcell_built.h',
  command: [
    py,
    '-c',
    'from pathlib import Path; Path("@OUTPUT@").write_text("#pragma once\\n")',
  ],
  build_by_default: true,
)
_link_idx = get_option('default_library') == 'static' ? 1 : 0
linkcell_dep = declare_dependency(
  link_args: [linkcell_cargo[_link_idx].full_path()],
  sources: linkcell_built_h,
  include_directories: include_directories('include'),
  dependencies: rust_sys,
  version: meson.project_version(),
)
meson.override_dependency('linkcell', linkcell_dep)

if _install
  install_headers('include/linkcell.h', 'include/linkcell.hpp')
  pkg = import('pkgconfig')
  pkg.generate(
    name: 'linkcell',
    filebase: 'linkcell',
    description: 'Periodic linked-cell k-nearest neighbour search',
    version: meson.project_version(),
    url: 'https://github.com/d-SEAMS/linkcell',
    # Raw `-lfoo` does not get `-L${libdir}`; wrap consumers then fail
    # at link with "cannot find -llinkcell".
    libraries: ['-L${libdir}', '-llinkcell'],
    libraries_private: rust_sys,
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
