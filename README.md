# cocoa

CoCoA - Collaborative Coupling Applications

## how to install

### spack installation

Currently, `ProXPDE` and `FEMuS` are not yet available in mainstream `spack`, clone our fork

```
git clone https://github.com/capitalaslash/spack.git -b 0.22_femus
```

### setup `spack` environment

* copy `spack-env/spack-template.yaml` in `spack/spack.yaml`, edit packages and install them

```bash
cp spack-env/spack-template.yaml spack-env/spack.yaml
edit spack-env/spack.yaml
spack env activate spack-env
spack concretize -f
spack install
```

### setup `OpenFOAM` (optional)

```bash
source <FOAM_DIR>/etc/bashrc
```

### setup `cocoa`

* copy `run_cmake-template.bash` in `run_cmake.bash` edit options and run it

```bash
cp run_cmake-template.bash run_cmake.bash
edit run_cmake.bash
bash ./run_cmake.bash
```

* compile library and app

```bash
cd build
ninja
```

## how to run

* set up environment

```bash
spack env activate spack-env
```

* copy/link configuration files in build directory

```bash
cd build
ln -s ../data/fd1d_heat.dat
ln -s ../data/fd1d_hc.dat
```

* run app

```bash
cd build
./main
```

* plot results of fd1d

```bash
./fd1dplot.py build/fd1d_heat.dat
```
