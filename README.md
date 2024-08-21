# cocoa

CoCoA - Collaborative Coupling Applications

## how to install

* copy `spack-env/spack-template.yaml` in `spack/spack.yaml`, edit packages and install them

```bash
cp spack-env/spack-template.yaml spack-env/spack.yaml
edit spack-env/spack.yaml
spack env activate spack-env
spack concretize -f
spack install
```

* setup OpenFOAM (optional)

```bash
source <FOAM_DIR>/etc/bashrc
```

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
ln -s ../data/fd1d1.dat
ln -s ../data/fd1d2.dat
```

* run app

```bash
cd build
./main
```

* plot results of fd1d

```bash
./fd1dplot.py build/fd1d1.dat
```
