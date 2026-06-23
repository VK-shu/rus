# RUS -- Repeat-Until-Success quantum circuit synthesis

## BUILD

### Requirements

- **C++17** compiler (GCC or Clang)
- **CMake** 3.16 or newer
- **OpenMP**
- **Boost** 1.70 or newer
  - `program_options`
  - `chrono`
  - `timer`
  - `system`
- **GMP** (`libgmp`, `libgmpxx`)
- **MPFR** (`libmpfr`)
- **PARI/GP** (`libpari`) — required by the norm solver; often must be built from source if not available from your package manager

Information about program use is available through the `--help` option.

### Install dependencies (Ubuntu / Debian / WSL)

```bash
sudo apt-get update
sudo apt-get install -y \
  build-essential cmake \
  libboost-all-dev \
  libgmp-dev libmpfr-dev
```

Install PARI if `libpari-dev` is available:

```bash
sudo apt-get install -y pari-gp libpari-dev
```

If PARI is not available or linking fails, build from source:

```bash
wget http://pari.math.u-bordeaux.fr/pub/pari/unix/pari-2.15.5.tar.gz
tar -xvf pari-2.15.5.tar.gz
cd pari-2.15.5
./Configure
make
sudo make install
cp misc/gprc.dft ~/.gprc
cd ..
sudo ldconfig
```

### Compile

From the `external/rus` directory, using an out-of-source build (recommended):

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

The binary is produced at `build/rusSyn`.

In-source build is also supported:

```bash
cmake -DCMAKE_BUILD_TYPE=Release .
make -j$(nproc)
```

The binary is produced at `./rusSyn`.

### Verify

```bash
./rusSyn --about
```

## ABOUT 
The program code based on results of https://arxiv.org/abs/1409.3552. 
In addition to Boost, The GNU Multiple Precision Arithmetic Library, The GNU MPFR Library the library 
mpfr::real by Christian Schneider <software(at)chschneider(dot)eu> is used for high precision
In addition, a significant portion of the code leverages from [SQCT](https://github.com/vadym-kl/sqct).

## DIRECTORY STRUCTURE

See **[docs/STRUCTURE.md](docs/STRUCTURE.md)** for the full file layout.

```
apps/     Entry points (rusSyn → main_rus.cpp)
core/     Shared library (rings, matrices, circuit I/O, utilities)
rus/      RUS algorithm (PSLQ → Normalization → synthesis)
appr/     Single-qubit unitary approximation
es/       Exact Clifford+T synthesis
sk/       Solovay–Kitaev (not linked into rusSyn)
theory/   Paper numerical verification (not linked into rusSyn)
ttmath/   Third-party big-integer library
docs/     Project documentation
output/   QASM output from rusSyn
tests/    Mathematica notebooks
```

## USAGE

### Build Single Gate

```bash
./rusSyn -T pi/16 -O out.qasm -F 5 -E 1e-5 -C g-count
```

Output is written to `output/out.qasm` by default (relative `-O` paths go under `output/`).

### Build Gate Database
Add `--database`

For example, 

```
./rusSyn -O output_folder -F 5 -E 1e-5 -C g-count -P pi/32 --qubit-name q[0] --ancil-name a[0] --database
```

will generate Rz gate from $+\pi$ to $-\pi$ per $\pi/32$ (Centered at 0). The output files name will be `${output_folder}/out${i}_pi|32_1e-5`, where $i$ is integer (which means $\theta=i*\pi/32$).

(The program may die for high effort and precision currently)
