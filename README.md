# Compilation instructions

Have SDPA downloaded and compiled as described  in https://github.com/beautifulv0id/SDPA. 

Then:

```bash
cd GRET-SDP-DEMO
mkdir build && cd build
cmake .. -DSDPA_ROOT_DIR=path/to/SDPA -DMOSEK_ROOT_DIR=path/to/MOSEK
make
```
