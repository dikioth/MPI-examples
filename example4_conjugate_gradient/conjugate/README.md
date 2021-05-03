# Conjugate gradient method with stencil-based matrix-vector multiplication


- the binary program is `conjugate` and will be compiled in the `bin` directory. Usage is: `./bin/conjugate n` where n is the number of intervall splits along each axis. 
- To make the binary file, run: `make`
- To test the algorithm (with n=256): run `make test`
- To perform the strong and weak scale tests on the cluster computer, run the bash scripts located in the `test` directory. 
