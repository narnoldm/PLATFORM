Example2                 {#example2}
=============



This example case uses the Tecio metadata extension to read in tecplot data set files and generate the basis modes using the parallelized SVD.

```bash 
cd examples/example2
mpirun -np <num processes> ../../build/bin/POD POD_tec.inp
```
This is the out of the box tecplot POD example. This executes the method of snapshots on the same dataset files used for the test suite and can be visualized using Tecplot. 

Depending on your installation you may need to run the case in intermediate steps. To achieve this open the input file POD_tec.inp

MOSstep should be set to 0 by default. However you can run the various intermediate steps independently if needed by setting MOSstep to 1,2, and 3 on successive runs. 




This example are based on the umich deepblue data set and documentation are availible [here](https://deepblue.lib.umich.edu/data/concern/data_sets/6w924c14h?locale=en)
Both example cases just use a small subset of this dataset to show the functionality of the code. 