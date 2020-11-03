Create a directory for the simulation:

```
mkdir Simulation
cd Simulation
```

To run simulator for a certain coverage level, create a directory for that experiment:

```
mkdir 20x
cd 20x
```

Run the simulator:

```
simulate_random.sh 20
```

This will simulate a trio, error-correct each member with ntEdit, create and FMD index using ropebwt and finally get the child-specific strings and analysze them. ntHits and ntEdit should be in path. The `scripts` directory for simulation should also be in path.
