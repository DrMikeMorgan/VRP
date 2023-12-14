# VRP

A (very) simple python facade for the vehicle routing problem solver in OR Tools

# Installation

Requires dependencies: ortools networkx matplotlib

```bash
pip install ortools networkx matplotlib
```
or
```bash
python -m pip install ortools networkx matplotlib
```
# Usage

Import as normal and check documentation for full details, example below

```Python
from VRP import VRP
help(VRP)              #for documentation
tsp = VRP(100)         #one vehicle (TSP) with 100 nodes randomly positioned in a square
vrp = VRP(100,5)       #100 nodes, max 5 vehicles
cvrp = VRP(100,5,capacities=[50]*5, demands=[2]*100)  #capacitated VRP with homogeneous demands and capacities
#Further options for e.g. time windows available, see docs

#run guided local search on the problem for 60 seconds, also options for tabu search and simulated annealing
solution = cvrp.guidedLocalSearch(60)

#display solution
print(solution.ObjectiveValue())
cvrp.render(solution)
cvrp.show()
``` 
