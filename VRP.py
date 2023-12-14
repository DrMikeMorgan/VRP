from ortools.constraint_solver import routing_enums_pb2
from ortools.constraint_solver import pywrapcp
import sys
import networkx as nx
import matplotlib.pyplot as plt
import random

class VRP:
	def __init__(self,n,v=1,max_dist=(2**31),scale=1000,pos=[], capacities=[],demands=[],windows=[],speed=1000,manhattan=0.0):
		"""Constructor, builds a VRP problem
		:param:	n:	The number of nodes 
		:param:	v:	The number of vehicles, default 1 (TSP)
		:param:	max_dist:	The maximum distance a single vehicle may travel, default 2^31
		:param:	scale:	The size of the square on which nodes will be uniformly distributed. 
						Used if positions have not been provided to the constructor, default 1000
		:param:	pos:	A list of tuples with (x,y) coordinates for the locations (optional). 
						Default randomises x and y between -scale/2 and scale/2 with the depot at (0,0)
		:param:	capacities:	The capacities for each vehicle, must also supply a vector of demands if this is used
		:param:	demands:	The demands for each location
		:param:	speed:	The average speed of a vehicle in units of scale per hour. 
						Default 1000 (the width/height of the default square), only relevant if time windows are supplied
		:param:	windows:	A list of tuples providing the time window for each location, expressed in hours from time 0
		:param: manhattan:	A factor which decides how distances will be computed. 
							0 gives euclidean distance (default), 1 gives manhattan distance, 
							any other value gives a weighted average between the two. Suggested values between 0 and 1 
		"""
		self.colours = [(random.random()*0.8,random.random()*0.8,random.random()*0.8,0.5) for x in range(v)]
		self.n=n
		if len(pos)>0:
			self.pos=pos
		else:
			self.pos = {i: (random.random()*scale-scale/2, random.random()*scale-scale/2) for i in range(n)}
			self.pos[0]=(0,0)
		self.G = nx.random_geometric_graph(n,0,pos=self.pos)
		euclid_distance = [ [ ((self.pos[i][0]-self.pos[j][0])**2 + (self.pos[i][1]-self.pos[j][1])**2) ** 0.5 for j in range(n)] for i in range(n) ]
		manhattan_distance = [ [ abs(self.pos[i][0]-self.pos[j][0]) + abs(self.pos[i][1]-self.pos[j][1]) for j in range(n)] for i in range(n) ]
		self.d = [ [ manhattan*manhattan_distance[i][j] + (1.0-manhattan)*euclid_distance[i][j] for j in range(n)] for i in range(n) ]
		self.manager = pywrapcp.RoutingIndexManager(n,v,0)
		self.routing = pywrapcp.RoutingModel(self.manager)
		def distance_callback(from_index, to_index):
			from_node = self.manager.IndexToNode(from_index)
			to_node = self.manager.IndexToNode(to_index)
			return self.d[from_node][to_node]
		self.cost = self.routing.RegisterTransitCallback(distance_callback)
		self.routing.SetArcCostEvaluatorOfAllVehicles(self.cost)
		self.capacities = capacities
		self.demands = demands
		self.windows = [(int(x*speed),int(y*speed)) for (x,y) in windows]
		self.vehicles = v
		if v > 1:
			if windows == []:
				dimension_name = "Distance"
				self.routing.AddDimension(self.cost,0,max_dist,True,dimension_name)
			else:
				dimension_name = "Time"
				self.routing.AddDimension(self.cost,scale*scale,max_dist,False,dimension_name)
				time_dimension = self.routing.GetDimensionOrDie(dimension_name)
				for location_idx, time_window in enumerate(self.windows):
					if location_idx == 0:
						continue
					index = self.manager.NodeToIndex(location_idx)
					time_dimension.CumulVar(index).SetRange(time_window[0], time_window[1])
				for i in range(self.vehicles):
					self.routing.AddVariableMinimizedByFinalizer(time_dimension.CumulVar(self.routing.Start(i)))
					self.routing.AddVariableMinimizedByFinalizer(time_dimension.CumulVar(self.routing.End(i)))
		
		if len(capacities) > 0:
			def demand_callback(from_index):
				from_node = self.manager.IndexToNode(from_index)
				return self.demands[from_node]
			demand_callback_index = self.routing.RegisterUnaryTransitCallback(demand_callback)
			self.routing.AddDimensionWithVehicleCapacity(demand_callback_index,0, self.capacities, True, "Capacity")

    
	def render(self, solution):
		"""creates internal graph structure for solution rendering, call show() to display on screen
		:param:	solution:	A solution object for display
		"""
		self.G.remove_edges_from(self.G.edges())
		for vehicle_id in range(self.vehicles):
			index = self.routing.Start(vehicle_id)
			while not self.routing.IsEnd(index):
				previous_index = index
				index = solution.Value(self.routing.NextVar(index))
				self.G.add_edge(self.manager.IndexToNode(previous_index), self.manager.IndexToNode(index),weight=self.d[self.manager.IndexToNode(previous_index)][self.manager.IndexToNode(index)],color=self.colours[vehicle_id%len(self.colours)])
	
	def show(self,name="Results"):
		"""shows edges assigned to solution passed to the most recent render() call"""
		nodeSizes = [5]*self.n
		nodeSizes[0] = 15
		nx.draw(self.G,self.pos,edge_color=nx.get_edge_attributes(self.G,'color').values(),node_size=nodeSizes)
		plt.get_current_fig_manager().set_window_title(name)
		plt.show()
		
	def solve(self):
		"""solves using OR Tools' default engine"""
		return self.routing.Solve()

	def guidedLocalSearch(self, t=0):
		"""Solves using Guided Local Search
		:param:	t:	the maximum time allowed for the process
		:return:	the solution object or None if no solution was found
		"""
		search_parameters = pywrapcp.DefaultRoutingSearchParameters()
		search_parameters.local_search_metaheuristic = (routing_enums_pb2.LocalSearchMetaheuristic.GUIDED_LOCAL_SEARCH)
		search_parameters.time_limit.seconds = t
		return self.routing.SolveWithParameters(search_parameters)
	
	def tabuSearch(self, t=0):
		"""Solves using Tabu Search
		:param:	t:	the maximum time allowed for the process
		:return:	the solution object or None if no solution was found
		"""
		search_parameters = pywrapcp.DefaultRoutingSearchParameters()
		search_parameters.local_search_metaheuristic = (routing_enums_pb2.LocalSearchMetaheuristic.TABU_SEARCH)
		search_parameters.time_limit.seconds = t
		return self.routing.SolveWithParameters(search_parameters)
	
	def simulatedAnnealing(self, t=0):
		"""Solves using Simulated Annealing
		:param:	t:	the maximum time allowed for the process
		:return:	the solution object or None if no solution was found
		"""
		search_parameters = pywrapcp.DefaultRoutingSearchParameters()
		search_parameters.local_search_metaheuristic = (routing_enums_pb2.LocalSearchMetaheuristic.SIMULATED_ANNEALING)
		search_parameters.time_limit.seconds = t
		return self.routing.SolveWithParameters(search_parameters)
	
	def distanceSaving(self, sol):
		"""Returns the distance saved compared to individual journeys from each location to the depot, as a percentage
		:param:	sol:	The solution to be compared
		:return:		The percentage reduction in cost, e.g. 20% reduction means the solution covers 80% of the individual journey distance 
		"""
		return 100 * (sum(self.d[0])*2 - sol.ObjectiveValue()) / (sum(self.d[0])*2) 
