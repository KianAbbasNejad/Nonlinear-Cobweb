# Nonlinear Dynamics in a Cobweb Model
This script calculates the following for a simple non-linear cobweb model: 
1. Iteration Loop
2. Time Series 
3. Bifurcation Diagram
4. Parameter Basins of Attraction - creates a meshgrid of parameter space, given an initial condition creates a time series of that initial condition on each point in the meshgrid, and calculates its period given a tolerance level. 

The code is inspired by Hommes (2013) and his lectures at University of
Amsterdam as well as the UvA CeNDEF software '*E&FChaos*'.

The code doesn't have GUI, but I've tried my best to make it easy to change the equations/parameters and play around with the graphs. :) 

__Sample Figures__

![basin_attraction](https://user-images.githubusercontent.com/45733935/81515600-287c3580-9335-11ea-8463-4ddc97b1398b.png)
![Bifurcation](https://user-images.githubusercontent.com/45733935/79876587-87602600-83eb-11ea-9672-bc0d607f4631.png)
![Time_Series](https://user-images.githubusercontent.com/45733935/79876593-88915300-83eb-11ea-9111-a7ebe9be03d6.png)
![Iteration](https://user-images.githubusercontent.com/45733935/79876595-8929e980-83eb-11ea-9de7-9b44b8166706.png)
