# Status of Jaguar in Corcovado National Park, Costa Rica

This repo contains code to run capture-recapture model to determine the status of jaguars in Corcovado National Park, Costa Rica. 

We implemented a (spatially-explicit) Jolly-Seber model with habitat selection in Stan. (See the `stan` folder). Since Stan does not allow discrete parameters, the state of the individual was manually marginalized during posterior sampling via a forward algorithm and the sample back using forward filtering backward sampling (FFBS) algorithm. One can use `make` to run part of all of the work flow

![density map](https://raw.githubusercontent.com/YunyiShen/Costa-Rica-Jaguar/master/res/Figs/js_null_den_est.png)
![population estimate](https://raw.githubusercontent.com/YunyiShen/Costa-Rica-Jaguar/master/res/Figs/js_null_pop_est.png)
