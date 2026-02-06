#include <iostream>
#include <fstream>
#include <vector>
#include <pybind11/stl.h>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

// Function to take a step in the SIR model
// state: vector of S, I, R
// beta: infection rate
// gamma: recovery rate
// dt: time step
std::vector<double> take_step(std::vector<double> state, double beta, double gamma, double dt, double N, double t){
    std::vector<double> new_state(4);
    //todo: implement the SIR model
    new_state[0] = state[0]-(beta*state[1]*state[0]/N)*dt;
    new_state[1] = state[1]+(beta*state[1]*state[0]/N-gamma*state[1])*dt;
    new_state[2] = state[2]+(gamma*state[1])*dt;
    new_state[3] = t; // time
    return new_state;
}
// Function simulating num_steps of the SIR model, saving the state every return_every steps and returning the results
// S0: initial number of susceptible individuals
// I0: initial number of infected individuals
// R0: initial number of recovered individuals
// beta: infection rate
// gamma: recovery rate
// dt: time step
// num_steps: number of steps to simulate
// return_every: save the state every return_every steps
pybind11::array integrate_system(double S0, double I0, double R0, double beta, double gamma, double dt, int tot_time){
    // TODO: implement the SIR model
    std::vector<double> state = {S0, I0, R0, 0};
    std::vector<std::vector<double>> results;
    std::cout << "Running simulation with dt = " << dt << std::endl;
    // TODO: implement the SIR model
    double N = S0 + I0 + R0;
    for (double t = dt; t <= tot_time; t += dt) {

        std::vector<double> new_state = take_step(state, beta, gamma, dt, N, t);
        results.push_back(new_state);
        state = new_state;
    }
    return pybind11::cast(results);
}

PYBIND11_MODULE(SIR_python, m) {
    m.doc() = "This is a Python binding for the SIR model";

    m.def("integrate_system", &integrate_system);
}
