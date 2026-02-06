#include <iostream>
#include <fstream>
#include <vector>
using namespace std;
// Function to take a step in the SIR model
// state: vector of S, I, R
// beta: infection rate
// gamma: recovery rate
// dt: time step
std::vector<double> take_step(std::vector<double> state, double beta, double gamma, double dt){
    std::vector<double> new_state;
    //todo: implement the SIR model
    double S = state[0];
    double I = state[1];
    double R = state[2];
    double N = S+I+R;

    //Differential change
    double dS = -beta*I*S/N;
    double dI = beta*I*S/N -gamma*I;
    double dR = gamma*I;


    new_state.assign({S + dS*dt, I + dI*dt,R + dR*dt });
    
    return new_state;
}

void simulate(const string filename,const double dt, const double max_time, double S, double I, double R, const double beta, const double gamma){    
    cout << "Running sim with dt = " <<  dt << "\n";
    std::vector<double> state = {S,I,R}; // Initial condition

    double time = 0.0;
    int t_end = max_time / dt;

    std::ofstream myFile(filename); 
    if (myFile.is_open()){
        myFile << "S I R time \n" ;
        for (double value : state){
            myFile << value << " ";
        }
        myFile << time << "\n";
    }

    for (int i = 0; i < t_end; i++){
        state = take_step(state,beta,gamma,dt);
        double old_time = time;
        time += dt;

        if (int( time +1e9 ) > int (old_time + 1e9) ){
            for (double value : state){
                myFile << value << " ";
            }
            myFile << time << "\n";
        }        
    }
}

//======================================================================================================
//======================== Main function ===============================================================
//======================================================================================================


int main(int argc, char* argv[]){

    double S = 999.0;
    double I = 1.0;
    double R = 0.0;

    double beta = 0.2;
    double gamma = 0.1;
    double max_time = 200.0;

    vector<string> filenames{"data/data1.txt","data/data2.txt","data/data3.txt","data/data4.txt","data/data5.txt"};
    vector<double> stepsizes{1,0.1,0.01, 0.001, 0.0001};
    vector<string> filenames_half{"data/data1_half.txt","data/data2_half.txt","data/data3_half.txt","data/data4_half.txt","data/data5_half.txt"};
    vector<double> stepsizes_half{0.5,0.05,0.005,0.0005,0.00005};

    for (int i = 0; i < filenames.size(); i++){
        simulate(filenames[i],stepsizes[i],max_time,S,I,R,beta,gamma);
        simulate(filenames_half[i],stepsizes_half[i],max_time,S,I,R,beta,gamma);
    }
    
    simulate("data/bigstep.txt",10.0,max_time, 999.0, 1.0, 0.0, 0.2, 0.1);
    simulate("data/data_beta0.txt",0.01, max_time , 999.0, 1.0, 0.0, 0.0, 0.1);
    simulate("data/data_gamma0.txt", 0.01, max_time ,999.0, 1.0, 0.0, 0.2 , 0.0);
    simulate("data/data_I0.txt", 0.01, max_time ,1000.0, 0.0, 0.0, 0.2 , 0.1);

}