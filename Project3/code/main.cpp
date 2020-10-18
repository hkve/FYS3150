#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <functional>
#include <math.h>
#include "time.h"
using namespace std;

double G_ = 6.67408e-11;
double AU = 1.496e+11;
double ME = 5.972e+24; 
double G = G_/pow(AU,3)*pow(3600*24,2)*ME;
double c = 299792458/AU*3600*24;


struct Body {
  int UUID;
  double *pos, *vel, m;  
};




class System{
    public:
    Body *bodies;
    int bodyCount, N, method;
    double ***pos, ***vel, dt, beta;
    bool GR;

    //void (*updateAcceleration)(double**, int, Body*, double);
    System(string initfile){
        readData(initfile);
        
    }
 

    void solve(double _dt = 0.1, int _N = 10000, int _method = 1, double _beta=2, bool _GR = false){
        /*
        Solves the solar system positions for time resolution dt and N steps.
        Method 0 is Forward Euler, Method 1 is Velocity Verlet (default)
        beta is the varying scalar in project 3e),
        if GR is true, the general relativity equation in 3i) will be used
        */

        method = _method;
        dt = _dt;
        N = _N;
        beta = _beta;
        GR = _GR;
        

        pos = new double**[N];
        vel = new double**[N];
        for(int i=0; i < N; i++){
            pos[i] = new double*[bodyCount]; 
            vel[i] = new double*[bodyCount]; 
        }

        double *a[bodyCount]; // array to be filled with the acceleration of each body. Its value is updated each time step
        for(int i = 0; i < bodyCount; i ++){
            a[i] = new double[3];
        }
        storePosVel(0); // store the initial pos/vel of bodies
        updateAcceleration(a); // calculate initial acceleration

        for(int t = 1; t < N; t ++){
            
            if(method == 0){ 
                FwdEulerStep(a, dt);
                updateAcceleration(a); // calculates the acceleration on each body (i.e fills array a)
                 }
            else if(method == 1) { 
                // No need to update acceleration each step here, as it is already dynamically done within the Velocity verlet integrtaion loop
                VelVerStep(a, dt); 
                }
           
            storePosVel(t);  // store positions / velocity updates

        }
        
    }
    

    void write(string filename, int dpts){
        /*
        Writes solved values + simulation specs to the file data/filename for each of the planets in the following way:
        
        UUID,dt,N,method,fpy
        x0,x1,x2,...,xN,
        y0,y1,y2,...,yN,
        z0,z1,z2,...,zN,
        vx0,vx1,vx2,...,vxN,
        vy0,vy1,vy2,...,vyN,
        vz0,vz1,vz2,...,vzN,
        *
        
        the * is a separator for easier parsing in python

        */
       
        int dptsDist = (int)(N/(dpts));
        if (dptsDist == 0){
            dptsDist= 1;
        }
        
        ofstream dataout;
        dataout.open("data/"+filename);
        for(int i = 0; i < bodyCount; i ++){
            dataout << bodies[i].UUID << "," << dt << "," << N << "," << method << "," <<dpts+1<< endl;
            for(int j = 0; j < 3; j++){
                for(int k = 0; k < dpts; k ++){
                    dataout << pos[k*dptsDist][i][j] << ",";
                }
                dataout << pos[N-1][i][j] <<endl;
            }  
            for(int j = 0; j < 3; j++){
                for(int k = 0; k < dpts; k ++){
                    dataout << vel[k*dptsDist][i][j] << ",";
                }
                dataout << pos[N-1][i][j] << endl;
            }  
            dataout << "*" << endl;
        }
        
        dataout.close();
    }

    private:


    void readData(string initfile){
        /* 
        Reads the data from given initfile and fills the class variable bodies with 
        instances of the planetary bodies in the initfile. Also sets the value of bodyCount (no. of bodies in system)
        */


        // first count number of bodies in body file to properly init array lengths
        bodyCount = 0;
        string dummy;
        ifstream lineCounterFile(initfile);

        while (getline(lineCounterFile, dummy))
            bodyCount++;
      
        bodies = new Body[bodyCount];
        
        ifstream data(initfile);
        int k = 0;
        string UUID, x,y,z,vx,vy,vz,m;
        while( data.good()){
            try{
                getline(data, UUID, data.widen(','));
                getline(data, x, data.widen(','));
                getline(data, y, data.widen(','));
                getline(data, z, data.widen(','));
                getline(data, vx, data.widen(','));
                getline(data, vy, data.widen(','));
                getline(data, vz, data.widen(','));
                getline(data, m, data.widen('\n'));
                
                double *pos = new double[3] {stod(x), stod(y), stod(z)};
                double *vel = new double[3] {stod(vx), stod(vy), stod(vz)};
            
                bodies[k].UUID = stoi(UUID);
                bodies[k].pos = pos;
                bodies[k].vel = vel;
                bodies[k].m = stod(m);
                k ++;
                } 
            catch (const std::invalid_argument)     {
                // Avoids parsing trouble if initfile has newline at eof
                bodyCount--;
                break;
            }

        }
        data.close();
        
    }

    void FwdEulerStep(double **a, double dt){
        for(int i = 0; i < bodyCount; i ++){
            for (int j = 0; j < 3; j ++){
                bodies[i].pos[j] = bodies[i].vel[j]*dt + bodies[i].pos[j];
                bodies[i].vel[j] = a[i][j]*dt + bodies[i].vel[j];
                //a[i][j] = 0;
            }
        }

    }
    
    void VelVerStep(double **a, double dt){
        /*
        Uses velocity Verlet to integrate
        */
        for(int i = 0; i < bodyCount; i ++){
            for (int j = 0; j < 3; j ++){
                bodies[i].pos[j] = bodies[i].pos[j] +bodies[i].vel[j]*dt +0.5*a[i][j]*dt*dt;
            }
        }
        // velocity verlet needs the acceleration at the next time step, so we update the acceleration array and store the old one
        double *a_old[bodyCount];
        for(int i = 0; i < bodyCount; i++){
            a_old[i] = new double[3];
            for(int j = 0; j < 3; j++){
                a_old[i][j] = a[i][j];
            }
        }
        updateAcceleration(a);  


        for(int i = 0; i < bodyCount; i ++){
            for (int j = 0; j < 3; j ++){
                bodies[i].vel[j] = bodies[i].vel[j] + (a_old[i][j]+ a[i][j])/2*dt;
            }
        }

        delete[] *a_old; // old one is no longer needed. Makes sure to remove it from heap

    }

    void storePosVel(int t){
        for(int i = 0 ; i < bodyCount; i ++){
                pos[t][i] = new double[3] {bodies[i].pos[0],bodies[i].pos[1], bodies[i].pos[2]};
                vel[t][i] = new double[3] {bodies[i].vel[0],bodies[i].vel[1], bodies[i].vel[2]};   
            }
    }

    double getGRterm(Body b1, Body b2, double dist){
        // calculates the GR term used in 3i)
        double relpos[3] = {b1.pos[0]-b2.pos[0],b1.pos[1]-b2.pos[1], b1.pos[2]-b2.pos[2]};
        double relvel[3] = {b1.vel[0]-b2.vel[0],b1.vel[1]-b2.vel[1], b1.vel[2]-b2.vel[2]};

        double l = (pow(relpos[0],2) + pow(relpos[1],2) + pow(relpos[2],2))
                 * (pow(relvel[0],2) + pow(relvel[1],2) + pow(relvel[2],2))
                 - pow(relpos[0]*relvel[0] + relpos[1]*relvel[1] + relpos[2]*relvel[2], 2);  
        
        return 3*l/(dist *c*c);
    }

    void updateAcceleration(double **a){
        double GMm, dist, GRterm;
   
         for(int i = 0; i < bodyCount; i ++){
            for (int j = 0; j < 3; j ++){
                a[i][j] = 0;
            }
         }
        // Calculate a on bodies
        for(int i = 0; i <bodyCount-1; i++){
            dist = 0;
            
            for(int j = i+1; j<bodyCount; j++){
                for(int k = 0; k < 3; k++){
                    dist += pow(bodies[i].pos[k]- bodies[j].pos[k],2);
                }
                
                GRterm = 0;
                if(GR){
                    // add some terms if the system is to be solved with general relativity approximation
                    GRterm = getGRterm(bodies[i], bodies[j], dist);
                }


                dist = pow(dist,(beta +1)/2);
                GMm = G*bodies[i].m*bodies[j].m;
                for(int k = 0; k < 3; k++){
                    a[i][k] -= GMm/dist*(bodies[i].pos[k]- bodies[j].pos[k])*(1+GRterm);
                    a[j][k] -= a[i][k]; // lookup is faster than doing the math
                }
            }
        }
        // Convert a into accelerations (still under name a, though)
        for(int i = 0; i < bodyCount; i++){
            for(int j = 0; j < 3; j++){
                a[i][j] /= bodies[i].m;
            }
        }

    }

    
};


int main(int argc, char** argv){
    /*
    Arguments:
        4 required arguments in the following order:

        systemInit (string): name of the file where initial planet data is read from
        systemOut (string): name of the file (in the data/ dir) where the data will be stored
        dpts (int): Numbers of datapoints per to be stored (default is N)
        dt (float): time step size (in days) of integration loop
        N (int): number of integration steps
        method (int): 0 or 1. Which method of integration to use
            - 0: Forward Euler
            - 1: Velocity Verlet
        beta (float): number between 2 and 3. Used for 3e). beta = 2 results in normal Newtonian gravity
        GR (bool/int): 0 or 1. Decides if the general relativity correction term should be included
        timeFormat (string): "days" or "years". Which time format to use. "days" means all times are in days, while "years" means all times are in years
        q (bool/int): 0 or 1. Run program quietly
    */
    clock_t start = clock();

    string initfile = (string)argv[1];
    string outfile = (string)argv[2];
    int dpts = atoi(argv[3]);
    double dt = atof(argv[4]);
    int N = atoi(argv[5]);
    int method = atoi(argv[6]);
    double beta = atof(argv[7]);
    bool GR = (bool)atoi(argv[8]);
    string timeFormat = (string)argv[9];
    bool q = (bool)atoi(argv[10]);

    
     
    
    if (timeFormat == "years") {
        G = G_/pow(AU,3)*pow(3600*24*365,2)*ME;
        c = 299792458/AU*3600*24*365;
    }
    else if( timeFormat != "days"){
        cout << "Error, unrecognized time format '" << timeFormat << "'. Proceeding with day-format" << endl; 
    }

    System simple(initfile);
    simple.solve(dt, N, method, beta, GR);
    clock_t stop= clock();
    if(not q){
    cout << "Solving done in " << ((stop - start) / (double)CLOCKS_PER_SEC) << "s " << endl;
    }
    start = clock();
    if(not q){
    cout << "Writing..." << endl;
    }
    simple.write(outfile, dpts);
    stop = clock();
    if(not q){
    cout << "Writing done in " << ((stop - start) / (double)CLOCKS_PER_SEC) << "s " << endl;
    }
    return 0;
}