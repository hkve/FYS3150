#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <functional>
#include <math.h>
using namespace std;

double G_ = 6.67408e-11;
double AU = 1.496e+11;
double ME = 5.972e+24; 
double G = G_/pow(AU,3)*pow(3600*24,2)*ME;
struct Body {
  int UUID;
  double *pos, *pos0, *vel, *vel0, m;  
};


class System{
    public:
    Body *bodies;
    string bodyfile;
    int bodyCount, N;
    double ***pos, ***vel;
    System(string _bodyfile){
        bodyfile = _bodyfile;
        readData();
        solve();
        write();
    }
 

    void readData(){
        // first count number of bodies in body file to properly init array lengths
        bodyCount = 0;

        ifstream datacount;
        datacount.open(bodyfile);
        while(datacount.good()){
            string foo;
            getline(datacount, foo);
            bodyCount ++;
        }
        datacount.close();

        // then init bodies array
        bodies = new Body[bodyCount];
        ifstream data;
        data.open(bodyfile);
        int k = 0;
        while( data.good()){
            string UUID, x,y,z,vx,vy,vz,m;
            
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
            double *pos0 = new double[3] {stod(x), stod(y), stod(z)};
            double *vel0 = new double[3] {stod(vx), stod(vy), stod(vz)};
           
            bodies[k].UUID = stoi(UUID);
            bodies[k].pos = pos;
            bodies[k].pos0 = pos0; // might not work
            bodies[k].vel = vel;
            bodies[k].vel0 = vel0; // same as above
            bodies[k].m = stod(m);
            k ++;
        }
        data.close();
    }

    void FwdEulerStep(double *forces[3], double dt){
        for(int i = 0; i < bodyCount; i ++){
            for (int j = 0; j < 3; j ++){
                bodies[i].vel[j] = forces[i][j]*dt + bodies[i].vel[j];
                bodies[i].pos[j] = bodies[i].vel[j]*dt + bodies[i].pos[j];
                forces[i][j] = 0;
            }
        }

    };
    void VelVerStep(double *forces[3]){};


    void solve(double _dt = 0.1, int _N = 10000, int _method= 0){
        double dt = _dt;
        N = _N;

        pos = new double**[N];
        vel = new double**[N];
        for(int i=0; i < N; i++){
            pos[i] = new double*[bodyCount]; 
            vel[i] = new double*[bodyCount]; 
        }
   

        double *forces[bodyCount]; 
        for(int i = 0; i < bodyCount; i ++){
            forces[i] = new double[3];
        }
        double C, d; //force const. , distance


        // void(System::*step)(double*[3]);
        // if (_method == 0){
        //     step = &System::FwdEulerStep;
        // }
        // else{
        //     step = &System::VelVerStep;
        // }

        for(int t = 0; t < N; t ++){
            // Calculate forces on bodies
            for(int i = 0; i <bodyCount-1; i++){
                d = 0;
                
                for(int j = i+1; j<bodyCount; j++){
                    for(int k = 0; k < 3; k++){
                        d += pow(bodies[i].pos[k]- bodies[j].pos[k],2);
                    }
                    d = pow(d,3/2);
                    C = G*bodies[i].m*bodies[j].m;
                    for(int k = 0; k < 3; k++){
                        forces[i][k] -= C/d*(bodies[i].pos[k]- bodies[j].pos[k]);
                        forces[j][k] -= forces[i][k];
                    }
                }
            }
            // Convert forces into accelerations (still under name forces, though)
            for(int i = 0; i < bodyCount; i++){
                for(int j = 0; j < 3; j++){
                    forces[i][j] /= bodies[i].m;
                }
            }
            FwdEulerStep(forces, dt=dt);
            
            // store positions / velocity updates
            for(int i = 0 ; i < bodyCount; i ++){
                pos[t][i] = new double[3] {bodies[i].pos[0],bodies[i].pos[1], bodies[i].pos[2]};
                vel[t][i] = new double[3] {bodies[i].vel[0],bodies[i].vel[1], bodies[i].vel[2]};
                
            }
            

        }
    }

    void write(){
        ofstream dataout;
        dataout.open("data/TEMPNAME.dat");
        for(int i =0; i < bodyCount; i ++){
            dataout << bodies[i].UUID << ","<< bodies[i].pos0[0] << ","<< bodies[i].pos0[1] << ","<< bodies[i].pos0[2] <<
             ","<< bodies[i].vel0[0] <<  ","<< bodies[i].vel0[1] <<  ","<< bodies[i].vel0[2] << ","<< bodies[i].m << endl;

            for(int j = 0; j < 3; j++){
                for(int k = 0; k < N; k ++){
                    dataout << pos[k][i][j] << ",";
                }
                dataout << endl;
            }  
            for(int j = 0; j < 3; j++){
                for(int k = 0; k < N; k ++){
                    dataout << vel[k][i][j] << ",";
                }
                dataout << endl;
            }  
            dataout << "|" << endl;
        }
        
        dataout.close();
    }
};


int main(int argc, char** argv){
    /*
    Arguments:
        
    */
    System A("sys1_mod.txt") ;
    return 0;
}