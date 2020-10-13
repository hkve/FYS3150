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
    System(string _bodyfile, double dt, int N, int method){
        bodyfile = _bodyfile;
        readData();
        
        
        solve(dt, N, method);
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
            catch (const std::invalid_argument)     {
                bodyCount--;
                break;
            }

        }
        

        data.close();
        
    }

    void FwdEulerStep(double *a[3], double dt){
        for(int i = 0; i < bodyCount; i ++){
            for (int j = 0; j < 3; j ++){
                bodies[i].vel[j] = a[i][j]*dt + bodies[i].vel[j];
                bodies[i].pos[j] = bodies[i].vel[j]*dt + bodies[i].pos[j];
                //a[i][j] = 0;
            }
        }

    }
    void VelVerStep(double **a, double dt){
        
        for(int i = 0; i < bodyCount; i ++){
            for (int j = 0; j < 3; j ++){
                bodies[i].pos[j] = bodies[i].pos[j] +bodies[i].vel[j]*dt +0.5*a[i][j]*dt*dt;
            }
        }
        double *a_next[bodyCount]; 
        for(int i = 0; i < bodyCount; i ++){a_next[i] = new double[3];}

        updateAcceleration(a_next);  
        for(int i = 0; i < bodyCount; i ++){
            for (int j = 0; j < 3; j ++){
                bodies[i].vel[j] = bodies[i].vel[j] + (a[i][j]+ a_next[i][j])/2*dt;
            }
        }
        for(int i = 0; i < bodyCount; i++){
            for(int j = 0; j < 3; j++){
                a[i][j] = a_next[i][j];
            }
        }
      

    // stuff
        delete[] *a_next;
        
       


    }

    void updateAcceleration(double **a){
        double C,d;
         for(int i = 0; i < bodyCount; i ++){
            for (int j = 0; j < 3; j ++){
                a[i][j] = 0;
            }
         }
        // Calculate a on bodies
        for(int i = 0; i <bodyCount-1; i++){
            d = 0;
            
            for(int j = i+1; j<bodyCount; j++){
                for(int k = 0; k < 3; k++){
                    d += pow(bodies[i].pos[k]- bodies[j].pos[k],2);
                }
                d = pow(d,3/2);
                C = G*bodies[i].m*bodies[j].m;
                for(int k = 0; k < 3; k++){
                    a[i][k] -= C/d*(bodies[i].pos[k]- bodies[j].pos[k]);
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

    void solve(double _dt = 0.1, int _N = 10000, int _method = 1){
        int method = _method;
        double dt = _dt;
        N = _N;
        
        pos = new double**[N];
        vel = new double**[N];
        for(int i=0; i < N; i++){
            pos[i] = new double*[bodyCount]; 
            vel[i] = new double*[bodyCount]; 
        }
        double *a[bodyCount]; 
        for(int i = 0; i < bodyCount; i ++){a[i] = new double[3];}
   
        for(int t = 0; t < N; t ++){
            updateAcceleration(a);
            if(method == 0){
                FwdEulerStep(a, dt);
            }
            else if(method == 1) {

                VelVerStep(a, dt);
            }
           
          
            
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
    double dt = atof(argv[1]);
    int N = atoi(argv[2]);
    int method = atoi(argv[3]);
    
    System A("sys1_mod.txt", dt, N, method) ;
    return 0;
}