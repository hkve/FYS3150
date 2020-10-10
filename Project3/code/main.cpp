#include <iostream>
#include <string>
#include <fstream>
#include <vector>
using namespace std;

class Body{
    public:
    int UUID;
    double *pos, *vel, m;
    Body(int _UUID, double* _pos, double* _vel, double _m){
        UUID = _UUID;
        pos = _pos;
        vel = _vel;
        m = _m;
    }
    Body(){};

};

class System{
    public:
    Body *bodies;
    int bodyCount;
    System(string bodyfile){
        // read from file bodyfile
       
        readData(bodyfile);

    }
 

    void readData(string bodyfile){
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
        int i = 0;
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
            bodies[i] = *(new Body(stoi(UUID), pos, vel, stod(m)));
            cout << stod(m) << endl;
            i ++;
        }
        data.close();
    }


    void solve(double _dt = 0.1, int _N = 1000, int _method= 0){
        double dt = _dt;
        int N = _N;
        int method = _method;

    }


};


int main(int argc, char** argv){
    /*
    Arguments:
        
    */
    System A("sys1.txt") ;
    cout << A.bodies[0].m << " " << A.bodies[1].m << endl;
    return 0;
}