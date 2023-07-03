#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
#include <string>

const double t = 1.0/10000.0;


class Vector {
    
    public:

        double x;
        double y;
        double z;

    Vector(double input_x, double input_y, double input_z) {
        x = input_x;
        y = input_y;
        z = input_z;
    }

    double GetX() {
        return x;
    }

    double GetY() {
        return y;
    }

    double GetZ() {
        return z;
    }

    void SetX(double input_x) {
        x = input_x;
    }

    void SetY(double input_y) {
        y = input_y;
    }

    void SetZ(double input_z) {
        z = input_z;
    }

    void AddVector(Vector vec){
        x += vec.GetX();
        y += vec.GetY();
        z += vec.GetZ();
    }

    Vector AddVectorReturn (Vector vec){
        Vector result(x + vec.GetX(), y + vec.GetY(), z + vec.GetZ());
        return result;
    }

    void MultiOnConst (double k) {
        x *= k;
        y *= k;
        z *= k;
    }

    Vector MultiOnConstReturn (double k) {
        Vector result(x*k, y*k, z*k);
        return result;
    }
};

class Atom {


    private:
    //coordinats
        Vector coordinates;

    //speed
        Vector speed;

    //force
        Vector force;

    public:

        Atom(double x = 0, double y = 0, double z = 0, double v_x = 0, double v_y = 0, double v_z = 0, double f_x = 0, double f_y = 0, double f_z = 0): coordinates(0,0,0), speed(0,0,0), force(0,0,0) {
            
            coordinates.SetX(x);
            coordinates.SetY(y);
            coordinates.SetZ(z);
            
            speed.SetX(v_x);
            speed.SetY(v_y);
            speed.SetZ(v_z);

            force.SetX(f_x);
            force.SetY(f_y);
            force.SetZ(f_z);

        }

        double GetX() {
            return coordinates.GetX();
        }

        double GetY() {
            return coordinates.GetY();
        }

        double GetZ() {
            return coordinates.GetZ();
        }

        void SetX(double x) {
            coordinates.SetX(x);
        }

        void SetY(double y) {
            coordinates.SetY(y);
        }

        void SetZ(double z) {
            coordinates.SetZ(z);
        }

        double GetXSpeed() {
            return speed.GetX();
        }

        double GetYSpeed() {
            return speed.GetY();
        }

        double GetZSpeed() {
            return speed.GetZ();
        }

        void SetXSpedd(double v_x) {
            speed.SetX(v_x);
        }

        void SetYSpeed(double v_y) {
            speed.SetY(v_y);
        }

        void SetZSpeed(double v_z) {
            speed.SetZ(v_z);
        }

        double GetXForce() {
            return force.GetX();
        }

        double GetYForce() {
            return force.GetY();
        }

        double GetZForce() {
            return force.GetZ();
        }

        void null_force() {
            force.SetX(0);
            force.SetY(0);
            force.SetZ(0);
        }

        void calculate_force(Atom atom) {

            double delta_x = atom.GetX() - coordinates.GetX();
            double delta_y = atom.GetY() - coordinates.GetY();
            double delta_z = atom.GetZ() - coordinates.GetZ();

            /*
            double r = std::sqrt(delta_x*delta_x + delta_y*delta_y + delta_z*delta_z);
            double norm_r_x = delta_x / r;
            double norm_r_y = delta_y / r;
            double norm_r_z = delta_z / r;

            double force_to_atom = -24*( std::pow((1/r),7) - 2*std::pow((1/r), 13));
            

            // можно не извлекать корень а работать с квадратом
            

            Vector force_to_atom_vec(force_to_atom*norm_r_x,
                            force_to_atom*norm_r_y,
                            force_to_atom*norm_r_z);

            */

            double r_squard = delta_x*delta_x + delta_y*delta_y + delta_z*delta_z;

            double force_to_atom = (1/r_squard);

            Vector force_to_atom_vec(force_to_atom*delta_x,
                            force_to_atom*delta_y,
                            force_to_atom*delta_z);

            force.AddVector(force_to_atom_vec);
        }

        void move() {
            speed.AddVector(force.MultiOnConstReturn(t));
            coordinates.AddVector(speed.MultiOnConstReturn(t));
        }

        double p2pPotEnerge (Atom atom) {

            double r_squard = (atom.GetX() - coordinates.GetX())*(atom.GetX() - coordinates.GetX())
                              + (atom.GetY() - coordinates.GetY())*(atom.GetY() - coordinates.GetY())
                              + (atom.GetZ() - coordinates.GetZ())*(atom.GetZ() - coordinates.GetZ());

            return -std::pow((1/r_squard), 0.5);

        }

        double GetKinEnerge() {
            return (speed.x*speed.x + speed.y*speed.y + speed.z*speed.z)/2;
        }
};

int StdOutputWrite(Atom atoms_arr[], int size) {

    for (int i = 0; i < size; i++) {
        std::cout 
            << atoms_arr[i].GetX() << ' ' 
            << atoms_arr[i].GetY() << ' '
            << atoms_arr[i].GetZ() << ' '
            << atoms_arr[i].GetXSpeed() << ' '
            << atoms_arr[i].GetYSpeed() << ' '
            << atoms_arr[i].GetZSpeed() << ' '
            << atoms_arr[i].GetXForce() << ' '
            << atoms_arr[i].GetYForce() << ' '
            << atoms_arr[i].GetZForce() << ' '
            << '\n';
    }

    return 0;

}

int CreateDrop (Atom atoms_arr[], double drops_radius, Vector drops_center, double lattice_step, double atoms_max_speed_proection, Vector drops_speed) {
    
    int n = 2*drops_radius/lattice_step;

    unsigned couter = 0;

    unsigned seed = 1000;
    std::default_random_engine rnd(seed);
    std::uniform_real_distribution<double> dstr(-atoms_max_speed_proection, atoms_max_speed_proection);


    for (int i = -n/2; i < n/2; ++i){
        for (int j = -n/2; j < n/2; ++j){
            for (int k = -n/2; k < n/2; ++k){
                if ((std::pow((i*lattice_step), 2) + std::pow((j*lattice_step), 2) + std::pow((k*lattice_step), 2)) < std::pow(drops_radius, 2)) {
                    atoms_arr[couter] = Atom(i*lattice_step, j*lattice_step, k*lattice_step, dstr(rnd),dstr(rnd),dstr(rnd));
                    ++couter;
                }
            }
        }
    }

    return ++couter;
}

int CreateCube (Atom atoms_arr[], double lattice_step, int size, double atoms_max_speed_proection) {

    unsigned couter = 0;

    unsigned seed = 1000;
    std::default_random_engine rnd(seed);
    std::uniform_real_distribution<double> dstr(-atoms_max_speed_proection, atoms_max_speed_proection);

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            for (int k = 0; k < size; k++) {
                atoms_arr[couter] = Atom(i*lattice_step, j*lattice_step, k*lattice_step, dstr(rnd),dstr(rnd),dstr(rnd));
                ++couter;
            }
        }
    }

    return couter;

}

int ForseCulc (Atom atoms_arr[], int size) {

    for (int i = 0; i < size; i++) {
        atoms_arr[i].null_force();
    }
    
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (i != j) atoms_arr[i].calculate_force(atoms_arr[j]);
        }
    }

    return 0;
}

int TimeStep (Atom atoms_arr[], int size) {

    ForseCulc(atoms_arr, size);

    for (int i = 0; i < size; i++) {
        atoms_arr[i].move();
    }

    return 0;
}

int WriteCoorToFile (Atom atoms_arr[], int size, unsigned step) {

    std::string file_name = "./data/atoms_step_" + std::to_string(step) + ".xyz";

    std::ofstream fout(file_name);

    fout << size << '\n' << "some atoms" << '\n';

    for (int i = 0; i < size; i++) {
        fout 
            << 'H' << ' '
            << std::to_string(atoms_arr[i].GetX()) << ' ' 
            << std::to_string(atoms_arr[i].GetY()) << ' '
            << std::to_string(atoms_arr[i].GetZ()) << ' '
            << '\n'; 
    }

    fout.close();

    return 0;
}

int EnegeCulc (double stepsEnerge[] ,Atom atoms_arr[], int size, unsigned step) {
    double energe = 0;
    for (int i = 0; i < size; i++) {
        stepsEnerge[step] += atoms_arr[i].GetKinEnerge();

        for (int j = i + 1; j < size; j++) {
            stepsEnerge[step] += atoms_arr[i].p2pPotEnerge(atoms_arr[j]);
        }
        
    }

    return 0;
}

double WriteEnergyToFile (double stepsEnerge[], unsigned steps) {

    std::string file_name = "./data/energe.txt";

    std::ofstream fout(file_name);

    for (int i = 0; i < steps - 1; i++) {

            fout << stepsEnerge[i] << '\n';
    }

    fout << stepsEnerge[steps - 1];

    fout.close();

    return 0;
}


int main () {

    Atom array[1024];

    unsigned size;

    int cube_size = 2;

    double init_lattice_step = 10, atoms_speed = 1;

    unsigned const steps = 1000000;

    double energe[steps];

    size = CreateCube(array, init_lattice_step, cube_size , atoms_speed );
    
    if (steps) {
        for (int i = 0; i < steps; i++) {
            TimeStep(array, size);
            EnegeCulc(energe, array, size, i);
            if (i % 1000 == 0) {
                WriteCoorToFile(array, size, i);
            }
        }
    }

    WriteEnergyToFile(energe, steps);

    return 0;
}