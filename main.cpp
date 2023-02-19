#include <iostream>
#include <cmath>

const double sigma = 1, epsilon = 1, mass = 1;
const double tau = std::sqrt(std::pow(sigma, 2)*mass/epsilon);
const double t = tau/1000;

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
        x += vec.x;
        y += vec.y;
        z += vec.z;
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

        Atom(double x, double y, double z, double v_x, double v_y, double v_z, double f_x = 0, double f_y = 0, double f_z = 0): coordinates(0,0,0), speed(0,0,0), force(0,0,0) {
            
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

            double r = std::sqrt(coordinates.GetX()*coordinates.GetX()
                                 + coordinates.GetY()*coordinates.GetY()
                                 + coordinates.GetZ()*coordinates.GetZ());
            double norm_r_x = (atom.GetX() - coordinates.GetX()) / r;
            double norm_r_y = (atom.GetY() - coordinates.GetY()) / r;
            double norm_r_z = (atom.GetZ() - coordinates.GetZ()) / r;

            double force_to_atom = -24*( std::pow((1/r),7) - 2*std::pow((1/r), 13));

            Vector force_to_atom_vec(force_to_atom*norm_r_x,
                            force_to_atom*norm_r_y,
                            force_to_atom*norm_r_z);

            force.AddVector(force_to_atom_vec);
        }

        void move() {
            
            speed.AddVector(force.MultiOnConstReturn(t));
            coordinates.AddVector(speed.MultiOnConstReturn(t));
        }
};

int CreateDrop (double drops_radius, Vector drops_center, double lattice_step, double atoms_max_speed_proection, Vector drops_speed) {
    return 0;
}

int main () {
    return 0;
}