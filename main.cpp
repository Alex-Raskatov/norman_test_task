#include <iostream>
#include <cmath>

const double sigma = 1, epsilon = 1, mass = 1;
const double tau = std::sqrt(std::pow(sigma, 2)*mass/epsilon);
const double t = tau/1000;

class Atom {

    private:
    //coordinats
        double x_coor;
        double y_coor;
        double z_coor;

    //speed
        double x_speed;
        double y_speed;
        double z_speed;

    //force
        double x_force;
        double y_force;
        double z_force;

    public:

        Atom(double x, double y, double z, double v_x, double v_y, double v_z){

            x_coor = x;
            y_coor = y;
            z_coor = z;

            x_speed = v_x;
            y_speed = v_y;
            z_speed = v_z;
        }

        double get_x_coor() {
            return x_coor;
        }

        double get_y_coor() {
            return y_coor;
        }

        double get_z_coor() {
            return z_coor;
        }

        void set_x_coor(double x) {
            x_coor = x;
        }

        void set_y_coor(double y) {
            y_coor = y;
        }

        void set_z_coor(double z) {
            z_coor = z;
        }

        double get_x_speed() {
            return x_speed;
        }

        double get_y_speed() {
            return y_speed;
        }

        double get_z_speed() {
            return z_speed;
        }

        void set_x_speed(double v_x) {
            x_speed = v_x;
        }

        void set_y_speed(double v_y) {
            y_speed = v_y;
        }

        void set_z_speed(double v_z) {
            z_speed = v_z;
        }

        double get_x_force() {
            return x_force;
        }

        double get_y_force() {
            return y_force;
        }

        double get_z_force() {
            return z_force;
        }

        void null_force() {
            x_force = 0;
            y_force = 0;
            z_force = 0;
        }

        void calculate_force(Atom atom) {

            double r = std::sqrt(x_coor*x_coor + y_coor*y_coor + z_coor*z_coor);
            double norm_r_x = (atom.get_x_coor() - x_coor) / r;
            double norm_r_y = (atom.get_y_coor() - y_coor) / r;
            double norm_r_z = (atom.get_z_coor() - z_coor) / r;

            double force = -24*( std::pow((1/r),7) - 2*std::pow((1/r), 13));

            x_force += force*norm_r_x;
            y_force += force*norm_r_y;
            z_force += force*norm_r_z;
        }

        void move() {
            
            x_speed += x_force*t;
            y_speed += y_force*t;
            z_speed += z_force*t;

            x_coor += x_speed*t;
            y_coor += y_speed*t;
            z_coor += z_speed*t;
        }
};

