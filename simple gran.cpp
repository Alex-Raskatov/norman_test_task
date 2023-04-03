#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
#include <string>
#include <vector>

void some_function() {};

const double t = 1.0/1000.0;
//const double L = 50;

void write_coords_to_file(double atom_x[], double atom_y[], double atom_z[], unsigned size, unsigned step) {
    
    std::string file_name = "./data_gran/atoms_step_" + std::to_string(step) + ".xyz";

    std::ofstream fout(file_name);

    fout << size << '\n' << "some atoms" << '\n';

    for (int i = 0; i < size; i++) {
        fout 
            << 'H' << ' '
            << std::to_string(atom_x[i]) << ' ' 
            << std::to_string(atom_y[i]) << ' '
            << std::to_string(atom_z[i]) << ' '
            << '\n'; 
    }

    fout.close();
}

void write_energy_to_file(std::vector<double>& full_energe, std::vector<double>& kin_energe, std::vector<double>& pot_energe, unsigned steps) {

    std::ofstream fout_full("./data_gran/full_energe.txt");
    std::ofstream fout_kin("./data_gran/kin_energe.txt");
    std::ofstream fout_pot("./data_gran/pot_energe.txt");

    for (int i = 0; i < steps; i++) {
            fout_full << full_energe[i] << '\n';
            fout_kin << kin_energe[i] << '\n';
            fout_pot << pot_energe[i] << '\n';
    }

    fout_full.close();
    fout_kin.close();
    fout_pot.close();
}

void write_maxwell_to_file(std::vector<double>& maxwell_abs, std::vector<double>& maxwell_x, std::vector<double>& maxwell_y, std::vector<double>& maxwell_z,  unsigned maxwell_size) {

    std::ofstream m_x("./data_gran/maxwell_x.txt");
    std::ofstream m_y("./data_gran/maxwell_y.txt");
    std::ofstream m_z("./data_gran/maxwell_z.txt");
    std::ofstream m_a("./data_gran/maxwell_a.txt");

    for (int i = 0; i < maxwell_size; i++) {
            m_a << maxwell_abs[i] << '\n';
            m_x << maxwell_x[i] << '\n';
            m_y << maxwell_y[i] << '\n';
            m_z << maxwell_z[i] << '\n';
    }

    m_a.close();
    m_x.close();
    m_y.close();
    m_z.close();
}

unsigned create_cube (double atom_x[], double atom_y[], double atom_z[], double lattice_step, int cube_size) {

    unsigned counter = 0;

    for (int i = -cube_size/2; i < cube_size/2; i++) {
        for (int j = -cube_size/2; j < cube_size/2; j++) {
            for (int k = -cube_size/2; k < cube_size/2; k++) {
                atom_x[counter] = lattice_step*i;
                atom_y[counter] = lattice_step*j;
                atom_z[counter] = lattice_step*k;
                counter++;
            }
        }
    }

    return counter;
}

int main() {

    unsigned const step = 300000; //количество шагов интегрирования

    unsigned const init_size = 300; // количество частиц

    // объявление массивов частиц

    double atoms_x[init_size] = {0}, atoms_y[init_size] = {0}, atoms_z[init_size] = {0};
    double atoms_speed_x[init_size] = {0}, atoms_speed_y[init_size] = {0}, atoms_speed_z[init_size] = {0};
    double atoms_force_x[init_size] = {0}, atoms_force_y[init_size] = {0}, atoms_force_z[init_size] = {0};

    // векторы для максвелла

    int maxwell_size = 100;

    std::vector<double> maxwell_x(maxwell_size), maxwell_y(maxwell_size), maxwell_z(maxwell_size), maxwell_abs(maxwell_size);

    for (int i = 0; i < maxwell_size; i++) {
        maxwell_abs[i] = 0;
        maxwell_x[i] = 0;
        maxwell_y[i] = 0;
        maxwell_z[i] = 0;
    }

    // векторы для энергий

    std::vector<double> kin_energe(step), pot_energe(step), full_energe(step);

    for (unsigned i = 0; i < step; i++) {
        kin_energe[i] = 0;
        pot_energe[i] = 0;
        full_energe[i] = 0;
    }

    // объявление вспомогательных переменных

    double r_square = 0, delta_x = 0, delta_y = 0, delta_z = 0, force_to_atom = 0, init_lattice_step = 5, atoms_speed = 0.001; 
    int cube_size = 5;

    // размеры ячейки

    const double L = init_lattice_step*(cube_size+2);

    // создание капли или другой конфигурации

    unsigned size = create_cube(atoms_x, atoms_y, atoms_z, init_lattice_step, cube_size);
    std::cout << size;

    // главный цикл интегрирования

    for (unsigned i = 0; i < step; i++) {

        //обнуляем все силы

        for(unsigned j = 0; j < size; j++) {
            atoms_force_x[j] = 0;
            atoms_force_y[j] = 0;
            atoms_force_z[j] = 0;
        }

        // считаем все силы и потенциальные энергии

        for (unsigned j = 0; j < size; j++) {
            for (unsigned k = j + 1; k < size; k++) {
                
                delta_x = atoms_x[j] - atoms_x[k];
                if (delta_x > L/2) delta_x = -(L - delta_x);
                else if (delta_x < - L/2) delta_x = L + delta_x;
                delta_y = atoms_y[j] - atoms_y[k];
                if (delta_y > L/2) delta_y = -(L - delta_y);
                else if (delta_y < - L/2) delta_y = L + delta_y;
                delta_z = atoms_z[j] - atoms_z[k];
                if (delta_z > L/2) delta_z = -(L - delta_z);
                else if (delta_z < - L/2) delta_z = L + delta_z;
                
                r_square = delta_x*delta_x + delta_y*delta_y + delta_z*delta_z; // квадрат расстояния между частими

                force_to_atom = -24*( std::pow((1/r_square),4) - 2*std::pow((1/r_square), 7)); // сила взаимодействия делить на r
                
                // добавляем силу к атому из главного цикла
                atoms_force_x[j] += force_to_atom*delta_x;
                atoms_force_y[j] += force_to_atom*delta_y;
                atoms_force_z[j] += force_to_atom*delta_z;

                // и к атому из вторичного цикла для того чтобы не считать 2 раза
                atoms_force_x[k] -= force_to_atom*delta_x;
                atoms_force_y[k] -= force_to_atom*delta_y;
                atoms_force_z[k] -= force_to_atom*delta_z;

                // потенциальная энергия от взаимодействия j и k частицы

                pot_energe[i] += 4*(std::pow((1/r_square),6) - std::pow((1/r_square), 3));
            }
        }

        // считаем скорости, смещение 

        for (unsigned j = 0; j < size; j++) {

            // считаем скорости

            atoms_speed_x[j] += atoms_force_x[j]*t;
            atoms_speed_y[j] += atoms_force_y[j]*t;
            atoms_speed_z[j] += atoms_force_z[j]*t;

            // считаем кинетическую энергию частицы

            kin_energe[i] += (atoms_speed_x[j]*atoms_speed_x[j] + atoms_speed_y[j]*atoms_speed_y[j] + atoms_speed_z[j]*atoms_speed_z[j])/2;

            // считаем смещения

            atoms_x[j] += atoms_speed_x[j]*t;
            if (atoms_x[j] > L/2) atoms_x[j] -= L;
            else if (atoms_x[j] < -L/2) atoms_x[j] += L;
            atoms_y[j] += atoms_speed_y[j]*t;
            if (atoms_y[j] > L/2) atoms_y[j] -= L;
            else if (atoms_y[j] < -L/2) atoms_y[j] += L;
            atoms_z[j] += atoms_speed_z[j]*t;
            if (atoms_z[j] > L/2) atoms_z[j] -= L;
            else if (atoms_z[j] < -L/2) atoms_z[j] += L;
        }

        // считаем полную энергию

        full_energe[i] = kin_energe[i] + pot_energe[i];

        int scale = 10;

        //if (i % 1000 == 0) write_coords_to_file(atoms_x, atoms_y, atoms_z, size, i);

        if (i % 100 == 0 and i >= 100000) {

            double max_v = 0;
            for (int j = 0; j < size; j++) {
                if (std::abs(std::floor(scale*atoms_speed_x[j])) < maxwell_size)
                    maxwell_x[std::abs(std::floor(scale*atoms_speed_x[j]))]++;
                double v = scale*(atoms_speed_x[j]*atoms_speed_x[j] + atoms_speed_y[j]*atoms_speed_y[j] + atoms_speed_z[j]*atoms_speed_z[j]);
                if (max_v < v) { 
                    max_v = v;
                    std::cout << i << " " << max_v << std::endl;
                }
                if (v < maxwell_size)
                    maxwell_abs[std::floor(v)] += 1;
            }
        }


    }

    //write_energy_to_file(full_energe, kin_energe, pot_energe, step);
    write_maxwell_to_file(maxwell_abs, maxwell_x, maxwell_y,maxwell_z, maxwell_size);

    return 0;
}