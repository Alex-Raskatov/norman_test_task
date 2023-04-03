#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
#include <string>
#include <vector>

void some_function() {};

const double t = 1.0/1000.0;

void write_coords_to_file(double atom_x[], double atom_y[], double atom_z[], unsigned size, unsigned step) {
    
    std::string file_name = "./data/atoms_step_" + std::to_string(step) + ".xyz";

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

    std::ofstream fout_full("./data/full_energe.txt");
    std::ofstream fout_kin("./data/kin_energe.txt");
    std::ofstream fout_pot("./data/pot_energe.txt");

    for (int i = 0; i < steps; i++) {
            fout_full << full_energe[i] << '\n';
            fout_kin << kin_energe[i] << '\n';
            fout_pot << pot_energe[i] << '\n';
    }

    fout_full.close();
    fout_kin.close();
    fout_pot.close();
}

unsigned create_cube (double atom_x[], double atom_y[], double atom_z[], double lattice_step, unsigned cube_size) {

    unsigned counter = 0;

    for (unsigned i = 0; i < cube_size; i++) {
        for (unsigned j = 0; j < cube_size; j++) {
            for (unsigned k = 0; k < cube_size; k++) {
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

    unsigned const step = 1000000; //количество шагов интегрирования

    unsigned const init_size = 65; // количество частиц

    // объявление массивов частиц

    double atoms_x[init_size] = {0}, atoms_y[init_size] = {0}, atoms_z[init_size] = {0};
    double atoms_speed_x[init_size] = {0}, atoms_speed_y[init_size] = {0}, atoms_speed_z[init_size] = {0};
    double atoms_force_x[init_size] = {0}, atoms_force_y[init_size] = {0}, atoms_force_z[init_size] = {0};

    // double kin_energe[step] = {0}, pot_energe[step] = {0}, full_energe[step] = {0};

    std::vector<double> kin_energe(step), pot_energe(step), full_energe(step);

    for (unsigned i = 0; i < step; i++) {
        kin_energe[i] = 0;
        pot_energe[i] = 0;
        full_energe[i] = 0;
    }

    // объявление вспомогательных переменных

    double r_square = 0, delta_x = 0, delta_y = 0, delta_z = 0, force_to_atom = 0, init_lattice_step = 5, atoms_speed = 0.001; 

    // создание капли или другой конфигурации

    unsigned size = create_cube(atoms_x, atoms_y, atoms_z, init_lattice_step, 4);

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
                delta_y = atoms_y[j] - atoms_y[k];
                delta_z = atoms_z[j] - atoms_z[k];
                
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
            atoms_y[j] += atoms_speed_y[j]*t;
            atoms_z[j] += atoms_speed_z[j]*t;
        }

        // считаем полную энергию

        full_energe[i] = kin_energe[i] + pot_energe[i];

        if (i % 5000 == 0) {
            write_coords_to_file(atoms_x, atoms_y, atoms_z, size, i);
        }
    }

    write_energy_to_file(full_energe, kin_energe, pot_energe, step);

    return 0;
}