#ifndef SPHERE_H
#define SPHERE_H

#include <iostream>

class sphere
{
private:
    int x, y, z;        // Центр электрона
    int R = 0, G = 0, B = 0, alpha = 255; // Цвет электрона
    int r;              // Радиус электрона
    float vx = 0, vy = 0, vz = 0; // Скорости электрона
    float ax = 0, ay = 0, az = 0; // Ускорения электрона
public:
    sphere(int x, int y, int z, int r, int R = 0, int G = 0, int B = 0, int alpha = 255);
    int get_x();
    int get_y();
    int get_z();
    int get_r();
    int get_vx();
    int get_vy();
    int get_vz();

    void set_va(int vx, int vy, int vz, int ax = 0, int ay = 0, int az = 0);
    void move(int t);

    void set_color(int red, int green, int blue);
    int get_r_color();
    int get_g_color();
    int get_b_color();
    int get_alpha_color();

};

#endif // SPHERE_H
