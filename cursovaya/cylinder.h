#ifndef CYLINDER_H
#define CYLINDER_H

#include <iostream>

class cylinder{
private:
    int x, y, z;          // Центр цилиндра
    int R = 0, G = 0, B = 0, alpha = 255; // Цвет цилиндра
    int r, r_o = 0, h;          // Радиус, радиус отверстия, высота
    float vx = 0, vy = 0, vz = 0; // Скорости по xyz
    float ax = 0, ay = 0, az = 0; // Ускорения по xyz
public:
    cylinder(int x, int y, int z, int r, int h, int r_o = 0, int R = 0, int G = 0, int B = 0, int alpha = 255);
    int get_x();
    int get_y();
    int get_z();
    int get_r();
    int get_r_o();
    int get_h();

    void set_va(int vx, int vy, int vz, int ax = 0, int ay = 0, int az = 0);
    void move(int t);

    void set_color(int red, int green, int blue);
    int get_r_color();
    int get_g_color();
    int get_b_color();
    int get_alpha_color();
};

#endif // CYLINDER_H
