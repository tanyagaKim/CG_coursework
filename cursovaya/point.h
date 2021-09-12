#ifndef POINT_H
#define POINT_H


class point
{
private:
    int x, y, z;        // Центр электрона
    int r;
    int R = 0, G = 0, B = 0; // Цвет электрона
public:
    point(int x, int y, int z, int r, int R = 0, int G = 0, int B = 0);
    int get_x();
    int get_y();
    int get_z();
    int get_r();

    int get_r_color();
    int get_g_color();
    int get_b_color();
};

#endif // POINT_H
