#include "point.h"

point::point(int x, int y, int z, int r, int R, int G, int B)
{
    this->x = x;
    this->y = y;
    this->z = z;
    this->r = r;
    this->R = R;
    this->G = G;
    this->B = B;
}

int point::get_x(){
    return this->x;
}

int point::get_y(){
    return this->y;
}

int point::get_z(){
    return this->z;
}

int point::get_r(){
    return this->r;
}

int point::get_r_color(){
    return this->R;
}

int point::get_g_color(){
    return this->G;
}

int point::get_b_color(){
    return this->B;
}
