#include "cylinder.h"

cylinder::cylinder(int x, int y, int z, int r, int h, int r_o, int R, int G, int B, int alpha)
{
    this->x = x;
    this->y = y;
    this->z = z;
    this->r = r;
    this->r_o = r_o;
    this->h = h;
    this->R = R;
    this->G = G;
    this->B = B;
    this->alpha = alpha;
}

int cylinder::get_x(){
    return this->x;
}

int cylinder::get_y(){
    return this->y;
}

int cylinder::get_z(){
    return this->z;
}

int cylinder::get_r(){
    return this->r;
}

int cylinder::get_r_o(){
    return this->r_o;
}

int cylinder::get_h(){
    return this->h;
}

void cylinder::set_va(int vx, int vy, int vz, int ax, int ay, int az){

    this->vx = vx;
    this->vy = vy;
    this->vz = vz;
    this->ax = ax;
    this->ay = ay;
    this->az = az;
}

void cylinder::set_color(int red, int green, int blue){

    this->R = red;
    this->G = green;
    this->B = blue;
}

int cylinder::get_r_color(){
    return this->R;
}

int cylinder::get_g_color(){
    return this->G;
}

int cylinder::get_b_color(){
    return this->B;
}

int cylinder::get_alpha_color(){
    return this->alpha;
}

void cylinder::move(int t){

    this->x = this->x + vx * t + ax * t * t / 2;
    this->y = this->y + vy * t + ay * t * t / 2;
    this->z = this->z + vz * t + az * t * t / 2;

    this->vx += this->ax * t;
    this->vy += this->ay * t;
    this->vz += this->az * t;
}

