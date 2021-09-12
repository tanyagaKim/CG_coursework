#include "sphere.h"

sphere::sphere(int x, int y, int z, int r, int R, int G, int B, int alpha)
{
    this->x = x;
    this->y = y;
    this->z = z;
    this->r = r;
    this->R = R;
    this->G = G;
    this->B = B;
    this->alpha = alpha;
}

int sphere::get_x(){
    return this->x;
}

int sphere::get_y(){
    return this->y;
}

int sphere::get_z(){
    return this->z;
}

int sphere::get_r(){
    return this->r;
}

int sphere::get_vx(){
    return this->vx;
}

int sphere::get_vy(){
    return this->vy;
}

int sphere::get_vz(){
    return this->vz;
}

void sphere::set_va(int vx, int vy, int vz, int ax, int ay, int az){

    this->vx = vx;
    this->vy = vy;
    this->vz = vz;
    this->ax = ax;
    this->ay = ay;
    this->az = az;
}

void sphere::set_color(int red, int green, int blue){

    this->R = red;
    this->G = green;
    this->B = blue;
}

int sphere::get_r_color(){
    return this->R;
}

int sphere::get_g_color(){
    return this->G;
}

int sphere::get_b_color(){
    return this->B;
}

int sphere::get_alpha_color(){
    return this->B;
}

void sphere::move(int t){

    this->x = this->x + vx * t + ax * t * t / 2;
    this->y = this->y + vy * t + ay * t * t / 2;
    this->z = this->z + vz * t + az * t * t / 2;

    this->vx += this->ax * t;
    this->vy += this->ay * t;
    this->vz += this->az * t;
}
