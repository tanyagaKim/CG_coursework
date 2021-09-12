#include "prisma.h"

prisma::prisma(int x, int y, int z, int l, int w, int h, int R, int G, int B, int alpha)
{
    this->x = x;
    this->y = y;
    this->z = z;
    this->l = l;
    this->w = w;
    this->h = h;
    this->R = R;
    this->G = G;
    this->B = B;
    this->alpha = alpha;
}

int prisma::get_x(){
    return this->x;
}

int prisma::get_y(){
    return this->y;
}

int prisma::get_z(){
    return this->z;
}

int prisma::get_l(){
    return this->l;
}

int prisma::get_w(){
    return this->w;
}

int prisma::get_h(){
    return this->h;
}

void prisma::set_va(int vx, int vy, int vz, int ax, int ay, int az){

    this->vx = vx;
    this->vy = vy;
    this->vz = vz;
    this->ax = ax;
    this->ay = ay;
    this->az = az;
}

void prisma::set_color(int red, int green, int blue){

    this->R = red;
    this->G = green;
    this->B = blue;
}

int prisma::get_r_color(){
    return this->R;
}

int prisma::get_g_color(){
    return this->G;
}

int prisma::get_b_color(){
    return this->B;
}

int prisma::get_alpha_color(){
    return this->alpha;
}

void prisma::move(int t){

    this->x = this->x + vx * t + ax * t * t / 2;
    this->y = this->y + vy * t + ay * t * t / 2;
    this->z = this->z + vz * t + az * t * t / 2;

    this->vx += this->ax * t;
    this->vy += this->ay * t;
    this->vz += this->az * t;
}

