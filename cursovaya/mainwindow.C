#include "mainwindow.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    // Добавление сцены
    QImage imag (width, height, QImage::Format_RGB32);
    im = imag;

    // Создание массивов
    zbuf = new int[width * height];
    cadr = new int[width * height];

    // Заполнение фоновым цветом сцены
    QColor *c = new QColor(Rs, Gs, Bs);
    im.fill(*c);

    ui->label->setPixmap(QPixmap::fromImage(im));
    ui->label->show();

    // Создание цилиндров-катодов в элт
    c_list.push_back(*(new cylinder(80, 240, 0, 30, 70, 0, 100, 100, 200, 10))); // нить накаливания
    c_list.push_back(*(new cylinder(200, 240, 0, 90, 40, 70, 100, 200, 100, 100))); // фокусирующий анод
    c_list.push_back(*(new cylinder(300, 240, 0, 90, 40, 70, 100, 200, 100, 200))); // ускоряющий анод    <- прозрачная фигня, которая долго отрисовывается :(
    elt = new cylinder(250, 240, 0, 150, 500, 0, 220, 250, 250, 10); // корпус ЭЛТ
    screen = new prisma(480, 240, 0, 10, 180, 180, 220, 255, 255, 10); // экран
    dif = new prisma(390, 240, 0, 10, 160, 140, 200, 200, 250);     // диффракционная решетка
    int len = 100, line = ui->lineEdit->text().split(" ")[0].toInt();
    int num = len / line + 1, dl = dif->get_h() / num;
    for (int i = 0; i < num; i++)
        if (i % 2 == 1) dif_list.push_back(*(new prisma(dif->get_x(), dif->get_y() - dl * num / 2 + dl * i, dif->get_z(),
                                                        dif->get_l(), dif->get_w(), dl, dif->get_r_color(), dif->get_g_color(),
                                                        dif->get_b_color())));

    sph = new sphere(0, 0, 0, 150, 200, 255, 255);
    pr = new prisma(0, 0, 0, 100,100, 100, 200, 255, 255);
    cyl = new cylinder(0, 0, 0, 100, 300, 0, 200, 255, 255);
    // Заполнение кадра цветами элт и z-буфера
    draw_elt();

    // Создание таймера для анимации
    timer = new QTimer(this);
    timer->setInterval(50);
    connect(timer, SIGNAL(timeout()), this, SLOT(animation()));
}

MainWindow::~MainWindow()
{
    /*for (iter = c_list.begin(); iter != c_list.end(); iter++){
        { cylinder *zx = &(*iter); delete zx; }
    }*/
    delete timer;
    delete ui;
    free(zbuf);
    free(cadr);
}

float convert_pi(int a){
    return a * 3.14 / 180;
}

void multiply_matrix(float C[], float A[4][4], float B[4])
{
    for (int i = 0; i < 4; ++i)
    {
        C[i] = 0;
        for (int k = 0; k < 4; ++k)
            C[i] += A[i][k] * B[k];
    }
}

// Далее алгоритмы отрисовки - заполнение z-буфера и матрицы цветов

// Отрисовка сферы
void MainWindow::standartdrawcircle(sphere s){

    int r = s.get_r();
    int xc = s.get_x();
    int yc = s.get_y();
    int zc = s.get_z();
    int Rm = s.get_r_color();
    int Gm = s.get_g_color();
    int Bm = s.get_b_color();
    int alpha = s.get_alpha_color();

    float xyz[] = {xl, yl, zl, 1};
    float xyz01[] = {0, 0, 0, 1};
    multiply_matrix(xyz01, matrix, xyz);

    int xll = xyz01[0], yll = xyz01[1], zll = xyz01[2];

    float xyz0[] = {xc, yc, zc, 1};
    float xyz00[] = {0, 0, 0, 1};
    multiply_matrix(xyz00, matrix, xyz0);

    int xcc = xyz00[0], ycc = xyz00[1], zcc = xyz00[2];

    for (int x = -r; x <= r; x++)
        for(int y = -r; y <= r; y++){

            int z_2 = r * r - x * x - y * y;
            if (z_2 < 0) continue;

            float z = -sqrt(z_2);

            int xm = x + xcc + dx;
            int ym = y + ycc + dy;
            int zm = z + zcc;

            float xlll = xll - xcc - x, ylll = yll - y, zlll = zll - zcc - z;
            float k = 0.5 + 0.5 * (xlll * x + ylll * y + zlll * z) / sqrt(x * x + y * y + z * z) / sqrt(xlll * xlll + ylll * ylll + zlll * zlll);

            int R = k * Rm, G = k * Gm, B = k * Bm;

            int ind = xm + ym * width;
            if (xm < 0 || ym < 0 || xm >= width || ym >= height) continue;
            if (cadr[ind] > zm && (cadr[ind] != 1000 || cadr[ind + 1] == 1000 || cadr[ind - 1] == 1000)){

                cadr[ind] = zm;
                elect.setPixel(xm, ym, qRgba(R, G, B, 255));
            }
        }
}

// Отрисовка каркаса элт
void MainWindow::standartdrawcarcascylinder(cylinder c){

    int h = c.get_h();
    int r = c.get_r();
    int xc = c.get_x();
    int yc = c.get_y();
    int zc = c.get_z();
    int rc = c.get_r_o();
    int Rm = c.get_r_color();
    int Gm = c.get_g_color();
    int Bm = c.get_r_color();
    int alpha = c.get_alpha_color();

    float step = (float)1/r;

    for(float t = 0; t <= 2 * M_PI; t += step){

        float x1 = xc - h/2, y1 = yc + r * sin(t), z1 = zc + r * cos(t);
        float x3 = xc + h/2, y3 = y1, z3 = z1;

        float xyz[] = {x1, y1, z1, 1};
        float xyz01[] = {0, 0, 0, 1};
        multiply_matrix(xyz01, matrix, xyz);
        x1 = xyz01[0] + dx, y1 = xyz01[1] + dy, z1 = xyz01[2];

        float xyz2[] = {x3, y3, z3, 1};
        float xyz02[] = {0, 0, 0, 1};
        multiply_matrix(xyz02, matrix, xyz2);
        x3 = xyz02[0] + dx, y3 = xyz02[1] + dy, z3 = xyz02[2];

        if (x1 < 0 || y1 < 0 || x1 >= width || y1 >= height) continue;

        int ind = x1 + y1 * width;
        if (zbuf[ind] > z1){

            zbuf[ind] = z1;
            im.setPixel(x1, y1, qRgba(10, 10, 10, 200));
        }
        if (x3 < 0 || y3 < 0 || x3 >= width || y3 >= height) continue;

        ind = x3 + y3 * width;
        if (zbuf[ind] > z3){

            zbuf[ind] = z3;
            im.setPixel(x3, y3, qRgba(10, 10, 10, 200));
        }
    }

    int y1 = yc + r, y2 = yc - r;
    int z1 = zc + r, z2 = zc - r;
    for (int xt = -h/2; xt < h/2; xt++){

        float xyz[] = {xc + xt, y1, zc, 1};
        float xyz01[] = {0, 0, 0, 1};
        multiply_matrix(xyz01, matrix, xyz);
        int xp = xyz01[0] + dx, yp = xyz01[1] + dy, zp = xyz01[2];

        if (xp < 0 || yp < 0 || xp >= width || yp >= height) continue;

        int ind = xp + yp * width;
        if (zbuf[ind] > zp){

            zbuf[ind] = zp;
            im.setPixel(xp, yp, qRgba(10, 10, 10, 200));
        }

        float xyz2[] = {xc + xt, y2, zc, 1};
        float xyz02[] = {0, 0, 0, 1};
        multiply_matrix(xyz02, matrix, xyz2);
        xp = xyz02[0] + dx, yp = xyz02[1] + dy, zp = xyz02[2];

        if (xp < 0 || yp < 0 || xp >= width || yp >= height) continue;

        ind = xp + yp * width;
        if (zbuf[ind] > zp){

            zbuf[ind] = zp;
            im.setPixel(xp, yp, qRgba(10, 10, 10, 200));
        }

        float xyz3[] = {xc + xt, yc, z1, 1};
        float xyz03[] = {0, 0, 0, 1};
        multiply_matrix(xyz03, matrix, xyz3);
        xp = xyz03[0] + dx, yp = xyz03[1] + dy, zp = xyz03[2];

        if (xp < 0 || yp < 0 || xp >= width || yp >= height) continue;

        ind = xp + yp * width;
        if (zbuf[ind] > zp){

            zbuf[ind] = zp;
            im.setPixel(xp, yp, qRgba(10, 10, 10, 200));
        }

        float xyz4[] = {xc + xt, yc, z2, 1};
        float xyz04[] = {0, 0, 0, 1};
        multiply_matrix(xyz04, matrix, xyz4);
        xp = xyz04[0] + dx, yp = xyz04[1] + dy, zp = xyz04[2];

        if (xp < 0 || yp < 0 || xp >= width || yp >= height) continue;

        ind = xp + yp * width;
        if (zbuf[ind] > zp){

            zbuf[ind] = zp;
            im.setPixel(xp, yp, qRgba(10, 10, 10, 200));
        }
    }
}

// Отрисовка призмы
void MainWindow::standartdrawprisma(prisma p){

    int h = p.get_h();
    int w = p.get_w();
    int l = p.get_l();
    int xc = p.get_x();
    int yc = p.get_y();
    int zc = p.get_z();
    int Rm = p.get_r_color();
    int Gm = p.get_g_color();
    int Bm = p.get_r_color();
    int alpha = p.get_alpha_color();

    // нижний и верхний горизонтальные прямоугольник
    float xt = 0, yt = -h/2, zt = 0;
    float xll = xl - xc, yll = yl - yc, zll = zl - zc + h/2;
    float k1 = 0.5 + 0.5 * (xll * xt + yll * yt + zll * zt) / sqrt(xt * xt + yt * yt + zt * zt) / sqrt(xll * xll + yll * yll + zll * zll);

    xt = 0, yt = h/2, zt = 0;
    xll = xl - xc, yll = yl - yc, zll = zl - zc - h/2;
    float k2 = 0.5 + 0.5 * (xll * xt + yll * yt + zll * zt) / sqrt(xt * xt + yt * yt + zt * zt) / sqrt(xll * xll + yll * yll + zll * zll);

    float y1 = yc - h/2, y2 = yc + h/2;

    for (int xt = -l/2; xt <= l/2; xt++)
        for (int zt = -w/2; zt <= w/2; zt++){

        float x = xc + xt, z = zc + zt;

        float xyz[] = {x, y1, z, 1};
        float xyz01[] = {0, 0, 0, 1};
        multiply_matrix(xyz01, matrix, xyz);
        int xm = xyz01[0] + dx, ym = xyz01[1] + dy, zm = xyz01[2];

        int ind = xm + ym * width;
        if (xm < 0 || ym < 0 || xm >= width || ym >= height) continue;
        if (zbuf[ind] > zm){
            zbuf[ind] = zm;
            im.setPixel(xm, ym, qRgba(Rm * k1, Gm * k1, Bm * k1, alpha));
        }
        float xyz2[] = {x, y2, z, 1};
        float xyz02[] = {0, 0, 0, 1};
        multiply_matrix(xyz02, matrix, xyz2);
        xm = xyz02[0] + dx, ym = xyz02[1] + dy, zm = xyz02[2];

        ind = xm + ym * width;
        if (xm < 0 || ym < 0 || xm >= width || ym >= height) continue;
        if (zbuf[ind] > zm){
            zbuf[ind] = zm;
            im.setPixel(xm, ym, qRgba(Rm * k2, Gm * k2, Bm * k2, alpha));
        }
    }

    xt = 0, yt = 0, zt = -w/2;
    xll = xl - xc, yll = yl - yc, zll = zl - zc - h/2;
    float k3 = 0.5 + 0.5 * (xll * xt + yll * yt + zll * zt) / sqrt(xt * xt + yt * yt + zt * zt) / sqrt(xll * xll + yll * yll + zll * zll);

    xt = 0, yt = 0, zt = w/2;
    xll = xl - xc, yll = yl - yc, zll = zl - zc - h/2;
    float k4 = 0.5 + 0.5 * (xll * xt + yll * yt + zll * zt) / sqrt(xt * xt + yt * yt + zt * zt) / sqrt(xll * xll + yll * yll + zll * zll);

    float z1 = zc - w/2, z2 = zc + w/2;

    for (int xt = -l/2; xt <= l/2; xt++)
        for (int yt = -h/2; yt <= h/2; yt++){

        float x = xc + xt, y = yc + yt;

        float xyz[] = {x, y, z1, 1};
        float xyz01[] = {0, 0, 0, 1};
        multiply_matrix(xyz01, matrix, xyz);
        int xm = xyz01[0] + dx, ym = xyz01[1] + dy, zm = xyz01[2];

        int ind = xm + ym * width;
        if (xm < 0 || ym < 0 || xm >= width || ym >= height) continue;
        if (zbuf[ind] > zm){
            zbuf[ind] = zm;
            im.setPixel(xm, ym, qRgba(Rm * k3, Gm * k3, Bm * k3, alpha));
        }
        float xyz2[] = {x, y, z2, 1};
        float xyz02[] = {0, 0, 0, 1};
        multiply_matrix(xyz02, matrix, xyz2);
        xm = xyz02[0] + dx, ym = xyz02[1] + dy, zm = xyz02[2];

        ind = xm + ym * width;
        if (xm < 0 || ym < 0 || xm >= width || ym >= height) continue;
        if (zbuf[ind] > zm){
            zbuf[ind] = zm;
            im.setPixel(xm, ym, qRgba(Rm * k4, Gm * k4, Bm * k4, alpha));
        }
    }

    xt = -l/2, yt = 0, zt = 0;
    xll = xl - xc, yll = yl - yc, zll = zl - zc - h/2;
    float k5 = 0.5 + 0.5 * (xll * xt + yll * yt + zll * zt) / sqrt(xt * xt + yt * yt + zt * zt) / sqrt(xll * xll + yll * yll + zll * zll);

    xt = l/2, yt = 0, zt = 0;
    xll = xl - xc, yll = yl - yc, zll = zl - zc - h/2;
    float k6 = 0.5 + 0.5 * (xll * xt + yll * yt + zll * zt) / sqrt(xt * xt + yt * yt + zt * zt) / sqrt(xll * xll + yll * yll + zll * zll);

    float x1 = xc - l/2, x2 = xc + l/2;

    for (int yt = -h/2; yt <= h/2; yt++)
        for (int zt = -w/2; zt <= w/2; zt++){

        float z = zc + zt, y = yc + yt;

        float xyz[] = {x1, y, z, 1};
        float xyz01[] = {0, 0, 0, 1};
        multiply_matrix(xyz01, matrix, xyz);
        int xm = xyz01[0] + dx, ym = xyz01[1] + dy, zm = xyz01[2];

        int ind = xm + ym * width;
        if (xm < 0 || ym < 0 || xm >= width || ym >= height) continue;
        if (zbuf[ind] > zm){
            zbuf[ind] = zm;
            im.setPixel(xm, ym, qRgba(Rm * k5, Gm * k5, Bm * k5, alpha));
        }
        float xyz2[] = {x2, y, z, 1};
        float xyz02[] = {0, 0, 0, 1};
        multiply_matrix(xyz02, matrix, xyz2);
        xm = xyz02[0] + dx, ym = xyz02[1] + dy, zm = xyz02[2];

        ind = xm + ym * width;
        if (xm < 0 || ym < 0 || xm >= width || ym >= height) continue;
        if (zbuf[ind] > zm){
            zbuf[ind] = zm;
            im.setPixel(xm, ym, qRgba(Rm * k6, Gm * k6, Bm * k6, alpha));
        }
    }
}

// Отрисовка цилиндра
void MainWindow::standartdrawcylinder(cylinder c){

    int h = c.get_h();
    int r = c.get_r();
    int xc = c.get_x();
    int yc = c.get_y();
    int zc = c.get_z();
    int rc = c.get_r_o();
    int Rm = c.get_r_color();
    int Gm = c.get_g_color();
    int Bm = c.get_r_color();
    int alpha = c.get_alpha_color();

    float step = (float)1/r;

    for(float t = 0; t <= 2 * M_PI; t += step){

        // горизонтально-расположенный цилиндр

        float x1 = xc - h/2, y1 = yc + r * sin(t), z1 = zc + r * cos(t);
        float x3 = xc + h/2, y3 = y1, z3 = z1;

        float xt = 0, yt = y1 - yc, zt = z1 - zc;
        float xll = xl - xc, yll = yl - y1, zll = zl - z1;
        float k = 0.5 + 0.5 * (xll * xt + yll * yt + zll * zt) / sqrt(xt * xt + yt * yt + zt * zt) / sqrt(xll * xll + yll * yll + zll * zll);

        float xyz[] = {x1, y1, z1, 1};
        float xyz01[] = {0, 0, 0, 1};
        multiply_matrix(xyz01, matrix, xyz);
        x1 = xyz01[0], y1 = xyz01[1], z1 = xyz01[2];

        float xyz2[] = {x3, y3, z3, 1};
        float xyz02[] = {0, 0, 0, 1};
        multiply_matrix(xyz02, matrix, xyz2);
        x3 = xyz02[0], y3 = xyz02[1], z3 = xyz02[2];

        int R = k * Rm, G = k * Gm, B = k * Bm;        
        QColor c(R, G, B, alpha);

        int Dx = abs(x3 - x1);
        int Dy = abs(y3 - y1);

        int sx = (x3 >= x1) ? 1 : -1;
        int sy = (y3 >= y1) ? 1 : -1;

        if (Dy <= Dx){

            int d = Dy * 2 - Dx;
            int d1 = Dy * 2;
            int d2 = (Dy - Dx) * 2;

            int x = x1;
            int y = y1;
            int z = z1;

            for (int i = 0; i <= Dx; i++, x += sx){

                if (d > 0){
                    d += d2;
                    y += sy;
                }
                else d += d1;
                int ind = (x + dx) + (y + dy) * width;
                if (x + dx - 1 < 0 || y + dy - 1 < 0 || x + dx >= width || y + dy >= height) continue;
                if (y3 != y1) z = (y - y1) * (float)(z3 - z1) / (y3 - y1) + z1;
                if (x3 != x1) z = (x - x1) * (float)(z3 - z1) / (x3 - x1) + z1;
                //if (zbuf[ind] > z && (zbuf[ind] != 1000 || zbuf[ind + 1] == 1000 || zbuf[ind - 1] == 1000)){
                if (zbuf[ind] > z){

                    zbuf[ind] = z;
                    im.setPixel(x + dx, y + dy, qRgba(R, G, B, 200));

                    zbuf[ind + 1] = z;
                    im.setPixel(x + dx + 1, y + dy, qRgba(R, G, B, 200));

                    zbuf[ind - 1] = z;
                    im.setPixel(x + dx - 1, y + dy, qRgba(R, G, B, 200));

                    im.setPixel(x + dx, y + dy + 1, qRgba(R, G, B, 200));

                    im.setPixel(x + dx, y + dy - 1, qRgba(R, G, B, 200));
                }
            }
        }
        else{
            int d = Dx * 2 - Dy;
            int d1 = Dx * 2;
            int d2 = (Dx - Dy) * 2;

            int x = x1;
            int y = y1;
            int z = z1;

            for (int i = 0; i <= Dy; i++, y += sy){
                if (d > 0){
                    d += d2;
                    x += sx;
                }
                else d += d1;
                int ind = (x + dx) + (y + dy) * width;
                if (x + dx - 1 < 0 || y + dy - 1 < 0 || x + dx >= width || y + dy >= height) continue;
                if (y3 != y1) z = (y - y1) * (float)(z3 - z1) / (y3 - y1) + z1;
                if (x3 != x1) z = (x - x1) * (float)(z3 - z1) / (x3 - x1) + z1;
                //if (zbuf[ind] > z && (zbuf[ind] != 1000 || zbuf[ind + 1] == 1000 || zbuf[ind - 1] == 1000)){
                if (zbuf[ind] > z){

                    zbuf[ind] = z;
                    im.setPixel(x + dx, y + dy, qRgba(R, G, B, 10));

                    zbuf[ind + 1] = z;
                    im.setPixel(x + dx + 1, y + dy, qRgba(R, G, B, 200));

                    zbuf[ind - 1] = z;
                    im.setPixel(x + dx - 1, y + dy, qRgba(R, G, B, 200));

                    im.setPixel(x + dx, y + dy + 1, qRgba(R, G, B, 200));

                    im.setPixel(x + dx, y + dy - 1, qRgba(R, G, B, 200));
                }
            }
        }
    }

    standartdrawroot(rc, r, h/2, xc , yc, zc, Rm, Gm, Bm, alpha);
    standartdrawroot(rc, r, -h/2, xc, yc, zc, Rm, Gm, Bm, alpha);
}

// Отрисовка крышки цилиндра
void MainWindow::standartdrawroot(int r1, int r2, int h, int xc, int yc, int zc, int Rm, int Gm, int Bm, int alpha){

    float k = 0.5 + 0.5 * (xl * h) / sqrt(h * h) / sqrt(xl * xl + yl * yl + zl * zl);

    int R = k * Rm, G = k * Gm, B = k * Bm;

    QColor c(R, G, B, alpha);

    for (int yt = -r2; yt <= r2; yt++)
        for (int zt = -r2; zt <= r2; zt++){

            int rt2 = zt * zt + yt * yt;
            if (rt2 < r1 * r1) continue;
            if (rt2 > r2 * r2) continue;

            float xyz[] = {float(xc + h), float(yt + yc), float(zt + zc), 1};
            float xyz01[] = {0, 0, 0, 1};
            multiply_matrix(xyz01, matrix, xyz);
            int x = xyz01[0] + dx, y = xyz01[1] + dy, z = xyz01[2];

            int ind = x + y * width;
            if (x < 0 || y < 0 || x >= width || y >= height) continue;
            if (zbuf[ind] > z){

                zbuf[ind] = z;
                im.setPixel(x, y, c.rgba());
            }
        }
}

// Отрисовка добавленного электрона
void MainWindow::standartdrawpoint(point p){

    int r = p.get_r();
    int xc = p.get_x();
    int yc = p.get_y();
    int zc = p.get_z();
    int Rm = p.get_r_color();
    int Gm = p.get_g_color();
    int Bm = p.get_b_color();

    for (int z = -r; z <= r; z++)
        for(int y = -r; y <= r; y++){

            int r_2 = z * z + y * y;
            if (r_2 > r * r) continue;

            float k = r_2 / r * r;

            float xyz0[] = {xc, yc + y, zc + z, 1};
            float xyz00[] = {0, 0, 0, 1};
            multiply_matrix(xyz00, matrix, xyz0);

            int xm = xyz00[0] + dx;
            int ym = xyz00[1] + dy;
            int zm = xyz00[2];

            int R = k * Rm, G = k * Gm, B = k * Bm;

            int ind = xm + ym * width;
            if (xm < 0 || ym < 0 || xm >= width || ym >= height) continue;
            if (cadr[ind] > zm){

                cadr[ind] = zm;
                elect.setPixel(xm, ym, qRgba(R, G, B, 255));
            }
        }
}

// Передвижение электронов
void MainWindow::animation(){

    if (time > max_time) { time = 0; timer->stop(); timer->start(); }

    if (!time)
        for (int i = 0; i < num; i++){
            sp_list.push_back(*(new sphere(x_st, y_st, z_st, 10)));
            sp_list.back().set_va( max_v/2 - rand() % max_v, max_v/2 - rand() % max_v, max_v/2 - rand() % max_v); // задание скоростей без ускорений
            sp_list.back().set_color(200, 200, 100);       // задание цвета
        }

    draw();

    for (it = sp_list.begin(); it != sp_list.end(); it++)
        (*it).move(1);
    time++;
}

// Заполнение кадра фоном и изображением элт (заполнение z-буфера и заполнение карты цветов)
void MainWindow::draw_elt(){

    im.fill(qRgb(Rs, Gs, Bs));

    // заполнение z-буфера максимальными значениями
    for (int x = 0; x < width; x++)
        for (int y = 0; y < height; y++)
            zbuf[x + y * width] = 1000;

    // отрисовка элт на кадре
    for (iter = c_list.begin(); iter != c_list.end(); iter++)
        standartdrawcylinder(*iter);

    // отрисовка диффракционной решетки на кадре
    for (dif_it = dif_list.begin(); dif_it != dif_list.end(); dif_it++)
        standartdrawprisma(*dif_it);
    standartdrawprisma(*screen);

    standartdrawcarcascylinder(*elt);
}

// Проверка электрона на существование
int MainWindow::check_electron(sphere s, int t){

    // текущее положение электрона
    int x1 = s.get_x();
    int y1 = s.get_y();
    int z1 = s.get_z();

    // следующее положение электрона
    int x2 = x1 + s.get_vx() * t;
    int y2 = y1 + s.get_vy() * t;
    int z2 = z1 + s.get_vz() * t;

    // Удаление электрона, если он вышел за пределы сцены
    {
        // проецирование электрона на плоскость
        //std::cout << x1 << " " << y1 << " " << z1 << "\n";
        float xyz[] = {float(x1), float(y1), float(z1), 1};
        float xyz01[] = {0, 0, 0, 1};
        multiply_matrix(xyz01, matrix, xyz);

        //std::cout << xyz01[1] + dy << " " <<  float(height) << "\n";

        if (int(xyz01[0] + dx) < 0 || int(xyz01[0] + dx) >= width || int(xyz01[1] + dy) < 0 || int(xyz01[1] + dy) >= height) return 1;
    }

    //Удаление электрона,если он пересекся с корпусом элт
    {
        int xc = elt->get_x();  // Параметры цилиндра элт
        int yc = elt->get_y();
        int zc = elt->get_z();
        int r = elt->get_r();
        int r_o = elt->get_r_o(); // Радиус отверстия
        int h = elt->get_h();
        int r_2 = r * r;
        int r_o_2 = r_o * r_o;

        // В проекции на оуz расстояние от центра цилиндра до двух точек прохождения электрона

        int r1 = (z1 - zc) * (z1 - zc) + (y1 - yc) * (y1 - yc);
        int r2 = (z2 - zc) * (z2 - zc) + (y2 - yc) * (y2 - yc);

        if ((r1 > r_2 && r2 <= r_2 && xc - h/2 < x2 && xc + h/2 > x2) || // если электрон вылетел через боковую поверхность цилиндра
                (r1 <= r_2 && r2 > r_2 && xc - h/2 < x1 && xc + h/2 > x1) || // если электрон влетел в боковую поверхность цилиндра
                ((xc - h/2 >= x1 && xc - h/2 <= x2) || (xc + h/2 >= x1 && xc + h/2 <= x2) ||  // если электрон попал в один из торцов цилиндра
                (xc - h/2 >= x2 && xc - h/2 <= x1) || (xc + h/2 >= x2 && xc + h/2 <= x1)))
            return 1;
    }

    // Удаление электрона, если он пересекается с поверхностью какого-либо цилиндра
    iter = c_list.begin(); // первый цилиндр - нить накаливания - испускает электроны
    for (++iter; iter != c_list.end(); iter++){

        int xc = iter->get_x();  // Параметры цилиндра
        int yc = iter->get_y();
        int zc = iter->get_z();
        int r = iter->get_r();
        int r_o = iter->get_r_o(); // Радиус отверстия
        int h = iter->get_h();
        int r_2 = r * r;
        int r_o_2 = r_o * r_o;

        // В проекции на оуz расстояние от центра цилиндра до двух точек прохождения электрона

        int r1 = (z1 - zc) * (z1 - zc) + (y1 - yc) * (y1 - yc);
        int r2 = (z2 - zc) * (z2 - zc) + (y2 - yc) * (y2 - yc);

        if ((r1 >= r_2 && r2 <= r_2 && xc - h/2  <= x2 && xc + h/2 >= x2) || // если электрон вылетел через боковую поверхность цилиндра
                (r1 <= r_2 && r2 >= r_2 && xc - h/2 <= x1 && xc + h/2 >= x1) || // если электрон влетел в боковую поверхность цилиндра
                ((x1 <= xc - h/2 && xc - h/2 <= x2) || (x1 <= xc + h/2 && xc + h/2 <= x2))&&
                 //r2 <= r_2 && r2 >= r_o_2)
                 !(r1 <= r_o_2 && r2 <= r_o_2 || r1 >= r_2 && r2 >= r_2)) // если электрон попал в один из торцов цилиндра
            return 1;
    }

    // Удаление электронов на диффракционной решетке
    for (dif_it = dif_list.begin(); dif_it != dif_list.end(); dif_it++){

        int xc = (*dif_it).get_x();  // Параметры призмы
        int yc = (*dif_it).get_y();
        int zc = (*dif_it).get_z();
        int l = (*dif_it).get_l();
        int w = (*dif_it).get_w();
        int h = (*dif_it).get_h();

        if (x1 <= xc - l/2 && xc - l/2 <= x2 && yc - h/2 <= y2 && y2 <= yc + h/2 &&
                zc - w/2 <= z2 && z2 <= zc + w/2) return 1;
    }

    // Проявление на экране
    {
        int xc = screen->get_x();  // Параметры призмы
        int yc = screen->get_y();
        int zc = screen->get_z();
        int l = screen->get_l();
        int w = screen->get_w();
        int h = screen->get_h();

        if (x1 <= xc - l/2 && xc - l/2 <= x2 && yc - h/2 <= y2 && y2 <= yc + h/2 &&
                zc - w/2 <= z2 && z2 <= zc + w/2) return -1;
    }

    return 0;
}

// Отрисовка на экране кадра (элт + электроны)
void MainWindow::draw(){

    // Очистка сцены и заполнение её фоновым цветом
    elect = QImage(im);

    // копирование z-буфера элт
    for (int x = 0; x < width; x++)
        for (int y = 0; y < height; y++)
            cadr[x + y * width] = zbuf[x + y * width];

    // отрисовка электронов на кадр (при их существовании)
    for (it = sp_list.begin(); it != sp_list.end(); it++){

        switch (check_electron((*it), 1)) {
            case 1:{ sp_list.erase(it); continue; }
            case -1:{
                point_list.push_back(*(new point((*it).get_x(), (*it).get_y(), (*it).get_z(), (*it).get_r(),
                                                 (*it).get_r_color(), (*it).get_g_color(), (*it).get_b_color())));
                sp_list.erase(it); continue; }
        }
        standartdrawcircle(*it);
    }

    // отрисовка добавленных электронов
    for (point_it = point_list.begin(); point_it != point_list.end(); point_it++)
       standartdrawpoint(*point_it);

    // Отрисовка корпуса ЭЛТ
    //standartdrawcylinder(*elt, cadr);

    // отрисовка сцены
    ui->label->setPixmap(QPixmap::fromImage(elect));
    ui->label->show();
}

void MainWindow::show_sphere(){

    if (time > max_time) { time = 0; timer->stop(); timer->start(); }

    im.fill(qRgb(Rs, Gs, Bs));

    // заполнение z-буфера максимальными значениями
    for (int x = 0; x < width; x++)
        for (int y = 0; y < height; y++)
            zbuf[x + y * width] = 1000;

    //standartdrawcircle(*sph);
    //standartdrawprisma(*pr);
    standartdrawcylinder(*cyl);

    ui->label->setPixmap(QPixmap::fromImage(im));
    ui->label->show();

    a += 1;
    if (b > 240) a -= 2;
    b += 5;

    if (a > 270 || b > 360)  { a = 0; b = 0; timer->stop(); }

    float pi_a = convert_pi(a), pi_b = convert_pi(b);

    float cosa = cos(pi_a), sina = sin(pi_a), cosb = cos(pi_b), sinb = sin(pi_b);

    matrix [0][0] = cosa, matrix [0][1] = -sina;
    matrix [1][0] = sina * cosb, matrix [1][1] = cosa * cosb;
    matrix [2][0] = sina * sinb, matrix [2][1] = cosa * sinb;
    matrix [1][2] = -sinb, matrix [2][2] = cosb;

    time++;
}

// кнопка запуска анимации
void MainWindow::on_pushButton_clicked()
{    
    dif_list.clear();

    // количество электронов
    num = ui->textEdit_4->toPlainText().split(" ")[0].toInt();

    int len = 100, line = ui->lineEdit->text().split(" ")[0].toInt();
    int num = len / line + 1, dl = dif->get_h() / num;
    for (int i = 0; i < num; i++)
        if (i % 2 == 1) dif_list.push_back(*(new prisma(dif->get_x(), dif->get_y() - dl * num / 2 + dl * i, dif->get_z(),
                                                        dif->get_l(), dif->get_w(), dl, dif->get_r_color(), dif->get_g_color(),
                                                        dif->get_b_color())));

    int val = ui->lineEdit_2->text().split(" ")[0].toInt();
    max_v = val / 2;

    draw_elt();

    timer->start(); // запуск анимации
}

// Прекращение анимации (остановка и сброс таймера)
void MainWindow::on_pushButton_9_clicked()
{
    sp_list.clear(); // очистка списка электронов и заполнение новыми значениями
    point_list.clear(); // очистка списка электронов на экране
    time = 0;
    timer->stop();
    draw_elt();
    draw();
}

// Дальше какие-то повороты и перемещения
void MainWindow::on_pushButton_5_clicked()
{
    a += 10;
    if (a >= 0) a -= 360;

    float pi_a = convert_pi(a), pi_b = convert_pi(b);
    float cosa = cos(pi_a), sina = sin(pi_a), cosb = cos(pi_b), sinb = sin(pi_b);

    matrix [0][0] = cosa, matrix [0][1] = -sina;
    matrix [1][0] = sina * cosb, matrix [1][1] = cosa * cosb;
    matrix [2][0] = sina * sinb, matrix [2][1] = cosa * sinb;

    draw_elt();
    draw();
}

void MainWindow::on_pushButton_4_clicked()
{
    a -= 10;
    if (a <= -360) a += 360;

    float pi_a = convert_pi(a), pi_b = convert_pi(b);
    float cosa = cos(pi_a), sina = sin(pi_a), cosb = cos(pi_b), sinb = sin(pi_b);

    matrix [0][0] = cosa, matrix [0][1] = -sina;
    matrix [1][0] = sina * cosb, matrix [1][1] = cosa * cosb;
    matrix [2][0] = sina * sinb, matrix [2][1] = cosa * sinb;

    draw_elt();
    draw();
}

void MainWindow::on_pushButton_6_clicked()
{
    b += 10;
    if (b >= 360) b -= 360;

    float pi_a = convert_pi(a), pi_b = convert_pi(b);
    float cosa = cos(pi_a), sina = sin(pi_a), cosb = cos(pi_b), sinb = sin(pi_b);

    matrix [1][2] = -sinb, matrix [2][2] = cosb;
    matrix [1][0] = sina * cosb, matrix [1][1] = cosa * cosb;
    matrix [2][0] = sina * sinb, matrix [2][1] = cosa * sinb;

    draw_elt();
    draw();
}

void MainWindow::on_pushButton_7_clicked()
{
    b -= 10;
    if (b < 0) b += 360;

    float pi_a = convert_pi(a), pi_b = convert_pi(b);
    float cosa = cos(pi_a), sina = sin(pi_a), cosb = cos(pi_b), sinb = sin(pi_b);

    matrix [1][2] = -sinb, matrix [2][2] = cosb;
    matrix [1][0] = sina * cosb, matrix [1][1] = cosa * cosb;
    matrix [2][0] = sina * sinb, matrix [2][1] = cosa * sinb;

    draw_elt();
    draw();
}

void MainWindow::on_pushButton_11_clicked()
{
    dx -= 30;

    draw_elt();
    draw();
}

void MainWindow::on_pushButton_10_clicked()
{
    dx += 30;

    draw_elt();
    draw();
}

void MainWindow::on_pushButton_12_clicked()
{
    dy -= 30;

    draw_elt();
    draw();
}

void MainWindow::on_pushButton_13_clicked()
{
    dy += 30;

    draw_elt();
    draw();
}

void MainWindow::on_pushButton_17_clicked()
{
    yl -= 40;

    draw_elt();
    draw();
}

void MainWindow::on_pushButton_14_clicked()
{
    yl += 40;

    draw_elt();
    draw();
}

void MainWindow::on_pushButton_15_clicked()
{
    zl -= 40;

    draw_elt();
    draw();
}

void MainWindow::on_pushButton_16_clicked()
{
    zl += 40;

    draw_elt();
    draw();
}

void MainWindow::on_pushButton_18_clicked()
{
    time = 0;
    timer->stop();
}

void MainWindow::on_pushButton_8_clicked()
{
    timer->stop();
    dif_list.clear();

    num = ui->textEdit_4->toPlainText().split(" ")[0].toInt();

    int len = 100, line = ui->lineEdit->text().split(" ")[0].toInt();
    int num = len / line + 1, dl = dif->get_h() / num;
    for (int i = 0; i < num; i++)
        if (i % 2 == 1) dif_list.push_back(*(new prisma(dif->get_x(), dif->get_y() - dl * num / 2 + dl * i, dif->get_z(),
                                                        dif->get_l(), dif->get_w(), dl, dif->get_r_color(), dif->get_g_color(),
                                                        dif->get_b_color())));

    int val = ui->lineEdit_2->text().split(" ")[0].toInt();
    max_v = val / 2;

    draw_elt();

    timer->start();

}
