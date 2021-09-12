#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QGraphicsScene>
#include <QGraphicsView>
#include <QColor>
#include <QPen>
#include <QBrush>
#include <QTimer>
#include <QPixmap>
#include <QImage>
#include <ui_mainwindow.h>
#include <memory.h>
#include <iostream>
#include <math.h>
#include <cstdlib>
#include <time.h>
#include <windows.h>
#include <list>
#include <thread>
#include <chrono>
#include "sphere.h"
#include "cylinder.h"
#include "prisma.h"
#include "point.h"

using namespace std;

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

 private:
    Ui::MainWindow *ui;
    QGraphicsScene *scene;
    QTimer *timer;
    QImage im, elect;    // массивы цветов для кадра с элт и с итоговым изображением

    int max_v = 10;
    sphere *sph;
    prisma *pr;
    cylinder *cyl;

    //int Rs  = 0, Gs = 0, Bs = 0; // Цвет сцены
    int Rs  = 250, Gs = 250, Bs = 150; // Цвет сцены
    int width = 650, height = 600;  //Ширина и высота сцены

    int *zbuf = nullptr, *cadr = nullptr;   // массивы z-буферов для элт и итогового изображения

    int a = 0, b = 0;      // уголы поворота камеры (вокруг z и вокруг x)

    int dx = 250, dy = 200;     // смещения на экране

    float xl = 0, yl = 150, zl = -250; // координаты источника

    std::list <sphere> sp_list;          // список электронов
    std::list <sphere> :: iterator it;   // итератор для списка электронов
    int num = 0;                    // необходимое количество генерируемых электронов
    int x_st = 80, y_st = 240, z_st = 0;    // точка генерации электронов

    std::list <cylinder> c_list;             // список катодов (цилиндров)
    std::list <cylinder> :: iterator iter;   // итератор на список цилиндров
    cylinder *elt;

    prisma *screen;
    prisma *dif;

    std::list <prisma> dif_list;             // список полосок решетки
    std::list <prisma> :: iterator dif_it;   // итератор на список диффракционной решетки

    std::list <point> point_list;             // список добавленных электронов
    std::list <point> :: iterator point_it;   // итератор на список добавленных электронов

    int time = 0, max_time = 50;   // текущее время анимации и максимальное время

    float matrix[4][4] = { {1, 0, 0 , 0},       // матрица для поворота камеры относительно текущей системы координат
                           {0, 1, 0, 0},
                           {0, 0, 1, 0},
                           {0, 0, 0, 1} };

private slots:
    void on_pushButton_clicked();   // функция анимации
    void standartdrawcircle(sphere s);  // функция отрисовки сферы
    void standartdrawprisma(prisma p);
    void standartdrawcarcascylinder(cylinder c);
    void standartdrawpoint(point p);
    void standartdrawcylinder(cylinder c);  // функция отрисовки цилиндра
    void standartdrawroot(int r1, int r2, int h, int xc, int yc, int zc, int Rm, int Gm, int Bm, int alpha); // функция отрисовки крышки цилиндра
    void draw();        // отрисовка итогового кадра
    void draw_elt();    // отрисовка кадра с фоном и элт
    //void drawpixel(int xc, int yc, int zc, int *z_b, QColor **col, int R, int G, int B, int alpha, int xr);
    //void stdrawcylinder(cylinder c, int *z_b, QColor **col);
    void animation();   // анимация изображения
    int check_electron(sphere s, int t);    // функция проверки электронов на существование
    void on_pushButton_5_clicked(); // какие-то повороты
    void on_pushButton_4_clicked();
    void on_pushButton_6_clicked();
    void on_pushButton_7_clicked();
    void on_pushButton_9_clicked();
    void on_pushButton_11_clicked();
    void on_pushButton_10_clicked();
    void on_pushButton_12_clicked();
    void on_pushButton_13_clicked();
    void on_pushButton_17_clicked();
    void on_pushButton_14_clicked();
    void on_pushButton_15_clicked();
    void on_pushButton_16_clicked();
    void show_sphere();
    void on_pushButton_18_clicked();
    void on_pushButton_8_clicked();
};

#endif // MAINWINDOW_H
