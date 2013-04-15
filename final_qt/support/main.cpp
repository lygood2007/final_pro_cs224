/** main.cpp
 ** Brief: This is the source file of the main function, which defines the entrance of program.
 ** Project: large-scale fluids
 ** Date: 04/10/2013
 ** Member: Scott, Hobarts, Yan Li
 **/

#include <QtGui/QApplication>
#include "mainwindow.h"

#include <iostream>

using namespace std;

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;
    w.show();

//    w.setWindowState(w.windowState() | Qt::WindowFullScreen);
    
    return a.exec();

}
