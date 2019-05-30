#include "mainwindow.h"
#include <QApplication>
#include <QLocale>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
	QLocale::setDefault(QLocale::c());	
	QLocale::system().toDouble();
    MainWindow w;
    w.show();

    return a.exec();
}
