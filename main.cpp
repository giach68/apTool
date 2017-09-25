#include "aptool.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    apTool w;
    w.show();

    return a.exec();
}
