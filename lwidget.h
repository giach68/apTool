#ifndef LWIDGET_H
#define LWIDGET_H
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <QWidget>
#include <QObject>
#include <QLabel>

class LWidget : public QLabel
{
    Q_OBJECT

public:
    LWidget(QWidget *parent = 0);

protected:
    void mousePressEvent(QMouseEvent *event);

signals:
    void clicked(int lx, int ly);

};

#endif // LWIDGET_H
