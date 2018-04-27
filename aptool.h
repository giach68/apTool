#ifndef APTOOL_H
#define APTOOL_H

#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include "imageview.h"

#include <QMainWindow>

namespace Ui {
class apTool;
}

class apTool : public QMainWindow
{
    Q_OBJECT

public:
    explicit apTool(QWidget *parent = 0);
    ~apTool();

    ImageView* iw;

private slots:
    void on_pushButton_clicked();

    void on_processButton_clicked();

    void on_showButton_clicked();

    void on_pushButton_2_clicked();

    void on_clearButton_clicked();

    void on_lWidget_clicked(int lx, int ly);

    void on_pushButton_4_clicked();

    void on_lxSpinBox_valueChanged(double i);

    void on_lySpinBox_valueChanged(double i);

    void on_pushButton_3_clicked();

    void on_pushButton_5_clicked();

    void showAP(cv::String windowName);

    void on_clearButton_2_clicked();

    void on_pushButton_6_clicked();

    void on_pushButton_7_clicked();

    void on_fullTestBut_clicked();

private:
    Ui::apTool *ui;

    QPixmap drawLWidget();
};

#endif // APTOOL_H
