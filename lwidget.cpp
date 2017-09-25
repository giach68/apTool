#include "lwidget.h"
#include <QMouseEvent>

LWidget::LWidget(QWidget *parent) : QLabel(parent) {}

void LWidget::mousePressEvent(QMouseEvent *event) {
    if(event->button()==Qt::LeftButton)
        emit clicked(event->x(), event->y());
}
